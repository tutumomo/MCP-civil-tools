import pyproj
import math
import logging  # 新增 logging 套件

# 設定 logging，讓錯誤訊息能正確輸出到伺服器日誌
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# --- 計算 UTM Zone ---
def get_utm_zone(longitude):
    return math.floor((longitude + 180) / 6) + 1

# --- 經緯度轉平面座標（根據不同 datum 判斷投影方式） ---
def latlon_to_projected(lat, lon, datum_key):
    try:
        if not (-90 <= lat <= 90):
            raise ValueError("緯度必須在 -90 到 90 度之間")
        if not (-180 <= lon <= 180):
            raise ValueError("經度必須在 -180 到 180 度之間")

        datum_epsg = DATUMS[datum_key]
        source_crs = pyproj.CRS(datum_epsg)

        if datum_key == 'TWD97':
            target_crs = pyproj.CRS("EPSG:3826")
            zone = "TM2-121"
            is_south = False
        else:
            zone = get_utm_zone(lon)
            is_south = lat < 0
            utm_epsg = 32700 + zone if is_south else 32600 + zone
            target_crs = pyproj.CRS.from_epsg(utm_epsg)

        transformer = pyproj.Transformer.from_crs(source_crs, target_crs, always_xy=True)
        x, y = transformer.transform(lon, lat)

        return x, y, zone, is_south

    except Exception as e:
        # 將錯誤訊息輸出到日誌，方便伺服器維運
        logger.error(f"轉換錯誤: {e}")
        return None, None, None, None


# --- 平面座標轉經緯度 ---
def projected_to_latlon(x, y, zone, is_south, datum_key):
    try:
        datum_epsg = DATUMS[datum_key]
        target_crs = pyproj.CRS(datum_epsg)

        if datum_key == 'TWD97':
            source_crs = pyproj.CRS("EPSG:3826")
        elif datum_key == 'WGS 84':
            if not isinstance(zone, int) or not (1 <= zone <= 60):
                raise ValueError("UTM Zone 必須是 1 到 60 之間的整數")
            utm_epsg = 32700 + zone if is_south else 32600 + zone
            source_crs = pyproj.CRS.from_epsg(utm_epsg)
        else:
            raise ValueError(f"不支援的大地基準: {datum_key}")

        transformer = pyproj.Transformer.from_crs(source_crs, target_crs, always_xy=True)
        lon, lat = transformer.transform(x, y)

        return lon, lat

    except Exception as e:
        # 將錯誤訊息輸出到日誌，方便伺服器維運
        logger.error(f"反向轉換錯誤: {e}")
        return None, None

# --- 定義常用大地基準 ---
DATUMS = {
    "WGS 84": "EPSG:4326",
    "TWD97": "EPSG:3824",
}

# --- 曼寧係數查詢（簡化範例，實際可擴充）---
MANNING_N_TABLE = {
    "混凝土光滑": 0.012,
    "混凝土粗糙": 0.017,
    "土渠整修": 0.022,
    "天然土渠": 0.030,
    "碎石河道": 0.035,
}

def query_manning_n(material: str):
    n = MANNING_N_TABLE.get(material)
    if n is not None:
        return material, n
    else:
        return None, None

# --- 土壓力係數計算 ---
def active_earth_pressure_coefficient(phi_deg: float) -> float:
    """主動土壓力係數 Ka (庫倫公式)"""
    import math
    phi = math.radians(phi_deg)
    return round((1 - math.sin(phi)) / (1 + math.sin(phi)), 4)

def passive_earth_pressure_coefficient(phi_deg: float) -> float:
    """被動土壓力係數 Kp (庫倫公式)"""
    import math
    phi = math.radians(phi_deg)
    return round((1 + math.sin(phi)) / (1 - math.sin(phi)), 4)

# --- 排水溝流速計算（曼寧公式）---
def channel_flow_velocity(n: float, r: float, s: float) -> float:
    """曼寧公式計算流速 v = (1/n) * R^(2/3) * S^(1/2)"""
    import math
    return (1/n) * (r ** (2/3)) * (s ** 0.5)

def channel_flow_discharge(v: float, a: float) -> float:
    """流量 Q = v * A"""
    return v * a
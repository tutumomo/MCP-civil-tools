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

# --- 曼寧係數查詢（依據附件 ALL.md）---
MANNING_N_TABLE = {
    "純細砂": 0.025,
    "不緻密之細砂": 0.027,
    "粗砂及細砂": 0.028,
    "平常砂土": 0.030,
    "砂質壤土": 0.032,
    "堅壤土及黏質壤土": 0.035,
    "平常礫土": 0.035,
    "全面密草生": 0.040,
    "粗礫、石礫及砂礫": 0.040,
    "礫岩、硬頁岩、軟岩、水成岩": 0.045,
    "硬岩": 0.050,
    "混凝土": 0.012,
}

# 最大容許流速表（依據附件 ALL.md）
MAX_VELOCITY_TABLE = {
    "純細砂": (0.23, 0.30),
    "不緻密之細砂": (0.30, 0.46),
    "粗砂及細砂": (0.46, 0.61),
    "平常砂土": (0.61, 0.76),
    "砂質壤土": (0.76, 0.84),
    "堅壤土及黏質壤土": (0.91, 1.14),
    "平常礫土": (1.23, 1.52),
    "全面密草生": (1.50, 2.50),
    "粗礫、石礫及砂礫": (1.52, 1.83),
    "礫岩、硬頁岩、軟岩、水成岩": (1.83, 2.44),
    "硬岩": (3.05, 4.57),
    "混凝土": (4.57, 6.10),
}

def query_manning_n(material: str):
    """
    回傳材料名稱、曼寧係數n、最大容許流速範圍。
    """
    n = MANNING_N_TABLE.get(material)
    v_range = MAX_VELOCITY_TABLE.get(material)
    if n is not None and v_range is not None:
        return material, n, v_range
    elif n is not None:
        return material, n, None
    else:
        return None, None, None

def get_max_velocity(material: str):
    """
    查詢最大容許流速範圍。
    """
    return MAX_VELOCITY_TABLE.get(material)

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

# --- 排水溝流速計算（曼寧公式，含流速檢核）---
def channel_flow_velocity(n: float, r: float, s: float, material: str = None):
    """
    曼寧公式計算流速 v = (1/n) * R^(2/3) * S^(1/2)
    若有指定材料名稱，則自動檢核最大容許流速。
    """
    import math
    v = (1/n) * (r ** (2/3)) * (s ** 0.5)
    warning = None
    if material:
        v_range = get_max_velocity(material)
        if v_range:
            v_min, v_max = v_range
            if v > v_max:
                warning = f"警告：計算流速 v = {v:.4f} m/s 已超過『{material}』最大容許流速 {v_max} m/s！"
    return v, warning

def channel_flow_discharge(v: float, a: float) -> float:
    """流量 Q = v * A"""
    return v * a

def slope_stability_safety_factor(slope: float, unit_weight: float = None, friction_angle: float = None, cohesion: float = None, water_table: float = 0, method: str = "簡化法", soil_type: str = None):
    """
    邊坡穩定安全係數計算（簡化法，支援土壤參數查表）。
    slope: 坡度（百分比）
    unit_weight: 土壤單位重（kN/m3）
    friction_angle: 摩擦角（°）
    cohesion: 凝聚力（kPa）
    water_table: 地下水位高（m）
    method: 計算方法（預設簡化法）
    soil_type: 土壤類型（可選，若有則自動查表）
    """
    # 1. 若有土壤類型，查表自動補齊參數
    if soil_type:
        key, uw, coh, fa, _, _ = material_parameter_query(soil_type)
        if unit_weight is None and uw: unit_weight = uw
        if cohesion is None and coh is not None: cohesion = coh
        if friction_angle is None and fa is not None: friction_angle = fa
    # 預設值
    if unit_weight is None: unit_weight = 18.0
    if friction_angle is None: friction_angle = 30.0
    if cohesion is None: cohesion = 0.0
    # 2. 計算坡角
    import math
    beta = math.atan(slope / 100)  # 坡度百分比轉坡角（弧度）
    # 3. 簡化法（平面滑動）安全係數公式
    # FS = (c + (γ * h * cosβ) * tanφ) / (γ * h * sinβ)
    h = 5.0  # 假設滑動面深度5m（可擴充為參數）
    γ = unit_weight
    φ = math.radians(friction_angle)
    c = cohesion
    # 地下水位修正（若有，則有效單位重減半）
    γ_eff = γ if water_table == 0 else γ * 0.5
    numerator = c + (γ_eff * h * math.cos(beta)) * math.tan(φ)
    denominator = γ_eff * h * math.sin(beta)
    FS = numerator / denominator if denominator > 0 else float('inf')
    is_pass = FS >= 1.5
    # 4. 說明
    msg = (
        f"簡化法FS={FS:.2f}（≧1.5 {'合格' if is_pass else '不合格'}）\n"
        f"參數：坡度={slope}%，單位重={γ}kN/m3，摩擦角={math.degrees(φ):.1f}°，凝聚力={c}kPa，滑動面深度={h}m，地下水位={water_table}m\n"
        f"公式：FS = (c + γh cosβ tanφ) / (γh sinβ)\n依據《水土保持技術規範》與常用設計手冊，建議安全係數FS≧1.5。"
    )
    return {
        'fs': FS,
        'is_pass': is_pass,
        'slope': slope,
        'unit_weight': γ,
        'friction_angle': math.degrees(φ),
        'cohesion': c,
        'water_table': water_table,
        'soil_type': soil_type,
        'method': method,
        'h': h,
        'message': msg
    }

# --- USLE 參數查表 ---
R_TABLE = {
    '台北市': 350, '新北市': 340, '桃園市': 320, '台中市': 370, '南投縣': 400, '高雄市': 300, '屏東縣': 320, '花蓮縣': 420, '台東縣': 410,
}
K_TABLE = {
    '黏土': 0.20, '砂土': 0.30, '壤土': 0.25, '礫石': 0.15, '一般黏土': 0.20, '砂壤土': 0.28,
}
C_TABLE = {
    '裸地': 1.0, '草地': 0.05, '農地': 0.3, '林地': 0.01, '水田': 0.02,
}
P_TABLE = {
    '無措施': 1.0, '等高耕作': 0.5, '梯田': 0.2, '覆蓋作物': 0.3,
}

def get_r_factor(region, rainfall=None):
    if region and region in R_TABLE:
        return R_TABLE[region], f"依據地區查表({region})"
    elif rainfall:
        return rainfall * 0.3, "依據年降雨量簡化推估"
    else:
        return 300, "預設值"

def get_k_factor(soil_type, soil_factor=None):
    if soil_type and soil_type in K_TABLE:
        return K_TABLE[soil_type], f"依據土壤類型查表({soil_type})"
    elif soil_factor is not None:
        return soil_factor, "由使用者輸入"
    else:
        return 0.25, "預設值"

def get_c_factor(land_use, cover_factor=None):
    if land_use and land_use in C_TABLE:
        return C_TABLE[land_use], f"依據土地利用查表({land_use})"
    elif cover_factor is not None:
        return cover_factor, "由使用者輸入"
    else:
        return 0.05, "預設值"

def get_p_factor(practice, practice_factor=None):
    if practice and practice in P_TABLE:
        return P_TABLE[practice], f"依據措施查表({practice})"
    elif practice_factor is not None:
        return practice_factor, "由使用者輸入"
    else:
        return 0.5, "預設值"

def soil_erosion_modulus(length, slope, rainfall, region=None, soil_type=None, land_use=None, practice=None, soil_factor=None, cover_factor=None, practice_factor=None, method='USLE'):
    """
    土壤侵蝕模數/流失量計算（USLE公式，支援地區、土壤、土地利用、措施查表）。
    """
    # 1. R 因子
    R, R_src = get_r_factor(region, rainfall)
    # 2. K 因子
    K, K_src = get_k_factor(soil_type, soil_factor)
    # 3. L 因子
    m = 0.5 if slope >= 5 else 0.2
    L = (length / 22.13) ** m
    # 4. S 因子
    import math
    theta = math.atan(slope / 100)
    S = 10.8 * math.sin(theta) + 0.03
    # 5. C 因子
    C, C_src = get_c_factor(land_use, cover_factor)
    # 6. P 因子
    P, P_src = get_p_factor(practice, practice_factor)
    # 7. USLE 計算
    A = R * K * L * S * C * P
    msg = (
        f"USLE公式：A = R*K*L*S*C*P\n"
        f"R(降雨侵蝕力)={R:.2f}（{R_src}），K(土壤可蝕性)={K:.3f}（{K_src}），L(坡長)={L:.3f}，S(坡度)={S:.3f}，C(覆蓋)={C:.3f}（{C_src}），P(水保)={P:.3f}（{P_src}）\n"
        f"依據美國USDA/NRCS與台灣水保局公開資料，參數可依地區/作物/措施調整。"
    )
    # 回傳所有參數明細
    return {
        'erosion_modulus': A,
        'soil_loss': A,
        'method': method,
        'region': region,
        'soil_type': soil_type,
        'land_use': land_use,
        'practice': practice,
        'R': R, 'R_src': R_src,
        'K': K, 'K_src': K_src,
        'L': L,
        'S': S,
        'C': C, 'C_src': C_src,
        'P': P, 'P_src': P_src,
        'message': msg
    }

# Rational 公式逕流係數查表
RUNOFF_COEFF_TABLE = {
    '都市': 0.8, '住宅區': 0.6, '工業區': 0.7, '道路': 0.85, '農地': 0.4, '林地': 0.2, '裸地': 0.5, '草地': 0.3,
}
# 地區降雨強度（mm/hr）簡化查表（假設10年重現期、60分鐘歷時）
RAINFALL_INTENSITY_TABLE = {
    '台北市': 110, '新北市': 105, '桃園市': 100, '台中市': 95, '南投縣': 120, '高雄市': 90, '屏東縣': 100, '花蓮縣': 130, '台東縣': 125,
}

def get_runoff_coeff(land_use, runoff_coeff=None):
    if land_use and land_use in RUNOFF_COEFF_TABLE:
        return RUNOFF_COEFF_TABLE[land_use], f"依據土地利用查表({land_use})"
    elif runoff_coeff is not None:
        return runoff_coeff, "由使用者輸入"
    else:
        return 0.5, "預設值"

def calc_idf_coeffs(P):
    """
    依據規範公式，計算A、B、C、G、H係數
    """
    A = (P / (189.96 + 0.31 * P)) ** 2
    B = 55
    C = (P / (381.71 + 1.45 * P)) ** 2
    G = (P / (42.89 + 1.33 * P)) ** 2
    H = (P / (65.33 + 1.836 * P)) ** 2
    return A, B, C, G, H

def calc_rainfall_intensity_by_spec(P, T, t):
    """
    依據規範公式，計算降雨強度I值
    P: 年平均降雨量(mm)
    T: 重現期(年)
    t: 集流時間/降雨延時(分鐘)
    """
    import math
    A, B, C, G, H = calc_idf_coeffs(P)
    I = (G + H * math.log(T)) ** (A / ((t + B) ** C))
    return I, {'A': A, 'B': B, 'C': C, 'G': G, 'H': H}

def get_rainfall_intensity(region=None, rainfall_intensity=None, return_period=None, duration=None, P=None):
    """
    若有P值，優先依規範公式推估I值，否則查表或用經驗公式。
    """
    if P is not None and return_period is not None and duration is not None:
        I, coeffs = calc_rainfall_intensity_by_spec(P, return_period, duration)
        src = f"依據規範公式(P={P}, T={return_period}, t={duration})，A={coeffs['A']:.3f},B={coeffs['B']},C={coeffs['C']:.3f},G={coeffs['G']:.3f},H={coeffs['H']:.3f}"
        return I, src
    if region and region in RAINFALL_INTENSITY_TABLE:
        return RAINFALL_INTENSITY_TABLE[region], f"依據地區查表({region})"
    elif rainfall_intensity is not None:
        return rainfall_intensity, "由使用者輸入"
    else:
        return 100, "預設值"

def catchment_peak_runoff(area, rainfall_intensity=None, runoff_coeff=None, land_use=None, region=None, return_period=None, duration=None, method="Rational", P=None):
    """
    集水區最大逕流量計算（Rational公式，支援土地利用、地區查表，並可依P自動推估I值）。
    """
    # 1. 逕流係數 C
    C, C_src = get_runoff_coeff(land_use, runoff_coeff)
    # 2. 降雨強度 I
    I, I_src = get_rainfall_intensity(region, rainfall_intensity, return_period, duration, P)
    # 3. 面積 A
    A = area  # ha
    # 4. Rational 公式 Q = C * I * A / 360
    Q = C * I * A / 360
    msg = (
        f"Rational公式：Q = C * I * A / 360\n"
        f"C(逕流係數)={C:.2f}（{C_src}），I(降雨強度)={I:.1f}mm/hr（{I_src}），A(面積)={A:.2f}ha\n"
        f"依據《水土保持技術規範》第17條與附件公式，參數可依土地利用、地區、重現期、歷時、年雨量調整。"
    )
    return {
        'peak_runoff': Q,
        'runoff_coeff': C,
        'runoff_coeff_src': C_src,
        'rainfall_intensity': I,
        'rainfall_intensity_src': I_src,
        'area': A,
        'method': method,
        'message': msg
    }

def retaining_wall_stability_check(height: float, thickness: float, unit_weight: float = None, friction_angle: float = None, cohesion: float = None, backfill_slope: float = 0, water_table: float = 0, soil_type: str = None):
    """
    護岸/擋土牆穩定檢核（滑動、傾倒、承載力，支援土壤參數查表）。
    height: 牆高(m)
    thickness: 牆底厚度(m)
    unit_weight: 回填土單位重(kN/m3)
    friction_angle: 回填土摩擦角(°)
    cohesion: 回填土凝聚力(kPa)
    backfill_slope: 回填坡度(°)
    water_table: 地下水位高(m)
    soil_type: 土壤類型（可選，若有則自動查表）
    """
    # 1. 若有土壤類型，查表自動補齊參數
    if soil_type:
        key, uw, coh, fa, _, _ = material_parameter_query(soil_type)
        if unit_weight is None and uw: unit_weight = uw
        if cohesion is None and coh is not None: cohesion = coh
        if friction_angle is None and fa is not None: friction_angle = fa
    # 預設值
    if unit_weight is None: unit_weight = 18.0
    if friction_angle is None: friction_angle = 30.0
    if cohesion is None: cohesion = 0.0
    # 幾何參數
    H = height
    B = thickness
    γ = unit_weight
    φ = friction_angle
    c = cohesion
    β = backfill_slope
    h_w = water_table
    # 2. 土壓力計算（簡化庫倫法，假設無黏聚力）
    import math
    φ_rad = math.radians(φ)
    β_rad = math.radians(β)
    K_a = (math.cos(β_rad) - math.sqrt(math.cos(β_rad)**2 - math.cos(φ_rad)**2)) / (math.cos(β_rad) + math.sqrt(math.cos(β_rad)**2 - math.cos(φ_rad)**2)) if β > 0 else (1 - math.sin(φ_rad)) / (1 + math.sin(φ_rad))
    P_a = 0.5 * γ * H**2 * K_a  # 單位長度主動土壓力
    # 3. 滑動檢核
    W = γ * B * H  # 牆重（假設單位長度，僅牆體）
    F_slide = (W * math.tan(φ_rad) + c * B) / P_a if P_a > 0 else float('inf')
    # 4. 傾倒檢核
    e = B/2 - (P_a * H/3) / W  # 牆底壓力中心偏移
    M_r = W * (B/2 - e)  # 抗傾倒力矩
    M_o = P_a * H/3  # 傾倒力矩
    F_overturn = M_r / M_o if M_o > 0 else float('inf')
    # 5. 承載力檢核
    q_max = W / B + P_a * H/3 / B  # 牆趾最大壓力
    q_allow = 200.0  # 假設基礎承載力標準值（kPa），可查表
    F_bearing = q_allow / q_max if q_max > 0 else float('inf')
    # 6. 合格判斷
    pass_slide = F_slide >= 1.5
    pass_overturn = F_overturn >= 2.0
    pass_bearing = F_bearing >= 2.5
    is_pass = pass_slide and pass_overturn and pass_bearing
    # 7. 說明
    msg = (
        f"滑動SF={F_slide:.2f}（≧1.5 {'合格' if pass_slide else '不合格'}），傾倒SF={F_overturn:.2f}（≧2.0 {'合格' if pass_overturn else '不合格'}），承載SF={F_bearing:.2f}（≧2.5 {'合格' if pass_bearing else '不合格'}）\n"
        f"土壤參數：單位重={γ}kN/m3，摩擦角={φ}°，凝聚力={c}kPa，Kₐ={K_a:.3f}\n"
        f"依據《水土保持技術規範》與常用設計手冊，建議安全係數：滑動≧1.5，傾倒≧2.0，承載≧2.5。"
    )
    return {
        'sf_slide': F_slide,
        'sf_overturn': F_overturn,
        'sf_bearing': F_bearing,
        'pass_slide': pass_slide,
        'pass_overturn': pass_overturn,
        'pass_bearing': pass_bearing,
        'is_pass': is_pass,
        'unit_weight': γ,
        'friction_angle': φ,
        'cohesion': c,
        'soil_type': soil_type,
        'height': H,
        'thickness': B,
        'backfill_slope': β,
        'water_table': h_w,
        'message': msg
    }

def vegetation_slope_suggestion(slope: float, soil_type: str, climate: str):
    """
    植生護坡設計建議。
    依據《水土保持技術規範》及附件，依坡度、土壤、氣候自動查表建議。
    """
    # 坡度分級
    if slope <= 5:
        grade = 1
    elif slope <= 15:
        grade = 2
    elif slope <= 30:
        grade = 3
    elif slope <= 40:
        grade = 4
    elif slope <= 55:
        grade = 5
    elif slope <= 100:
        grade = 6
    else:
        grade = 7

    # 建議工法與草種
    if grade <= 3:
        method = "噴播草皮或直接播種，分區分期施工"
        species = "百慕達草、狗牙根、地毯草等耐旱耐沖蝕草種"
        coverage = 90.0
        extra = ""
    elif grade <= 5:
        method = "噴播草皮+格框或土工網，分區分期施工"
        species = "百慕達草、狗牙根、地毯草等耐旱耐沖蝕草種"
        coverage = 90.0
        extra = "坡度較大時應加強排水設施與坡面穩定措施"
    else:
        method = "格框+噴播草皮或植生袋，分區分期施工"
        species = "百慕達草、狗牙根、地毯草等耐旱耐沖蝕草種，並可輔以灌木或喬木"
        coverage = 90.0
        extra = "坡度極大時須加強排水、坡面穩定與緩衝帶設置"

    # 註明依據條文
    regulation = (
        "依據《水土保持技術規範》第8、167、172條及附件坡度分級、樣區面積、覆蓋率等規定，"
        "建議坡面植生覆蓋率達90%，草種以百慕達草、狗牙根等耐旱耐沖蝕草種為主，"
        "坡度大於30%時應加強排水與坡面穩定措施，並分區分期施工以減少裸露面積。"
    )
    msg = f"{extra}\n{regulation}"
    return slope, soil_type, climate, method, species, coverage, msg

def material_parameter_query(material: str):
    """
    常用材料設計參數查詢。
    """
    # 常用材料參數表（依據土木工程常用標準值與規範彙整）
    MATERIAL_PARAMETER_TABLE = {
        # 土壤類
        "一般黏土": {"unit_weight": 18.0, "cohesion": 20.0, "friction_angle": 25.0, "strength": 200.0, "message": "依據常用土壤工程手冊"},
        "砂土": {"unit_weight": 17.5, "cohesion": 0.0, "friction_angle": 32.0, "strength": 150.0, "message": "依據常用土壤工程手冊"},
        "礫石": {"unit_weight": 20.0, "cohesion": 0.0, "friction_angle": 38.0, "strength": 250.0, "message": "依據常用土壤工程手冊"},
        # 混凝土類
        "混凝土": {"unit_weight": 24.0, "cohesion": None, "friction_angle": None, "strength": 21000.0, "message": "f'c=210kg/cm²，依據混凝土設計規範"},
        # 石材類
        "花崗岩": {"unit_weight": 26.0, "cohesion": None, "friction_angle": None, "strength": 100000.0, "message": "依據岩石工程手冊"},
        "石灰岩": {"unit_weight": 25.0, "cohesion": None, "friction_angle": None, "strength": 80000.0, "message": "依據岩石工程手冊"},
        # 其他常見材料
        "頁岩": {"unit_weight": 23.0, "cohesion": None, "friction_angle": None, "strength": 40000.0, "message": "依據岩石工程手冊"},
        "砂岩": {"unit_weight": 22.0, "cohesion": None, "friction_angle": None, "strength": 35000.0, "message": "依據岩石工程手冊"},
    }
    # 支援模糊查詢
    key = None
    for k in MATERIAL_PARAMETER_TABLE.keys():
        if k in material or material in k:
            key = k
            break
    if key:
        param = MATERIAL_PARAMETER_TABLE[key]
        return (key, param["unit_weight"], param["cohesion"], param["friction_angle"], param["strength"], param["message"])
    else:
        return (material, None, None, None, None, "查無此材料，請確認名稱或參考規範。")

# 坡面保護工法建議查表
SLOPE_PROTECTION_TABLE = [
    # (坡度上限, 建議工法, 適用說明)
    (5,   "草皮或直接播種", "坡度平緩，適合草皮、噴播或直接播種。"),
    (15,  "噴播草皮+分區施工", "坡度較緩，建議分區分期噴播草皮，必要時加覆蓋物。"),
    (30,  "噴播草皮+格框/土工網", "坡度中等，建議噴播草皮並加格框或土工網。"),
    (40,  "格框+噴播草皮/土工網", "坡度較大，建議格框、土工網加強，並輔以噴播草皮。"),
    (55,  "格框+混凝土/石籠+草皮", "坡度陡峭，建議格框、混凝土或石籠加草皮，並加強排水。"),
    (100, "混凝土/石籠+格框+植生袋", "坡度極大，建議混凝土、石籠、格框與植生袋組合，並設緩衝帶。"),
    (999, "結構型護坡（混凝土、石籠）", "超陡坡，需採結構型護坡並加強排水與穩定。"),
]

# 土壤鬆散、降雨量大時加強建議
def slope_protection_suggestion(slope: float, soil_type: str = None, rainfall: float = None, region: str = None):
    """
    坡面保護工法建議（依坡度分級、土壤、降雨/地區自動查表）。
    """
    # 坡度分級查表
    for upper, method, desc in SLOPE_PROTECTION_TABLE:
        if slope <= upper:
            suggested_method = method
            method_desc = desc
            break
    else:
        suggested_method = "結構型護坡"
        method_desc = "坡度極大，需採結構型護坡。"

    # 土壤鬆散、降雨量大時加強建議
    extra = ""
    loose_soil = soil_type and ("砂" in soil_type or "鬆" in soil_type)
    high_rain = False
    if rainfall and rainfall > 2500:
        high_rain = True
    if region:
        # 若有地區，嘗試查表取得年雨量
        R_TABLE = {
            '台北市': 2400, '新北市': 2500, '桃園市': 2100, '台中市': 2100, '南投縣': 2800, '高雄市': 2200, '屏東縣': 2300, '花蓮縣': 3000, '台東縣': 2900,
        }
        rain = R_TABLE.get(region)
        if rain and rain > 2500:
            high_rain = True
    if loose_soil and high_rain:
        extra = "土壤鬆散且降雨量大，建議加強排水與結構型工法。"
    elif loose_soil:
        extra = "土壤鬆散，建議加強格框、土工網或結構型工法。"
    elif high_rain:
        extra = "降雨量大，建議加強排水與穩定措施。"

    # 依據說明
    regulation = (
        "依據《水土保持技術規範》及附件坡度分級、常用工法建議，坡度愈大、土壤愈鬆、降雨愈多，應採更強化之坡面保護工法，並加強排水與分區分期施工。"
    )
    msg = f"{method_desc} {extra}\n{regulation}"
    return slope, soil_type, rainfall, suggested_method, msg

# 土壤滲透係數查表（m/hr）
SOIL_K_TABLE = {
    '砂土': 1.0,
    '壤土': 0.2,
    '黏土': 0.01,
    '礫石': 5.0,
    '砂壤土': 0.5,
    '一般黏土': 0.01,
}

def get_soil_k(soil_type, k=None):
    if soil_type and soil_type in SOIL_K_TABLE:
        return SOIL_K_TABLE[soil_type], f"依據土壤類型查表({soil_type})"
    elif k is not None:
        return k, "由使用者輸入"
    else:
        return 0.2, "預設值"

def infiltration_facility_design(facility_type: str, k: float = None, area: float = None, rainfall: float = None, soil_type: str = None):
    """
    滲水設施設計（支援多型式、土壤查表、流量自動計算與尺寸建議）。
    facility_type: 設施型式（滲井、滲渠、滲溝、滲透池）
    k: 土壤滲透係數（m/hr）
    area: 集水面積（m2）
    rainfall: 設計降雨量（mm）
    soil_type: 土壤類型（可選，若有則自動查表）
    """
    # 1. 土壤滲透係數查表
    K, K_src = get_soil_k(soil_type, k)
    # 2. 設計流量 Q = A * P / 3600 (m3/hr)，P=降雨量(mm)
    if area and rainfall:
        Q = area * rainfall / 1000 / 1  # m3/次（假設一次降雨）
        Q_hr = area * rainfall / 1000  # m3/次（同上，1小時內排除）
        Q_msg = f"依集水面積{area}m2與降雨量{rainfall}mm計算設計流量Q={Q_hr:.2f}m3/次"
    else:
        Q_hr = 1.0
        Q_msg = "預設設計流量1.0m3/次"
    # 3. 設施尺寸建議（簡化公式，依型式）
    if facility_type == '滲井':
        # 滲井設計：Q = π*D*H*K，假設D=1m，H=1.5m
        D = 1.0
        H = 1.5
        q_infil = 3.14 * D * H * K  # m3/hr
        n = max(1, int(Q_hr / q_infil) + 1)
        size = f"直徑{D}m,深{H}m,數量{n}口"
        size_msg = f"每口滲井可處理{q_infil:.2f}m3/hr，建議設置{n}口。"
    elif facility_type == '滲渠':
        # 滲渠設計：Q = 寬*深*長*K，假設寬0.5m,深1m,長5m
        W = 0.5
        H = 1.0
        L = 5.0
        q_infil = W * H * L * K  # m3/hr
        n = max(1, int(Q_hr / q_infil) + 1)
        size = f"寬{W}m,深{H}m,長{L}m,數量{n}條"
        size_msg = f"每條滲渠可處理{q_infil:.2f}m3/hr，建議設置{n}條。"
    elif facility_type == '滲溝':
        # 滲溝設計：Q = 寬*深*長*K，假設寬0.3m,深0.5m,長10m
        W = 0.3
        H = 0.5
        L = 10.0
        q_infil = W * H * L * K
        n = max(1, int(Q_hr / q_infil) + 1)
        size = f"寬{W}m,深{H}m,長{L}m,數量{n}條"
        size_msg = f"每條滲溝可處理{q_infil:.2f}m3/hr，建議設置{n}條。"
    elif facility_type == '滲透池':
        # 滲透池設計：Q = 面積*深*K，假設面積10m2,深1.2m
        A = 10.0
        H = 1.2
        q_infil = A * H * K
        n = max(1, int(Q_hr / q_infil) + 1)
        size = f"面積{A}m2,深{H}m,數量{n}座"
        size_msg = f"每座滲透池可處理{q_infil:.2f}m3/hr，建議設置{n}座。"
    else:
        size = "型式不明，請輸入滲井、滲渠、滲溝或滲透池"
        size_msg = ""
    # 4. 說明
    msg = (
        f"{Q_msg}，{size_msg}\n土壤滲透係數K={K}m/hr（{K_src}）\n依據《水土保持技術規範》與常用設計手冊，建議設施型式、尺寸及數量可依實際需求調整。"
    )
    return {
        'facility_type': facility_type,
        'design_flow': Q_hr,
        'soil_k': K,
        'soil_k_src': K_src,
        'suggested_size': size,
        'area': area,
        'rainfall': rainfall,
        'soil_type': soil_type,
        'message': msg
    }

# IDF曲線查表（部分地區、重現期、歷時，mm/hr）
IDF_TABLE = {
    '台北市': {
        10: {10: 200, 30: 140, 60: 110, 120: 80, 180: 60},
        25: {10: 240, 30: 160, 60: 130, 120: 95, 180: 70},
        50: {10: 270, 30: 180, 60: 150, 120: 110, 180: 80},
    },
    '新北市': {
        10: {10: 190, 30: 135, 60: 105, 120: 75, 180: 55},
        25: {10: 230, 30: 155, 60: 125, 120: 90, 180: 65},
        50: {10: 260, 30: 170, 60: 140, 120: 105, 180: 75},
    },
    '台中市': {
        10: {10: 170, 30: 120, 60: 95, 120: 70, 180: 50},
        25: {10: 210, 30: 140, 60: 120, 120: 85, 180: 60},
        50: {10: 240, 30: 160, 60: 135, 120: 100, 180: 70},
    },
}
# 若無查表資料，採經驗公式 I = a/(t+b)^n
IDF_COEFF = {
    'default': {'a': 800, 'b': 15, 'n': 0.7},
    '台北市': {'a': 900, 'b': 12, 'n': 0.72},
    '新北市': {'a': 850, 'b': 13, 'n': 0.71},
    '台中市': {'a': 700, 'b': 15, 'n': 0.68},
}

def idf_curve_query(location: str, return_period: float, duration: float):
    """
    降雨強度-歷時-頻率（IDF）曲線查詢（支援地區/重現期/歷時查表與推估）。
    location: 地區名稱
    return_period: 重現期（年）
    duration: 歷時（分鐘）
    """
    loc = location
    rp = int(round(return_period))
    dur = int(round(duration))
    # 1. 查表
    intensity = None
    src = ""
    if loc in IDF_TABLE and rp in IDF_TABLE[loc]:
        # 找最接近的歷時
        durations = sorted(IDF_TABLE[loc][rp].keys())
        closest = min(durations, key=lambda x: abs(x-dur))
        intensity = IDF_TABLE[loc][rp][closest]
        src = f"依據IDF查表({loc}, {rp}年, {closest}分)"
    # 2. 若無查表，採經驗公式
    if intensity is None:
        coeff = IDF_COEFF.get(loc, IDF_COEFF['default'])
        a, b, n = coeff['a'], coeff['b'], coeff['n']
        intensity = a / ((dur + b) ** n)
        src = f"依據經驗公式I=a/(t+b)^n，a={a},b={b},n={n}"
    msg = (
        f"IDF曲線查詢：地區={loc}，重現期={rp}年，歷時={dur}分鐘\n"
        f"降雨強度={intensity:.1f} mm/hr（{src}）\n"
        f"依據《水土保持技術規範》與各地IDF曲線，若無查表資料則採經驗公式推估。"
    )
    return {
        'location': loc,
        'return_period': rp,
        'duration': dur,
        'rainfall_intensity': intensity,
        'source': src,
        'message': msg
    }
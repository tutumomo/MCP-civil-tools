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

def get_manning_materials_list():
    """
    回傳所有支援的曼寧係數材料名稱清單。
    """
    return list(MANNING_N_TABLE.keys())

# 最大容許流速表（依據附件 ALL.md 與規範）
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

def get_max_velocity_materials_list():
    """
    回傳所有支援的最大容許流速材料名稱清單。
    """
    return list(MAX_VELOCITY_TABLE.keys())

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
def active_earth_pressure_coefficient(phi_deg: float):
    """主動土壓力係數 Ka (庫倫公式)"""
    import math
    phi = math.radians(phi_deg)
    ka = round((1 - math.sin(phi)) / (1 + math.sin(phi)), 4)
    formula = "Ka = (1 - sinφ) / (1 + sinφ)"
    calculation_steps = f"Ka = (1 - sin{phi_deg:.2f}°) / (1 + sin{phi_deg:.2f}°) = {ka}"
    regulation = "依據《水土保持技術規範》第117、118、164條"
    explanation = "主動土壓力常用庫倫公式計算，設計時應依規範選用正確內摩擦角。"
    input_params = {'phi': phi_deg}
    return {
        'Ka': ka,
        'regulation': regulation,
        'explanation': explanation,
        'input_params': input_params,
        'formula': formula,
        'calculation_steps': calculation_steps
    }

def passive_earth_pressure_coefficient(phi_deg: float):
    """被動土壓力係數 Kp (庫倫公式)"""
    import math
    phi = math.radians(phi_deg)
    kp = round((1 + math.sin(phi)) / (1 - math.sin(phi)), 4)
    formula = "Kp = (1 + sinφ) / (1 - sinφ)"
    calculation_steps = f"Kp = (1 + sin{phi_deg:.2f}°) / (1 - sin{phi_deg:.2f}°) = {kp}"
    regulation = "依據《水土保持技術規範》第117、118、164條"
    explanation = "被動土壓力常用庫倫公式計算，設計時應依規範選用正確內摩擦角。"
    input_params = {'phi': phi_deg}
    return {
        'Kp': kp,
        'regulation': regulation,
        'explanation': explanation,
        'input_params': input_params,
        'formula': formula,
        'calculation_steps': calculation_steps
    }

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
    formula = "v = (1/n) × R^(2/3) × S^(1/2)"
    calculation_steps = f"v = (1/{n:.3f}) × {r:.3f}^(2/3) × {s:.4f}^(1/2) = {v:.4f} m/s"
    regulation = "依據《水土保持技術規範》第19、85條及附件曼寧公式、最大容許流速表"
    explanation = "明渠流速常用曼寧公式計算，並應檢核材料最大容許流速以符合法規。"
    input_params = {
        'n': n,
        'r': r,
        's': s,
        'material': material
    }
    return {
        'velocity': v,
        'warning': warning,
        'regulation': regulation,
        'explanation': explanation,
        'input_params': input_params,
        'formula': formula,
        'calculation_steps': calculation_steps
    }

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
    formula = "FS = (c + γh cosβ tanφ) / (γh sinβ)"
    calculation_steps = (
        f"FS = ({c:.2f} + ({γ_eff:.2f} × {h:.2f} × cos{math.degrees(beta):.2f}° × tan{math.degrees(φ):.2f}°)) / "
        f"({γ_eff:.2f} × {h:.2f} × sin{math.degrees(beta):.2f}°) = {FS:.2f}"
    )
    regulation = "依據《水土保持技術規範》第152、154條及附件最小安全係數表"
    explanation = "坡地穩定設計須符合規範建議安全係數，常用簡化法計算平面滑動安全係數。"
    input_params = {
        'slope': slope,
        'unit_weight': unit_weight,
        'friction_angle': friction_angle,
        'cohesion': cohesion,
        'water_table': water_table,
        'method': method,
        'soil_type': soil_type
    }
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
        'message': msg,
        'regulation': regulation,
        'explanation': explanation,
        'input_params': input_params,
        'formula': formula,
        'calculation_steps': calculation_steps
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

def get_supported_regions():
    """
    回傳所有支援的地區名稱清單（R 因子/降雨/IDF/年雨量等）。
    """
    return list(R_TABLE.keys())

def get_supported_soil_types():
    """
    回傳所有支援的土壤類型清單（K 因子/土壤參數/滲透係數等）。
    """
    return list(K_TABLE.keys())

def get_supported_land_uses():
    """
    回傳所有支援的土地利用型態清單（C 因子/逕流係數等）。
    """
    return list(C_TABLE.keys())

def get_supported_practices():
    """
    回傳所有支援的水保措施清單（P 因子）。
    """
    return list(P_TABLE.keys())

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
    formula = "A = R × K × L × S × C × P"
    calculation_steps = (
        f"A = {R:.2f} × {K:.3f} × {L:.3f} × {S:.3f} × {C:.3f} × {P:.3f} = {A:.2f}"
    )
    regulation = "依據《水土保持技術規範》第35、92條及附件USLE公式、參數表"
    explanation = "USLE公式為國際通用土壤流失量推估方法，台灣水保規範明定可用於坡地土壤侵蝕評估。"
    input_params = {
        'length': length,
        'slope': slope,
        'rainfall': rainfall,
        'region': region,
        'soil_type': soil_type,
        'land_use': land_use,
        'practice': practice,
        'soil_factor': soil_factor,
        'cover_factor': cover_factor,
        'practice_factor': practice_factor,
        'method': method
    }
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
        'message': msg,
        'regulation': regulation,
        'explanation': explanation,
        'input_params': input_params,
        'formula': formula,
        'calculation_steps': calculation_steps
    }

# 逕流係數表（依據水土保持技術規範第18條）
RUNOFF_COEFF_TABLE = {
    '陡峻山地': {
        'min': 0.75,
        'max': 0.90,
        '開發中': 0.95
    },
    '山麓區': {
        'min': 0.70,
        'max': 0.80,
        '開發中': 0.90
    },
    '丘陵地或森林': {
        'min': 0.50,
        'max': 0.75,
        '開發中': 0.90
    },
    '平坦耕地': {
        'min': 0.45,
        'max': 0.60,
        '開發中': 0.85
    },
    '非農業使用': {
        'min': 0.75,
        'max': 0.95,
        '開發中': 1.00
    }
}

def get_runoff_coeff(land_use, runoff_coeff=None, is_developing=False):
    """
    查詢逕流係數C值
    land_use: 土地利用類型或集水區狀況，必須是以下之一：
        - 陡峻山地
        - 山麓區
        - 丘陵地或森林
        - 平坦耕地
        - 非農業使用
    runoff_coeff: 直接指定的逕流係數值
    is_developing: 是否為開發中狀態
    """
    # 驗證land_use是否為有效值
    valid_land_uses = ['陡峻山地', '山麓區', '丘陵地或森林', '平坦耕地', '非農業使用']
    if land_use not in valid_land_uses:
        return None, f"土地利用類型必須是以下之一：{', '.join(valid_land_uses)}"
    
    # 如果直接指定逕流係數，直接返回
    if runoff_coeff is not None:
        return runoff_coeff, "由使用者輸入"
    
    # 開發中狀態直接返回1.0
    if is_developing:
        return 1.0, "開發中狀態，逕流係數以1.0計算"
    
    # 查表並計算平均值
    if land_use in RUNOFF_COEFF_TABLE:
        if is_developing:
            return RUNOFF_COEFF_TABLE[land_use]['開發中'], f"開發中狀態，{land_use}"
        else:
            # 計算範圍內的平均值
            min_val = RUNOFF_COEFF_TABLE[land_use]['min']
            max_val = RUNOFF_COEFF_TABLE[land_use]['max']
            avg = (min_val + max_val) / 2
            return avg, f"{land_use}逕流係數範圍：{min_val}~{max_val}，採用平均值{avg:.2f}"
    
    return None, f"找不到土地利用類型：{land_use}"

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
    elif rainfall_intensity is not None:
        return rainfall_intensity, "由使用者輸入"
    elif region and region in RAINFALL_INTENSITY_TABLE:
        return RAINFALL_INTENSITY_TABLE[region], f"依據地區查表({region})"
    else:
        # 使用經驗公式 I = a/(t+b)^n
        coeff = IDF_COEFF.get(region, IDF_COEFF['default'])
        a, b, n = coeff['a'], coeff['b'], coeff['n']
        I = a / ((duration + b) ** n)
        src = f"依據經驗公式I=a/(t+b)^n，a={a},b={b},n={n}"
        return I, src

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
    formula = "Q = C × I × A / 360"
    calculation_steps = f"Q = {C:.2f} × {I:.1f} × {A:.2f} / 360 = {Q:.4f} cms"
    regulation = "依據《水土保持技術規範》第16、17、18條及附件降雨強度、逕流係數表"
    explanation = "集水區面積小於1000公頃時，無實測資料可採Rational公式計算洪峰流量。"
    input_params = {
        'area': area,
        'rainfall_intensity': rainfall_intensity,
        'runoff_coeff': runoff_coeff,
        'land_use': land_use,
        'region': region,
        'return_period': return_period,
        'duration': duration,
        'method': method,
        'P': P
    }
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
        'message': msg,
        'regulation': regulation,
        'explanation': explanation,
        'input_params': input_params,
        'formula': formula,
        'calculation_steps': calculation_steps
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
    # 4. 傾倒檢核（修正）
    M_r = W * (B/2)  # 牆重對趾部的抗傾倒力矩
    M_o = P_a * H/3  # 主動土壓力對趾部的傾倒力矩
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
    formula = "滑動: SF = (W tanφ + cB) / Pa\n傾倒: SF = Mr / Mo\n承載: SF = q_allow / q_max"
    calculation_steps = (
        f"滑動: SF = ({W:.2f} × tan{φ:.2f}° + {c:.2f} × {B:.2f}) / {P_a:.2f} = {F_slide:.2f}\n"
        f"傾倒: SF = {M_r:.2f} / {M_o:.2f} = {F_overturn:.2f}\n"
        f"承載: SF = {q_allow:.2f} / {q_max:.2f} = {F_bearing:.2f}"
    )
    regulation = "依據《水土保持技術規範》第117、118、164條及附件最小安全係數表"
    explanation = "擋土牆設計須同時檢核滑動、傾倒、承載三項安全係數，並符合規範建議標準。"
    input_params = {
        'height': height,
        'thickness': thickness,
        'unit_weight': unit_weight,
        'friction_angle': friction_angle,
        'cohesion': cohesion,
        'backfill_slope': backfill_slope,
        'water_table': water_table,
        'soil_type': soil_type
    }
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
        'message': msg,
        'regulation': regulation,
        'explanation': explanation,
        'input_params': input_params,
        'formula': formula,
        'calculation_steps': calculation_steps
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
        return (material, None, None, None, None, f"查無此材料，支援查詢的材料有：{', '.join(MATERIAL_PARAMETER_TABLE.keys())}")

def get_supported_materials():
    """
    回傳所有支援的常用材料設計參數材料名稱清單。
    """
    MATERIAL_PARAMETER_TABLE = {
        # ... existing code ...
    }
    return list(MATERIAL_PARAMETER_TABLE.keys())

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

def get_supported_slope_protection_methods():
    """
    回傳所有支援的坡面保護工法建議清單。
    """
    return list(set([row[1] for row in SLOPE_PROTECTION_TABLE]))

# 土壤滲透係數查表（m/hr）
SOIL_K_TABLE = {
    '砂土': 1.0,
    '壤土': 0.2,
    '黏土': 0.01,
    '礫石': 5.0,
    '砂壤土': 0.5,
    '一般黏土': 0.01,
}

def get_supported_soil_k_types():
    """
    回傳所有支援的土壤滲透係數土壤類型清單。
    """
    return list(SOIL_K_TABLE.keys())

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

def get_supported_idf_locations():
    """
    回傳所有支援的IDF曲線地點清單。
    """
    return list(IDF_TABLE.keys())

# 若無查表資料，採經驗公式 I = a/(t+b)^n
IDF_COEFF = {
    'default': {'a': 800, 'b': 15, 'n': 0.7},
    '台北市': {'a': 900, 'b': 12, 'n': 0.72},
    '新北市': {'a': 850, 'b': 13, 'n': 0.71},
    '台中市': {'a': 700, 'b': 15, 'n': 0.68},
}

def idf_curve_query(location: str, return_period: float, duration: float, P: float = None):
    """
    降雨強度-歷時-頻率（IDF）曲線查詢（支援地區/重現期/歷時查表與推估，並可依P自動推估A、B、C、G、H）。
    location: 地區名稱
    return_period: 重現期（年）
    duration: 歷時（分鐘）
    P: 年平均降雨量（mm，可選，若有則依規範公式自動推估I值）
    """
    loc = location
    rp = int(round(return_period))
    dur = int(round(duration))
    # 1. 查表
    intensity = None
    src = ""
    detail = ""
    if loc in IDF_TABLE and rp in IDF_TABLE[loc]:
        # 找最接近的歷時
        durations = sorted(IDF_TABLE[loc][rp].keys())
        closest = min(durations, key=lambda x: abs(x-dur))
        intensity = IDF_TABLE[loc][rp][closest]
        src = f"依據IDF查表({loc}, {rp}年, {closest}分)"
    # 2. 若無查表，優先依P用規範公式
    if intensity is None and P is not None:
        I, coeffs = calc_rainfall_intensity_by_spec(P, return_period, duration)
        intensity = I
        src = f"依據規範公式(P={P}, T={return_period}, t={duration})"
        detail = (
            f"A={coeffs['A']:.3f}, B={coeffs['B']}, C={coeffs['C']:.3f}, G={coeffs['G']:.3f}, H={coeffs['H']:.3f}\n"
            f"I = (G + H * ln(T)) ^ (A / (t + B)^C)"
        )
    # 3. 若無P則用經驗公式
    if intensity is None:
        coeff = IDF_COEFF.get(loc, IDF_COEFF['default'])
        a, b, n = coeff['a'], coeff['b'], coeff['n']
        intensity = a / ((dur + b) ** n)
        src = f"依據經驗公式I=a/(t+b)^n，a={a},b={b},n={n}"
        detail = "I = a / (t + b)^n\nt 為歷時（分鐘）"
    msg = (
        f"IDF曲線查詢：地區={loc}，重現期={rp}年，歷時={dur}分鐘\n"
        f"降雨強度={intensity:.1f} mm/hr（{src}）\n"
        f"{detail}"
        f"依據《水土保持技術規範》與各地IDF曲線，若無查表資料則採規範公式或經驗公式推估。"
    )
    return {
        'location': loc,
        'return_period': rp,
        'duration': dur,
        'rainfall_intensity': intensity,
        'source': src,
        'message': msg
    }

def slope_protection_suggestion(slope: float, soil_type: str = None, rainfall: float = None, region: str = None):
    """
    坡面保護工法建議查詢（依坡度分級、土壤、降雨/地區自動查表）。
    參數：坡度、土壤、降雨、地區
    回傳：建議工法、說明
    """
    # 根據坡度查表
    method = ""
    for limit, m, desc in SLOPE_PROTECTION_TABLE:
        if slope <= limit:
            method = m
            break
    if not method:
        method = "結構型護坡（混凝土、石籠）"
        desc = "超陡坡，需採結構型護坡並加強排水與穩定。"
    # 組合說明
    msg = f"坡度={slope}%，土壤={soil_type or '-'}，降雨={rainfall or '-'}mm，建議工法：{method}。{desc}"
    regulation = "依據《水土保持技術規範》第8、167、172條及附件坡度分級、樣區面積、覆蓋率等規定"
    return slope, soil_type, rainfall, method, f"{msg}\n{regulation}"

def gabion_stability_check(
    height: float,
    width: float,
    wall_weight: float,
    phi: float,
    delta: float,
    theta: float,
    i: float,
    gamma: float,
    friction_coef: float,
    pressure_mode: str = "active"
) -> dict:
    """
    土石籠擋土牆穩定分析（主動/被動土壓力）。
    參數：
      - height: 土石籠高度 (m)
      - width: 土石籠寬度 (m)
      - wall_weight: 擋土牆總重 (kN/m)
      - phi: 土壤內摩擦角 (°)
      - delta: 牆背摩擦角 (°)
      - theta: 牆背傾斜角 (°)
      - i: 地表傾斜角 (°)
      - gamma: 土壤飽和單位重 (kN/m³)
      - friction_coef: 摩擦係數
      - pressure_mode: 土壓力模式 ("active" 或 "passive")
    回傳：dict，含 success, data, message, report（繁體中文完整報告書）
    """
    import math
    from datetime import datetime
    
    try:
        # 1. 計算土壓力係數
        phi_rad = math.radians(phi)
        theta_rad = math.radians(theta)
        delta_rad = math.radians(delta)
        i_rad = math.radians(i)
        
        if i == 0:
            # 簡化公式 (當地表水平時)
            if pressure_mode == "active":
                k = (1 - math.sin(phi_rad)) / (1 + math.sin(phi_rad))
                formula = "Ka = (1 - sinφ) / (1 + sinφ)"
            else:
                k = (1 + math.sin(phi_rad)) / (1 - math.sin(phi_rad))
                formula = "Kp = (1 + sinφ) / (1 - sinφ)"
        else:
            # 完整庫倫公式
            if pressure_mode == "active":
                numerator = math.cos(phi_rad - theta_rad)**2
                denominator = math.cos(theta_rad)**2 * math.cos(delta_rad + theta_rad)
                sqrt_part = math.sqrt(
                    (math.sin(phi_rad + delta_rad) * math.sin(phi_rad - i_rad)) / 
                    (math.cos(delta_rad + theta_rad) * math.cos(theta_rad - i_rad))
                )
                k = numerator / (denominator * (1 + sqrt_part)**2)
                formula = "Ka = cos²(φ-θ) / [cos²θ·cos(δ+θ)·(1+√Q)²], Q = [sin(φ+δ)·sin(φ-i)] / [cos(δ+θ)·cos(θ-i)]"
            else:
                numerator = math.cos(phi_rad + theta_rad)**2
                denominator = math.cos(theta_rad)**2 * math.cos(delta_rad - theta_rad)
                sqrt_part = math.sqrt(
                    (math.sin(phi_rad + delta_rad) * math.sin(phi_rad + i_rad)) / 
                    (math.cos(delta_rad - theta_rad) * math.cos(theta_rad - i_rad))
                )
                k = numerator / (denominator * (1 - sqrt_part)**2)
                formula = "Kp = cos²(φ+θ) / [cos²θ·cos(δ-θ)·(1-√Q)²], Q = [sin(φ+δ)·sin(φ+i)] / [cos(δ-θ)·cos(θ-i)]"
        
        # 2. 計算土壓力總和
        pa = 0.5 * gamma * height**2 * k
        
        # 3. 計算分力
        pv = pa * math.sin(delta_rad + theta_rad)
        ph = pa * math.cos(delta_rad + theta_rad)
        
        # 4. 計算力矩
        mr = wall_weight * width / 2
        mf = ph * height / 3
        
        # 5. 計算安全係數
        fs_overturning = mr / mf if mf != 0 else float('inf')
        fs_sliding = ((wall_weight + pv) * friction_coef) / ph if ph != 0 else float('inf')
        
        # 穩定性評估
        assessment = []
        if fs_overturning >= 1.5:
            assessment.append("抗傾覆: 安全")
        else:
            assessment.append("抗傾覆: 不安全!")
        
        if fs_sliding >= 1.5:
            assessment.append("抗滑動: 安全")
        else:
            assessment.append("抗滑動: 不安全!")
        
        # 報告書
        report = (
            f"【土石籠擋土牆{pressure_mode}土壓力分析計算書】\n"
            f"生成時間: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n"
            f"## 輸入參數\n"
            f"- 分析模式: {pressure_mode}土壓力分析\n"
            f"- 土壤內摩擦角 φ = {phi}°\n"
            f"- 牆背摩擦角 δ = {delta}°\n"
            f"- 牆背傾斜角 θ = {theta}°\n"
            f"- 地表傾斜角 i = {i}°\n"
            f"- 土壤飽和單位重 γ = {gamma} kN/m³\n"
            f"- 擋土牆總重 W = {wall_weight} kN/m\n"
            f"- 土石籠高度 H = {height} m\n"
            f"- 土石籠寬度 B = {width} m\n"
            f"- 摩擦係數 μ = {friction_coef}\n\n"
            f"## 計算公式\n{formula}\n\n"
            f"## 計算結果\n"
            f"- 土壓力係數 K{'a' if pressure_mode == 'active' else 'p'} = {k:.4f}\n"
            f"- 土壓力總和 P{'a' if pressure_mode == 'active' else 'p'} = {pa:.2f} kN/m\n"
            f"- 垂直分力 Pv = {pv:.2f} kN/m\n"
            f"- 水平分力 Ph = {ph:.2f} kN/m\n"
            f"- 扶正力矩 Mr = {mr:.2f} kN·m/m\n"
            f"- 傾倒力矩 Mf = {mf:.2f} kN·m/m\n"
            f"- **穩定性評估**: {' | '.join(assessment)}\n"
        )
        
        msg = f"土壓力係數: {k:.4f}, 土壓力總和: {pa:.2f} kN/m, 安全係數: 抗傾覆={fs_overturning:.2f}, 抗滑動={fs_sliding:.2f}"
        
        return {
            "success": True,
            "data": {
                "pressure_coef": k,
                "total_pressure": pa,
                "vertical_force": pv,
                "horizontal_force": ph,
                "restoring_moment": mr,
                "overturning_moment": mf,
                "fs_overturning": fs_overturning,
                "fs_sliding": fs_sliding,
                "assessment": assessment
            },
            "message": msg,
            "report": report
        }
        
    except Exception as e:
        return {
            "success": False,
            "message": f"計算錯誤: {str(e)}",
            "report": ""
        }

def u_channel_rebar_calculation(
    height: float,
    wall_slope: float,
    soil_slope: float,
    soil_angle: float,
    effective_depth: float,
    soil_weight: float = 18.0
) -> dict:
    """
    U型溝鋼筋量計算
    參數：
      - height: 溝高 (m)
      - wall_slope: 溝壁傾角 (m)
      - soil_slope: 土方傾角 (°)
      - soil_angle: 安息角 (°)
      - effective_depth: 有效厚度 (m)
      - soil_weight: 土重 (kN/m³)，預設 18.0
    回傳：dict，含計算結果與報告書
    """
    import math
    
    # 角度轉換
    wall_slope_rad = math.atan(wall_slope)
    soil_slope_rad = math.radians(soil_slope)
    soil_angle_rad = math.radians(soil_angle)
    
    # 計算土壓力係數 Ka
    numerator = math.sin(soil_angle_rad) * math.sin(soil_angle_rad - soil_slope_rad)
    denominator = math.cos(wall_slope_rad + soil_slope_rad) * math.cos(wall_slope_rad)
    
    if numerator / denominator < 0:
        return {
            "success": False,
            "message": "根據公式此資料不能計算",
            "report": "【U型溝鋼筋量計算報告】\n\n輸入參數：\n- 溝高 H = {:.3f} m\n- 溝壁傾角 m = {:.3f}\n- 土方傾角 i = {:.2f}°\n- 安息角 ψ = {:.2f}°\n- 有效厚度 d = {:.3f} m\n- 土重 γ = {:.1f} kN/m³\n\n計算結果：\n根據公式此資料不能計算，請檢查輸入參數是否合理。"
        }
    
    ka = (math.cos(soil_angle_rad + wall_slope_rad))**2 / (
        (math.cos(wall_slope_rad))**2 * 
        (1 + math.sqrt(numerator / denominator))**2
    )
    
    # 計算土壓力 P
    p = soil_weight * height**2 * ka / (2 * math.cos(wall_slope_rad))
    
    # 計算彎矩 M
    m = soil_weight * height**3 * ka / (6 * math.cos(wall_slope_rad))
    
    # 計算鋼筋量 As (使用容許應力 1400 kgf/cm² * 0.875)
    jfs = 1400 * 0.875  # kgf/cm²
    as_rebar = m / (jfs * effective_depth) * 1000000 / 1000  # cm²
    
    # 生成報告書
    report = (
        "【U型溝鋼筋量計算報告】\n\n"
        "輸入參數：\n"
        f"- 溝高 H = {height:.3f} m\n"
        f"- 溝壁傾角 m = {wall_slope:.3f}\n"
        f"- 土方傾角 i = {soil_slope:.2f}°\n"
        f"- 安息角 ψ = {soil_angle:.2f}°\n"
        f"- 有效厚度 d = {effective_depth:.3f} m\n"
        f"- 土重 γ = {soil_weight:.1f} kN/m³\n\n"
        "計算公式：\n"
        "1. 土壓力係數 Ka = cos²(ψ+m) / [cos²m·(1+√Q)²]\n"
        "   其中 Q = [sinψ·sin(ψ-i)] / [cos(m+i)·cosm]\n"
        "2. 土壓力 P = γ·H²·Ka / (2·cosm)\n"
        "3. 彎矩 M = γ·H³·Ka / (6·cosm)\n"
        "4. 鋼筋量 As = M / (fs·d) × 10⁶ / 1000\n\n"
        "計算結果：\n"
        f"- 土壓力係數 Ka = {ka:.4f}\n"
        f"- 土壓力 P = {p:.3f} kN/m\n"
        f"- 彎矩 M = {m:.3f} kN·m/m\n"
        f"- 鋼筋量 As = {as_rebar:.3f} cm²/m"
    )
    
    return {
        "success": True,
        "data": {
            "earth_pressure_coef": ka,
            "earth_pressure": p,
            "moment": m,
            "rebar_area": as_rebar
        },
        "message": f"土壓力係數 Ka = {ka:.4f}, 土壓力 P = {p:.3f} kN/m, 彎矩 M = {m:.3f} kN·m/m, 鋼筋量 As = {as_rebar:.3f} cm²/m",
        "report": report
    }

# 鋼筋規格資料表
REBAR_TABLE = {
    "#3": {"diameter": 9.53, "area": 0.71, "weight": 0.559, "perimeter": 29.9, "name": "10"},
    "#4": {"diameter": 12.7, "area": 1.27, "weight": 0.994, "perimeter": 39.9, "name": "13"},
    "#5": {"diameter": 15.9, "area": 1.98, "weight": 1.55, "perimeter": 49.9, "name": "16"},
    "#6": {"diameter": 19.1, "area": 2.85, "weight": 2.24, "perimeter": 60.0, "name": "19"},
    "#7": {"diameter": 22.2, "area": 3.88, "weight": 3.05, "perimeter": 69.7, "name": "22"},
    "#8": {"diameter": 25.4, "area": 5.07, "weight": 3.98, "perimeter": 79.8, "name": "25"},
    "#9": {"diameter": 28.7, "area": 6.45, "weight": 5.06, "perimeter": 90.2, "name": "29"},
    "#10": {"diameter": 32.2, "area": 8.17, "weight": 6.41, "perimeter": 101.2, "name": "32"},
    "#11": {"diameter": 35.8, "area": 10.08, "weight": 7.91, "perimeter": 112.5, "name": "36"},
    "#14": {"diameter": 43.0, "area": 14.52, "weight": 11.4, "perimeter": 135.1, "name": "43"},
    "#18": {"diameter": 57.3, "area": 25.81, "weight": 20.3, "perimeter": 180.0, "name": "57"}
}

def get_rebar_info(rebar_number: str) -> dict:
    """
    查詢鋼筋規格資料
    參數：
      - rebar_number: 鋼筋編號（如 "#3"）
    回傳：dict，含直徑、截面積、單位重量、周長
    """
    return REBAR_TABLE.get(rebar_number)

def get_all_rebar_numbers() -> list:
    """
    回傳所有可用的鋼筋編號
    """
    return list(REBAR_TABLE.keys())

def calculate_rebar_weight(rebar_number: str, length: float) -> float:
    """
    計算鋼筋重量
    參數：
      - rebar_number: 鋼筋編號（如 "#3"）
      - length: 長度（m）
    回傳：重量（kg）
    """
    rebar_info = get_rebar_info(rebar_number)
    if rebar_info:
        return rebar_info["weight"] * length
    return None
"""
MCP Civil Tools 伺服器

本伺服器提供座標轉換與多種常用土木工程計算工具，供 LLM 或 MCP 客戶端調用。
"""
from mcp.server.fastmcp import FastMCP
from util import (
    latlon_to_projected, projected_to_latlon, query_manning_n, 
    active_earth_pressure_coefficient, passive_earth_pressure_coefficient, 
    channel_flow_velocity, channel_flow_discharge, slope_stability_safety_factor, 
    soil_erosion_modulus, catchment_peak_runoff, retaining_wall_stability_check, 
    vegetation_slope_suggestion, material_parameter_query, idf_curve_query, 
    SLOPE_PROTECTION_TABLE, u_channel_rebar_calculation, get_rebar_info, 
    get_all_rebar_numbers, calculate_rebar_weight
)
from utm_types import UTMResult, LatLonResult, ErrorResponse, ManningNResult, EarthPressureResult, ChannelFlowResult, VegetationSlopeSuggestion, SoilErosionResult  # 導入自訂型別
from fastapi import Query
# 新增支援清單查詢函式
from util import get_manning_materials_list, get_max_velocity_materials_list, get_supported_materials, get_supported_regions, get_supported_soil_types, get_supported_land_uses, get_supported_practices, get_supported_runoff_land_uses, get_supported_slope_protection_methods, get_supported_soil_k_types, get_supported_idf_locations

# 建立 MCP 伺服器
mcp = FastMCP("MCP Civil Tools")

@mcp.tool()
def latlon_to_utm(latitude: float, longitude: float, datum: str = None) -> dict:
    """
    將經緯度轉換為 UTM 座標，只回傳 X,Y 字串（小數點下4位）。
    """
    try:
        if datum is None:
            datum = "TWD97"
        x, y, zone, is_south = latlon_to_projected(latitude, longitude, datum)
        if datum == "TWD97":
            zone = "TM2-121"
            is_south = False
        if x is None or y is None:
            return {"success": False, "data": None, "message": "錯誤: 座標轉換失敗，請檢查輸入值。"}
        # 四捨五入到小數點下4位，並以字串 "X,Y" 格式回傳
        xy_str = f"{x:.4f},{y:.4f}"
        return {"success": True, "data": xy_str, "message": ""}
    except Exception as e:
        return {"success": False, "data": None, "message": f"錯誤: {str(e)}"}

@mcp.tool()
def utm_to_latlon(easting: float, northing: float, zone: int = None, south: bool = None, datum: str = None) -> dict:
    """
    將 UTM 座標轉換為經緯度，只回傳 "緯度,經度" 字串（小數點下15位）。
    """
    try:
        if datum is None:
            datum = "TWD97"
        if zone is None:
            zone = "TM2-121" if datum == "TWD97" else 51
        if south is None:
            south = False
        lon, lat = projected_to_latlon(easting, northing, zone, south, datum)
        if lat is None or lon is None:
            return {"success": False, "data": None, "message": "錯誤: 座標轉換失敗，請檢查輸入值。"}
        latlon_str = f"{lat:.15f},{lon:.15f}"
        return {"success": True, "data": latlon_str, "message": ""}
    except Exception as e:
        return {"success": False, "data": None, "message": f"錯誤: {str(e)}"}

# --- 新增土木工程工具 ---
@mcp.tool()
def get_manning_n(material: str) -> dict:
    """
    查詢水土保持技術規範常見材料的曼寧係數與最大容許流速。
    """
    desc, n, v_range = query_manning_n(material)
    if n is not None:
        msg = f"{desc} 的曼寧係數 n = {n}"
        if v_range:
            msg += f"，最大容許流速範圍：{v_range[0]}~{v_range[1]} m/s"
        return {"success": True, "data": {"material": desc, "n": n, "v_range": v_range}, "message": msg}
    else:
        supported = get_manning_materials_list()
        return {"success": False, "data": supported, "message": "查無此材料，支援查詢的材料如下"}

@mcp.tool()
def calc_active_earth_pressure(phi: float) -> dict:
    """
    計算主動土壓力係數 Ka（輸入內摩擦角，單位：度）。
    """
    result = active_earth_pressure_coefficient(phi)
    msg = (
        f"主動土壓力係數 Ka = {result['Ka']}\n"
        f"依據: {result['regulation']}\n說明: {result['explanation']}\n"
        f"輸入參數: {result['input_params']}\n公式: {result['formula']}\n計算過程: {result['calculation_steps']}"
    )
    return {"success": True, "data": result, "message": msg}

@mcp.tool()
def calc_passive_earth_pressure(phi: float) -> dict:
    """
    計算被動土壓力係數 Kp（輸入內摩擦角，單位：度）。
    """
    result = passive_earth_pressure_coefficient(phi)
    msg = (
        f"被動土壓力係數 Kp = {result['Kp']}\n"
        f"依據: {result['regulation']}\n說明: {result['explanation']}\n"
        f"輸入參數: {result['input_params']}\n公式: {result['formula']}\n計算過程: {result['calculation_steps']}"
    )
    return {"success": True, "data": result, "message": msg}

@mcp.tool()
def calc_slope_stability(slope: float, unit_weight: float, friction_angle: float, cohesion: float, water_table: float = 0, method: str = "簡化法") -> dict:
    """
    邊坡穩定安全係數計算。
    """
    result = slope_stability_safety_factor(slope, unit_weight, friction_angle, cohesion, water_table, method)
    msg = (
        f"安全係數 = {result['fs']:.2f}，方法：{result['method']}，合格：{result['is_pass']}\n"
        f"依據: {result['regulation']}\n說明: {result['explanation']}\n"
        f"輸入參數: {result['input_params']}\n公式: {result['formula']}\n計算過程: {result['calculation_steps']}"
    )
    return {"success": True, "data": result, "message": msg}

@mcp.tool()
def calc_soil_erosion(
    length: float,
    slope: float,
    rainfall: float,
    region: str = None,
    soil_type: str = None,
    land_use: str = None,
    practice: str = None,
    soil_factor: float = None,
    cover_factor: float = None,
    practice_factor: float = None,
    method: str = "USLE"
) -> dict:
    """
    土壤侵蝕模數/流失量計算（支援地區、土壤、土地利用、措施查表）。
    """
    result = soil_erosion_modulus(
        length, slope, rainfall,
        region=region,
        soil_type=soil_type,
        land_use=land_use,
        practice=practice,
        soil_factor=soil_factor,
        cover_factor=cover_factor,
        practice_factor=practice_factor,
        method=method
    )
    msg = (
        f"土壤侵蝕模數/流失量 = {result['erosion_modulus']:.2f}\n"
        f"依據: {result['regulation']}\n說明: {result['explanation']}\n"
        f"輸入參數: {result['input_params']}\n公式: {result['formula']}\n計算過程: {result['calculation_steps']}"
    )
    return {"success": True, "data": result, "message": msg}

@mcp.tool()
def calc_catchment_runoff(
    area: float,
    rainfall_intensity: float = None,
    runoff_coeff: float = None,
    land_use: str = None,
    region: str = None,
    return_period: float = None,
    duration: float = None,
    method: str = "Rational",
    P: float = None
) -> dict:
    """
    集水區最大逕流量計算（支援土地利用、地區、重現期、歷時查表，及年平均降雨量P自動推估I值）。
    """
    result = catchment_peak_runoff(
        area,
        rainfall_intensity=rainfall_intensity,
        runoff_coeff=runoff_coeff,
        land_use=land_use,
        region=region,
        return_period=return_period,
        duration=duration,
        method=method,
        P=P
    )
    msg = (
        f"最大逕流量 Q = {result['peak_runoff']:.4f} cms\n"
        f"依據: {result['regulation']}\n說明: {result['explanation']}\n"
        f"輸入參數: {result['input_params']}\n公式: {result['formula']}\n計算過程: {result['calculation_steps']}"
    )
    return {"success": True, "data": result, "message": msg}

@mcp.tool()
def check_retaining_wall(
    height: float,
    thickness: float,
    unit_weight: float = None,
    friction_angle: float = None,
    cohesion: float = None,
    backfill_slope: float = 0,
    water_table: float = 0,
    soil_type: str = None
) -> dict:
    """
    護岸/擋土牆穩定檢核（滑動、傾倒、承載力，支援土壤參數查表）。
    """
    result = retaining_wall_stability_check(
        height,
        thickness,
        unit_weight=unit_weight,
        friction_angle=friction_angle,
        cohesion=cohesion,
        backfill_slope=backfill_slope,
        water_table=water_table,
        soil_type=soil_type
    )
    msg = (
        f"滑動SF={result['sf_slide']:.2f}，傾倒SF={result['sf_overturn']:.2f}，承載SF={result['sf_bearing']:.2f}，合格：{result['is_pass']}\n"
        f"依據: {result['regulation']}\n說明: {result['explanation']}\n"
        f"輸入參數: {result['input_params']}\n公式: {result['formula']}\n計算過程: {result['calculation_steps']}"
    )
    return {"success": True, "data": result, "message": msg}

@mcp.tool()
def suggest_vegetation_slope(slope: float, soil_type: str, climate: str) -> dict:
    """
    植生護坡設計建議。
    """
    s, st, c, method, species, coverage, msg = vegetation_slope_suggestion(slope, soil_type, climate)
    return {
        "success": True,
        "data": {
            "slope": s,
            "soil_type": st,
            "climate": c,
            "suggested_method": method,
            "suggested_species": species,
            "coverage": coverage
        },
        "message": msg
    }

@mcp.tool()
def query_material_parameter(material: str) -> dict:
    """
    查詢常用材料的設計參數（如單位重、凝聚力、摩擦角、強度等）。
    範例：
    - 查詢：壤土的材料參數？
    - 查詢：請給我一般黏土的設計參數
    回傳：材料名稱、單位重、凝聚力、摩擦角、強度等資訊。
    """
    m, uw, coh, fa, strength, msg = material_parameter_query(material)
    if uw is None and coh is None and fa is None and strength is None:
        return {"success": False, "data": None, "message": msg}
    return {
        "success": True,
        "data": {"material": m, "unit_weight": uw, "cohesion": coh, "friction_angle": fa, "strength": strength},
        "message": msg
    }

@mcp.tool()
def suggest_slope_protection(slope: float, soil_type: str = None, rainfall: float = None, region: str = None) -> dict:
    """
    坡面保護工法建議查詢（依坡度分級、土壤、降雨/地區自動查表）。
    """
    # 取得 method, desc
    method = ""
    desc = ""
    for limit, m, d in SLOPE_PROTECTION_TABLE:
        if slope <= limit:
            method = m
            desc = d
            break
    if not method:
        method = "結構型護坡（混凝土、石籠）"
        desc = "超陡坡，需採結構型護坡並加強排水與穩定。"
    msg = f"坡度={slope}%，土壤={soil_type or '-'}，降雨={rainfall or '-'}mm，建議工法：{method}。{desc}"
    regulation = "依據《水土保持技術規範》第8、167、172條及附件坡度分級、樣區面積、覆蓋率等規定"
    return {
        "success": True,
        "data": {
            "slope": slope,
            "soil_type": soil_type,
            "rainfall": rainfall,
            "suggested_method": method,
            "description": desc,
            "regulation": regulation
        },
        "message": f"{msg}\n{regulation}"
    }

@mcp.tool()
def design_infiltration_facility(
    facility_type: str,
    k: float = None,
    area: float = None,
    rainfall: float = None,
    soil_type: str = None
) -> dict:
    """
    滲水設施設計（支援多型式、土壤查表、流量自動計算與尺寸建議）。
    """
    result = infiltration_facility_design(
        facility_type,
        k=k,
        area=area,
        rainfall=rainfall,
        soil_type=soil_type
    )
    return {"success": True, "data": result, "message": ""}

@mcp.tool()
def query_idf_curve(
    location: str,
    return_period: float,
    duration: float
) -> dict:
    """
    降雨強度-歷時-頻率（IDF）曲線查詢（支援地區/重現期/歷時查表與推估）。
    """
    result = idf_curve_query(location, return_period, duration)
    return {"success": True, "data": result, "message": ""}

# --- 支援清單查詢 API ---
@mcp.tool()
def list_supported_materials() -> dict:
    """
    查詢所有支援的常用材料設計參數材料名稱清單。
    範例：
    - 查詢：有哪些常用材料？
    - 請列出可以查詢的材料種類
    回傳：材料名稱清單（如 ['一般黏土', '砂土', ...]）
    """
    return {"success": True, "data": get_supported_materials(), "message": ""}

@mcp.tool()
def list_supported_manning_materials() -> dict:
    """
    查詢目前支援的所有水溝鋪面/渠道材料的曼寧係數種類清單。
    範例：
    - 查詢：可以查詢哪些水溝鋪面的曼寧係數？
    - 有哪些渠道材料的曼寧係數可以查？
    回傳：材料名稱清單（如 ['純細砂', '混凝土', ...]）
    """
    return {"success": True, "data": get_manning_materials_list(), "message": ""}

@mcp.tool()
def list_supported_max_velocity_materials() -> dict:
    """
    查詢所有支援的最大容許流速材料名稱清單。
    範例：
    - 查詢：有哪些材料有最大容許流速？
    - 請列出可以查最大流速的材料
    回傳：材料名稱清單（如 ['混凝土', '全面密草生', ...]）
    """
    return {"success": True, "data": get_max_velocity_materials_list(), "message": ""}

@mcp.tool()
def list_supported_regions() -> dict:
    """
    查詢所有支援的地區名稱清單（R 因子/降雨/IDF/年雨量等）。
    範例：
    - 查詢：有哪些地區可以查詢？
    - 請列出支援的地區
    回傳：地區名稱清單（如 ['台北市', '新北市', ...]）
    """
    return {"success": True, "data": get_supported_regions(), "message": ""}

@mcp.tool()
def list_supported_soil_types() -> dict:
    """
    查詢所有支援的土壤類型清單（K 因子/土壤參數/滲透係數等）。
    範例：
    - 查詢：有哪些土壤類型？
    - 請列出可以查詢的土壤種類
    回傳：土壤類型清單（如 ['壤土', '砂土', ...]）
    """
    return {"success": True, "data": get_supported_soil_types(), "message": ""}

@mcp.tool()
def list_supported_land_uses() -> dict:
    """
    查詢所有支援的土地利用型態清單（C 因子/逕流係數等）。
    範例：
    - 查詢：有哪些土地利用型態？
    - 請列出可以查詢的土地利用
    回傳：土地利用型態清單（如 ['農業區', '都市區', ...]）
    """
    return {"success": True, "data": get_supported_land_uses(), "message": ""}

@mcp.tool()
def list_supported_practices() -> dict:
    """
    查詢所有支援的水保措施清單（P 因子）。
    範例：
    - 查詢：有哪些水保措施？
    - 請列出可以查詢的水保措施
    回傳：水保措施清單（如 ['等高耕作', '覆蓋作物', ...]）
    """
    return {"success": True, "data": get_supported_practices(), "message": ""}

@mcp.tool()
def list_supported_runoff_land_uses() -> dict:
    """
    查詢所有支援的土地利用型態清單（逕流係數）。
    範例：
    - 查詢：有哪些土地利用型態有逕流係數？
    - 請列出可以查逕流係數的土地利用
    回傳：土地利用型態清單（如 ['農業區', '都市區', ...]）
    """
    return {"success": True, "data": get_supported_runoff_land_uses(), "message": ""}

@mcp.tool()
def list_supported_slope_protection_methods() -> dict:
    """
    查詢所有支援的坡面保護工法建議清單。
    範例：
    - 查詢：有哪些坡面保護工法？
    - 請列出可以查詢的坡面保護工法
    回傳：工法名稱清單（如 ['草皮或直接播種', '噴播草皮+格框/土工網', ...]）
    """
    return {"success": True, "data": get_supported_slope_protection_methods(), "message": ""}

@mcp.tool()
def list_supported_soil_k_types() -> dict:
    """
    查詢所有支援的土壤滲透係數土壤類型清單。
    範例：
    - 查詢：有哪些土壤有滲透係數？
    - 請列出可以查滲透係數的土壤種類
    回傳：土壤類型清單（如 ['壤土', '砂土', ...]）
    """
    return {"success": True, "data": get_supported_soil_k_types(), "message": ""}

@mcp.tool()
def list_supported_idf_locations() -> dict:
    """
    查詢所有支援的IDF曲線地點清單。
    範例：
    - 查詢：有哪些IDF地點？
    - 請列出可以查詢的IDF曲線地點
    回傳：地點名稱清單（如 ['台北市', '新北市', ...]）
    """
    return {"success": True, "data": get_supported_idf_locations(), "message": ""}

# 新增：多斷面排水流速/流深/流量計算工具
@mcp.tool()
def calc_channel_section_flow(
    cross_section_type: str,
    flow: float,
    slope: float,
    manning_n: float,
    material: str = None,
    diameter: float = None,
    rect_width: float = None,
    rect_height: float = None,
    trap_bottom: float = None,
    trap_top: float = None,
    trap_height: float = None
) -> dict:
    """
    排水斷面流速/流深/流量計算（支援圓形、矩形、梯形，含流速檢核與完整報告）。
    參數：
      - cross_section_type: 斷面型式（圓形、矩形、梯形）
      - flow: 流量 Q (cms)
      - slope: 坡度 (%)
      - manning_n: 曼寧係數
      - material: 渠道材質（可選，檢核流速）
      - diameter: 管徑 (cm, 圓形用)
      - rect_width: 底寬 (cm, 矩形/梯形用)
      - rect_height: 高度 (cm, 矩形用)
      - trap_bottom: 底寬 (cm, 梯形用)
      - trap_top: 頂寬 (cm, 梯形用)
      - trap_height: 高度 (cm, 梯形用)
    回傳：dict，含 success, data, message, report（繁體中文完整報告書）
    """
    import math
    from math import sin, acos, sqrt
    # 單位轉換
    S = slope / 100
    n = manning_n
    Q = flow
    V = y = area = perimeter = R = None
    formula = ""
    calc_steps = ""
    section_desc = ""
    # 計算
    try:
        if cross_section_type == "圓形":
            if not diameter:
                return {"success": False, "message": "請輸入管徑 (cm)", "report": ""}
            D = diameter / 100
            epsilon = 1e-6
            y_low, y_high = 0.0001, D - epsilon
            y = (y_low + y_high) / 2
            for _ in range(100):
                theta = 2 * acos(1 - 2 * y / D)
                area = (D ** 2 / 8) * (theta - sin(theta))
                perimeter = (D / 2) * theta
                R = area / perimeter
                V = (1 / n) * (R ** (2 / 3)) * (S ** 0.5)
                Q_cal = area * V
                if abs(Q_cal - Q) < 1e-7:
                    break
                if Q_cal < Q:
                    y_low = y
                else:
                    y_high = y
                y = (y_low + y_high) / 2
            section_desc = f"圓形管徑 D={D*100:.1f}cm"
            formula = "Q = A × V, V = (1/n) × R^(2/3) × S^(1/2)"
            calc_steps = f"θ = 2×arccos(1-2y/D), A = (D²/8)(θ-sinθ), P = (D/2)θ, R = A/P, V = (1/n)R^(2/3)S^(1/2), Q = A×V"
        elif cross_section_type == "矩形":
            if not rect_width or not rect_height:
                return {"success": False, "message": "請輸入底寬與高度 (cm)", "report": ""}
            b = rect_width / 100
            h = rect_height / 100
            epsilon = 1e-6
            y_low, y_high = 0.0001, h - epsilon
            y = (y_low + y_high) / 2
            for _ in range(100):
                area = b * y
                perimeter = b + 2 * y
                R = area / perimeter
                V = (1 / n) * (R ** (2 / 3)) * (S ** 0.5)
                Q_cal = area * V
                if abs(Q_cal - Q) < 1e-7:
                    break
                if Q_cal < Q:
                    y_low = y
                else:
                    y_high = y
                y = (y_low + y_high) / 2
            section_desc = f"矩形底寬 b={b*100:.1f}cm, 高度 h={h*100:.1f}cm"
            formula = "Q = A × V, V = (1/n) × R^(2/3) × S^(1/2)"
            calc_steps = f"A = b×y, P = b+2y, R = A/P, V = (1/n)R^(2/3)S^(1/2), Q = A×V"
        elif cross_section_type == "梯形":
            if not trap_bottom or not trap_top or not trap_height:
                return {"success": False, "message": "請輸入底寬、頂寬、高度 (cm)", "report": ""}
            b1 = trap_bottom / 100
            b2 = trap_top / 100
            h = trap_height / 100
            epsilon = 1e-6
            y_low, y_high = 0.0001, h - epsilon
            y = (y_low + y_high) / 2
            for _ in range(100):
                area = (b1 + b2) * y / 2
                side = sqrt(y**2 + ((b2 - b1) / 2) ** 2)
                perimeter = b2 + 2 * side
                R = area / perimeter
                V = (1 / n) * (R ** (2 / 3)) * (S ** 0.5)
                Q_cal = area * V
                if abs(Q_cal - Q) < 1e-7:
                    break
                if Q_cal < Q:
                    y_low = y
                else:
                    y_high = y
                y = (y_low + y_high) / 2
            section_desc = f"梯形底寬 b1={b1*100:.1f}cm, 頂寬 b2={b2*100:.1f}cm, 高度 h={h*100:.1f}cm"
            formula = "Q = A × V, V = (1/n) × R^(2/3) × S^(1/2)"
            calc_steps = f"A = (b1+b2)×y/2, P = b2+2×√(y²+((b2-b1)/2)²), R = A/P, V = (1/n)R^(2/3)S^(1/2), Q = A×V"
        else:
            return {"success": False, "message": "不支援的斷面型式，請選擇圓形、矩形或梯形。", "report": ""}
        # 流速限制檢核
        check_msg = ""
        min_safe = max_safe = None
        if material:
            from util import get_max_velocity
            v_range = get_max_velocity(material)
            if v_range:
                min_safe, max_safe = v_range
                if V < min_safe:
                    check_msg = f"【檢核警告】計算流速 {V:.3f} m/s 低於安全下限 {min_safe:.2f} m/s，可能導致泥砂淤積。"
                elif V > max_safe:
                    check_msg = f"【檢核警告】計算流速 {V:.3f} m/s 超過安全上限 {max_safe:.2f} m/s，請考慮設置消能設施。"
                else:
                    check_msg = "計算結果符合安全流速規範。"
            else:
                check_msg = "無法取得該材質對應的流速限制。"
        # 滿流警告
        full_flow_warning = ""
        rel_tol = 1e-3
        if cross_section_type == "圓形":
            if (D - y) / D < rel_tol:
                full_flow_warning = "【檢核警告】計算流深已接近管徑，可能表示管道已滿流。"
        elif cross_section_type == "矩形":
            if (h - y) / h < rel_tol:
                full_flow_warning = "【檢核警告】計算流深已接近通道設計高度，可能表示通道已滿流。"
        elif cross_section_type == "梯形":
            if (h - y) / h < rel_tol:
                full_flow_warning = "【檢核警告】計算流深已接近通道設計高度，可能表示通道已滿流。"
        # 報告書
        report = (
            f"【排水斷面流速/流深/流量計算報告】\n"
            f"斷面型式：{cross_section_type}\n{section_desc}\n"
            f"流量 Q = {Q:.3f} cms\n坡度 S = {slope:.3f}%\n曼寧係數 n = {n:.3f}\n"
            f"{f'渠道材質：{material}' if material else ''}\n"
            f"\n【計算公式】\n{formula}\n【計算步驟】\n{calc_steps}\n"
            f"\n【計算結果】\n流速 V = {V:.3f} m/s\n流深 y = {y:.3f} m\n斷面積 A = {area:.4f} m²\n水力半徑 R = {R:.4f} m\n周長 P = {perimeter:.4f} m\n"
            f"\n{check_msg}\n{full_flow_warning}"
        )
        msg = f"流速: {V:.3f} m/s，流深: {y:.3f} m。{check_msg}"
        return {
            "success": True,
            "data": {
                "velocity": V,
                "flow_depth": y,
                "area": area,
                "perimeter": perimeter,
                "hydraulic_radius": R,
                "section_desc": section_desc,
                "check_msg": check_msg,
                "full_flow_warning": full_flow_warning
            },
            "message": msg,
            "report": report
        }
    except Exception as e:
        return {"success": False, "message": f"計算錯誤: {str(e)}", "report": ""}

@mcp.tool()
def check_gabion_stability(
    height: float,
    width: float,
    wall_weight: float,
    phi: float,
    delta: float,
    theta: float = 0,
    i: float = 0,
    gamma: float = 18,
    friction_coef: float = 0.5,
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
      - theta: 牆背傾斜角 (°)，預設 0
      - i: 地表傾斜角 (°)，預設 0
      - gamma: 土壤飽和單位重 (kN/m³)，預設 18
      - friction_coef: 摩擦係數，預設 0.5
      - pressure_mode: 土壓力模式 ("active" 或 "passive")，預設 "active"
    回傳：dict，含 success, data, message, report（繁體中文完整報告書）
    """
    result = gabion_stability_check(
        height=height,
        width=width,
        wall_weight=wall_weight,
        phi=phi,
        delta=delta,
        theta=theta,
        i=i,
        gamma=gamma,
        friction_coef=friction_coef,
        pressure_mode=pressure_mode
    )
    return result

@mcp.tool()
def calc_u_channel_rebar(
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
    try:
        result = u_channel_rebar_calculation(
            height=height,
            wall_slope=wall_slope,
            soil_slope=soil_slope,
            soil_angle=soil_angle,
            effective_depth=effective_depth,
            soil_weight=soil_weight
        )
        return result
    except Exception as e:
        return {
            "success": False,
            "message": f"計算過程中發生錯誤: {str(e)}",
            "report": f"【U型溝鋼筋量計算報告】\n\n計算過程中發生錯誤: {str(e)}\n請檢查輸入參數是否正確。"
        }

@mcp.tool()
def list_rebar_numbers() -> dict:
    """
    列出所有可用的鋼筋編號及其規格資料
    回傳：dict，含 success, data, message
    """
    try:
        rebar_numbers = get_all_rebar_numbers()
        rebar_data = {}
        for number in rebar_numbers:
            info = get_rebar_info(number)
            if info:
                rebar_data[number] = info
        return {
            "success": True,
            "data": rebar_data,
            "message": f"可用的鋼筋編號：{', '.join(rebar_numbers)}"
        }
    except Exception as e:
        return {
            "success": False,
            "message": f"查詢失敗：{str(e)}"
        }

@mcp.tool()
def get_rebar_specs(rebar_number: str) -> dict:
    """
    查詢鋼筋規格資料
    參數：
      - rebar_number: 鋼筋編號（如 "#3"）
    回傳：dict，含 success, data, message
    """
    try:
        rebar_info = get_rebar_info(rebar_number)
        if rebar_info:
            return {
                "success": True,
                "data": rebar_info,
                "message": f"鋼筋 {rebar_number} 規格：直徑 {rebar_info['diameter']}mm，截面積 {rebar_info['area']}cm²，單位重量 {rebar_info['weight']}kg/m，周長 {rebar_info['perimeter']}mm"
            }
        else:
            return {
                "success": False,
                "message": f"找不到鋼筋編號 {rebar_number}，可用的鋼筋編號：{', '.join(get_all_rebar_numbers())}"
            }
    except Exception as e:
        return {
            "success": False,
            "message": f"查詢失敗：{str(e)}"
        }

@mcp.tool()
def calculate_rebar_weight(rebar_number: str, length: float) -> dict:
    """
    計算鋼筋重量
    參數：
      - rebar_number: 鋼筋編號（如 "#3"）
      - length: 長度（m）
    回傳：dict，含 success, data, message
    """
    try:
        weight = calculate_rebar_weight(rebar_number, length)
        if weight is not None:
            return {
                "success": True,
                "data": {"weight": weight},
                "message": f"鋼筋 {rebar_number} 長度 {length}m 的重量為 {weight:.2f}kg"
            }
        else:
            return {
                "success": False,
                "message": f"找不到鋼筋編號 {rebar_number}，可用的鋼筋編號：{', '.join(get_all_rebar_numbers())}"
            }
    except Exception as e:
        return {
            "success": False,
            "message": f"計算失敗：{str(e)}"
        }

# 提供 ASGI 應用給 uvicorn 啟動 HTTP 服務
app = mcp.sse_app()  # 讓 uvicorn 可以直接啟動 HTTP 伺服器

# 若要用 CLI 方式啟動 stdio 服務，仍可保留以下程式
if __name__ == "__main__":
    mcp.run()
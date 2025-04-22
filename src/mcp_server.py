"""
MCP Civil Tools 伺服器

本伺服器提供座標轉換與多種常用土木工程計算工具，供 LLM 或 MCP 客戶端調用。
"""
from mcp.server.fastmcp import FastMCP
from util import latlon_to_projected, projected_to_latlon, query_manning_n, active_earth_pressure_coefficient, passive_earth_pressure_coefficient, channel_flow_velocity, channel_flow_discharge, slope_stability_safety_factor, soil_erosion_modulus, catchment_peak_runoff, retaining_wall_stability_check, vegetation_slope_suggestion, material_parameter_query, slope_protection_suggestion, infiltration_facility_design, idf_curve_query
from utm_types import UTMResult, LatLonResult, ErrorResponse, ManningNResult, EarthPressureResult, ChannelFlowResult, VegetationSlopeSuggestion, SoilErosionResult  # 導入自訂型別
from fastapi import Query

# 建立 MCP 伺服器
mcp = FastMCP("MCP Civil Tools")

@mcp.tool()
def latlon_to_utm(latitude: float, longitude: float, datum: str = None) -> str:
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
            return "錯誤: 座標轉換失敗，請檢查輸入值。"
        # 四捨五入到小數點下4位，並以字串 "X,Y" 格式回傳
        return f"{x:.4f},{y:.4f}"
    except Exception as e:
        return f"錯誤: {str(e)}"

@mcp.tool()
def utm_to_latlon(easting: float, northing: float, zone: int = None, south: bool = None, datum: str = None) -> str:
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
            return "錯誤: 座標轉換失敗，請檢查輸入值。"
        # 四捨五入到小數點下15位，並以字串 "緯度,經度" 格式回傳
        return f"{lat:.15f},{lon:.15f}"
    except Exception as e:
        return f"錯誤: {str(e)}"

# --- 新增土木工程工具 ---
@mcp.tool()
def get_manning_n(material: str) -> str:
    """
    查詢水土保持技術規範常見材料的曼寧係數與最大容許流速。
    """
    desc, n, v_range = query_manning_n(material)
    if n is not None:
        msg = f"{desc} 的曼寧係數 n = {n}"
        if v_range:
            msg += f"，最大容許流速範圍：{v_range[0]}~{v_range[1]} m/s"
        return msg
    else:
        return "查無此材料，請確認名稱。"

@mcp.tool()
def calc_active_earth_pressure(phi: float) -> str:
    """
    計算主動土壓力係數 Ka（輸入內摩擦角，單位：度）。
    """
    ka = active_earth_pressure_coefficient(phi)
    return f"主動土壓力係數 Ka = {ka}"

@mcp.tool()
def calc_passive_earth_pressure(phi: float) -> str:
    """
    計算被動土壓力係數 Kp（輸入內摩擦角，單位：度）。
    """
    kp = passive_earth_pressure_coefficient(phi)
    return f"被動土壓力係數 Kp = {kp}"

@mcp.tool()
def calc_channel_flow_velocity(n: float, r: float, s: float, material: str = None) -> str:
    """
    曼寧公式計算排水溝流速（n: 曼寧係數, r: 水力半徑(m), s: 坡度, material: 材料名稱，可自動檢核流速）。
    """
    v, warning = channel_flow_velocity(n, r, s, material)
    msg = f"流速 v = {v:.4f} m/s"
    if warning:
        msg += f"\n{warning}"
    return msg

@mcp.tool()
def calc_channel_flow_discharge(v: float, a: float) -> str:
    """
    計算排水溝流量（v: 流速(m/s), a: 斷面積(m2)）。
    """
    q = channel_flow_discharge(v, a)
    return f"流量 Q = {q:.4f} cms"

@mcp.tool()
def calc_slope_stability(slope: float, unit_weight: float, friction_angle: float, cohesion: float, water_table: float = 0, method: str = "簡化法") -> str:
    """
    邊坡穩定安全係數計算。
    """
    sf, m, is_pass, msg = slope_stability_safety_factor(slope, unit_weight, friction_angle, cohesion, water_table, method)
    return f"安全係數 = {sf:.2f}，方法：{m}，合格：{is_pass}\n{msg}"

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
) -> SoilErosionResult:
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
    return SoilErosionResult(
        erosion_modulus=result['erosion_modulus'],
        soil_loss=result['soil_loss'],
        method=result['method'],
        message=result['message']
    )

@mcp.tool()
def calc_catchment_runoff(
    area: float,
    rainfall_intensity: float = None,
    runoff_coeff: float = None,
    land_use: str = None,
    region: str = None,
    return_period: float = None,
    duration: float = None,
    method: str = "Rational"
) -> dict:
    """
    集水區最大逕流量計算（支援土地利用、地區、重現期、歷時查表）。
    """
    result = catchment_peak_runoff(
        area,
        rainfall_intensity=rainfall_intensity,
        runoff_coeff=runoff_coeff,
        land_use=land_use,
        region=region,
        return_period=return_period,
        duration=duration,
        method=method
    )
    return result

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
    return result

@mcp.tool()
def suggest_vegetation_slope(slope: float, soil_type: str, climate: str) -> VegetationSlopeSuggestion:
    """
    植生護坡設計建議。
    """
    s, st, c, method, species, coverage, msg = vegetation_slope_suggestion(slope, soil_type, climate)
    return VegetationSlopeSuggestion(
        slope=s,
        soil_type=st,
        climate=c,
        suggested_method=method,
        suggested_species=species,
        coverage=coverage,
        message=msg
    )

@mcp.tool()
def query_material_parameter(material: str) -> str:
    """
    常用材料設計參數查詢。
    """
    m, uw, coh, fa, strength, msg = material_parameter_query(material)
    return f"材料：{m}，單位重：{uw}kN/m3，凝聚力：{coh}kPa，摩擦角：{fa}°，強度：{strength}kPa\n{msg}"

@mcp.tool()
def suggest_slope_protection(slope: float, soil_type: str, rainfall: float) -> str:
    """
    坡面保護工法建議。
    """
    s, st, r, method, msg = slope_protection_suggestion(slope, soil_type, rainfall)
    return f"建議工法：{method}\n{msg}"

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
    return result

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
    return result

@mcp.post("/slope_protection_suggestion")
def api_slope_protection_suggestion(
    slope: float = Query(..., description="坡度 (%)"),
    soil_type: str = Query(None, description="土壤類型（如：壤土、砂土、黏土等，可選）"),
    rainfall: float = Query(None, description="年降雨量（mm，可選）"),
    region: str = Query(None, description="地區名稱（如：台北市，可選）")
):
    """
    坡面保護工法建議查詢（依坡度分級、土壤、降雨/地區自動查表）。
    - slope: 坡度（百分比）
    - soil_type: 土壤類型（可選）
    - rainfall: 年降雨量（mm，可選）
    - region: 地區名稱（可選，若有則自動查表年雨量）
    回傳：
    {
        "slope": 坡度,
        "soil_type": 土壤類型,
        "rainfall": 年降雨量,
        "suggested_method": 建議工法,
        "message": 詳細說明
    }
    """
    slope, soil_type, rainfall, method, msg = slope_protection_suggestion(slope, soil_type, rainfall, region)
    return {
        "slope": slope,
        "soil_type": soil_type,
        "rainfall": rainfall,
        "suggested_method": method,
        "message": msg
    }

# 提供 ASGI 應用給 uvicorn 啟動 HTTP 服務
app = mcp.sse_app()  # 讓 uvicorn 可以直接啟動 HTTP 伺服器

# 若要用 CLI 方式啟動 stdio 服務，仍可保留以下程式
if __name__ == "__main__":
    mcp.run()
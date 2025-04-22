"""
MCP Civil Tools 伺服器

本伺服器提供座標轉換與多種常用土木工程計算工具，供 LLM 或 MCP 客戶端調用。
"""
from mcp.server.fastmcp import FastMCP
from util import latlon_to_projected, projected_to_latlon, query_manning_n, active_earth_pressure_coefficient, passive_earth_pressure_coefficient, channel_flow_velocity, channel_flow_discharge
from utm_types import UTMResult, LatLonResult, ErrorResponse, ManningNResult, EarthPressureResult, ChannelFlowResult  # 導入自訂型別

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
    查詢水土保持技術規範常見材料的曼寧係數。
    """
    desc, n = query_manning_n(material)
    if n is not None:
        return f"{desc} 的曼寧係數 n = {n}"
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
def calc_channel_flow_velocity(n: float, r: float, s: float) -> str:
    """
    曼寧公式計算排水溝流速（n: 曼寧係數, r: 水力半徑(m), s: 坡度）。
    """
    v = channel_flow_velocity(n, r, s)
    return f"流速 v = {v:.4f} m/s"

@mcp.tool()
def calc_channel_flow_discharge(v: float, a: float) -> str:
    """
    計算排水溝流量（v: 流速(m/s), a: 斷面積(m2)）。
    """
    q = channel_flow_discharge(v, a)
    return f"流量 Q = {q:.4f} cms"

# 提供 ASGI 應用給 uvicorn 啟動 HTTP 服務
app = mcp.sse_app()  # 讓 uvicorn 可以直接啟動 HTTP 伺服器

# 若要用 CLI 方式啟動 stdio 服務，仍可保留以下程式
if __name__ == "__main__":
    mcp.run()
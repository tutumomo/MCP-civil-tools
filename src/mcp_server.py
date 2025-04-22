"""
MCP UTM 座標轉換伺服器

本伺服器將經緯度與UTM座標互轉功能，包裝為MCP工具，供LLM或MCP客戶端調用。
"""
from mcp.server.fastmcp import FastMCP
from utm_converter import latlon_to_projected, projected_to_latlon
from utm_types import UTMResult, LatLonResult, ErrorResponse  # 導入自訂型別

# 建立 MCP 伺服器
mcp = FastMCP("UTM Converter")

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

# 提供 ASGI 應用給 uvicorn 啟動 HTTP 服務
app = mcp.sse_app()  # 讓 uvicorn 可以直接啟動 HTTP 伺服器

# 若要用 CLI 方式啟動 stdio 服務，仍可保留以下程式
if __name__ == "__main__":
    mcp.run()
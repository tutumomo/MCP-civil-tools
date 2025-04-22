# This file is intentionally left blank.

from pydantic import BaseModel

# UTM 轉換結果資料結構
class UTMResult(BaseModel):
    easting: float  # 東移座標
    northing: float  # 北移座標
    zone: str  # UTM 分帶（數字或 TM2-121 字串）
    hemisphere: str  # 半球（North/South）

# 經緯度轉換結果資料結構
class LatLonResult(BaseModel):
    latitude: float  # 緯度
    longitude: float  # 經度

# 統一錯誤回傳格式
class ErrorResponse(BaseModel):
    error: str  # 錯誤訊息
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

# 曼寧係數查詢結果
class ManningNResult(BaseModel):
    description: str  # 材料描述
    n_value: float    # 曼寧係數

# 土壓力係數查詢結果
class EarthPressureResult(BaseModel):
    method: str       # 計算方法
    phi: float        # 內摩擦角
    k_value: float    # 土壓力係數

# 排水溝流速計算結果
class ChannelFlowResult(BaseModel):
    velocity: float   # 流速 (m/s)
    q: float          # 流量 (cms)
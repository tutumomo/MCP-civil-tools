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

# 邊坡穩定安全係數
class SlopeStabilityResult(BaseModel):
    safety_factor: float
    method: str
    is_pass: bool
    message: str

# 土壤侵蝕模數/流失量
class SoilErosionResult(BaseModel):
    erosion_modulus: float
    soil_loss: float
    method: str
    message: str

# 集水區逕流量
class RunoffResult(BaseModel):
    peak_runoff: float
    method: str
    message: str

# 護岸/擋土牆穩定檢核
class RetainingWallCheckResult(BaseModel):
    sliding_sf: float
    overturning_sf: float
    bearing_sf: float
    is_pass: bool
    message: str

# 植生護坡設計建議
class VegetationSlopeSuggestion(BaseModel):
    slope: float
    soil_type: str
    climate: str
    suggested_method: str
    suggested_species: str
    coverage: float
    message: str

# 常用材料設計參數查詢
class MaterialParameterResult(BaseModel):
    material: str
    unit_weight: float = None
    cohesion: float = None
    friction_angle: float = None
    strength: float = None
    message: str

# 坡面保護工法建議
class SlopeProtectionSuggestion(BaseModel):
    slope: float
    soil_type: str
    rainfall: float
    suggested_method: str
    message: str

# 滲水設施設計
class InfiltrationFacilityResult(BaseModel):
    facility_type: str
    design_flow: float
    suggested_size: str
    message: str

# IDF 曲線查詢
class IDFQueryResult(BaseModel):
    location: str
    return_period: float
    duration: float
    intensity: float
    message: str
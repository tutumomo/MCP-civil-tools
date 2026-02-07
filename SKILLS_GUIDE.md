# MCP Civil Tools 完整功能指南 (SKILLS_GUIDE)

本專案將《水土保持技術規範》與專業土木工程計算邏輯數位化，透過 MCP 協議提供 37 項核心工具。

## 1. 法規與知識檢索 (Regulations & Knowledge)
- **search_civil_regulations**: 智慧檢索《水土保持技術規範》。支援輸入條號（如「第18條」）或關鍵字搜尋。
- **list_supported_regions**: 列出支援 R 因子、IDF 或年雨量查詢的地區。
- **list_supported_materials**: 列出支援設計參數查詢的常用材料。
- **list_supported_manning_materials**: 列出支援曼寧係數查詢的材料清單。
- **list_supported_max_velocity_materials**: 列出支援最大容許流速檢核的材料。
- **list_supported_soil_types**: 列出支援 K 因子或物理參數查詢的土壤類型。
- **list_supported_land_uses**: 列出支援 C 因子或逕流係數查詢的土地利用型態。
- **list_supported_practices**: 列出支援 P 因子查詢的水保措施。
- **list_supported_slope_protection_methods**: 列出支援的坡面保護建議工法。
- **list_supported_soil_k_types**: 列出支援滲透係數查詢的土壤種類。
- **list_supported_idf_locations**: 列出支援 IDF 曲線查詢的地點清單。

## 2. 結構穩定分析 (Structural Stability)
- **check_retaining_wall**: 護岸/擋土牆穩定檢核。支援常時與地震時之滑動、傾倒、承載力檢核，並自動生成 Excalidraw 視覺化腳本。
- **check_gabion_stability**: 土石籠擋土牆穩定性分析，支援主動與被動土壓力模式。
- **calc_active_earth_pressure**: 計算主動土壓力係數 Ka (庫倫公式)。
- **calc_passive_earth_pressure**: 計算被動土壓力係數 Kp (庫倫公式)。
- **calc_slope_stability**: 邊坡穩定安全係數計算（平面滑動簡化法）。

## 3. 水理與水文計算 (Hydraulics & Hydrology)
- **calc_catchment_runoff**: 集水區最大逕流量計算（合理化公式 Rational Method）。
- **calc_channel_section_flow**: 排水斷面水理計算。支援矩形、圓形、梯形斷面，含出水高與流速規範檢核。
- **query_idf_curve**: 降雨強度-歷時-頻率 (IDF) 曲線查詢。
- **query_runoff_coeff**: 逕流係數 C 值查詢工具（依據規範第 18 條）。
- **design_infiltration_facility**: 滲水設施設計建議與設計流量計算。

## 4. 土壤侵蝕與護坡 (Soil Erosion & Protection)
- **calc_soil_erosion**: USLE 土壤流失量計算。
- **query_r_factor_tool**: 查詢降雨沖蝕指數 R 值（規範第 35 條）。
- **query_k_factor_tool**: 查詢土壤沖蝕指數 K 值（規範第 35 條）。
- **query_c_factor_tool**: 查詢覆蓋與管理因子 C 值（規範第 35 條）。
- **query_p_factor_tool**: 查詢水土保持處理因子 P 值（規範第 35 條）。
- **suggest_vegetation_slope**: 植生護坡設計建議、草種選取與覆蓋率檢核。
- **suggest_slope_protection**: 坡面保護工法建議查詢（依坡度分級）。

## 5. 基礎工具、材料與配筋 (Tools & Materials)
- **latlon_to_utm**: 經緯度座標轉換為 TWD97 (TM2) 或 UTM 座標。
- **utm_to_latlon**: TWD97 (TM2) 或 UTM 座標轉換為經緯度。
- **query_material_parameter**: 查詢常用材料之單位重、凝聚力、摩擦角與強度。
- **calc_u_channel_rebar**: U型溝鋼筋量智慧計算與配筋建議。
- **list_rebar_numbers**: 列出所有可用鋼筋編號 (#3 ~ #11)。
- **get_rebar_specs**: 查詢特定編號鋼筋之直徑、截面積、單位重與周長。
- **calculate_rebar_weight**: 根據編號與長度計算鋼筋總重。

---
*本指南由瑪麗蓮為 TOMO 總監特別整理。💋*

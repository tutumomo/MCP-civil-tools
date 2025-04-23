# MCP Civil Tools 伺服器

本專案是一個基於 MCP 協議的 Python 伺服器，提供經緯度與 UTM/TWD97 座標互轉，以及多種常用土木工程計算工具（目前主要針對水土保持技術規範提供如曼寧係數查詢、土壓力係數、排水溝流速、邊坡穩定、土壤侵蝕、逕流量、擋土牆檢核、植生建議、材料參數、坡面保護、滲水設施、IDF曲線等功能），適用於 LLM 工具、Claude Desktop 等 AI 應用整合。

---

## 目錄結構

```
MCP-civil-tools/
├── src/
│   ├── mcp_server.py         # MCP 伺服器主程式
│   ├── util.py              # 座標轉換與土木工程工具邏輯
│   └── utm_types/
│       └── __init__.py       # 型別定義
├── requirements.txt          # 依賴套件
├── README.md                 # 專案說明
└── .gitignore                # Git 忽略規則
```

---

## 安裝

1. 建議先建立虛擬環境：
   ```bash
   python -m venv .venv
   .venv\Scripts\activate  # Windows
   # 或 source .venv/bin/activate  # Linux/macOS
   ```
2. 安裝依賴：
   ```bash
   pip install -r requirements.txt
   ```

---

## 啟動方式

### 1. CLI (stdio) 模式

```bash
python src/mcp_server.py
```

### 2. HTTP 服務模式（建議用於 API 測試或 HTTP 整合）

```bash
uvicorn src.mcp_server:app --port 8000
```

---

## mcp.json 設定範例

若要讓 Claude Desktop 或其他 MCP 客戶端自動偵測並啟動本伺服器，MCP setting JSON 內容如下（請依實際路徑調整）：

```json
{
  "mcpServers": {
    "MCP-civil-tools": {
      "command": "path/.venv/Scripts/python.exe",
      "args": [
        "path/src/mcp_server.py"
      ]
    }
  }
}
或是
    "MCP-civil-tools": {
      "command": "C:/TOMO/MCP-civil-tools/.venv/Scripts/python.exe",
      "args": [
        "C:/TOMO/MCP-civil-tools/src/mcp_server.py"
      ],
      "disabled": true,
      "autoApprove": []
    },

```
- `command` 請填入你虛擬環境的 python.exe 絕對路徑。
- `args` 請填入 mcp_server.py 的絕對路徑。

以 Windowss 系統，Command args，輸入格式如下：
C:\TOMO\MCP-civil-tools\.venv\Scripts\python.exe C:\TOMO\MCP-civil-tools\src\mcp_server.py

以 MAC 系統，Command args，輸入格式如下：
/Users/tuchengshin/Documents/MCP/MCP-civil-tools/.venv/bin/python3 
/Users/tuchengshin/Documents/MCP/MCP-civil-tools/src/mcp_server.py

---

## 參數預設行為說明

- 使用者只輸入緯度、經度時，沒有輸入其他資訊時，預設的 UTM/TM2 Zone 就是 TM2-121，預設的半球就是北半球。
- 當使用者只輸入平面座標 X, Y 時，預設的 UTM Zone 是 TWD97，半球是北半球。

---

## 主要功能更新

- **所有查表型工具皆支援「支援清單查詢」API**：
  - 例如：`list_supported_materials`（常用材料）、`list_supported_manning_materials`（曼寧係數材料）、`list_supported_max_velocity_materials`（最大流速材料）、`list_supported_regions`（地區/IDF/年雨量）、`list_supported_soil_types`（土壤類型）、`list_supported_land_uses`（土地利用）、`list_supported_practices`（水保措施）、`list_supported_runoff_land_uses`（逕流係數土地利用）、`list_supported_slope_protection_methods`（坡面保護工法）、`list_supported_soil_k_types`（滲透係數土壤）、`list_supported_idf_locations`（IDF曲線地點）等。
  - 查詢時若輸入錯誤或查無資料，會自動提示所有可查詢的支援項目，提升使用體驗。

---

## 個別工具使用範例

### 經緯度轉 UTM
- 輸入：緯度、經度（可選 datum，預設 TWD97）
- 回傳：`"X,Y"` 字串，數值四捨五入到小數點下4位

#### 範例
```
輸入：24.125193616011536, 120.64098341751337
回傳：203650.6040,2670482.4250
```

### UTM 轉經緯度
- 輸入：X, Y（可選 zone, datum, south，預設 TWD97 北半球）
- 回傳：`"緯度,經度"` 字串，數值四捨五入到小數點下15位

#### 範例
```
輸入：203650.604, 2670482.425
回傳：24.125193616011536,120.64098341751337
```

### 曼寧係數查詢
- 輸入：材料名稱（如「混凝土」、「純細砂」、「全面密草生」等）
- 回傳：該材料的曼寧係數 n 及最大容許流速範圍

#### 範例
```
輸入：全面密草生的曼寧係數?
回傳：全面密草生 的曼寧係數 n = 0.040，最大容許流速範圍：1.5~2.5 m/s
混凝土的曼寧係數?
```

### 主動土壓力係數計算
- 輸入：內摩擦角 phi（度）時，主動土壓力係數 Ka 是多少?
- 回傳：主動土壓力係數 Ka

#### 範例
```
輸入：內摩擦角 phi（度）=30 時，主動土壓力係數 Ka 是多少?
回傳：主動土壓力係數 Ka = 0.3333
```

### 被動土壓力係數計算
- 輸入：內摩擦角 phi（度）=30 時，被動土壓力係數 Kp 是多少?
- 回傳：被動土壓力係數 Kp

#### 範例
```
輸入：內摩擦角 phi（度）=30 時，被動土壓力係數 Kp 是多少?
回傳：被動土壓力係數 Kp = 3.0
```

### 排水斷面流速/流深/流量計算（多斷面支援）
- 輸入：
  - cross_section_type（斷面型式：圓形、矩形、梯形）
  - flow（流量 Q, cms）
  - slope（坡度 %）
  - manning_n（曼寧係數）
  - material（渠道材質，可選，檢核流速）
  - diameter（管徑 cm，圓形用）
  - rect_width（底寬 cm，矩形/梯形用）
  - rect_height（高度 cm，矩形用）
  - trap_bottom（底寬 cm，梯形用）
  - trap_top（頂寬 cm，梯形用）
  - trap_height（高度 cm，梯形用）
- 回傳：dict，含 success, data, message, report（繁體中文完整報告書，含輸入參數、公式、步驟、檢核、結果）

#### 範例
```
輸入：cross_section_type="圓形", flow=0.5, slope=1.2, manning_n=0.013, material="混凝土", diameter=60
回傳：
{
  "success": true,
  "data": {
    "velocity": 2.12,
    "flow_depth": 0.28,
    ...
  },
  "message": "流速: 2.12 m/s，流深: 0.28 m。計算結果符合安全流速規範。",
  "report": "【排水斷面流速/流深/流量計算報告】\n斷面型式：圓形\n圓形管徑 D=60.0cm\n流量 Q = 0.500 cms...（略）"
}
```
- 報告內容包含：輸入參數、計算公式、步驟、流速檢核、滿流警告、詳細結果。

### 邊坡穩定安全係數計算
- 輸入：坡度、單位重、摩擦角、凝聚力、地下水位、方法
- 回傳：安全係數、方法、是否合格、說明

#### 範例
```
輸入：坡度=30, 單位重=18, 摩擦角=30, 凝聚力=10
回傳：安全係數 = 1.50，方法：簡化法，合格：True
預設安全係數1.5，僅為範例。
```

### 土壤侵蝕模數/流失量計算
- 輸入：坡長、坡度、降雨、土壤因子、植生因子、保護措施因子、方法
- 回傳：侵蝕模數、流失量、方法、說明

#### 範例
```
輸入：坡長=100, 坡度=10, 降雨=1200, 土壤因子=0.3, 植生因子=0.5, 保護因子=0.8
回傳：侵蝕模數 = 100.00，流失量 = 10.00，方法：USLE
預設侵蝕模數100，流失量10，僅為範例。
```

### 集水區最大逕流量計算
- 輸入：面積、降雨強度、逕流係數、方法
- 回傳：最大逕流量、方法、說明

#### 範例
```
輸入：面積=2, 強度=100, 逕流係數=0.6
回傳：最大逕流量 Q = 5.00 cms，方法：Rational
預設逕流量5cms，僅為範例。
```

### 護岸/擋土牆穩定檢核
- 輸入：牆高、厚度、單位重、摩擦角、凝聚力、背填坡度、地下水位
- 回傳：滑動、傾倒、承載安全係數、是否合格、說明

#### 範例
```
輸入：牆高=3, 厚=1, 單位重=18, 摩擦角=30, 凝聚力=10, 背填坡度=10
回傳：滑動SF=2.00，傾倒SF=2.50，承載SF=3.00，合格：True
預設安全係數均合格，僅為範例。
```

### 植生護坡設計建議
- 輸入：坡度、土壤類型、氣候
- 回傳：建議工法、草種、覆蓋率、說明

#### 範例
```
輸入：坡度=30, 土壤=壤土, 氣候=亞熱帶
回傳：建議工法：噴播草皮，草種：百慕達草，覆蓋率：90.0%
預設建議，僅為範例。
```

### 常用材料設計參數查詢
- 輸入：材料名稱
- 回傳：單位重、凝聚力、摩擦角、強度、說明

#### 範例
```
輸入：材料=壤土
回傳：材料：壤土，單位重：18.0kN/m3，凝聚力：10.0kPa，摩擦角：30.0°，強度：200.0kPa
預設參數，僅為範例。
```

### 坡面保護工法建議
- 輸入：坡度、土壤、降雨
- 回傳：建議工法、說明

#### 範例
```
輸入：坡度=30, 土壤=砂土, 降雨=1200
回傳：建議工法：格框+草皮
預設建議，僅為範例。
```

### 滲水設施設計
- 輸入：設施型式、土壤滲透係數、集水面積、降雨
- 回傳：設計流量、建議尺寸、說明

#### 範例
```
輸入：型式=滲水井, k=0.001, 面積=100, 降雨=1200
回傳：設施型式：滲水井，設計流量：1.0cms，建議尺寸：直徑1m,深1.5m
預設建議，僅為範例。
```

### IDF曲線查詢
- 輸入：地點、重現期、歷時
- 回傳：降雨強度、說明

#### 範例
```
輸入：地點=台中, 重現期=10, 歷時=60
回傳：地點：台中，重現期：10年，歷時：60分鐘，強度：100.0 mm/hr
預設強度100mm/hr，僅為範例。
```

### 查詢支援清單
- 你可以直接查詢有哪些可用的材料、地區、工法等：

```
查詢：有哪些常用材料？
回傳：['一般黏土', '砂土', '礫石', '混凝土', ...]

查詢：可以查詢哪些水溝鋪面的曼寧係數？
回傳：['純細砂', '混凝土', '全面密草生', ...]

查詢：有哪些坡面保護工法？
回傳：['草皮或直接播種', '噴播草皮+格框/土工網', ...]

查詢：有哪些IDF地點？
回傳：['台北市', '新北市', '台中市']
```

### 查詢失敗時自動提示
- 若查詢時輸入錯誤，會自動回傳所有支援查詢的項目：

```
輸入：查詢不存在的材料名稱
回傳：查無此材料，支援查詢的材料有：一般黏土, 砂土, 礫石, 混凝土, ...
```
## 土石籠擋土牆穩定分析

此功能用於分析土石籠擋土牆的穩定性，包括主動和被動土壓力計算。

### 使用方法

```python
result = check_gabion_stability(
    height=3.0,          # 土石籠高度 (m)
    width=2.0,           # 土石籠寬度 (m)
    wall_weight=100,     # 擋土牆總重 (kN/m)
    phi=30,              # 土壤內摩擦角 (°)
    delta=20,            # 牆背摩擦角 (°)，預設 0
    theta=0,             # 牆背傾斜角 (°)，預設 0
    i=0,                 # 地表傾斜角 (°)，預設 0
    gamma=18,            # 土壤飽和單位重 (kN/m³)，預設 18
    friction_coef=0.5,   # 摩擦係數，預設 0.5
    pressure_mode="active"  # 土壓力模式 ("active" 或 "passive")，預設 "active"
)
```

### 回傳結果

函數回傳一個字典，包含以下內容：

- `success`: 布林值，表示計算是否成功
- `data`: 計算結果數據，包含：
  - `earth_pressure_coef`: 土壓力係數
  - `total_pressure`: 總土壓力 (kN/m)
  - `vertical_force`: 垂直力分量 (kN/m)
  - `horizontal_force`: 水平力分量 (kN/m)
  - `restoring_moment`: 抗傾覆力矩 (kN·m/m)
  - `overturning_moment`: 傾覆力矩 (kN·m/m)
  - `overturning_safety_factor`: 抗傾覆安全係數
  - `sliding_safety_factor`: 抗滑動安全係數
- `message`: 計算結果摘要
- `report`: 完整的計算報告書（Markdown 格式）

### 計算報告書內容

報告書包含以下章節：

1. 輸入參數
2. 計算公式
3. 計算結果
4. 穩定評估

### U型溝鋼筋量計算
- 輸入：
  - height: 溝高 (m)
  - wall_slope: 溝壁傾角 (m)
  - soil_slope: 土方傾角 (°)
  - soil_angle: 安息角 (°)
  - effective_depth: 有效厚度 (m)
  - soil_weight: 土重 (kN/m³)，預設 18.0
- 回傳：dict，含計算結果與報告書

#### 範例
```
輸入：混凝土溝，height=1.5, wall_slope=0.5, soil_slope=15, soil_angle=30, effective_depth=0.2，計算水溝的配筋?
回傳：
{
  "success": true,
  "data": {
    "earth_pressure_coef": 0.3333,
    "earth_pressure": 6.750,
    "moment": 3.375,
    "rebar_area": 13.780
  },
  "message": "土壓力係數 Ka = 0.3333, 土壓力 P = 6.750 kN/m, 彎矩 M = 3.375 kN·m/m, 鋼筋量 As = 13.780 cm²/m",
  "report": "【U型溝鋼筋量計算報告】\n\n輸入參數：\n- 溝高 H = 1.500 m\n- 溝壁傾角 m = 0.500\n- 土方傾角 i = 15.00°\n- 安息角 ψ = 30.00°\n- 有效厚度 d = 0.200 m\n- 土重 γ = 18.0 kN/m³\n\n計算公式：\n1. 土壓力係數 Ka = cos²(ψ+m) / [cos²m·(1+√Q)²]\n   其中 Q = [sinψ·sin(ψ-i)] / [cos(m+i)·cosm]\n2. 土壓力 P = γ·H²·Ka / (2·cosm)\n3. 彎矩 M = γ·H³·Ka / (6·cosm)\n4. 鋼筋量 As = M / (fs·d) × 10⁶ / 1000\n\n計算結果：\n- 土壓力係數 Ka = 0.3333\n- 土壓力 P = 6.750 kN/m\n- 彎矩 M = 3.375 kN·m/m\n- 鋼筋量 As = 13.780 cm²/m"
}
```
## 鋼筋查詢功能
本工具提供以下鋼筋資料查詢功能：

1. 列出所有可用的鋼筋編號
   ```python
   list_rebar_numbers()
   ```

2. 查詢特定鋼筋編號的規格資料
   ```python
   get_rebar_specs(rebar_number="#3")
   ```

3. 計算鋼筋重量
   ```python
   calculate_rebar_weight(rebar_number="#3", length=10.0)
   ```

可用的鋼筋編號包括：#3、#4、#5、#6、#7、#8、#9、#10、#11，每個編號對應的規格資料包括：
- 直徑（mm）
- 截面積（cm²）
- 單位重量（kg/m）
- 周長（mm）

### 鋼筋規格查詢
輸入「鋼筋規格 #3」或「鋼筋資料 #3」可查詢特定鋼筋的詳細資料，包括：
- 鋼筋編號
- 直徑 (mm)
- 截面積 (cm²)
- 單位重量 (kg/m)
- 周長 (mm)

### 鋼筋重量計算
輸入「鋼筋重量 長度 6 #3」可計算指定長度的鋼筋重量，例如：
- 長度：6m
- 鋼筋：#3
- 重量：2.04 kg

### 鋼筋截面積查詢
輸入「鋼筋截面積 #3」可查詢特定鋼筋的截面積，例如：
- #3 鋼筋截面積：0.71 cm²

### 所有鋼筋列表
輸入「所有鋼筋」可列出所有可用的鋼筋編號，包括：
- #3 至 #11 鋼筋
- 各鋼筋的基本規格

### U 型溝配筋建議
輸入「U型溝配筋 面積 10」可查詢建議的配筋方式，系統會：
1. 根據輸入的鋼筋斷面積 (cm²/m)
2. 自動計算並建議主筋與副筋的配筋方式
3. 間距會取整到最接近的 5cm
4. 提供完整的配筋建議報告

### 逕流係數查詢
逕流係數查詢工具可依據水土保持技術規範第18條，查詢不同土地利用類型或集水區狀況的逕流係數C值。

使用方式：
```python
# 查詢特定土地利用類型的逕流係數
result = query_runoff_coeff(land_use="農業區")

# 查詢開發中狀態的逕流係數
result = query_runoff_coeff(land_use="農業區", is_developing=True)

# 直接指定逕流係數值
result = query_runoff_coeff(runoff_coeff=0.5)
```

查詢結果包含：
- 逕流係數值
- 來源說明
- 規範依據
- 土地利用類型
- 開發狀態
- 數值範圍（如適用）


## 使用範例(功能導向，讓大模型自動調用相關工具來求解)

- 有一集水區面積約5ha，農業區，新北市，重現期50年，降雨延時60min，調用工具取得該地區的降雨強度值，計算該集水區最大逕流量?
- 有一集水區面積約5ha，農業區，新北市，降雨延時60min，調用工具取得該地區的降雨強度值，請分別計算重現期25年、50年該集水區最大逕流量? 並輸出完整計算式的報告書。
- 依據這個逕流量，設計一條寬50cm，深70cm的混凝土溝，設計坡度容許範圍? 
- 依據這個逕流量，設計一條寬50cm，溝深60cm的混凝土溝，坡度2.5%，檢核該設計是否OK，並出具一份完整的檢核報告書。 
- 有一個破碎岩盤的坡面，角度約60度，位於熱帶多雨地區，請提出坡面保護建議。
- 無基礎的混凝土重力式擋土牆，全高2.3M，牆頂寬50cm，牆底寬100cm，土壤單位重18kN/m3，摩擦角30度，凝聚力10kPa，地下水位在牆頂以下2M，請調用MCP工具，檢核其穩定性，並輸出完整計算報告書。
- 混凝土溝，height=1.5, wall_slope=0.5, soil_slope=15, soil_angle=30, effective_depth=0.2，調用工具，計算水溝所需的鋼筋量及建議配筋?

---

## 授權

本專案採用 MIT License 授權，歡迎自由使用與貢獻。
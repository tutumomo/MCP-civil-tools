# MCP Civil Tools 伺服器

本專案是一個基於 MCP 協議的 Python 伺服器，提供經緯度與 UTM/TWD97 座標互轉，以及多種常用土木工程計算工具（如曼寧係數查詢、土壓力係數、排水溝流速等），適用於 LLM 工具、Claude Desktop 等 AI 應用整合。

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
```
- `command` 請填入你虛擬環境的 python.exe 絕對路徑。
- `args` 請填入 mcp_server.py 的絕對路徑。

以 Windowss 系統，Command args，輸入格式如下：
c:\Users\tutumomo\Documents\MCP\MCP-civil-tools\.venv\Scripts\python.exe c:\Users\tutumomo\Documents\MCP\MCP-civil-tools\src\mcp_server.py

以 MAC 系統，Command args，輸入格式如下：
/Users/tuchengshin/Documents/MCP/MCP-civil-tools/.venv/bin/python3 
/Users/tuchengshin/Documents/MCP/MCP-civil-tools/src/mcp_server.py

---

## 參數預設行為說明

- 使用者只輸入緯度、經度時，沒有輸入其他資訊時，預設的 UTM/TM2 Zone 就是 TM2-121，預設的半球就是北半球。
- 當使用者只輸入平面座標 X, Y 時，預設的 UTM Zone 是 TWD97，半球是北半球。

---

## API 使用範例

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
- 輸入：材料名稱（如「混凝土光滑」、「天然土渠」等）
- 回傳：該材料的曼寧係數 n

#### 範例
```
輸入：混凝土光滑
回傳：混凝土光滑 的曼寧係數 n = 0.012
```

### 主動土壓力係數計算
- 輸入：內摩擦角 phi（度）
- 回傳：主動土壓力係數 Ka

#### 範例
```
輸入：30
回傳：主動土壓力係數 Ka = 0.3333
```

### 被動土壓力係數計算
- 輸入：內摩擦角 phi（度）
- 回傳：被動土壓力係數 Kp

#### 範例
```
輸入：30
回傳：被動土壓力係數 Kp = 3.0
```

### 排水溝流速計算（曼寧公式）
- 輸入：n（曼寧係數）、r（水力半徑 m）、s（坡度）
- 回傳：流速 v (m/s)

#### 範例
```
輸入：n=0.012, r=0.5, s=0.01
回傳：流速 v = 3.4925 m/s
```

### 排水溝流量計算
- 輸入：v（流速 m/s）、a（斷面積 m2）
- 回傳：流量 Q (cms)

#### 範例
```
輸入：v=2.0, a=0.8
回傳：流量 Q = 1.6000 cms
```

---

## 授權

本專案採用 MIT License 授權，歡迎自由使用與貢獻。
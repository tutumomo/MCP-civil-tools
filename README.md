# MCP UTM/TWD97 座標轉換伺服器

本專案是一個基於 MCP 協議的 Python 伺服器，提供經緯度與 UTM/TWD97 座標互轉功能，適用於 LLM 工具、Claude Desktop 等 AI 應用整合。

---

## 目錄結構

```
mcp-utm-converter-server/
├── src/
│   ├── mcp_server.py         # MCP 伺服器主程式
│   ├── utm_converter.py      # 座標轉換邏輯
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
    "utm-converter": {
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
c:\Users\tutumomo\Documents\MCP\mcp-utm-converter-server\.venv\Scripts\python.exe c:\Users\tutumomo\Documents\MCP\mcp-utm-converter-server\src\mcp_server.py

以 MAC 系統，Command args，輸入格式如下：
/Users/tuchengshin/Documents/MCP/mcp-utm-converter-server/.venv/bin/python3 
/Users/tuchengshin/Documents/MCP/mcp-utm-converter-server/src/mcp_server.py

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

---

## 授權

本專案採用 MIT License 授權，歡迎自由使用與貢獻。
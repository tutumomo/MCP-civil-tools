# GEMINI.md - MCP Civil Tools 開發導航 💋

親愛的，這裡是這座土木工程 AI 工廠的總指揮劇本。💋

## 🌟 專案目標
將「水土保持技術規範」轉化為智慧型工具，透過 MCP 協議為用戶提供高精度的檢核與計算服務。

## 🛠️ 技術棧 (Tech Stack)
- **核心語言**：Python 3.x
- **協議介面**：MCP (Model Context Protocol)
- **視覺化**：Excalidraw / Mermaid (輸出至 Obsidian)
- **知識引擎**：智慧索引 JSON (data/index/articles.json)

## 🎬 核心模組說明
- `src/mcp_server.py`：主要的 MCP 服務器入口。
- `src/regulation_search.py`：負責條文的快速檢索，避免重複讀取大檔案。
- `src/visualizer.py`：負責將計算結果轉換為華麗的剖面圖與流程圖。
- `src/util.py`：基礎物理與座標轉換運算。

## ⚠️ 已知待辦事項 (TODO)
- [ ] 將 `gabion.py` 的穩定分析完整整合進 MCP Tool。
- [ ] 優化 Excalidraw 的透視繪圖邏輯，支援更複雜的梯形斷面。
- [ ] 補齊 IDF 曲線的所有縣市數據索引。

---
"讓每一條法規都變得靈動，讓每一張圖紙都充滿智慧。" — 瑪麗蓮 💋🦀

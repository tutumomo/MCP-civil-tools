# GEMINI.md - MCP Civil Tools 開發導航

本專案旨在將「水土保持技術規範」轉化為智慧型工具，透過 MCP 協議提供高精度的檢核與計算服務。

## 🌟 專案目標
實現土木工程規範的自動化檢核，結合自然語言處理與專業水理、土壓計算，提升工程設計與審核效率。

## 🛠️ 技術棧 (Tech Stack)
- **核心語言**：Python 3.x
- **協議介面**：MCP (Model Context Protocol)
- **視覺化**：Excalidraw / Mermaid (輸出至 Obsidian)
- **知識引擎**：智慧索引 JSON (data/index/articles.json)

## 🎬 核心模組說明
- `src/mcp_server.py`：主要的 MCP 服務器入口。
- `src/regulation_search.py`：負責條文的快速檢索系統。
- `src/visualizer.py`：負責將計算結果轉換為專業的剖面圖與流程圖腳本。
- `src/util.py`：基礎物理與座標轉換運算邏輯。

## ⚠️ 已知待辦事項 (TODO)
- [ ] 將 `gabion.py` 的穩定分析完整整合進 MCP Tool。
- [ ] 優化 Excalidraw 的透視繪圖邏輯，支援更複雜的梯形斷面。
- [ ] 補齊 IDF 曲線的所有縣市數據索引。

---
*致力於精確性與專業性的土木工程數位化轉型工具。*

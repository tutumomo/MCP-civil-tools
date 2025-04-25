# 代碼風格和模式

請遵循以下代碼風格和常用模式：

- **組件命名**: React 組件應使用 **PascalCase** 命名 (例如，`VocabularyCard`, `WordInput`)。
- **函數命名**: 輔助函數應使用 **camelCase** 命名 (例如，`formatWord`, `generateStory`)。
- **變量命名**: 變量應使用 **camelCase** 命名 (例如，`userWords`, `generatedText`)。常量可以使用 **UPPER_SNAKE_CASE** (例如，`MAX_WORDS`)。
- **組件類型**: 除非有特殊理由，否則優先使用**函數式組件**和**箭頭函數**語法。
- **狀態管理**: 優先考慮使用 React 的 `useState` 和 `useEffect` 進行組件級別的狀態管理。如果需要更複雜的狀態管理，請根據具體情況考慮 `useContext` 或其他的輕量級狀態管理庫。
- **錯誤處理**: 在進行 API 調用和處理用戶輸入時，務必進行適當的錯誤處理。
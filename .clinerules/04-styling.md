# 樣式指南

- **Tailwind CSS 優先**: 在編寫組件樣式時，**始終優先使用 Tailwind CSS 的工具類**直接在 JSX 中進行樣式化。
- **自定義 CSS**: 如果需要自定義 CSS，請盡量通過在 `tailwind.config.js` 中擴展 Tailwind 的配置，或者使用 Tailwind 的 **`@apply` 指令**來組合 Tailwind 的工具類。避免編寫大量的完全自定義的 CSS。
- **Tailwind 配置**: 任何對 Tailwind CSS 的自定義（例如，添加新的顏色、字體、斷點）都應該在 `tailwind.config.js` 文件中進行。
- **響應式設計**: 在設計頁面和組件時，請始終考慮**響應式設計**，並利用 Tailwind CSS 的響應式前綴 (例如，`sm:`, `md:`, `lg:`) 來適配不同的屏幕尺寸。
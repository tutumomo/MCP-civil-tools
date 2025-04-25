# 重要文件位置（請勿修改）

為了項目的安全性和穩定性，請指示 Cline **不要讀取或修改** 以下文件或目錄：

- **環境變量文件**: `.env`，`.env.local`，`.env.development.local`，`.env.production.local` 等包含 API 密鑰或其他敏感環境變量的文件。
- **Ollama 相關配置**: 任何 Ollama 的配置文件或模型存儲目錄（如果 Cline 需要直接操作這些文件）。通常情況下，Cline 只需要通過 API 與 Ollama 交互，不需要直接修改其配置。
- **構建輸出目錄**: 通常是 `/.next/`，`/build/`，`/dist/` 等。
- **包管理器的鎖文件**: `package-lock.json` (npm) 或 `yarn.lock` (yarn)。

請確保 Cline 在執行任何操作時都遵守這些限制，以防止意外修改或暴露敏感信息。
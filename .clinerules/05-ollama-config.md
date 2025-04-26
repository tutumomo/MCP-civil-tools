# Ollama 配置

本項目使用本地運行的 Ollama 服務來生成短文。請注意以下配置和步驟：

- **Ollama 服務**: 確保 Ollama 服務已經在本地運行，並且可以通過默認端口 `http://localhost:11434/` 訪問.
- **模型**: 使用的模型是 **`gemma3:12b-it-qat`**。在與 Ollama 交互之前，請確保該模型已經通過 Ollama 下載 . 可以使用命令 `ollama run gemma3:3b-it-qat` （或者 `gemma3:7b-it-qat` 或 `gemma3:12b-it-qat`，根據實際下載的版本）來下載和運行模型。
- **API 交互**: 當需要生成包含用戶輸入單詞的短文時，請在 Next.js 的 API 路由中實現與 Ollama API 的交互。可以使用 `fetch` 或其他 HTTP 庫向 Ollama API 發送 POST 請求，請求生成文本。
- **API 端點**: Ollama 的文本生成 API 端點通常是 `/api/generate`。
- **請求體**: 發送給 Ollama API 的請求體應包含 `model` 參數（設置為 `gemma 3:12b-it-qat`）以及 `prompt` 參數，該 prompt 應該包含用戶輸入的單詞，並清晰地指示模型生成一個包含這些單詞的短文以幫助記憶。
- **示例 Prompt**: 一個示例 prompt 可能是：「請使用以下單詞創作一個簡短的故事：[用戶輸入的單詞列表，用逗號分隔]。」
- **錯誤處理**: 在與 Ollama API 交互時，請務必處理可能出現的錯誤，例如網絡連接錯誤或模型生成錯誤。
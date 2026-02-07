import json
import os

# Path to the indexed articles
INDEX_FILE = os.path.join(os.path.dirname(__file__), "../data/index/articles.json")

def query_regulation(query_text: str):
    """
    智慧檢索水土保持技術規範。
    輸入可以是「第18條」或關鍵字。
    """
    if not os.path.exists(INDEX_FILE):
        return "錯誤：尚未建立規範索引。請先執行 build_index.py。"

    with open(INDEX_FILE, 'r', encoding='utf-8') as f:
        articles = json.load(f)

    # 1. Exact match for "第 XX 條"
    import re
    article_match = re.search(r'第\s*(\d+)\s*條', query_text)
    if article_match:
        target = f"第 {article_match.group(1)} 條"
        if target in articles:
            return f"【{target}】\n\n{articles[target]}"
        
        # Try without space just in case
        target_no_space = f"第{article_match.group(1)}條"
        for key in articles:
            if target_no_space in key.replace(" ", ""):
                return f"【{key}】\n\n{articles[key]}"

    # 2. Keyword search
    results = []
    for title, content in articles.items():
        if query_text.lower() in content.lower() or query_text.lower() in title.lower():
            results.append(f"【{title}】\n{content[:200]}...")
            if len(results) >= 3: # Limit results
                break
    
    if results:
        return "\n\n---\n\n".join(results)
    
    return "查無相關條文內容。"

if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        print(query_regulation(sys.argv[1]))

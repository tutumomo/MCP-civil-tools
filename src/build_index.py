import os
import json
import re

DATA_PATH = "/Users/tuchengshin/MCP-civil-tools/水土保持技術規範(含附件).txt"
INDEX_DIR = "/Users/tuchengshin/MCP-civil-tools/data/index"

def build_index():
    if not os.path.exists(DATA_PATH):
        print("Data file not found.")
        return

    with open(DATA_PATH, 'r', encoding='utf-8') as f:
        content = f.read()

    # Simple logic: split by "第 XX 條"
    articles = re.split(r'(第\s*\d+\s*條)', content)
    
    index = {}
    if len(articles) > 1:
        for i in range(1, len(articles), 2):
            title = articles[i].strip()
            # The content is the next part
            text = articles[i+1].strip() if i+1 < len(articles) else ""
            index[title] = text
            
    # Save index to a JSON file
    index_file = os.path.join(INDEX_DIR, "articles.json")
    with open(index_file, 'w', encoding='utf-8') as f:
        json.dump(index, f, ensure_ascii=False, indent=2)
    
    print(f"Index built with {len(index)} articles.")

if __name__ == "__main__":
    os.makedirs(INDEX_DIR, exist_ok=True)
    build_index()

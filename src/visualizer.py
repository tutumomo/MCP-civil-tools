import os
import json
import uuid

def generate_retaining_wall_excalidraw(height, top_width, bottom_width, sf_sliding, sf_overturning):
    """
    生成擋土牆剖面圖的 Excalidraw JSON。
    """
    # Base coordinates
    base_x = 500
    base_y = 500
    scale = 100 # 1m = 100px
    
    h_px = height * scale
    tw_px = top_width * scale
    bw_px = bottom_width * scale
    
    # Points for the wall (polygon)
    # 0: top-left, 1: top-right, 2: bottom-right, 3: bottom-left
    # Assume rectangular/trapezoidal wall facing right
    p0 = {"x": base_x, "y": base_y - h_px}
    p1 = {"x": base_x + tw_px, "y": base_y - h_px}
    p2 = {"x": base_x + bw_px, "y": base_y}
    p3 = {"x": base_x, "y": base_y}
    
    elements = [
        # The Wall Body
        {
            "type": "rectangle",
            "x": p3["x"], "y": p0["y"],
            "width": bw_px, "height": h_px,
            "backgroundColor": "#e5e7eb",
            "strokeColor": "#374151",
            "id": "wall_body"
        },
        # Title
        {
            "type": "text",
            "x": base_x, "y": base_y - h_px - 40,
            "text": f"擋土牆剖面圖 (H={height}m)",
            "fontSize": 20,
            "fontFamily": 5
        },
        # SF Info
        {
            "type": "text",
            "x": base_x + bw_px + 20, "y": base_y - h_px,
            "text": f"安全係數:\n滑動: {sf_sliding}\n傾倒: {sf_overturning}",
            "fontSize": 14,
            "strokeColor": "#1e40af" if sf_sliding >= 1.5 else "#ef4444"
        }
    ]
    
    # Simple Excalidraw wrapper
    data = {
        "type": "excalidraw",
        "version": 2,
        "source": "https://openclaw.ai",
        "elements": elements,
        "appState": {"viewBackgroundColor": "#ffffff"},
        "files": {}
    }
    
    return json.dumps(data, ensure_ascii=False)

def wrap_in_obsidian_excalidraw(json_data):
    return f"""---
excalidraw-plugin: parsed
tags: [excalidraw, civil-tools]
---
# Excalidraw Data
## Text Elements
%%
## Drawing
```json
{json_data}
```
%%"""

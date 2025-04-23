import tkinter as tk
from tkinter import ttk
from math import sin, acos, sqrt

# 曼寧係數範圍 (根據水土保持技術規範)
MANNING_RANGES = {
    "平整混凝土面": (0.014, 0.016),
    "未粉飾混凝土面": (0.017, 0.020),
    "噴漿面": (0.020, 0.025),
    "混凝土砌卵石側牆 (混凝土打底)": (0.025, 0.030),
    "混凝土砌卵石側牆 (天然渠底)": (0.027, 0.030),
    "乾砌卵石": (0.030, 0.035),
    "乾砌塊石": (0.032, 0.035),
    "水泥砂漿砌磚": (0.015, 0.018),
    "平滑瀝青面": (0.014, 0.015),
    "粗糙瀝青面": (0.017, 0.018),
    "清潔土渠": (0.020, 0.025),
    "短草土渠": (0.027, 0.033),
    "密生水草高及水面土渠": (0.080, 0.120),
    "密生灌木高及水面土渠": (0.100, 0.140),
    "均勻岩面": (0.030, 0.040),
    "不規則岩石面": (0.040, 0.050),
    "內面光滑混凝土管": (0.013, 0.014),
    "排水瓦管": (0.017, 0.018)
}

# 渠道材質流速限制 (單位：m/s)
CHANNEL_VELOCITY_LIMITS = {
    "純細砂": (0.23, 0.30),
    "不緻密之細砂": (0.30, 0.46),
    "粗石及細砂土": (0.46, 0.61),
    "平常砂土": (0.61, 0.76),
    "砂質壤土": (0.76, 0.84),
    "堅壤土及粘質壤土": (0.91, 1.14),
    "平常礫土": (1.23, 1.52),
    "全面密草生": (1.50, 2.50),
    "粗礫、石礫及砂礫": (1.52, 1.83),
    "礫岩、硬土層、軟質水成岩": (1.83, 2.44),
    "硬岩": (3.05, 4.57),
    "混凝土": (4.57, 6.10)
}

# 產生曼寧選項顯示字串與對應範圍的對照字典
manning_options = []
manning_display_to_range = {}
for name, (low, high) in MANNING_RANGES.items():
    avg = (low + high) / 2
    display = f"{name} ({low:.3f}~{high:.3f}, 平均{avg:.3f})"
    manning_options.append(display)
    manning_display_to_range[display] = (low, high)

# 產生渠道材質選項顯示字串與流速限制的對照字典
channel_options = []
channel_display_to_limits = {}
for name, (v_low, v_high) in CHANNEL_VELOCITY_LIMITS.items():
    display = f"{name} ({v_low:.2f}~{v_high:.2f} m/s)"
    channel_options.append(display)
    channel_display_to_limits[display] = (v_low, v_high)

def update_inputs(event=None):
    # 隱藏所有斷面相關輸入欄位
    diameter_label.grid_remove()
    diameter_entry.grid_remove()
    rect_width_label.grid_remove()
    rect_width_entry.grid_remove()
    rect_height_label.grid_remove()
    rect_height_entry.grid_remove()
    trap_bottom_label.grid_remove()
    trap_bottom_entry.grid_remove()
    trap_top_label.grid_remove()
    trap_top_entry.grid_remove()
    trap_height_label.grid_remove()
    trap_height_entry.grid_remove()

    cs = cross_section_combobox.get()
    if cs == "圓形 (RCP 管)":
        diameter_label.grid(row=5, column=0, padx=5, pady=5)
        diameter_entry.grid(row=5, column=1, padx=5, pady=5)
    elif cs == "矩形溝":
        rect_width_label.grid(row=5, column=0, padx=5, pady=5)
        rect_width_entry.grid(row=5, column=1, padx=5, pady=5)
        rect_height_label.grid(row=6, column=0, padx=5, pady=5)
        rect_height_entry.grid(row=6, column=1, padx=5, pady=5)
    elif cs == "梯形溝":
        trap_bottom_label.grid(row=5, column=0, padx=5, pady=5)
        trap_bottom_entry.grid(row=5, column=1, padx=5, pady=5)
        trap_top_label.grid(row=6, column=0, padx=5, pady=5)
        trap_top_entry.grid(row=6, column=1, padx=5, pady=5)
        trap_height_label.grid(row=7, column=0, padx=5, pady=5)
        trap_height_entry.grid(row=7, column=1, padx=5, pady=5)

def update_manning_value(event=None):
    """ 當曼寧係數下拉選單改變時，自動更新右側的輸入欄位預設值為該選項的平均值 """
    selection = manning_combobox.get()
    if selection in manning_display_to_range:
        low, high = manning_display_to_range[selection]
        avg = (low + high) / 2
        manning_entry.delete(0, tk.END)
        manning_entry.insert(0, f"{avg:.3f}")

def compute_rect_flow_depth(Q, n, S, b, h):
    epsilon = 1e-6
    y_low, y_high = 0.0001, h - epsilon
    y = (y_low + y_high) / 2
    for _ in range(100):
        area = b * y
        perimeter = b + 2 * y
        R = area / perimeter
        V = (1 / n) * (R ** (2 / 3)) * (S ** 0.5)
        Q_cal = area * V
        if abs(Q_cal - Q) < 1e-7:
            break
        if Q_cal < Q:
            y_low = y
        else:
            y_high = y
        y = (y_low + y_high) / 2
    return V, y

def compute_trap_flow_depth(Q, n, S, b1, b2, h):
    epsilon = 1e-6
    y_low, y_high = 0.0001, h - epsilon
    y = (y_low + y_high) / 2
    for _ in range(100):
        area = (b1 + b2) * y / 2
        side = sqrt(y**2 + ((b2 - b1) / 2) ** 2)
        perimeter = b2 + 2 * side
        R = area / perimeter
        V = (1 / n) * (R ** (2 / 3)) * (S ** 0.5)
        Q_cal = area * V
        if abs(Q_cal - Q) < 1e-7:
            break
        if Q_cal < Q:
            y_low = y
        else:
            y_high = y
        y = (y_low + y_high) / 2
    return V, y

def compute_circular_flow_depth(D, Q, n, S):
    epsilon = 1e-6
    y_low, y_high = 0.0001, D - epsilon
    y = (y_low + y_high) / 2
    for _ in range(100):
        theta = 2 * acos(1 - 2 * y / D)
        A = (D ** 2 / 8) * (theta - sin(theta))
        P = (D / 2) * theta
        R = A / P
        V = (1 / n) * (R ** (2 / 3)) * (S ** 0.5)
        Q_cal = A * V
        if abs(Q_cal - Q) < 1e-7:
            break
        if Q_cal < Q:
            y_low = y
        else:
            y_high = y
        y = (y_low + y_high) / 2
    return V, y

def calculate():
    try:
        Q = float(flow_entry.get())  # 流量 (cms)
        slope_percent = float(slope_entry.get())
        S = slope_percent / 100  # 坡度 S = 百分比/100

        # 取得曼寧係數，由使用者可修改的文字框取得
        n = float(manning_entry.get())
        
        cs = cross_section_combobox.get()
        if cs == "圓形 (RCP 管)":
            D = float(diameter_entry.get()) / 100  # 單位轉換成公尺
            V, y = compute_circular_flow_depth(D, Q, n, S)
        elif cs == "矩形溝":
            b = float(rect_width_entry.get()) / 100
            h = float(rect_height_entry.get()) / 100
            V, y = compute_rect_flow_depth(Q, n, S, b, h)
        elif cs == "梯形溝":
            b1 = float(trap_bottom_entry.get()) / 100
            b2 = float(trap_top_entry.get()) / 100
            h = float(trap_height_entry.get()) / 100
            V, y = compute_trap_flow_depth(Q, n, S, b1, b2, h)
        else:
            result_label.config(text="請選擇有效的排水斷面類型")
            return

        # 根據所選渠道材質取得安全流速範圍
        material_disp = material_combobox.get()
        if material_disp in channel_display_to_limits:
            min_safe, max_safe = channel_display_to_limits[material_disp]
            if V < min_safe:
                check_msg = f"【檢核警告】計算流速 {V:.3f} m/s 低於安全下限 {min_safe:.2f} m/s，可能導致泥砂淤積。"
            elif V > max_safe:
                check_msg = f"【檢核警告】計算流速 {V:.3f} m/s 超過安全上限 {max_safe:.2f} m/s，請考慮設置消能設施。"
            else:
                check_msg = "計算結果符合安全流速規範。"
        else:
            check_msg = "無法取得該材質對應的流速限制。"

        # 改用相對誤差檢查滿流狀態
        full_flow_warning = ""
        rel_tol = 1e-3  # 當差值比例小於 0.1% 時，認定為滿流
        if cs == "圓形 (RCP 管)":
            if (D - y) / D < rel_tol:
                full_flow_warning = "【檢核警告】計算流深已接近管徑，可能表示管道已滿流。"
        elif cs in ["矩形溝", "梯形溝"]:
            if (h - y) / h < rel_tol:
                full_flow_warning = "【檢核警告】計算流深已接近通道設計高度，可能表示通道已滿流。"

        final_message = f"流速: {V:.3f} m/s\n流深: {y:.3f} m\n{check_msg}"
        if full_flow_warning:
            final_message += "\n" + full_flow_warning

        result_label.config(text=final_message)
    except ValueError:
        result_label.config(text="請輸入有效的數字")


# GUI 建立
window = tk.Tk()
window.title("排水斷面計算")

# 排水斷面類型下拉選單
ttk.Label(window, text="排水斷面類型:").grid(row=0, column=0, padx=5, pady=5)
cross_section_combobox = ttk.Combobox(window, values=["圓形 (RCP 管)", "矩形溝", "梯形溝"], state="readonly")
cross_section_combobox.grid(row=0, column=1, padx=5, pady=5)
cross_section_combobox.current(0)
cross_section_combobox.bind("<<ComboboxSelected>>", update_inputs)

# 流量與坡度輸入
ttk.Label(window, text="流量 Q (cms):").grid(row=1, column=0, padx=5, pady=5)
flow_entry = ttk.Entry(window)
flow_entry.grid(row=1, column=1, padx=5, pady=5)

ttk.Label(window, text="坡度 (%):").grid(row=2, column=0, padx=5, pady=5)
slope_entry = ttk.Entry(window)
slope_entry.grid(row=2, column=1, padx=5, pady=5)

# 曼寧係數下拉選單與可修改的數值輸入
ttk.Label(window, text="曼寧係數範圍:").grid(row=3, column=0, padx=5, pady=5)
manning_combobox = ttk.Combobox(window, values=manning_options, state="readonly")
manning_combobox.grid(row=3, column=1, padx=5, pady=5)
manning_combobox.current(0)
manning_combobox.bind("<<ComboboxSelected>>", update_manning_value)

ttk.Label(window, text="使用曼寧係數 (可修改):").grid(row=3, column=2, padx=5, pady=5)
manning_entry = ttk.Entry(window)
# 預設填入第一個選項的平均值
default_low, default_high = manning_display_to_range[manning_options[0]]
default_avg = (default_low + default_high) / 2
manning_entry.insert(0, f"{default_avg:.3f}")
manning_entry.grid(row=3, column=3, padx=5, pady=5)

# 新增：渠道材質 (流速限制) 下拉選單，預設設為混凝土
ttk.Label(window, text="渠道材質 (流速限制):").grid(row=4, column=0, padx=5, pady=5)
material_combobox = ttk.Combobox(window, values=channel_options, state="readonly")
material_combobox.grid(row=4, column=1, padx=5, pady=5)
# 設定預設值為混凝土
for i, option in enumerate(channel_options):
    if option.startswith("混凝土"):
        material_combobox.current(i)
        break

# 斷面輸入欄位
diameter_label = ttk.Label(window, text="管徑 (cm):")
diameter_entry = ttk.Entry(window)

rect_width_label = ttk.Label(window, text="底寬 (cm):")
rect_width_entry = ttk.Entry(window)
rect_height_label = ttk.Label(window, text="高度 (cm):")
rect_height_entry = ttk.Entry(window)

trap_bottom_label = ttk.Label(window, text="底寬 (cm):")
trap_bottom_entry = ttk.Entry(window)
trap_top_label = ttk.Label(window, text="頂寬 (cm):")
trap_top_entry = ttk.Entry(window)
trap_height_label = ttk.Label(window, text="高度 (cm):")
trap_height_entry = ttk.Entry(window)

update_inputs()

ttk.Button(window, text="計算", command=calculate).grid(row=8, column=0, columnspan=4, pady=10)
result_label = ttk.Label(window, text="")
result_label.grid(row=9, column=0, columnspan=4, pady=5)

window.mainloop()

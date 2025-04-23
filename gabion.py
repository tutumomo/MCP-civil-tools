'''
取自 土石籠穩定分析.xlsx，經由 Python 轉換而成。
使用 chatwise + MCP工具(excel、filesystem 讀取 google 雲端硬碟內檔案後，透過 DeepSeek V3 0324 編寫
'''

import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import math
from datetime import datetime

class GabionStabilityApp:
    def __init__(self, root):
        self.root = root
        self.root.title("土石籠擋土牆穩定分析程式 (完整公式版)")
        self.root.geometry("1100x850")
        
        # 主框架容器
        main_frame = ttk.Frame(root)
        main_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # 頂部框架（模式選擇+輸入+結果）
        top_frame = ttk.Frame(main_frame)
        top_frame.pack(fill=tk.X, pady=5)
        
        # 計算模式選擇
        self.pressure_mode = tk.StringVar(value="active")
        self.create_mode_selector(top_frame)
        
        # 中間框架（輸入+結果）
        middle_frame = ttk.Frame(top_frame)
        middle_frame.pack(fill=tk.X, pady=5)
        
        self.create_input_section(middle_frame)
        self.create_result_section(middle_frame)
        
        # 計算過程日誌
        self.create_calculation_log(main_frame)
        
        # 底部按鈕框架（確保在最下方）
        button_frame = ttk.Frame(main_frame)
        button_frame.pack(side=tk.BOTTOM, fill=tk.X, padx=10, pady=10)
        self.create_buttons(button_frame)
        
        self.calculation_steps = []
    
    def create_mode_selector(self, parent):
        """創建計算模式選擇區"""
        mode_frame = ttk.LabelFrame(parent, text="分析模式", padding=10)
        mode_frame.pack(fill=tk.X, padx=5, pady=5)
        
        ttk.Radiobutton(mode_frame, text="主動土壓力 (Ka)", variable=self.pressure_mode, 
                       value="active").pack(side=tk.LEFT, padx=10)
        ttk.Radiobutton(mode_frame, text="被動土壓力 (Kp)", variable=self.pressure_mode, 
                       value="passive").pack(side=tk.LEFT, padx=10)
    
    def create_input_section(self, parent):
        """創建輸入參數區"""
        input_frame = ttk.LabelFrame(parent, text="設計參數輸入", padding=10)
        input_frame.pack(side=tk.LEFT, fill=tk.Y, padx=5, pady=5)
        
        # 土壤參數
        params = [
            ("土壤內摩擦角 φ (°)", "phi", 30),
            ("牆背摩擦角 δ (°)", "delta", 20),
            ("牆背傾斜角 θ (°)", "theta", 10),
            ("地表傾斜角 i (°)", "i", 0),
            ("土壤飽和單位重 γ (kN/m³)", "gamma", 18),
            ("擋土牆總重 W (kN/m)", "wall_weight", 50),
            ("土石籠高度 H (m)", "height", 3),
            ("土石籠寬度 B (m)", "width", 2),
            ("摩擦係數 μ", "friction_coef", 0.5)
        ]
        
        for i, (label, name, default) in enumerate(params):
            ttk.Label(input_frame, text=label).grid(row=i, column=0, sticky=tk.W, pady=2)
            setattr(self, name, tk.DoubleVar(value=default))
            ttk.Entry(input_frame, textvariable=getattr(self, name), width=15).grid(row=i, column=1, pady=2)
    
    def create_result_section(self, parent):
        """創建結果顯示區"""
        result_frame = ttk.LabelFrame(parent, text="分析結果", padding=10)
        result_frame.pack(side=tk.LEFT, fill=tk.Y, padx=5, pady=5, expand=True)
        
        # 土壓力結果
        ttk.Label(result_frame, text="土壓力係數:").grid(row=0, column=0, sticky=tk.W, pady=2)
        self.ka_result = ttk.Label(result_frame, text="", width=20)
        self.ka_result.grid(row=0, column=1, sticky=tk.W, pady=2)
        
        ttk.Label(result_frame, text="土壓力總和 (kN/m):").grid(row=1, column=0, sticky=tk.W, pady=2)
        self.pa_result = ttk.Label(result_frame, text="", width=20)
        self.pa_result.grid(row=1, column=1, sticky=tk.W, pady=2)
        
        ttk.Label(result_frame, text="垂直分力 (kN/m):").grid(row=2, column=0, sticky=tk.W, pady=2)
        self.pv_result = ttk.Label(result_frame, text="", width=20)
        self.pv_result.grid(row=2, column=1, sticky=tk.W, pady=2)
        
        ttk.Label(result_frame, text="水平分力 (kN/m):").grid(row=3, column=0, sticky=tk.W, pady=2)
        self.ph_result = ttk.Label(result_frame, text="", width=20)
        self.ph_result.grid(row=3, column=1, sticky=tk.W, pady=2)
        
        # 穩定性分析
        ttk.Label(result_frame, text="扶正力矩 (kN·m/m):").grid(row=4, column=0, sticky=tk.W, pady=2)
        self.mr_result = ttk.Label(result_frame, text="", width=20)
        self.mr_result.grid(row=4, column=1, sticky=tk.W, pady=2)
        
        ttk.Label(result_frame, text="傾倒力矩 (kN·m/m):").grid(row=5, column=0, sticky=tk.W, pady=2)
        self.mf_result = ttk.Label(result_frame, text="", width=20)
        self.mf_result.grid(row=5, column=1, sticky=tk.W, pady=2)
        
        ttk.Label(result_frame, text="抗傾覆安全係數:").grid(row=6, column=0, sticky=tk.W, pady=2)
        self.fs_overturning = ttk.Label(result_frame, text="", width=20)
        self.fs_overturning.grid(row=6, column=1, sticky=tk.W, pady=2)
        
        ttk.Label(result_frame, text="抗滑動安全係數:").grid(row=7, column=0, sticky=tk.W, pady=2)
        self.fs_sliding = ttk.Label(result_frame, text="", width=20)
        self.fs_sliding.grid(row=7, column=1, sticky=tk.W, pady=2)
        
        # 安全評估
        ttk.Label(result_frame, text="整體穩定性評估:").grid(row=8, column=0, sticky=tk.W, pady=2)
        self.stability_assessment = ttk.Label(result_frame, text="", width=30, font=('Arial', 10, 'bold'))
        self.stability_assessment.grid(row=8, column=1, columnspan=2, sticky=tk.W, pady=2)
    
    def create_calculation_log(self, parent):
        """創建計算過程日誌區"""
        log_frame = ttk.LabelFrame(parent, text="詳細計算過程", padding=10)
        log_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)
        
        self.log_text = tk.Text(log_frame, wrap=tk.WORD, height=20, font=('Courier New', 10))
        scrollbar = ttk.Scrollbar(log_frame, orient="vertical", command=self.log_text.yview)
        self.log_text.configure(yscrollcommand=scrollbar.set)
        
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        self.log_text.pack(fill=tk.BOTH, expand=True)
    
    def create_buttons(self, parent):
        """創建功能按鈕區 - 修正可見性問題"""
        # 按鈕樣式設定
        style = ttk.Style()
        style.configure('TButton', font=('Arial', 10), padding=5)
        
        # 計算按鈕
        calc_btn = ttk.Button(parent, text="計算", command=self.calculate, style='TButton')
        calc_btn.pack(side=tk.LEFT, padx=10, ipadx=10)
        
        # 清除按鈕
        clear_btn = ttk.Button(parent, text="清除", command=self.clear, style='TButton')
        clear_btn.pack(side=tk.LEFT, padx=10, ipadx=10)
        
        # 輸出按鈕
        export_btn = ttk.Button(parent, text="輸出計算書", command=self.export_report, style='TButton')
        export_btn.pack(side=tk.LEFT, padx=10, ipadx=10)
        
        # 退出按鈕
        quit_btn = ttk.Button(parent, text="退出", command=self.root.quit, style='TButton')
        quit_btn.pack(side=tk.RIGHT, padx=10, ipadx=10)
        
        # 確保按鈕可見
        parent.tkraise()    
    def calculate_earth_pressure_coefficient(self):
        """計算土壓力係數 (Ka或Kp)"""
        phi = self.phi.get()
        theta = self.theta.get()
        delta = self.delta.get()
        i = self.i.get()
        mode = self.pressure_mode.get()
        
        self.log_step("\n1. 計算土壓力係數 " + ("Ka (主動)" if mode == "active" else "Kp (被動)"))
        self.log_step(f"   輸入參數: φ={phi}°, θ={theta}°, δ={delta}°, i={i}°")
        
        phi_rad = math.radians(phi)
        theta_rad = math.radians(theta)
        delta_rad = math.radians(delta)
        i_rad = math.radians(i)
        
        if i == 0:
            # 簡化公式 (當地表水平時)
            if mode == "active":
                k = (1 - math.sin(phi_rad)) / (1 + math.sin(phi_rad))
                self.log_step("   使用簡化主動土壓力公式:")
                self.log_step("   Ka = (1 - sinφ) / (1 + sinφ)")
            else:
                k = (1 + math.sin(phi_rad)) / (1 - math.sin(phi_rad))
                self.log_step("   使用簡化被動土壓力公式:")
                self.log_step("   Kp = (1 + sinφ) / (1 - sinφ)")
            
            self.log_step(f"   = (1 - sin{phi}°) / (1 + sin{phi}°)" if mode == "active" else 
                         f"   = (1 + sin{phi}°) / (1 - sin{phi}°)")
            self.log_step(f"   = {k:.4f}")
        else:
            # 完整庫倫公式
            if mode == "active":
                numerator = math.cos(phi_rad - theta_rad)**2
                denominator = math.cos(theta_rad)**2 * math.cos(delta_rad + theta_rad)
                sqrt_part = math.sqrt(
                    (math.sin(phi_rad + delta_rad) * math.sin(phi_rad - i_rad)) / 
                    (math.cos(delta_rad + theta_rad) * math.cos(theta_rad - i_rad))
                )
                k = numerator / (denominator * (1 + sqrt_part)**2)
                
                self.log_step("   使用完整庫倫主動土壓力公式:")
                self.log_step("   Ka = cos²(φ-θ) / [cos²θ·cos(δ+θ)·(1+√Q)²]")
                self.log_step("   其中 Q = [sin(φ+δ)·sin(φ-i)] / [cos(δ+θ)·cos(θ-i)]")
            else:
                numerator = math.cos(phi_rad + theta_rad)**2
                denominator = math.cos(theta_rad)**2 * math.cos(delta_rad - theta_rad)
                sqrt_part = math.sqrt(
                    (math.sin(phi_rad + delta_rad) * math.sin(phi_rad + i_rad)) / 
                    (math.cos(delta_rad - theta_rad) * math.cos(theta_rad - i_rad))
                )
                k = numerator / (denominator * (1 - sqrt_part)**2)
                
                self.log_step("   使用完整庫倫被動土壓力公式:")
                self.log_step("   Kp = cos²(φ+θ) / [cos²θ·cos(δ-θ)·(1-√Q)²]")
                self.log_step("   其中 Q = [sin(φ+δ)·sin(φ+i)] / [cos(δ-θ)·cos(θ-i)]")
            
            self.log_step(f"\n   逐步計算:")
            self.log_step(f"   - 分子 cos²({phi}{'-' if mode == 'active' else '+'}{theta}) = {numerator:.4f}")
            self.log_step(f"   - 分母 cos²{theta}·cos({delta}{'+' if mode == 'active' else '-'}{theta}) = {denominator:.4f}")
            self.log_step(f"   - Q = {sqrt_part**2:.4f} → √Q = {sqrt_part:.4f}")
            self.log_step(f"   - 最終 K{'a' if mode == 'active' else 'p'} = {k:.4f}")
        
        return k
    
    def calculate(self):
        """執行完整計算流程"""
        try:
            self.log_text.delete(1.0, tk.END)
            self.calculation_steps = []
            
            mode = "主動" if self.pressure_mode.get() == "active" else "被動"
            self.log_step(f"=== 土石籠擋土牆{mode}土壓力分析 ===")
            self.log_step(f"分析時間: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            
            # 1. 計算土壓力係數
            k = self.calculate_earth_pressure_coefficient()
            self.ka_result.config(text=f"{k:.4f} ({'Ka' if self.pressure_mode.get() == 'active' else 'Kp'})")
            
            # 2. 計算土壓力總和
            pa = 0.5 * self.gamma.get() * self.height.get()**2 * k
            self.pa_result.config(text=f"{pa:.2f}")
            
            self.log_step("\n2. 計算土壓力總和")
            self.log_step(f"   P{'a' if self.pressure_mode.get() == 'active' else 'p'} = 0.5 × γ × H² × K{'a' if self.pressure_mode.get() == 'active' else 'p'}")
            self.log_step(f"      = 0.5 × {self.gamma.get()} × {self.height.get()}² × {k:.4f}")
            self.log_step(f"      = {0.5 * self.gamma.get()} × {self.height.get()**2} × {k:.4f}")
            self.log_step(f"      = {pa:.2f} kN/m")
            
            # 3. 計算分力
            delta_rad = math.radians(self.delta.get())
            theta_rad = math.radians(self.theta.get())
            pv = pa * math.sin(delta_rad + theta_rad)
            ph = pa * math.cos(delta_rad + theta_rad)
            self.pv_result.config(text=f"{pv:.2f}")
            self.ph_result.config(text=f"{ph:.2f}")
            
            self.log_step("\n3. 計算土壓力分力")
            self.log_step(f"   Pv = P{'a' if self.pressure_mode.get() == 'active' else 'p'} × sin(δ+θ)")
            self.log_step(f"      = {pa:.2f} × sin({self.delta.get()}+{self.theta.get()})")
            self.log_step(f"      = {pa:.2f} × {math.sin(delta_rad + theta_rad):.4f}")
            self.log_step(f"      = {pv:.2f} kN/m (垂直分力)")
            
            self.log_step(f"\n   Ph = P{'a' if self.pressure_mode.get() == 'active' else 'p'} × cos(δ+θ)")
            self.log_step(f"      = {pa:.2f} × cos({self.delta.get()}+{self.theta.get()})")
            self.log_step(f"      = {pa:.2f} × {math.cos(delta_rad + theta_rad):.4f}")
            self.log_step(f"      = {ph:.2f} kN/m (水平分力)")
            
            # 4. 計算力矩
            mr = self.wall_weight.get() * self.width.get() / 2
            mf = ph * self.height.get() / 3
            self.mr_result.config(text=f"{mr:.2f}")
            self.mf_result.config(text=f"{mf:.2f}")
            
            self.log_step("\n4. 計算力矩")
            self.log_step("   扶正力矩 Mr = W × B/2")
            self.log_step(f"      = {self.wall_weight.get()} × {self.width.get()}/2")
            self.log_step(f"      = {mr:.2f} kN·m/m")
            
            self.log_step("\n   傾倒力矩 Mf = Ph × H/3")
            self.log_step(f"      = {ph:.2f} × {self.height.get()}/3")
            self.log_step(f"      = {mf:.2f} kN·m/m")
            
            # 5. 計算安全係數
            fs_overturning = mr / mf if mf != 0 else float('inf')
            fs_sliding = ((self.wall_weight.get() + pv) * self.friction_coef.get()) / ph if ph != 0 else float('inf')
            self.fs_overturning.config(text=f"{fs_overturning:.2f}")
            self.fs_sliding.config(text=f"{fs_sliding:.2f}")
            
            self.log_step("\n5. 計算安全係數")
            self.log_step("   抗傾覆安全係數 F = Mr / Mf")
            self.log_step(f"      = {mr:.2f} / {mf:.2f}")
            self.log_step(f"      = {fs_overturning:.2f}")
            
            self.log_step("\n   抗滑動安全係數 F = (W + Pv) × μ / Ph")
            self.log_step(f"      = ({self.wall_weight.get()} + {pv:.2f}) × {self.friction_coef.get()} / {ph:.2f}")
            self.log_step(f"      = {self.wall_weight.get() + pv:.2f} × {self.friction_coef.get()} / {ph:.2f}")
            self.log_step(f"      = {fs_sliding:.2f}")
            
            # 穩定性評估
            assessment = []
            if fs_overturning >= 1.5:
                assessment.append("抗傾覆: 安全")
            else:
                assessment.append("抗傾覆: 不安全!")
            
            if fs_sliding >= 1.5:
                assessment.append("抗滑動: 安全")
            else:
                assessment.append("抗滑動: 不安全!")
            
            color = "green" if (fs_overturning >= 1.5 and fs_sliding >= 1.5) else "red"
            self.stability_assessment.config(text=" | ".join(assessment), foreground=color)
            
            self.log_step("\n[穩定性評估]")
            self.log_step(" - " + " | ".join(assessment))
            self.log_step("\n=== 計算完成 ===")
            
        except Exception as e:
            messagebox.showerror("計算錯誤", f"發生錯誤: {str(e)}")
            self.log_step(f"\n[錯誤] 計算過程中發生錯誤: {str(e)}")
    
    def export_report(self):
        """導出Markdown格式計算書"""
        if not self.calculation_steps:
            messagebox.showwarning("警告", "請先執行計算再輸出報告")
            return
        
        file_path = filedialog.asksaveasfilename(
            defaultextension=".md",
            filetypes=[("Markdown文件", "*.md"), ("所有文件", "*.*")],
            title="保存計算書"
        )
        
        if not file_path:
            return
        
        try:
            with open(file_path, "w", encoding="utf-8") as f:
                mode = "主動" if self.pressure_mode.get() == "active" else "被動"
                f.write(f"# 土石籠擋土牆{mode}土壓力分析計算書\n\n")
                f.write(f"**生成時間**: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
                
                f.write("## 輸入參數\n")
                f.write(f"- 分析模式: {mode}土壓力分析\n")
                f.write(f"- 土壤內摩擦角 φ = {self.phi.get()}°\n")
                f.write(f"- 牆背摩擦角 δ = {self.delta.get()}°\n")
                f.write(f"- 牆背傾斜角 θ = {self.theta.get()}°\n")
                f.write(f"- 地表傾斜角 i = {self.i.get()}°\n")
                f.write(f"- 土壤飽和單位重 γ = {self.gamma.get()} kN/m³\n")
                f.write(f"- 擋土牆總重 W = {self.wall_weight.get()} kN/m\n")
                f.write(f"- 土石籠高度 H = {self.height.get()} m\n")
                f.write(f"- 土石籠寬度 B = {self.width.get()} m\n")
                f.write(f"- 摩擦係數 μ = {self.friction_coef.get()}\n\n")
                
                f.write("## 詳細計算過程\n```\n")
                for step in self.calculation_steps:
                    f.write(step + "\n")
                f.write("```\n\n")
                
                f.write("## 分析結果\n")
                f.write(f"- 土壓力係數 K{'a' if self.pressure_mode.get() == 'active' else 'p'} = {self.ka_result['text']}\n")
                f.write(f"- 土壓力總和 P{'a' if self.pressure_mode.get() == 'active' else 'p'} = {self.pa_result['text']} kN/m\n")
                f.write(f"- 垂直分力 Pv = {self.pv_result['text']} kN/m\n")
                f.write(f"- 水平分力 Ph = {self.ph_result['text']} kN/m\n")
                f.write(f"- 扶正力矩 Mr = {self.mr_result['text']} kN·m/m\n")
                f.write(f"- 傾倒力矩 Mf = {self.mf_result['text']} kN·m/m\n")
                f.write(f"- 抗傾覆安全係數 = {self.fs_overturning['text']}\n")
                f.write(f"- 抗滑動安全係數 = {self.fs_sliding['text']}\n")
                f.write(f"- **穩定性評估**: {self.stability_assessment['text']}\n")
                
            messagebox.showinfo("成功", f"計算書已成功導出至:\n{file_path}")
        except Exception as e:
            messagebox.showerror("導出錯誤", f"無法導出計算書:\n{str(e)}")
    
    def log_step(self, step):
        """記錄計算步驟"""
        self.calculation_steps.append(step)
        self.log_text.insert(tk.END, step + "\n")
        self.log_text.see(tk.END)
    
    def clear(self):
        """清除所有輸入和結果"""
        for var in [self.phi, self.delta, self.theta, self.i, self.gamma, 
                   self.wall_weight, self.height, self.width, self.friction_coef]:
            var.set(0)
        
        for label in [self.ka_result, self.pa_result, self.pv_result, self.ph_result,
                     self.mr_result, self.mf_result, self.fs_overturning, 
                     self.fs_sliding, self.stability_assessment]:
            label.config(text="")
        
        self.log_text.delete(1.0, tk.END)
        self.calculation_steps = []
        
if __name__ == "__main__":
    root = tk.Tk()
    app = GabionStabilityApp(root)
    root.mainloop()


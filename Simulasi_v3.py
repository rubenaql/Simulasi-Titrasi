import math
from dataclasses import dataclass

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

st.set_page_config(
    page_title="Simulator Titrasi",
    layout="wide",
    page_icon="⚗️",
    initial_sidebar_state="expanded",
)

# =========================
# FUNGSI KIMIA
# =========================
def hitung_kw(temp_c: float) -> float:
    return 10 ** (-14 + 0.031 * (temp_c - 25))

def pH_dari_H(H: float) -> float:
    H = max(H, 1e-14)
    return -math.log10(H)

def hitung_asam_kuat_basa_kuat(type_, c0, v0_liter, c_add, v_add_ml, Kw):
    v_add = v_add_ml / 1000
    Vt = v0_liter + v_add
    if Vt <= 0:
        return 7.0, "Volume total nol"
    eps_mol = 1e-12

    if type_ == "strongA_strongB":
        n_asam = c0 * v0_liter
        n_basa = c_add * v_add
        sisa_basa = n_basa - n_asam
        if abs(sisa_basa) < eps_mol:
            status = "Titik ekuivalen"
            H = math.sqrt(Kw)
        elif sisa_basa > 0:
            status = "Kelebihan basa"
            C_b = sisa_basa / Vt
            a = 1.0; b = C_b; c = -Kw
            H = (-b + math.sqrt(b*b - 4*a*c)) / (2*a)
            H = max(H, 1e-14)
        else:
            status = "Kelebihan asam"
            C_a = (-sisa_basa) / Vt
            a = 1.0; b = C_a; c = -Kw
            H = (-b + math.sqrt(b*b - 4*a*c)) / (2*a)
            H = max(H, 1e-14)
        return pH_dari_H(H), status

    if type_ == "strongB_strongA":
        n_basa = c0 * v0_liter
        n_asam = c_add * v_add
        sisa_asam = n_asam - n_basa
        if abs(sisa_asam) < eps_mol:
            status = "Titik ekuivalen"
            H = math.sqrt(Kw)
        elif sisa_asam > 0:
            status = "Kelebihan asam"
            C_a = sisa_asam / Vt
            a = 1.0; b = C_a; c = -Kw
            H = (-b + math.sqrt(b*b - 4*a*c)) / (2*a)
            H = max(H, 1e-14)
        else:
            status = "Kelebihan basa"
            C_b = (-sisa_asam) / Vt
            a = 1.0; b = C_b; c = -Kw
            H = (-b + math.sqrt(b*b - 4*a*c)) / (2*a)
            H = max(H, 1e-14)
        return pH_dari_H(H), status
    raise ValueError(f"Type tidak dikenal: {type_}")

def hitung_asam_lemah(c0, v0_liter, c_add, v_add_ml, pKa, Kw):
    Ka = 10 ** (-pKa)
    v_add = v_add_ml / 1000
    Vt = v0_liter + v_add
    nHA = c0 * v0_liter
    nOH = c_add * v_add
    eps = 1e-12

    if nOH < nHA - eps:
        sisa_HA = nHA - nOH
        terbentuk_A = nOH
        if terbentuk_A <= 0:
            Ca = nHA / Vt
            a = 1.0; b = Ka; c = -Ka * Ca
            H = (-b + math.sqrt(b*b - 4*a*c)) / (2*a)
            H = max(H, 1e-14)
            status = "Asam lemah (belum dititrasi)"
            return pH_dari_H(H), status
        else:
            pH = pKa + math.log10(terbentuk_A / sisa_HA)
            status = "Daerah buffer"
            return pH, status
    if abs(nOH - nHA) < eps:
        C_garam = nHA / Vt
        Kb = Kw / Ka
        OH = math.sqrt(Kb * C_garam)
        H = Kw / OH
        status = "Titik ekuivalen"
        return pH_dari_H(H), status
    kelebihan = nOH - nHA
    C_b = kelebihan / Vt
    a = 1.0; b = C_b; c = -Kw
    H = (-b + math.sqrt(b*b - 4*a*c)) / (2*a)
    H = max(H, 1e-14)
    status = "Kelebihan basa"
    return pH_dari_H(H), status

@dataclass
class Parameter:
    type_: str
    temp_c: float
    c0: float
    v0_ml: float
    c_add: float
    v_add_ml: float
    v_max: float
    pKa: float
    
    @property
    def v0_liter(self) -> float:
        return self.v0_ml / 1000.0

def hitung_ph(params: Parameter):
    Kw = hitung_kw(params.temp_c)
    v0_liter = params.v0_liter
    if params.type_ == "HCl_NaOH":
        pH, status = hitung_asam_kuat_basa_kuat("strongA_strongB", params.c0, v0_liter, params.c_add, params.v_add_ml, Kw)
    elif params.type_ == "NaOH_HCl":
        pH, status = hitung_asam_kuat_basa_kuat("strongB_strongA", params.c0, v0_liter, params.c_add, params.v_add_ml, Kw)
    else:
        pH, status = hitung_asam_lemah(params.c0, v0_liter, params.c_add, params.v_add_ml, params.pKa, Kw)
    pH = max(0.0, min(14.0, pH))
    return pH, status, Kw

# =========================
# FUNGSI WARNA INDIKATOR
# =========================
def get_indicator_color(pH, indicator):
    if indicator == "Phenolphthalein":
        if pH < 8.2:
            return "rgba(255, 255, 255, 0.1)"
        elif pH < 10.0:
            ratio = (pH - 8.2) / 1.8
            r = 255
            g = int(255 - 150 * ratio)
            b = int(255 - 75 * ratio)
            return f"#{r:02x}{g:02x}{b:02x}"
        else:
            return "#ff69b4"
    elif indicator == "Methyl Orange":
        if pH < 3.1:
            return "#ff0000"
        elif pH < 4.4:
            ratio = (pH - 3.1) / 1.3
            r = 255
            g = int(ratio * 255)
            b = 0
            return f"#{r:02x}{g:02x}{b:02x}"
        else:
            return "#ffff00"
    elif indicator == "Bromothymol Blue":
        if pH < 6.0:
            return "#ffff00"
        elif pH < 7.6:
            ratio = (pH - 6.0) / 1.6
            if ratio < 0.5:
                r = int(255 * (1 - 2 * ratio))
                g = 255
                b = 0
            else:
                r = 0
                g = int(255 * (2 - 2 * ratio))
                b = int(255 * (2 * ratio - 1))
            return f"#{r:02x}{g:02x}{b:02x}"
        else:
            return "#0000ff"
    return "#ffffff"

# =========================
# INISIALISASI SESSION STATE
# =========================
if 'v_add_slider' not in st.session_state:
    st.session_state.v_add_slider = 0
if 'v_add_number' not in st.session_state:
    st.session_state.v_add_number = 0.0

def update_number_from_slider():
    st.session_state.v_add_number = float(st.session_state.v_add_slider)

def update_slider_from_number():
    st.session_state.v_add_slider = int(round(st.session_state.v_add_number))

# =========================
# UI SIDEBAR
# =========================
st.title("Simulator Titrasi Interaktif")
st.markdown("Simulasi titrasi asam-basa dengan kurva pH dan perubahan warna indikator.")

with st.sidebar:
    st.header("Pengaturan Titrasi")
    type_ = st.selectbox(
        "Jenis Titrasi",
        options=["HCl_NaOH", "NaOH_HCl", "CH3COOH_NaOH"],
        format_func=lambda x: {
            "HCl_NaOH": "Asam Kuat (HCl) + Basa Kuat (NaOH)",
            "NaOH_HCl": "Basa Kuat (NaOH) + Asam Kuat (HCl)",
            "CH3COOH_NaOH": "Asam Lemah (CH3COOH) + Basa Kuat (NaOH)",
        }.get(x, x),
        key="jenis_titrasi"
    )
    st.markdown("---")
    st.subheader("Larutan Analit")
    c0 = st.number_input("Konsentrasi (M)", min_value=0.0, value=0.1, step=0.01, format="%.4f", key="c0")
    v0_ml = st.number_input("Volume (mL)", min_value=1.0, value=50.0, step=5.0, format="%.1f", key="v0_ml")
    
    st.subheader("Larutan Titran")
    c_add = st.number_input("Konsentrasi (M)", min_value=0.0, value=0.1, step=0.01, format="%.4f", key="c_add")
    v_max = st.slider("Volume Maksimum (mL)", min_value=10, max_value=100, value=100, key="v_max_slider")
    
    st.subheader("Parameter Tambahan")
    temp_c = st.slider("Suhu (°C)", min_value=20.0, max_value=30.0, value=25.0, step=0.5, key="temp_c_slider")
    
    # Slider volume ditambahkan
    st.slider(
        "Volume Ditambahkan (mL)",
        min_value=0,
        max_value=int(v_max),
        value=st.session_state.v_add_slider,
        step=1,
        key="v_add_slider",
        on_change=update_number_from_slider
    )
    st.number_input(
        "Volume Ditambahkan (mL) [Manual]",
        min_value=0.0,
        max_value=float(v_max),
        value=st.session_state.v_add_number,
        step=0.1,
        format="%.1f",
        key="v_add_number",
        on_change=update_slider_from_number
    )
    v_add_ml = float(st.session_state.v_add_slider)
    
    if type_ == "CH3COOH_NaOH":
        pKa = st.number_input("pKa Asam Lemah", min_value=0.0, value=4.76, step=0.1, format="%.2f", key="pKa")
    else:
        pKa = 4.76
    
    st.subheader("Indikator pH")
    indicator = st.selectbox(
        "Pilih indikator",
        options=["Phenolphthalein", "Methyl Orange", "Bromothymol Blue"],
        format_func=lambda x: {
            "Phenolphthalein": "Phenolphthalein (8.2-10)",
            "Methyl Orange": "Methyl Orange (3.1-4.4)",
            "Bromothymol Blue": "Bromothymol Blue (6.0-7.6)",
        }.get(x, x),
        key="indicator_select"
    )

# =========================
# PERHITUNGAN UTAMA
# =========================
params = Parameter(
    type_=type_, temp_c=temp_c, c0=c0, v0_ml=v0_ml,
    c_add=c_add, v_add_ml=v_add_ml, v_max=v_max, pKa=pKa
)

pH, status, Kw = hitung_ph(params)
Ve = (c0 * (v0_ml / 1000) / c_add) * 1000 if c_add > 0 else 0
solution_color = get_indicator_color(pH, indicator)

# Layout utama
left, right = st.columns([1, 2])
with left:
    st.subheader("Larutan")
    st.markdown(
        f"""
        <div style="width:220px; height:320px; border:2px solid #cccccc; border-radius:10px; margin:auto; 
                    background:{solution_color}; position:relative; overflow:hidden; box-shadow:0 4px 8px rgba(0,0,0,0.1);">
            <div style="position:absolute; bottom:10px; left:0; right:0; text-align:center; 
                        background:rgba(255,255,255,0.6); padding:5px; font-size:12px; font-weight:bold;">
                pH: {pH:.2f} | {indicator}<br>
                {status}
            </div>
        </div>
        """,
        unsafe_allow_html=True,
    )
    col1, col2 = st.columns(2)
    with col1:
        st.metric("pH", f"{pH:.2f}")
        st.metric("Volume Ekuivalen", f"{Ve:.2f} mL")
    with col2:
        st.metric("Status", status)
        st.metric("Volume Ditambahkan", f"{v_add_ml:.1f} mL")

# Kurva titrasi
vs = np.linspace(0, v_max, 250)
phs = []
for v in vs:
    tmp_params = Parameter(
        type_=type_, temp_c=temp_c, c0=c0, v0_ml=v0_ml,
        c_add=c_add, v_add_ml=float(v), v_max=v_max, pKa=pKa
    )
    curve_ph, _, _ = hitung_ph(tmp_params)
    phs.append(curve_ph)

fig = go.Figure()
fig.add_trace(go.Scatter(x=vs, y=phs, mode="lines", line=dict(width=4, color="#00ff99"), name="Kurva pH"))
if c_add > 0 and 0 <= Ve <= v_max:
    fig.add_vline(x=Ve, line_dash="dash", line_color="red", annotation_text="Titik Ekuivalen", annotation_position="top")
fig.add_trace(go.Scatter(x=[v_add_ml], y=[pH], mode="markers", marker=dict(size=15, color="red", symbol="circle"), name="Titik Saat Ini"))

if indicator == "Phenolphthalein":
    fig.add_hrect(y0=8.2, y1=10, line_width=0, fillcolor="pink", opacity=0.2, annotation_text="Range Phenolphthalein")
elif indicator == "Methyl Orange":
    fig.add_hrect(y0=3.1, y1=4.4, line_width=0, fillcolor="orange", opacity=0.2, annotation_text="Range Methyl Orange")
elif indicator == "Bromothymol Blue":
    fig.add_hrect(y0=6.0, y1=7.6, line_width=0, fillcolor="green", opacity=0.2, annotation_text="Range Bromothymol Blue")

fig.update_layout(
    template="plotly_white", height=500,
    xaxis_title="Volume Ditambahkan (mL)", yaxis_title="pH",
    yaxis=dict(range=[0, 14], gridcolor='lightgray'),
    xaxis=dict(gridcolor='lightgray'),
)
with right:
    st.subheader("Kurva Titrasi")
    st.plotly_chart(fig, use_container_width=True)

# Reaksi kimia sederhana
st.subheader("Reaksi Kimia")
if type_ in ["HCl_NaOH", "NaOH_HCl"]:
    st.latex(r"HCl + NaOH \rightarrow NaCl + H_2O")
else:
    st.latex(r"CH_3COOH + NaOH \rightarrow CH_3COONa + H_2O")
    st.latex(r"pH = pK_a + \log\frac{[A^-]}{[HA]}")

st.caption("Simulator Titrasi Interaktif - Kelompok 3 LPK")

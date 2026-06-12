# Simulator Titrasi Lanjutan berbasis Streamlit
import math
import time
from dataclasses import dataclass

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

st.set_page_config(
    page_title="Simulator Titrasi",
    layout="wide",
    page_icon="",
    initial_sidebar_state="expanded",
)

# =========================
# CSS untuk styling
# =========================
st.markdown(
    """
    <style>
    div[data-testid="stNumberInput"] div[data-baseweb="input"],
    div[data-testid="stNumberInput"] div[data-baseweb="input"] > div {
        border: none !important;
        border-color: transparent !important;
        box-shadow: none !important;
        background-color: transparent !important;
    }
    div[data-testid="stNumberInput"] div[data-baseweb="input"]:focus-within,
    div[data-testid="stNumberInput"] div[data-baseweb="input"] > div:focus-within {
        box-shadow: none !important;
        border: none !important;
    }
    div[data-testid="stNumberInput"] input {
        border: none !important;
        box-shadow: none !important;
        outline: none !important;
        background: transparent !important;
        padding: 0.25rem 0.5rem !important;
        width: auto !important;
        min-width: 80px !important;
    }
    div[data-testid="stNumberInput"] {
        width: auto !important;
    }
    
    .solution-container {
        width: 220px;
        height: 320px;
        border: 2px solid #cccccc;
        border-radius: 10px;
        margin: auto;
        position: relative;
        overflow: hidden;
        transition: background-color 0.3s ease;
        box-shadow: 0 4px 8px rgba(0,0,0,0.1);
    }
    .solution-label {
        position: absolute;
        bottom: 10px;
        left: 0;
        right: 0;
        text-align: center;
        background: rgba(255,255,255,0.6);
        padding: 5px;
        font-size: 12px;
        font-weight: bold;
        backdrop-filter: blur(4px);
    }
    </style>
    """,
    unsafe_allow_html=True,
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

def hitung_asam_diprotik(c0, v0_liter, c_add, v_add_ml, pKa1, pKa2, Kw):
    """
    Titrasi asam diprotik H2A dengan basa kuat NaOH.
    Mengembalikan pH dan status.
    """
    Ka1 = 10 ** (-pKa1)
    Ka2 = 10 ** (-pKa2)
    v_add = v_add_ml / 1000
    Vt = v0_liter + v_add
    n_H2A = c0 * v0_liter
    n_OH = c_add * v_add
    eps = 1e-12

    # Volume ekuivalen pertama dan kedua
    Ve1 = (n_H2A / c_add) * 1000  # mL
    Ve2 = 2 * Ve1

    # Tahap 1: Sebelum titik ekuivalen pertama (n_OH < n_H2A)
    if n_OH < n_H2A - eps:
        # Buffer H2A/HA- : H2A + OH- -> HA- + H2O
        sisa_H2A = n_H2A - n_OH
        terbentuk_HA = n_OH
        if terbentuk_HA <= 0:
            # Belum ada OH- sama sekali, hitung pH asam lemah diprotik (pendekatan H2A)
            # Gunakan rumus untuk asam lemah dengan Ka1 dominan
            Ca = n_H2A / Vt
            H = (-Ka1 + math.sqrt(Ka1*Ka1 + 4*Ka1*Ca)) / 2
            H = max(H, 1e-14)
            status = "Asam diprotik (belum dititrasi)"
            return pH_dari_H(H), status
        else:
            pH = pKa1 + math.log10(terbentuk_HA / sisa_H2A)
            status = "Buffer tahap 1 (H2A/HA-)"
            return pH, status

    # Tahap 2: Tepat titik ekuivalen pertama
    if abs(n_OH - n_H2A) < eps:
        # Garam NaHA, bersifat amfoter
        C_HA = n_H2A / Vt
        # pH untuk garam amfoter: [H+] = sqrt(Ka1*Ka2*C_HA / (Ka1 + C_HA)) tetapi pendekatan umum: [H+] = sqrt(Ka1*Ka2)
        # Untuk presisi lebih baik:
        if C_HA > 0:
            H = math.sqrt((Ka1 * Ka2 * C_HA + Ka1 * Kw) / (Ka1 + C_HA))
            H = max(H, 1e-14)
        else:
            H = math.sqrt(Ka1 * Ka2)
        status = "Titik ekuivalen pertama (NaHA)"
        return pH_dari_H(H), status

    # Tahap 3: Antara titik ekuivalen pertama dan kedua (n_H2A < n_OH < 2*n_H2A)
    if n_OH < 2 * n_H2A - eps:
        # Buffer HA-/A2- : HA- + OH- -> A2- + H2O
        sisa_HA = 2 * n_H2A - n_OH
        terbentuk_A2 = n_OH - n_H2A
        pH = pKa2 + math.log10(terbentuk_A2 / sisa_HA)
        status = "Buffer tahap 2 (HA-/A2-)"
        return pH, status

    # Tahap 4: Tepat titik ekuivalen kedua
    if abs(n_OH - 2 * n_H2A) < eps:
        # Garam Na2A, hidrolisis basa
        C_A2 = n_H2A / Vt
        Kb1 = Kw / Ka2
        OH = math.sqrt(Kb1 * C_A2)
        H = Kw / OH
        status = "Titik ekuivalen kedua (Na2A)"
        return pH_dari_H(H), status

    # Tahap 5: Kelebihan basa
    kelebihan_OH = n_OH - 2 * n_H2A
    C_b = kelebihan_OH / Vt
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
    pKa: float      # untuk monoprotik
    pKa1: float     # untuk diprotik
    pKa2: float
    
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
    elif params.type_ == "CH3COOH_NaOH":
        pH, status = hitung_asam_lemah(params.c0, v0_liter, params.c_add, params.v_add_ml, params.pKa, Kw)
    elif params.type_ == "AsamOksalat_NaOH":
        pH, status = hitung_asam_diprotik(params.c0, v0_liter, params.c_add, params.v_add_ml, params.pKa1, params.pKa2, Kw)
    else:
        raise ValueError("Titrasi tidak dikenal")
    pH = max(0.0, min(14.0, pH))
    return pH, status, Kw

# =========================
# FUNGSI WARNA INDIKATOR
# =========================
def get_indicator_color(pH, indicator):
    if indicator == "Phenolphthalein":
        if pH < 8.2:
            return "#ffffff"
        elif pH < 10.0:
            ratio = (pH - 8.2) / 1.8
            r = 255
            g = int(200 + (255 - 200) * ratio)
            b = int(200 + (255 - 200) * ratio)
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
            r = int(255 * (1 - ratio))
            g = 255
            b = int(255 * ratio)
            return f"#{r:02x}{g:02x}{b:02x}"
        else:
            return "#0000ff"
    return "#ffffff"

# =========================
# UI
# =========================
st.title("Simulator Titrasi Interaktif")
st.markdown("Simulasi titrasi asam-basa dengan kurva pH dan perubahan warna indikator.")

with st.sidebar:
    st.header("Pengaturan Titrasi")
    type_ = st.selectbox(
        "Jenis Titrasi",
        options=["HCl_NaOH", "NaOH_HCl", "CH3COOH_NaOH", "AsamOksalat_NaOH"],
        format_func=lambda x: {
            "HCl_NaOH": "Asam Kuat (HCl) + Basa Kuat (NaOH)",
            "NaOH_HCl": "Basa Kuat (NaOH) + Asam Kuat (HCl)",
            "CH3COOH_NaOH": "Asam Lemah (CH3COOH) + Basa Kuat (NaOH)",
            "AsamOksalat_NaOH": "Asam Oksalat (H2C2O4) + Basa Kuat (NaOH)",
        }.get(x, x),
        key="jenis_titrasi"
    )
    st.markdown("---")
    st.subheader("Larutan Analit")
    c0 = st.number_input("Konsentrasi (M)", min_value=0.0, value=0.1, step=0.01, format="%.4f", key="c0")
    v0_ml = st.number_input("Volume (mL)", min_value=1.0, value=50.0, step=5.0, format="%.1f", key="v0_ml")
    st.subheader("Larutan Titran")
    c_add = st.number_input("Konsentrasi (M)", min_value=0.0, value=0.1, step=0.01, format="%.4f", key="c_add")
    v_max = st.slider("Volume Maksimum (mL)", min_value=10, max_value=100, value=50, key="v_max_slider")
    st.subheader("Parameter Tambahan")
    temp_c = st.slider("Suhu (C)", min_value=20.0, max_value=30.0, value=25.0, step=0.5, key="temp_c_slider")
    v_add_ml = st.slider("Volume Ditambahkan (mL)", min_value=0, max_value=int(v_max), value=0, step=1, key="v_add_slider")
    
    if type_ == "CH3COOH_NaOH":
        pKa = st.number_input("pKa Asam Lemah", min_value=0.0, value=4.76, step=0.1, format="%.2f", key="pKa")
        pKa1, pKa2 = pKa, 0.0  # tidak dipakai
    elif type_ == "AsamOksalat_NaOH":
        pKa1 = st.number_input("pKa1 Asam Oksalat", min_value=0.0, value=1.25, step=0.1, format="%.2f", key="pKa1")
        pKa2 = st.number_input("pKa2 Asam Oksalat", min_value=0.0, value=4.14, step=0.1, format="%.2f", key="pKa2")
        pKa = 0.0
    else:
        pKa = 4.76
        pKa1, pKa2 = 0.0, 0.0
    
    st.markdown("---")
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
    st.markdown("---")
    st.subheader("Informasi Cepat")
    if c_add > 0:
        n_analit = c0 * (v0_ml / 1000)
        if type_ == "AsamOksalat_NaOH":
            ve1 = (n_analit / c_add) * 1000
            ve2 = 2 * ve1
            st.metric("Volume Ekuivalen 1", f"{ve1:.2f} mL")
            st.metric("Volume Ekuivalen 2", f"{ve2:.2f} mL")
        else:
            ve = (n_analit / c_add) * 1000
            st.metric("Volume Ekuivalen", f"{ve:.2f} mL")
    else:
        st.warning("Konsentrasi titran nol, tidak dapat menghitung volume ekuivalen.")
    
    # Rekomendasi indikator
    if type_ == "HCl_NaOH":
        rec_ind = "Bromothymol Blue atau Phenolphthalein"
        if indicator not in ["Bromothymol Blue", "Phenolphthalein"]:
            st.warning(f"Indikator {indicator} kurang sesuai.")
    elif type_ == "NaOH_HCl":
        rec_ind = "Methyl Orange atau Bromothymol Blue"
        if indicator not in ["Methyl Orange", "Bromothymol Blue"]:
            st.warning(f"Indikator {indicator} kurang sesuai.")
    elif type_ == "CH3COOH_NaOH":
        rec_ind = "Phenolphthalein"
        if indicator != "Phenolphthalein":
            st.warning("Indikator yang tepat adalah Phenolphthalein.")
    else:  # Asam Oksalat
        rec_ind = "Phenolphthalein (untuk titik ekuivalen kedua) atau Methyl Orange (untuk titik pertama)"
        st.info(f"Indikator yang disarankan: {rec_ind}")
    
    st.info(f"Indikator yang disarankan: {rec_ind}")
    st.caption("Geser slider dan amati perubahan pH serta warna larutan.")

# Parameter objek
params = Parameter(
    type_=type_, temp_c=temp_c, c0=c0, v0_ml=v0_ml,
    c_add=c_add, v_add_ml=v_add_ml, v_max=v_max,
    pKa=pKa, pKa1=pKa1, pKa2=pKa2
)

pH, status, Kw = hitung_ph(params)
# Hitung volume ekuivalen untuk ditampilkan di kurva
if c_add > 0:
    n_analit = c0 * (v0_ml / 1000)
    if type_ == "AsamOksalat_NaOH":
        Ve1 = (n_analit / c_add) * 1000
        Ve2 = 2 * Ve1
        Ve = Ve2  # untuk referensi
    else:
        Ve = (n_analit / c_add) * 1000
else:
    Ve = 0

solution_color = get_indicator_color(pH, indicator)

# Layout utama
left, right = st.columns([1, 2])
with left:
    st.subheader("Larutan")
    st.markdown(
        f"""
        <div class="solution-container" style="background: {solution_color};">
            <div class="solution-label">
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
        if type_ == "AsamOksalat_NaOH":
            st.metric("Volume Ekuivalen 1", f"{Ve1:.2f} mL")
            st.metric("Volume Ekuivalen 2", f"{Ve2:.2f} mL")
        else:
            st.metric("Volume Ekuivalen", f"{Ve:.2f} mL")
        st.metric("Suhu", f"{temp_c:.1f} C")
    with col2:
        st.metric("Status", status)
        st.metric("Volume Ditambahkan", f"{v_add_ml:.1f} mL")
        st.metric("Kw", f"{Kw:.2e}")
    
    with st.expander("Trayek Indikator"):
        if indicator == "Phenolphthalein":
            st.markdown("pH < 8.2 : tidak berwarna\n\npH 8.2-10 : pink\n\npH > 10 : merah muda")
        elif indicator == "Methyl Orange":
            st.markdown("pH < 3.1 : merah\n\npH 3.1-4.4 : jingga\n\npH > 4.4 : kuning")
        else:
            st.markdown("pH < 6.0 : kuning\n\npH 6.0-7.6 : hijau\n\npH > 7.6 : biru")

# Kurva titrasi
vs = np.linspace(0, v_max, 250)
phs = []
for v in vs:
    tmp_params = Parameter(
        type_=type_, temp_c=temp_c, c0=c0, v0_ml=v0_ml,
        c_add=c_add, v_add_ml=float(v), v_max=v_max,
        pKa=pKa, pKa1=pKa1, pKa2=pKa2
    )
    curve_ph, _, _ = hitung_ph(tmp_params)
    phs.append(curve_ph)

fig = go.Figure()
fig.add_trace(go.Scatter(x=vs, y=phs, mode="lines", line=dict(width=4, color="#00ff99"), name="Kurva pH"))

# Garis titik ekuivalen
if c_add > 0:
    if type_ == "AsamOksalat_NaOH":
        if 0 <= Ve1 <= v_max:
            fig.add_vline(x=Ve1, line_dash="dash", line_color="orange", annotation_text="VE1", annotation_position="top")
        if 0 <= Ve2 <= v_max:
            fig.add_vline(x=Ve2, line_dash="dash", line_color="red", annotation_text="VE2", annotation_position="top")
    else:
        if 0 <= Ve <= v_max:
            fig.add_vline(x=Ve, line_dash="dash", line_color="red", annotation_text="Titik Ekuivalen", annotation_position="top")

fig.add_trace(go.Scatter(x=[v_add_ml], y=[pH], mode="markers", marker=dict(size=15, color="red", symbol="circle"), name="Titik Saat Ini"))

# Area trayek indikator
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

# Reaksi dan perhitungan
st.subheader("Reaksi Kimia")
if type_ == "HCl_NaOH":
    st.latex(r"HCl + NaOH \rightarrow NaCl + H_2O")
    st.markdown("**Asam Kuat + Basa Kuat → Garam + Air**")
elif type_ == "NaOH_HCl":
    st.latex(r"NaOH + HCl \rightarrow NaCl + H_2O")
elif type_ == "CH3COOH_NaOH":
    st.latex(r"CH_3COOH + NaOH \rightarrow CH_3COONa + H_2O")
    st.latex(r"pH = pK_a + \log\frac{[A^-]}{[HA]}")
else:
    st.latex(r"H_2C_2O_4 + NaOH \rightarrow NaHC_2O_4 + H_2O")
    st.latex(r"NaHC_2O_4 + NaOH \rightarrow Na_2C_2O_4 + H_2O")
    st.latex(r"pH = pK_{a1} + \log\frac{[HA^-]}{[H_2A]}, \quad pH = pK_{a2} + \log\frac{[A^{2-}]}{[HA^-]}")

st.subheader("Perhitungan")
with st.expander("Lihat Langkah Perhitungan"):
    v0_liter = v0_ml / 1000
    mol_awal = c0 * v0_liter
    mol_titran = c_add * (v_add_ml / 1000)
    st.write(f"Volume analit = {v0_ml:.1f} mL = {v0_liter:.4f} L")
    st.write(f"Mol analit awal = {mol_awal:.6f} mol")
    st.write(f"Mol titran = {mol_titran:.6f} mol")
    if type_ == "AsamOksalat_NaOH":
        st.write(f"Titik ekuivalen 1: {Ve1:.2f} mL, titik ekuivalen 2: {Ve2:.2f} mL")
    else:
        st.write(f"Volume ekuivalen: {Ve:.2f} mL")
    st.write(f"pH terukur = {pH:.4f}")
    st.write(f"Status = {status}")

# Ekspor data
st.subheader("Ekspor Data")
export_df = pd.DataFrame({"Volume_mL": vs, "pH": phs})
with st.expander("Tabel Data Kurva"):
    st.dataframe(export_df.round(3), use_container_width=True, height=300)
csv = export_df.to_csv(index=False)
st.download_button("Unduh CSV", data=csv, file_name=f"titration_curve_{type_}.csv", mime="text/csv", key="download_csv")

# Titrasi otomatis
st.subheader("Titrasi Otomatis")
auto_speed = st.select_slider("Kecepatan simulasi", options=["Lambat", "Normal", "Cepat"], value="Normal", key="auto_speed")
delay = {"Lambat": 0.15, "Normal": 0.08, "Cepat": 0.03}[auto_speed]

if st.button("Mulai Simulasi Otomatis", type="primary", key="start_auto"):
    progress = st.progress(0)
    status_text = st.empty()
    chart_placeholder = st.empty()
    base_fig = go.Figure()
    base_fig.add_trace(go.Scatter(x=vs, y=phs, mode="lines", name="Kurva pH", line=dict(color="#00ff99", width=3)))
    base_fig.update_layout(template="plotly_white", height=450, yaxis=dict(range=[0, 14]), xaxis_title="Volume Ditambahkan (mL)", yaxis_title="pH")
    volumes = np.linspace(0, v_max, 50)
    for i, vol in enumerate(volumes):
        tmp_params = Parameter(
            type_=type_, temp_c=temp_c, c0=c0, v0_ml=v0_ml,
            c_add=c_add, v_add_ml=float(vol), v_max=v_max,
            pKa=pKa, pKa1=pKa1, pKa2=pKa2
        )
        temp_ph, temp_status, _ = hitung_ph(tmp_params)
        temp_fig = go.Figure(base_fig)
        temp_fig.add_trace(go.Scatter(x=[vol], y=[temp_ph], mode="markers", marker=dict(size=12, color="red"), name="Titik saat ini", showlegend=(i==0)))
        chart_placeholder.plotly_chart(temp_fig, use_container_width=True)
        progress.progress((i+1)/len(volumes))
        status_text.info(f"Volume: {vol:.1f} mL | pH: {temp_ph:.2f} | {temp_status}")
        time.sleep(delay)
    status_text.success("Simulasi selesai!")
    time.sleep(1)
    status_text.empty()

# Ringkasan
st.divider()
st.subheader("Ringkasan")
if type_ == "AsamOksalat_NaOH":
    ringkasan = pd.DataFrame({
        "Parameter": ["pH", "Status", "Volume Ekuivalen 1 (mL)", "Volume Ekuivalen 2 (mL)", "Indikator", "Suhu (C)", "Kw"],
        "Nilai": [round(pH,2), status, f"{Ve1:.2f}", f"{Ve2:.2f}", rec_ind, f"{temp_c:.1f}", f"{Kw:.2e}"]
    })
else:
    ringkasan = pd.DataFrame({
        "Parameter": ["pH", "Status", "Volume Ekuivalen (mL)", "Indikator", "Suhu (C)", "Kw"],
        "Nilai": [round(pH,2), status, f"{Ve:.2f}", rec_ind, f"{temp_c:.1f}", f"{Kw:.2e}"]
    })
st.table(ringkasan)
st.caption("Courtesy Of Kelompok 3 LPK | Simulator Titrasi Interaktif")

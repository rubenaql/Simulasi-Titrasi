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
# CSS KUSTOM (sama seperti sebelumnya)
# =========================
st.markdown(
    """
    <style>
    .status-box {
        background-color: rgba(128, 128, 128, 0.15);
        border-radius: 0.5rem;
        padding: 0.5rem;
        margin: 0.5rem 0;
    }
    .status-box p:first-child {
        margin: 0;
        font-size: 14px;
        font-weight: 600;
        color: inherit;
    }
    .status-box p:last-child {
        margin: 0;
        font-size: 22px;
        font-weight: 500;
        word-break: break-word;
        white-space: normal;
        color: inherit;
        line-height: 1.3;
    }
    .solution-info {
        background-color: rgba(0, 0, 0, 0.65) !important;
        color: white !important;
    }
    .recommendation-box {
        background-color: rgba(0, 128, 0, 0.1);
        border-left: 4px solid #00cc00;
        padding: 8px 12px;
        margin-top: 10px;
        border-radius: 5px;
        font-size: 13px;
    }
    .legend-box {
        background-color: rgba(100, 100, 100, 0.1);
        border: 1px solid #cccccc;
        border-radius: 5px;
        padding: 8px;
        margin-top: 8px;
        font-size: 11px;
    }
    </style>
    """,
    unsafe_allow_html=True,
)

# =========================
# FUNGSI KIMIA (sama seperti sebelumnya)
# =========================
def hitung_kw(temp_c: float) -> float:
    return 10 ** (-14 + 0.031 * (temp_c - 25))

def pH_dari_H(H: float) -> float:
    H = max(H, 1e-14)
    return -math.log10(H)

def hitung_asam_kuat_basa_kuat(type_, c0, v0_liter, c_add, v_add_ml, Kw):
    # ... (kode sama, tidak diubah)
    pass

def hitung_asam_oksalat(c0, v0_liter, c_add, v_add_ml, Kw):
    # ... (kode sama)
    pass

def hitung_boraks(c0, v0_liter, c_add, v_add_ml, Kw):
    # ... (kode sama)
    pass

def hitung_kompleksometri(c0, v0_liter, c_add, v_add_ml, Kf=10**10.7):
    # ... (kode sama)
    pass

def hitung_permanganometri(c0, v0_liter, c_add, v_add_ml, E0_Fe3_Fe2=0.77, E0_MnO4_Mn2=1.51):
    # ... (kode sama)
    pass

def get_indicator_color(pH, indicator):
    # ... (kode sama)
    pass

def get_kompleksometri_color(pCa, status):
    # ... (kode sama)
    pass

def get_permanganometri_color(E, status):
    # ... (kode sama)
    pass

def get_indicator_recommendation(jenis_titrasi):
    # ... (kode sama)
    pass

@dataclass
class Parameter:
    jenis_titrasi: str
    temp_c: float
    c0: float
    v0_ml: float
    c_add: float
    v_add_ml: float
    v_max: float
    pKa: float
    logKf: float = 10.7
    E0_Fe: float = 0.77
    E0_Mn: float = 1.51

@st.cache_data(ttl=3600, show_spinner=False)
def hitung_kurva(jenis_titrasi, temp_c, c0, v0_ml, c_add, v_max, pKa, logKf, E0_Fe, E0_Mn):
    vs = np.linspace(0, v_max, 250)
    ys = []
    for v in vs:
        tmp_params = Parameter(
            jenis_titrasi=jenis_titrasi,
            temp_c=temp_c,
            c0=c0,
            v0_ml=v0_ml,
            c_add=c_add,
            v_add_ml=float(v),
            v_max=v_max,
            pKa=pKa,
            logKf=logKf,
            E0_Fe=E0_Fe,
            E0_Mn=E0_Mn
        )
        y, _, _ = hitung_nilai(tmp_params)
        ys.append(y)
    return vs, ys

def hitung_nilai(params: Parameter):
    # ... (kode sama, tidak diubah)
    pass

# =========================
# INISIALISASI SESSION STATE DARI QUERY PARAMS
# =========================
# Nilai default
defaults = {
    "jenis_titrasi": "HCl_NaOH",
    "c0": 0.1,
    "v0_ml": 50.0,
    "c_add": 0.1,
    "v_max": 100,
    "temp_c": 25.0,
    "v_add": 0.0,
    "pKa": 4.76,
    "logKf": 10.7,
    "E0_Fe": 0.77,
    "E0_Mn": 1.51,
    "indicator_select": "Phenolphthalein"
}

# Baca query params
query = st.query_params
for key, default in defaults.items():
    if key in query:
        val = query[key]
        # Konversi tipe data
        if isinstance(default, float):
            val = float(val)
        elif isinstance(default, int):
            val = int(val)
        st.session_state[key] = val
    else:
        if key not in st.session_state:
            st.session_state[key] = default

# =========================
# UI SIDEBAR
# =========================
st.title("Simulator Titrasi Interaktif")
st.markdown("Simulasi berbagai jenis titrasi: asam-basa, kompleksometri, dan permanganometri.")

with st.sidebar:
    st.header("Pengaturan Titrasi")
    jenis_titrasi = st.selectbox(
        "Jenis Titrasi",
        options=[
            "HCl_NaOH",
            "NaOH_HCl",
            "CH3COOH_NaOH",
            "NaOH_AsamOksalat",
            "HCl_Boraks",
            "Kompleksometri_EDTA_Ca",
            "Permanganometri_Fe"
        ],
        format_func=lambda x: {
            "HCl_NaOH": "Asam Kuat (HCl) + Basa Kuat (NaOH)",
            "NaOH_HCl": "Basa Kuat (NaOH) + Asam Kuat (HCl)",
            "CH3COOH_NaOH": "Asam Lemah (CH3COOH) + Basa Kuat (NaOH)",
            "NaOH_AsamOksalat": "Basa Kuat (NaOH) + Asam Oksalat (H2C2O4)",
            "HCl_Boraks": "Asam Kuat (HCl) + Boraks (Na2B4O7)",
            "Kompleksometri_EDTA_Ca": "Kompleksometri: EDTA vs Ca2+ (pCa)",
            "Permanganometri_Fe": "Permanganometri: KMnO4 vs Fe2+ (potensial)"
        }.get(x, x),
        key="jenis_titrasi"
    )
    st.markdown("---")
    st.subheader("Larutan Analit")
    c0 = st.number_input(
        "Konsentrasi Analit (M)", min_value=0.0, value=st.session_state.c0, step=0.01, format="%.4f", key="c0",
        help="Konsentrasi zat yang akan dititrasi (dalam molar, M)"
    )
    v0_ml = st.number_input(
        "Volume Analit (mL)", min_value=1.0, value=st.session_state.v0_ml, step=5.0, format="%.1f", key="v0_ml",
        help="Volume larutan analit dalam mililiter"
    )

    st.subheader("Larutan Titran")
    c_add = st.number_input(
        "Konsentrasi Titran (M)", min_value=0.0, value=st.session_state.c_add, step=0.01, format="%.4f", key="c_add",
        help="Konsentrasi larutan penitrasi (dalam molar, M)"
    )
    v_max = st.number_input(
        "Volume Buret (mL)", min_value=10, max_value=500, value=st.session_state.v_max, step=10, key="v_max",
        help="Volume maksimum buret yang akan ditampilkan pada kurva"
    )

    st.subheader("Parameter Tambahan")
    temp_c = st.number_input(
        "Suhu Ruangan (°C)", min_value=0.0, max_value=50.0, value=st.session_state.temp_c, step=1.0, key="temp_c",
        help="Suhu mempengaruhi nilai Kw dan pH netral. Rentang 0-50°C."
    )
    
    if c_add <= 0:
        st.warning("⚠️ Konsentrasi titran harus > 0 untuk melakukan titrasi. Volume ditambahkan dinonaktifkan.")
        v_add_ml = 0.0
        disabled_add = True
    else:
        disabled_add = False
        v_add_ml = st.number_input(
            "Volume Ditambahkan (mL)", min_value=0.0, max_value=float(v_max), value=st.session_state.v_add, step=1.0, format="%.1f", key="v_add",
            help="Volume titran yang telah ditambahkan (dalam mL). Geser atau ketik langsung.",
            disabled=disabled_add
        )

    if jenis_titrasi == "CH3COOH_NaOH":
        pKa = st.number_input(
            "pKa Asam Lemah", min_value=0.0, value=st.session_state.pKa, step=0.1, format="%.2f", key="pKa",
            help="Nilai pKa asam asetat = 4,76. Untuk asam lain dapat disesuaikan."
        )
    else:
        pKa = 4.76
        st.session_state.pKa = pKa

    if jenis_titrasi == "Kompleksometri_EDTA_Ca":
        logKf = st.number_input(
            "log Kf (Ca-EDTA)", min_value=0.0, value=st.session_state.logKf, step=0.1, format="%.1f", key="logKf",
            help="Konstanta stabilitas kompleks Ca-EDTA (log Kf = 10,7 pada pH 10)"
        )
    else:
        logKf = 10.7
        st.session_state.logKf = logKf

    if jenis_titrasi == "Permanganometri_Fe":
        E0_Fe = st.number_input(
            "E0 Fe3+/Fe2+ (V)", min_value=0.0, value=st.session_state.E0_Fe, step=0.01, format="%.2f", key="E0_Fe",
            help="Potensial standar reduksi pasangan Fe3+/Fe2+ dalam suasana asam 1 M"
        )
        E0_Mn = st.number_input(
            "E0 MnO4-/Mn2+ (V)", min_value=0.0, value=st.session_state.E0_Mn, step=0.01, format="%.2f", key="E0_Mn",
            help="Potensial standar reduksi pasangan MnO4-/Mn2+ dalam suasana asam 1 M"
        )
    else:
        E0_Fe, E0_Mn = 0.77, 1.51
        st.session_state.E0_Fe = E0_Fe
        st.session_state.E0_Mn = E0_Mn

    if jenis_titrasi in ["HCl_NaOH", "NaOH_HCl", "CH3COOH_NaOH", "NaOH_AsamOksalat", "HCl_Boraks"]:
        st.subheader("Indikator pH")
        indicator = st.selectbox(
            "Pilih indikator",
            options=["Phenolphthalein", "Methyl Orange", "Bromothymol Blue"],
            format_func=lambda x: {
                "Phenolphthalein": "Phenolphthalein (8.2-10)",
                "Methyl Orange": "Methyl Orange (3.1-4.4)",
                "Bromothymol Blue": "Bromothymol Blue (6.0-7.6)",
            }.get(x, x),
            key="indicator_select",
            help="Pilih indikator untuk melihat perubahan warna sesuai trayek pH"
        )
        rec = get_indicator_recommendation(jenis_titrasi)
        st.markdown(f'<div class="recommendation-box">{rec}</div>', unsafe_allow_html=True)
    else:
        indicator = None
        rec = get_indicator_recommendation(jenis_titrasi)
        st.markdown(f'<div class="recommendation-box">{rec}</div>', unsafe_allow_html=True)

    # Tombol reset menggunakan query params
    if st.button("🔄 Reset ke Default", use_container_width=True):
        # Set query params ke nilai default
        for key, val in defaults.items():
            st.query_params[key] = str(val)
        st.rerun()

# =========================
# PERHITUNGAN UTAMA
# =========================
params = Parameter(
    jenis_titrasi=jenis_titrasi,
    temp_c=temp_c,
    c0=c0,
    v0_ml=v0_ml,
    c_add=c_add,
    v_add_ml=v_add_ml,
    v_max=v_max,
    pKa=pKa,
    logKf=logKf,
    E0_Fe=E0_Fe,
    E0_Mn=E0_Mn
)

nilai, status, satuan = hitung_nilai(params)
Ve = (c0 * (v0_ml / 1000) / c_add) * 1000 if c_add > 0 else 0

if c_add > 0 and Ve > v_max:
    st.sidebar.warning(f"⚠️ Volume ekuivalen teoritis ({Ve:.1f} mL) melebihi volume maksimum buret ({v_max} mL). Naikkan 'Volume Buret' untuk melihat titik ekuivalen pada kurva.")
elif c_add > 0:
    st.sidebar.success(f"📌 Volume ekuivalen teoritis: {Ve:.2f} mL")

# Tentukan warna larutan (sama seperti sebelumnya)
# ... (kode penentuan warna dan legend, tidak diubah)
# ... (layout kolom, kurva, ekspor data, reaksi, panduan)
# Pastikan semua bagian di bawah ini sama seperti kode sebelumnya, hanya bagian reset yang berbeda

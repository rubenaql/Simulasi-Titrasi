import math
from dataclasses import dataclass

import numpy as np
import plotly.graph_objects as go
import streamlit as st

st.set_page_config(
    page_title="Simulator Titrasi",
    layout="wide",
    page_icon="⚗️",
    initial_sidebar_state="expanded",
)

# =========================
# CSS KUSTOM UNTUK DARK MODE DAN LEGENDA
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
# FUNGSI KIMIA UMUM
# =========================
def hitung_kw(temp_c: float) -> float:
    return 10 ** (-14 + 0.031 * (temp_c - 25))

def pH_dari_H(H: float) -> float:
    H = max(H, 1e-14)
    return -math.log10(H)

# =========================
# TITRASI ASAM KUAT - BASA KUAT
# =========================
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
            a = 1.0
            b = C_b
            c = -Kw
            H = (-b + math.sqrt(b * b - 4 * a * c)) / (2 * a)
            H = max(H, 1e-14)
        else:
            status = "Kelebihan asam"
            C_a = (-sisa_basa) / Vt
            a = 1.0
            b = C_a
            c = -Kw
            H = (-b + math.sqrt(b * b - 4 * a * c)) / (2 * a)
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
            a = 1.0
            b = C_a
            c = -Kw
            H = (-b + math.sqrt(b * b - 4 * a * c)) / (2 * a)
            H = max(H, 1e-14)
        else:
            status = "Kelebihan basa"
            C_b = (-sisa_asam) / Vt
            a = 1.0
            b = C_b
            c = -Kw
            H = (-b + math.sqrt(b * b - 4 * a * c)) / (2 * a)
            H = max(H, 1e-14)
        return pH_dari_H(H), status
    raise ValueError(f"Type tidak dikenal: {type_}")

# =========================
# TITRASI BASA KUAT DENGAN ASAM DIPROTIK (ASAM OKSALAT)
# =========================
def hitung_asam_oksalat(c0, v0_liter, c_add, v_add_ml, Kw):
    Ka1 = 10 ** -1.25
    Ka2 = 10 ** -4.27
    v_add = v_add_ml / 1000
    Vt = v0_liter + v_add
    n_OH_awal = c0 * v0_liter
    n_H2A = c_add * v_add
    n_OH_sisa = n_OH_awal - 2 * n_H2A

    if abs(n_OH_sisa) < 1e-12:
        Kb = Kw / Ka2
        C_A2 = n_H2A / Vt
        OH = math.sqrt(Kb * C_A2)
        H = Kw / OH
        status = "Titik ekuivalen 2 (semua OH- habis)"
        return pH_dari_H(H), status
    elif n_OH_sisa > 0:
        C_OH = n_OH_sisa / Vt
        H = Kw / C_OH
        status = "Kelebihan basa (setelah titik ekuivalen)"
        return pH_dari_H(H), status
    else:
        n_H2A_awal = n_H2A
        if n_OH_awal <= n_H2A_awal - 1e-12:
            n_H2A_sisa = n_H2A_awal - n_OH_awal
            n_HA_terbentuk = n_OH_awal
            pH = -math.log10(Ka1) + math.log10(n_HA_terbentuk / n_H2A_sisa)
            status = "Daerah buffer (H2C2O4 / HC2O4-)"
            return pH, status
        else:
            n_HA_sisa = 2 * n_H2A_awal - n_OH_awal
            n_A2_terbentuk = n_OH_awal - n_H2A_awal
            pH = -math.log10(Ka2) + math.log10(n_A2_terbentuk / n_HA_sisa)
            status = "Daerah buffer (HC2O4- / C2O42-)"
            return pH, status

# =========================
# TITRASI HCl vs BORAKS
# =========================
def hitung_boraks(c0, v0_liter, c_add, v_add_ml, Kw):
    v_add = v_add_ml / 1000
    Vt = v0_liter + v_add
    n_boraks = c0 * v0_liter
    n_HCl = c_add * v_add
    n_HCl_sisa = n_HCl - 2 * n_boraks
    if abs(n_HCl_sisa) < 1e-12:
        C_H3BO3 = 4 * n_boraks / Vt
        Ka = 10 ** -9.24
        H = math.sqrt(Ka * C_H3BO3)
        H = max(H, 1e-14)
        status = "Titik ekuivalen (larutan H3BO3)"
        return pH_dari_H(H), status
    elif n_HCl_sisa > 0:
        C_HCl = n_HCl_sisa / Vt
        H = C_HCl
        status = "Kelebihan asam kuat"
        return pH_dari_H(H), status
    else:
        n_boraks_sisa = -n_HCl_sisa / 2
        C_OH = 2 * n_boraks_sisa / Vt
        H = Kw / C_OH
        status = "Kelebihan boraks (basa)"
        return pH_dari_H(H), status

# =========================
# TITRASI KOMPLEKSOMETRI (EDTA vs Ca2+)
# =========================
def hitung_kompleksometri(c0, v0_liter, c_add, v_add_ml, Kf=10**10.7):
    v_add = v_add_ml / 1000
    Vt = v0_liter + v_add
    n_Ca = c0 * v0_liter
    n_EDTA = c_add * v_add

    if n_EDTA >= n_Ca - 1e-12:
        if abs(n_EDTA - n_Ca) < 1e-12:
            C_Ca_total = n_Ca / Vt
            Ca = math.sqrt(C_Ca_total / Kf)
            status = "Titik ekuivalen"
        else:
            C_EDTA_lebih = (n_EDTA - n_Ca) / Vt
            Ca = (n_Ca / Vt) / (Kf * C_EDTA_lebih)
            status = "Kelebihan EDTA"
    else:
        sisa_Ca = (n_Ca - n_EDTA) / Vt
        Ca = sisa_Ca
        status = "Kelebihan Ca2+"
    Ca = max(Ca, 1e-14)
    pCa = -math.log10(Ca)
    return pCa, status

# =========================
# TITRASI PERMANGANOMETRI (KMnO4 vs Fe2+)
# =========================
def hitung_permanganometri(c0, v0_liter, c_add, v_add_ml, E0_Fe3_Fe2=0.77, E0_MnO4_Mn2=1.51):
    v_add = v_add_ml / 1000
    Vt = v0_liter + v_add
    n_Fe2 = c0 * v0_liter
    n_MnO4 = c_add * v_add
    n_Fe2_sisa = n_Fe2 - 5 * n_MnO4
    if abs(n_Fe2_sisa) < 1e-12:
        E = (5 * E0_MnO4_Mn2 + 1 * E0_Fe3_Fe2) / 6
        status = "Titik ekuivalen"
    elif n_Fe2_sisa > 0:
        n_Fe3_terbentuk = 5 * n_MnO4
        if n_Fe3_terbentuk <= 0:
            E = E0_Fe3_Fe2
        else:
            E = E0_Fe3_Fe2 + 0.0591 * math.log10(n_Fe3_terbentuk / n_Fe2_sisa)
        status = "Kelebihan Fe2+"
    else:
        n_MnO4_sisa = -n_Fe2_sisa / 5
        n_Mn2_terbentuk = n_Fe2 / 5
        E = E0_MnO4_Mn2 + (0.0591 / 5) * math.log10((n_MnO4_sisa / Vt) / (n_Mn2_terbentuk / Vt))
        status = "Kelebihan KMnO4"
    return E, status

# =========================
# FUNGSI WARNA INDIKATOR (asam-basa)
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
# FUNGSI WARNA KOMPLEKSOMETRI (EBT)
# =========================
def get_kompleksometri_color(pCa, status):
    if status == "Kelebihan Ca2+":
        return "#8B0000"      # merah anggur
    elif status == "Titik ekuivalen":
        return "#0000CD"      # biru medium
    else:
        return "#0000FF"      # biru terang

# =========================
# FUNGSI WARNA PERMANGANOMETRI
# =========================
def get_permanganometri_color(E, status):
    if status == "Kelebihan KMnO4":
        return "#CC00CC"      # ungu
    elif status == "Kelebihan Fe2+":
        return "#FFFFCC"      # kuning pucat
    else:
        return "#FFDDDD"      # merah muda pucat

# =========================
# REKOMENDASI INDIKATOR
# =========================
def get_indicator_recommendation(jenis_titrasi):
    recommendations = {
        "HCl_NaOH": "**Rekomendasi:** Bromothymol Blue (6.0-7.6) atau Phenolphthalein (8.2-10). Titik ekuivalen pH 7.",
        "NaOH_HCl": "**Rekomendasi:** Bromothymol Blue (6.0-7.6) atau Phenolphthalein (8.2-10). Titik ekuivalen pH 7.",
        "CH3COOH_NaOH": "**Rekomendasi:** Phenolphthalein (8.2-10). Titik ekuivalen pH > 7 (basa).",
        "NaOH_AsamOksalat": "**Rekomendasi:** Phenolphthalein (8.2-10). Titik ekuivalen kedua sekitar pH 8-9.",
        "HCl_Boraks": "**Rekomendasi:** Methyl Orange (3.1-4.4) atau Methyl Red. Titik ekuivalen pH sekitar 5.",
        "Kompleksometri_EDTA_Ca": "**Indikator khusus:** Eriochrome Black T (EBT). Perubahan warna: merah anggur → biru.",
        "Permanganometri_Fe": "**Autoindikator:** KMnO4 itu sendiri (ungu). Tidak perlu indikator tambahan."
    }
    return recommendations.get(jenis_titrasi, "")

# =========================
# STRUKTUR PARAMETER
# =========================
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

# =========================
# FUNGSI UTAMA PERHITUNGAN (dengan caching)
# =========================
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
    Kw = hitung_kw(params.temp_c)
    v0_liter = params.v0_ml / 1000.0
    jenis = params.jenis_titrasi

    if jenis == "HCl_NaOH":
        pH, status = hitung_asam_kuat_basa_kuat("strongA_strongB", params.c0, v0_liter, params.c_add, params.v_add_ml, Kw)
        return pH, status, "pH"
    elif jenis == "NaOH_HCl":
        pH, status = hitung_asam_kuat_basa_kuat("strongB_strongA", params.c0, v0_liter, params.c_add, params.v_add_ml, Kw)
        return pH, status, "pH"
    elif jenis == "CH3COOH_NaOH":
        Ka = 10 ** (-params.pKa)
        v_add = params.v_add_ml / 1000
        Vt = v0_liter + v_add
        nHA = params.c0 * v0_liter
        nOH = params.c_add * v_add
        eps = 1e-12
        if nOH < nHA - eps:
            sisa_HA = nHA - nOH
            terbentuk_A = nOH
            if terbentuk_A <= 0:
                Ca = nHA / Vt
                H = (-Ka + math.sqrt(Ka*Ka + 4*Ka*Ca)) / 2
                H = max(H, 1e-14)
                status = "Asam lemah (belum dititrasi)"
                return pH_dari_H(H), status, "pH"
            else:
                pH = params.pKa + math.log10(terbentuk_A / sisa_HA)
                status = "Daerah buffer"
                return pH, status, "pH"
        elif abs(nOH - nHA) < eps:
            C_garam = nHA / Vt
            Kb = Kw / Ka
            OH = math.sqrt(Kb * C_garam)
            H = Kw / OH
            status = "Titik ekuivalen"
            return pH_dari_H(H), status, "pH"
        else:
            kelebihan = nOH - nHA
            C_b = kelebihan / Vt
            H = Kw / C_b
            status = "Kelebihan basa"
            return pH_dari_H(H), status, "pH"
    elif jenis == "NaOH_AsamOksalat":
        pH, status = hitung_asam_oksalat(params.c0, v0_liter, params.c_add, params.v_add_ml, Kw)
        return pH, status, "pH"
    elif jenis == "HCl_Boraks":
        pH, status = hitung_boraks(params.c0, v0_liter, params.c_add, params.v_add_ml, Kw)
        return pH, status, "pH"
    elif jenis == "Kompleksometri_EDTA_Ca":
        pCa, status = hitung_kompleksometri(params.c0, v0_liter, params.c_add, params.v_add_ml, 10**params.logKf)
        return pCa, status, "pCa"
    elif jenis == "Permanganometri_Fe":
        E, status = hitung_permanganometri(params.c0, v0_liter, params.c_add, params.v_add_ml, params.E0_Fe, params.E0_Mn)
        return E, status, "E (V)"
    else:
        raise ValueError("Jenis titrasi tidak dikenal")

# =========================
# RESET SESSION STATE
# =========================
def reset_to_default():
    st.session_state.jenis_titrasi = "HCl_NaOH"
    st.session_state.c0 = 0.1
    st.session_state.v0_ml = 50.0
    st.session_state.c_add = 0.1
    st.session_state.v_max = 100
    st.session_state.temp_c = 25.0
    st.session_state.v_add = 0.0
    st.session_state.pKa = 4.76
    st.session_state.logKf = 10.7
    st.session_state.E0_Fe = 0.77
    st.session_state.E0_Mn = 1.51
    st.session_state.indicator_select = "Phenolphthalein"
    # Force rerun
    st.rerun()

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
        "Konsentrasi Analit (M)", min_value=0.0, value=0.1, step=0.01, format="%.4f", key="c0",
        help="Konsentrasi zat yang akan dititrasi (dalam molar, M)"
    )
    v0_ml = st.number_input(
        "Volume Analit (mL)", min_value=1.0, value=50.0, step=5.0, format="%.1f", key="v0_ml",
        help="Volume larutan analit dalam mililiter"
    )

    st.subheader("Larutan Titran")
    c_add = st.number_input(
        "Konsentrasi Titran (M)", min_value=0.0, value=0.1, step=0.01, format="%.4f", key="c_add",
        help="Konsentrasi larutan penitrasi (dalam molar, M)"
    )
    v_max = st.number_input(
        "Volume Buret (mL)", min_value=10, max_value=500, value=100, step=10, key="v_max",
        help="Volume maksimum buret yang akan ditampilkan pada kurva"
    )

    st.subheader("Parameter Tambahan")
    temp_c = st.number_input(
        "Suhu Ruangan (°C)", min_value=0.0, max_value=50.0, value=25.0, step=1.0, key="temp_c",
        help="Suhu mempengaruhi nilai Kw dan pH netral. Rentang 0-50°C."
    )
    
    # Peringatan jika konsentrasi titran nol
    if c_add <= 0:
        st.warning("⚠️ Konsentrasi titran harus > 0 untuk melakukan titrasi. Volume ditambahkan dinonaktifkan.")
        v_add_ml = 0.0
        disabled_add = True
    else:
        disabled_add = False
        v_add_ml = st.number_input(
            "Volume Ditambahkan (mL)", min_value=0.0, max_value=float(v_max), value=0.0, step=1.0, format="%.1f", key="v_add",
            help="Volume titran yang telah ditambahkan (dalam mL). Geser atau ketik langsung.",
            disabled=disabled_add
        )

    if jenis_titrasi == "CH3COOH_NaOH":
        pKa = st.number_input(
            "pKa Asam Lemah", min_value=0.0, value=4.76, step=0.1, format="%.2f", key="pKa",
            help="Nilai pKa asam asetat = 4,76. Untuk asam lain dapat disesuaikan."
        )
    else:
        pKa = 4.76

    if jenis_titrasi == "Kompleksometri_EDTA_Ca":
        logKf = st.number_input(
            "log Kf (Ca-EDTA)", min_value=0.0, value=10.7, step=0.1, format="%.1f", key="logKf",
            help="Konstanta stabilitas kompleks Ca-EDTA (log Kf = 10,7 pada pH 10)"
        )
    else:
        logKf = 10.7

    if jenis_titrasi == "Permanganometri_Fe":
        E0_Fe = st.number_input(
            "E0 Fe3+/Fe2+ (V)", min_value=0.0, value=0.77, step=0.01, format="%.2f", key="E0_Fe",
            help="Potensial standar reduksi pasangan Fe3+/Fe2+ dalam suasana asam 1 M"
        )
        E0_Mn = st.number_input(
            "E0 MnO4-/Mn2+ (V)", min_value=0.0, value=1.51, step=0.01, format="%.2f", key="E0_Mn",
            help="Potensial standar reduksi pasangan MnO4-/Mn2+ dalam suasana asam 1 M"
        )
    else:
        E0_Fe, E0_Mn = 0.77, 1.51

    # Indikator dan rekomendasi
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

    # Tombol reset
    if st.button("🔄 Reset ke Default", use_container_width=True):
        reset_to_default()

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

# Peringatan jika Ve > v_max
if c_add > 0 and Ve > v_max:
    st.sidebar.warning(f"⚠️ Volume ekuivalen teoritis ({Ve:.1f} mL) melebihi volume maksimum buret ({v_max} mL). Naikkan 'Volume Buret' untuk melihat titik ekuivalen pada kurva.")
elif c_add > 0:
    st.sidebar.success(f"📌 Volume ekuivalen teoritis: {Ve:.2f} mL")

# Tentukan warna larutan
if jenis_titrasi in ["HCl_NaOH", "NaOH_HCl", "CH3COOH_NaOH", "NaOH_AsamOksalat", "HCl_Boraks"] and indicator is not None:
    solution_color = get_indicator_color(nilai, indicator)
    info_indicator = f" | {indicator}"
elif jenis_titrasi == "Kompleksometri_EDTA_Ca":
    solution_color = get_kompleksometri_color(nilai, status)
    info_indicator = " | Indikator EBT"
    # Tambahkan legenda warna di sidebar
    with st.sidebar:
        st.markdown("""
        <div class="legend-box">
        <b>🎨 Legenda warna EBT:</b><br>
        🔴 Merah anggur → Kelebihan Ca²⁺<br>
    🔵 Biru medium → Titik ekuivalen<br>
    💙 Biru terang → Kelebihan EDTA
        </div>
        """, unsafe_allow_html=True)
elif jenis_titrasi == "Permanganometri_Fe":
    solution_color = get_permanganometri_color(nilai, status)
    info_indicator = " | Autoindikator KMnO₄"
    with st.sidebar:
        st.markdown("""
        <div class="legend-box">
        <b>🎨 Legenda warna KMnO₄:</b><br>
        🟣 Ungu → Kelebihan KMnO₄<br>
        🟡 Kuning pucat → Kelebihan Fe²⁺ (warna Fe³⁺)<br>
        🩷 Merah muda pucat → Titik ekuivalen
        </div>
        """, unsafe_allow_html=True)
else:
    solution_color = "#f0f0f0"
    info_indicator = ""

# Layout utama
left, right = st.columns([1, 2])
with left:
    st.subheader("Larutan")
    st.markdown(
        f"""
        <div style="width:220px; height:340px; border:2px solid #cccccc; border-radius:10px; margin:auto; 
                    background:{solution_color}; position:relative; overflow:hidden; box-shadow:0 4px 8px rgba(0,0,0,0.1);">
            <div class="solution-info" style="position:absolute; bottom:10px; left:0; right:0; text-align:center; 
                        background:rgba(0,0,0,0.65); padding:8px; font-size:11px; font-weight:bold;
                        word-wrap:break-word; white-space:normal; max-height:90px; overflow-y:auto; color:white;">
                {satuan.upper()}: {nilai:.3f}{info_indicator}<br>
                {status}
            </div>
        </div>
        """,
        unsafe_allow_html=True,
    )
    col1, col2 = st.columns(2)
    with col1:
        st.metric(satuan.upper(), f"{nilai:.3f} " + ("V" if satuan == "E (V)" else ""))
        st.metric("Volume Ekuivalen", f"{Ve:.2f} mL")
    with col2:
        st.metric("Volume Ditambahkan", f"{v_add_ml:.1f} mL")
        st.markdown(f"""
            <div class="status-box">
                <p>Status</p>
                <p>{status}</p>
            </div>
        """, unsafe_allow_html=True)

# Kurva titrasi (menggunakan cache)
if c_add > 0:
    vs, phs = hitung_kurva(jenis_titrasi, temp_c, c0, v0_ml, c_add, v_max, pKa, logKf, E0_Fe, E0_Mn)
else:
    vs = np.linspace(0, v_max, 250)
    phs = [0]*len(vs)

fig = go.Figure()
fig.add_trace(
    go.Scatter(
        x=vs, y=phs, mode="lines", line=dict(width=4, color="#00ff99"), name="Kurva"
    )
)
if c_add > 0 and 0 <= Ve <= v_max:
    fig.add_vline(
        x=Ve,
        line_dash="dash",
        line_color="red",
        annotation_text="Titik Ekuivalen",
        annotation_position="top",
    )
fig.add_trace(
    go.Scatter(
        x=[v_add_ml],
        y=[nilai],
        mode="markers",
        marker=dict(size=15, color="red", symbol="circle"),
        name="Titik Saat Ini",
    )
)

if satuan == "pH":
    y_title = "pH"
    y_range = [0, 14]
elif satuan == "pCa":
    y_title = "pCa"
    y_range = [0, 10]
else:
    y_title = "Potensial (V)"
    y_range = [0, 1.8]

fig.update_layout(
    template="plotly_white",
    height=500,
    xaxis_title="Volume Ditambahkan (mL)",
    yaxis_title=y_title,
    yaxis=dict(range=y_range, gridcolor="lightgray"),
    xaxis=dict(gridcolor="lightgray"),
)
with right:
    st.subheader("Kurva Titrasi")
    st.plotly_chart(fig, use_container_width=True)

# Reaksi kimia sederhana
st.subheader("Reaksi Kimia")
if jenis_titrasi == "HCl_NaOH":
    st.latex(r"HCl + NaOH \rightarrow NaCl + H_2O")
elif jenis_titrasi == "NaOH_HCl":
    st.latex(r"NaOH + HCl \rightarrow NaCl + H_2O")
elif jenis_titrasi == "CH3COOH_NaOH":
    st.latex(r"CH_3COOH + NaOH \rightarrow CH_3COONa + H_2O")
    st.latex(r"pH = pK_a + \log\frac{[A^-]}{[HA]}")
elif jenis_titrasi == "NaOH_AsamOksalat":
    st.latex(r"H_2C_2O_4 + 2NaOH \rightarrow Na_2C_2O_4 + 2H_2O")
    st.markdown("pKa₁ = 1,25; pKa₂ = 4,27")
elif jenis_titrasi == "HCl_Boraks":
    st.latex(r"Na_2B_4O_7 + 2HCl + 5H_2O \rightarrow 4H_3BO_3 + 2NaCl")
    st.markdown("Asam borat (H₃BO₃) pKa = 9,24")
elif jenis_titrasi == "Kompleksometri_EDTA_Ca":
    st.latex(r"Ca^{2+} + Y^{4-} \rightarrow CaY^{2-}")
    st.markdown(f"log Kf = {logKf} (pH 10, buffer amonia). Indikator EBT: merah anggur (Ca²⁺ bebas) → biru (kelebihan EDTA)")
elif jenis_titrasi == "Permanganometri_Fe":
    st.latex(r"MnO_4^- + 5Fe^{2+} + 8H^+ \rightarrow Mn^{2+} + 5Fe^{3+} + 4H_2O")
    st.markdown(f"E° Fe³⁺/Fe²⁺ = {E0_Fe} V, E° MnO₄⁻/Mn²⁺ = {E0_Mn} V. Autoindikator: ungu (kelebihan KMnO₄)")

# Panduan lengkap
with st.expander("📘 Panduan Lengkap & Cara Penggunaan", expanded=False):
    st.markdown("""
    ### Cara Menggunakan Simulator
    1. **Pilih jenis titrasi** di sidebar.
    2. **Atur konsentrasi dan volume** larutan analit (yang dititrasi) dan titran (penitrasi).
    3. **Tentukan volume maksimum buret** (tampilan kurva).
    4. **Masukkan volume titran yang ditambahkan** untuk melihat pH/pCa/potensial saat itu.
    5. **Pilih indikator** (untuk titrasi asam-basa) dan amati perubahan warna wadah larutan.
    6. **Gunakan tombol reset** untuk mengembalikan semua nilai ke default.
    
    ### Rumus yang Digunakan
    - **Asam kuat + basa kuat**: pH dihitung dari kelebihan H⁺ atau OH⁻, dengan koreksi autoprotolisis air.
    - **Asam lemah + basa kuat**: Persamaan Henderson-Hasselbalch untuk daerah buffer, dan hidrolisis garam pada titik ekuivalen.
    - **Asam oksalat + NaOH**: Pendekatan dua tahap dengan Ka1 dan Ka2.
    - **Boraks + HCl**: Reaksi menghasilkan H₃BO₃ (asam lemah).
    - **Kompleksometri**: Perhitungan [Ca²⁺] bebas berdasarkan konstanta stabilitas Kf.
    - **Permanganometri**: Potensial dihitung dengan persamaan Nernst untuk setiap daerah.
    
    ### Tips
    - Pastikan konsentrasi titran > 0.
    - Jika volume ekuivalen melebihi volume maksimum, naikkan nilai "Volume Buret".
    - Suhu 25°C memberikan nilai Kw = 1e-14 (pH netral 7). Suhu lain menggeser pH netral.
    - Indikator yang tepat akan berubah warna tepat di sekitar titik ekuivalen.
    """)

st.caption("Simulator Titrasi Lengkap - Kelompok 3 LPK")

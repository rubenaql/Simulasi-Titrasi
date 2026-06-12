import math
from dataclasses import dataclass

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

st.set_page_config(
    page_title="Simulator Titrasi - Belajar Kimia Interaktif",
    layout="wide",
    page_icon="🧪",
    initial_sidebar_state="expanded",
)

# =========================
# CSS KUSTOM
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
# FUNGSI KIMIA (disingkat untuk menghemat tempat, lengkap)
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
            status = "Titik ekuivalen (larutan garam netral)"
            H = math.sqrt(Kw)
        elif sisa_basa > 0:
            status = "Kelebihan basa (larutan bersifat basa)"
            C_b = sisa_basa / Vt
            a = 1.0; b = C_b; c = -Kw
            H = (-b + math.sqrt(b * b - 4 * a * c)) / (2 * a)
            H = max(H, 1e-14)
        else:
            status = "Kelebihan asam (larutan bersifat asam)"
            C_a = (-sisa_basa) / Vt
            a = 1.0; b = C_a; c = -Kw
            H = (-b + math.sqrt(b * b - 4 * a * c)) / (2 * a)
            H = max(H, 1e-14)
        return pH_dari_H(H), status
    if type_ == "strongB_strongA":
        n_basa = c0 * v0_liter
        n_asam = c_add * v_add
        sisa_asam = n_asam - n_basa
        if abs(sisa_asam) < eps_mol:
            status = "Titik ekuivalen (larutan garam netral)"
            H = math.sqrt(Kw)
        elif sisa_asam > 0:
            status = "Kelebihan asam (larutan bersifat asam)"
            C_a = sisa_asam / Vt
            a = 1.0; b = C_a; c = -Kw
            H = (-b + math.sqrt(b * b - 4 * a * c)) / (2 * a)
            H = max(H, 1e-14)
        else:
            status = "Kelebihan basa (larutan bersifat basa)"
            C_b = (-sisa_asam) / Vt
            a = 1.0; b = C_b; c = -Kw
            H = (-b + math.sqrt(b * b - 4 * a * c)) / (2 * a)
            H = max(H, 1e-14)
        return pH_dari_H(H), status
    raise ValueError(f"Type tidak dikenal: {type_}")

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
        status = "Titik ekuivalen (larutan natrium oksalat, bersifat basa)"
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
            status = "Daerah buffer (campuran asam oksalat dan oksalat asam)"
            return pH, status
        else:
            n_HA_sisa = 2 * n_H2A_awal - n_OH_awal
            n_A2_terbentuk = n_OH_awal - n_H2A_awal
            pH = -math.log10(Ka2) + math.log10(n_A2_terbentuk / n_HA_sisa)
            status = "Daerah buffer (campuran oksalat asam dan oksalat)"
            return pH, status

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
        status = "Titik ekuivalen (larutan asam borat, bersifat asam lemah)"
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

def hitung_kompleksometri(c0, v0_liter, c_add, v_add_ml, Kf=10**10.7):
    v_add = v_add_ml / 1000
    Vt = v0_liter + v_add
    n_Ca = c0 * v0_liter
    n_EDTA = c_add * v_add
    if n_EDTA >= n_Ca - 1e-12:
        if abs(n_EDTA - n_Ca) < 1e-12:
            C_Ca_total = n_Ca / Vt
            Ca = math.sqrt(C_Ca_total / Kf)
            status = "Titik ekuivalen (semua Ca²⁺ terikat EDTA)"
        else:
            C_EDTA_lebih = (n_EDTA - n_Ca) / Vt
            Ca = (n_Ca / Vt) / (Kf * C_EDTA_lebih)
            status = "Kelebihan EDTA (tidak ada ion Ca²⁺ bebas)"
    else:
        sisa_Ca = (n_Ca - n_EDTA) / Vt
        Ca = sisa_Ca
        status = "Kelebihan ion kalsium (Ca²⁺ bebas)"
    Ca = max(Ca, 1e-14)
    pCa = -math.log10(Ca)
    return pCa, status

def hitung_permanganometri(c0, v0_liter, c_add, v_add_ml, E0_Fe3_Fe2=0.77, E0_MnO4_Mn2=1.51):
    v_add = v_add_ml / 1000
    Vt = v0_liter + v_add
    n_Fe2 = c0 * v0_liter
    n_MnO4 = c_add * v_add
    n_Fe2_sisa = n_Fe2 - 5 * n_MnO4
    if abs(n_Fe2_sisa) < 1e-12:
        E = (5 * E0_MnO4_Mn2 + 1 * E0_Fe3_Fe2) / 6
        status = "Titik ekuivalen (semua Fe²⁺ teroksidasi)"
    elif n_Fe2_sisa > 0:
        n_Fe3_terbentuk = 5 * n_MnO4
        if n_Fe3_terbentuk <= 0:
            E = E0_Fe3_Fe2
        else:
            E = E0_Fe3_Fe2 + 0.0591 * math.log10(n_Fe3_terbentuk / n_Fe2_sisa)
        status = "Kelebihan besi(II) (Fe²⁺)"
    else:
        n_MnO4_sisa = -n_Fe2_sisa / 5
        n_Mn2_terbentuk = n_Fe2 / 5
        E = E0_MnO4_Mn2 + (0.0591 / 5) * math.log10((n_MnO4_sisa / Vt) / (n_Mn2_terbentuk / Vt))
        status = "Kelebihan permanganat (ungu)"
    return E, status

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

def get_kompleksometri_color(pCa, status):
    if "Kelebihan ion kalsium" in status:
        return "#8B0000"
    elif "Titik ekuivalen" in status:
        return "#0000CD"
    else:
        return "#0000FF"

def get_permanganometri_color(E, status):
    if "Kelebihan permanganat" in status:
        return "#CC00CC"
    elif "Kelebihan besi" in status:
        return "#FFFFCC"
    else:
        return "#FFDDDD"

def get_indicator_recommendation(jenis_titrasi):
    rec = {
        "HCl_NaOH": "✅ Rekomendasi: Phenolphthalein atau Bromothymol Blue.",
        "NaOH_HCl": "✅ Rekomendasi: Phenolphthalein atau Bromothymol Blue.",
        "CH3COOH_NaOH": "✅ Rekomendasi: Phenolphthalein.",
        "NaOH_AsamOksalat": "✅ Rekomendasi: Phenolphthalein.",
        "HCl_Boraks": "✅ Rekomendasi: Methyl Orange atau Methyl Red.",
        "Kompleksometri_EDTA_Ca": "🔬 Indikator khusus: Eriochrome Black T (EBT).",
        "Permanganometri_Fe": "🔬 Autoindikator: KMnO₄ sendiri (ungu)."
    }
    return rec.get(jenis_titrasi, "")

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
                status = "Daerah buffer (campuran asam asetat dan asetat)"
                return pH, status, "pH"
        elif abs(nOH - nHA) < eps:
            C_garam = nHA / Vt
            Kb = Kw / Ka
            OH = math.sqrt(Kb * C_garam)
            H = Kw / OH
            status = "Titik ekuivalen (larutan natrium asetat, bersifat basa)"
            return pH_dari_H(H), status, "pH"
        else:
            kelebihan = nOH - nHA
            C_b = kelebihan / Vt
            H = Kw / C_b
            status = "Kelebihan basa (larutan bersifat basa)"
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

# =========================
# INISIALISASI SESSION STATE DARI QUERY PARAMS
# =========================
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

query = st.query_params
for key, default_val in defaults.items():
    if key in query:
        val = query[key]
        if isinstance(default_val, float):
            val = float(val)
        elif isinstance(default_val, int):
            val = int(val)
        st.session_state[key] = val
    else:
        if key not in st.session_state:
            st.session_state[key] = default_val

# =========================
# UI SIDEBAR
# =========================
st.title("🧪 Simulator Titrasi Interaktif")
st.markdown("Belajar titrasi asam-basa, kompleksometri, dan permanganometri dengan mudah!")

with st.sidebar:
    st.header("⚙️ Pengaturan Percobaan")
    with st.expander("📖 Panduan Cepat", expanded=False):
        st.markdown("""
        **Apa itu titrasi?**  
        Titrasi adalah proses menambahkan larutan penitrasi ke larutan yang dititrasi hingga reaksi selesai (titik ekuivalen).  
        
        **Istilah penting:**  
        - **Larutan yang dititrasi** = zat yang akan diukur kadarnya.  
        - **Larutan penitrasi** = zat yang ditambahkan sedikit demi sedikit.  
        - **Titik ekuivalen** = saat jumlah penitrasi tepat habis bereaksi.  
        - **Indikator** = zat yang berubah warna di sekitar titik ekuivalen.  
        
        **Cara pakai:**  
        1. Pilih jenis titrasi.  
        2. Atur konsentrasi dan volume larutan.  
        3. Masukkan volume penitrasi yang ditambahkan.  
        4. Amati perubahan pH / warna / potensial.  
        5. Gunakan tombol **Reset** untuk memulai ulang.
        """)
    st.markdown("---")
    st.subheader("🧴 Larutan yang Dititrasi")
    st.number_input("Konsentrasi (M)", min_value=0.0, value=st.session_state.c0, step=0.01, format="%.4f", key="c0",
                    help="Semakin besar konsentrasi, semakin banyak zat yang perlu dititrasi.")
    st.number_input("Volume (mL)", min_value=1.0, value=st.session_state.v0_ml, step=5.0, format="%.1f", key="v0_ml",
                    help="Volume larutan yang akan dititrasi. Biasanya 50 mL.")
    st.subheader("💧 Larutan Penitrasi")
    st.number_input("Konsentrasi (M)", min_value=0.0, value=st.session_state.c_add, step=0.01, format="%.4f", key="c_add",
                    help="Larutan yang ditambahkan dari buret.")
    st.number_input("Volume Buret (mL)", min_value=10, max_value=500, value=st.session_state.v_max, step=10, key="v_max",
                    help="Kapasitas maksimum buret. Pastikan cukup untuk mencapai titik ekuivalen.")
    st.subheader("🌡️ Kondisi Percobaan")
    st.number_input("Suhu (°C)", min_value=0.0, max_value=50.0, value=st.session_state.temp_c, step=1.0, key="temp_c",
                    help="Suhu mempengaruhi nilai pH netral. Gunakan 25°C untuk kondisi standar.")
    if st.session_state.c_add <= 0:
        st.warning("⚠️ Konsentrasi penitrasi harus > 0 untuk memulai titrasi. Masukkan nilai positif.")
        st.session_state.v_add = 0.0
        disabled_add = True
    else:
        disabled_add = False
        st.number_input("Volume Ditambahkan (mL)", min_value=0.0, max_value=float(st.session_state.v_max), 
                        value=st.session_state.v_add, step=1.0, format="%.1f", key="v_add",
                        help="Jumlah penitrasi yang sudah ditambahkan. Geser atau ketik angka.", disabled=disabled_add)
    
    # Parameter khusus berdasarkan jenis titrasi
    # Kita gunakan st.session_state.jenis_titrasi langsung
    if st.session_state.jenis_titrasi == "CH3COOH_NaOH":
        st.number_input("pKa Asam Lemah", min_value=0.0, value=st.session_state.pKa, step=0.1, format="%.2f", key="pKa",
                        help="Untuk asam asetat, pKa = 4,76. Semakin kecil pKa, semakin kuat asam.")
    elif st.session_state.jenis_titrasi == "Kompleksometri_EDTA_Ca":
        st.number_input("log Kf (Ca-EDTA)", min_value=0.0, value=st.session_state.logKf, step=0.1, format="%.1f", key="logKf",
                        help="Nilai 10,7 untuk Ca-EDTA. Semakin besar, semakin stabil kompleks.")
    elif st.session_state.jenis_titrasi == "Permanganometri_Fe":
        st.number_input("E0 Fe³⁺/Fe²⁺ (V)", min_value=0.0, value=st.session_state.E0_Fe, step=0.01, format="%.2f", key="E0_Fe",
                        help="Potensial standar. Nilai 0,77 V untuk besi.")
        st.number_input("E0 MnO₄⁻/Mn²⁺ (V)", min_value=0.0, value=st.session_state.E0_Mn, step=0.01, format="%.2f", key="E0_Mn",
                        help="Potensial standar permanganat. Nilai 1,51 V.")
    
    # Indikator
    if st.session_state.jenis_titrasi in ["HCl_NaOH", "NaOH_HCl", "CH3COOH_NaOH", "NaOH_AsamOksalat", "HCl_Boraks"]:
        st.subheader("🌈 Indikator pH")
        st.selectbox("Pilih indikator",
                     options=["Phenolphthalein", "Methyl Orange", "Bromothymol Blue"],
                     format_func=lambda x: {
                         "Phenolphthalein": "Phenolphthalein (berubah pada pH 8,2-10 → pink)",
                         "Methyl Orange": "Methyl Orange (berubah pada pH 3,1-4,4 → merah ke kuning)",
                         "Bromothymol Blue": "Bromothymol Blue (berubah pada pH 6,0-7,6 → kuning ke biru)",
                     }.get(x, x),
                     key="indicator_select", help="Pilih indikator untuk melihat perubahan warna larutan.")
        rec = get_indicator_recommendation(st.session_state.jenis_titrasi)
        st.markdown(f'<div class="recommendation-box">{rec}</div>', unsafe_allow_html=True)
    else:
        rec = get_indicator_recommendation(st.session_state.jenis_titrasi)
        st.markdown(f'<div class="recommendation-box">{rec}</div>', unsafe_allow_html=True)
    
    if st.button("🔄 Reset ke Default", use_container_width=True):
        for key, val in defaults.items():
            st.query_params[key] = str(val)
        st.rerun()

# =========================
# AMBIL NILAI TERBARU DARI SESSION STATE
# =========================
jenis_titrasi = st.session_state.jenis_titrasi
c0 = st.session_state.c0
v0_ml = st.session_state.v0_ml
c_add = st.session_state.c_add
v_max = st.session_state.v_max
temp_c = st.session_state.temp_c
v_add_ml = st.session_state.v_add
pKa = st.session_state.pKa
logKf = st.session_state.logKf
E0_Fe = st.session_state.E0_Fe
E0_Mn = st.session_state.E0_Mn
indicator = st.session_state.indicator_select

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
    st.sidebar.warning(f"⚠️ Titik ekuivalen teoritis terjadi pada {Ve:.1f} mL, melebihi volume buret ({v_max} mL). Naikkan 'Volume Buret' agar kurva lengkap.")
elif c_add > 0:
    st.sidebar.success(f"📌 Titik ekuivalen diperkirakan pada {Ve:.2f} mL.")

# Tentukan warna larutan
if jenis_titrasi in ["HCl_NaOH", "NaOH_HCl", "CH3COOH_NaOH", "NaOH_AsamOksalat", "HCl_Boraks"] and indicator is not None:
    solution_color = get_indicator_color(nilai, indicator)
    info_indicator = f" | {indicator}"
elif jenis_titrasi == "Kompleksometri_EDTA_Ca":
    solution_color = get_kompleksometri_color(nilai, status)
    info_indicator = " | Indikator EBT"
    with st.sidebar:
        st.markdown("""
        <div class="legend-box">
        <b>🎨 Perubahan warna indikator EBT:</b><br>
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
        <b>🎨 Perubahan warna larutan:</b><br>
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
    st.subheader("🧪 Larutan dalam Labu")
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
        if satuan == "pH":
            st.metric("Tingkat Keasaman (pH)", f"{nilai:.2f}", help="pH < 7 = asam, pH = 7 netral, pH > 7 = basa")
        elif satuan == "pCa":
            st.metric("pCa", f"{nilai:.2f}", help="pCa = -log[Ca²⁺]. Semakin besar pCa, semakin sedikit ion kalsium bebas.")
        else:
            st.metric("Potensial (V)", f"{nilai:.3f} V", help="Potensial listrik larutan. Semakin tinggi, semakin kuat sifat oksidator.")
        st.metric("Titik Ekuivalen", f"{Ve:.2f} mL", help="Volume penitrasi yang dibutuhkan agar reaksi tepat habis.")
    with col2:
        st.metric("Volume Ditambahkan", f"{v_add_ml:.1f} mL", help="Jumlah penitrasi yang sudah dimasukkan.")
        st.markdown(f"""
            <div class="status-box">
                <p>Status Larutan</p>
                <p>{status}</p>
            </div>
        """, unsafe_allow_html=True)

# Kurva titrasi
if c_add > 0:
    vs, phs = hitung_kurva(jenis_titrasi, temp_c, c0, v0_ml, c_add, v_max, pKa, logKf, E0_Fe, E0_Mn)
else:
    vs = np.linspace(0, v_max, 250)
    phs = [0] * len(vs)

fig = go.Figure()
fig.add_trace(go.Scatter(x=vs, y=phs, mode="lines", line=dict(width=4, color="#00ff99"), name="Kurva"))
if c_add > 0 and 0 <= Ve <= v_max:
    fig.add_vline(x=Ve, line_dash="dash", line_color="red", annotation_text="Titik Ekuivalen", annotation_position="top")
fig.add_trace(go.Scatter(x=[v_add_ml], y=[nilai], mode="markers", marker=dict(size=15, color="red", symbol="circle"), name="Posisi saat ini"))

if satuan == "pH":
    y_title = "pH (tingkat keasaman)"
    y_range = [0, 14]
elif satuan == "pCa":
    y_title = "pCa (log konsentrasi kalsium bebas)"
    y_range = [0, 10]
else:
    y_title = "Potensial (V) - daya oksidasi"
    y_range = [0, 1.8]

fig.update_layout(
    template="plotly_white",
    height=500,
    xaxis_title="Volume Penitrasi Ditambahkan (mL)",
    yaxis_title=y_title,
    yaxis=dict(range=y_range, gridcolor="lightgray"),
    xaxis=dict(gridcolor="lightgray"),
)
with right:
    st.subheader("📈 Kurva Titrasi")
    st.plotly_chart(fig, use_container_width=True)
    st.caption("Grafik menunjukkan perubahan sifat larutan seiring penambahan penitrasi. Titik merah menunjukkan posisi saat ini.")

# Ekspor data
st.subheader("📥 Simpan Data Kurva")
export_df = pd.DataFrame({
    "Volume Penitrasi (mL)": vs,
    f"({satuan})": phs
})
st.dataframe(export_df.round(3), use_container_width=True, height=200)
csv = export_df.to_csv(index=False).encode('utf-8')
st.download_button(label="📥 Unduh CSV", data=csv, file_name=f"kurva_titrasi_{jenis_titrasi}.csv", mime="text/csv", help="Simpan data ke file CSV untuk analisis lanjutan.")

# Reaksi kimia
st.subheader("⚗️ Reaksi yang Terjadi")
if jenis_titrasi == "HCl_NaOH":
    st.latex(r"HCl + NaOH \rightarrow NaCl + H_2O")
    st.markdown("Asam klorida (HCl) bereaksi dengan natrium hidroksida (NaOH) menghasilkan garam dapur (NaCl) dan air. Titik ekuivalen terjadi pada pH 7.")
elif jenis_titrasi == "NaOH_HCl":
    st.latex(r"NaOH + HCl \rightarrow NaCl + H_2O")
    st.markdown("Sama seperti di atas, basa kuat (NaOH) dinetralkan oleh asam kuat (HCl).")
elif jenis_titrasi == "CH3COOH_NaOH":
    st.latex(r"CH_3COOH + NaOH \rightarrow CH_3COONa + H_2O")
    st.markdown("Asam asetat (cuka) bereaksi dengan NaOH membentuk natrium asetat dan air. Titik ekuivalen bersifat basa (pH > 7).")
elif jenis_titrasi == "NaOH_AsamOksalat":
    st.latex(r"H_2C_2O_4 + 2NaOH \rightarrow Na_2C_2O_4 + 2H_2O")
    st.markdown("Asam oksalat (asam diprotik) bereaksi dengan NaOH. Ada dua titik ekuivalen, yang kedua bersifat basa.")
elif jenis_titrasi == "HCl_Boraks":
    st.latex(r"Na_2B_4O_7 + 2HCl + 5H_2O \rightarrow 4H_3BO_3 + 2NaCl")
    st.markdown("Boraks (garam) bereaksi dengan HCl menghasilkan asam borat (asam lemah) dan garam. Titik ekuivalen bersifat asam lemah (pH sekitar 5).")
elif jenis_titrasi == "Kompleksometri_EDTA_Ca":
    st.latex(r"Ca^{2+} + Y^{4-} \rightarrow CaY^{2-}")
    st.markdown("Ion kalsium (Ca²⁺) membentuk kompleks stabil dengan EDTA. Indikator EBT berubah dari merah ke biru saat semua Ca²⁺ terikat.")
elif jenis_titrasi == "Permanganometri_Fe":
    st.latex(r"MnO_4^- + 5Fe^{2+} + 8H^+ \rightarrow Mn^{2+} + 5Fe^{3+} + 4H_2O")
    st.markdown("Permanganat (ungu) mengoksidasi besi(II) menjadi besi(III). Kelebihan permanganat memberikan warna ungu pada titik akhir.")

# Panduan lengkap
with st.expander("📘 Penjelasan Lengkap untuk Pemula", expanded=False):
    st.markdown("""
    ### Apa itu titrasi?
    Titrasi adalah teknik di laboratorium untuk menentukan kadar suatu zat. Caranya: meneteskan larutan yang sudah diketahui konsentrasinya (penitrasi) ke dalam larutan yang ingin diukur (dititrasi) sampai reaksi selesai.
    
    ### Komponen penting:
    - **Larutan yang dititrasi** : zat yang akan dicari kadarnya (misalnya cuka, air sadah).
    - **Larutan penitrasi** : zat yang konsentrasinya sudah diketahui (misalnya NaOH 0,1 M).
    - **Titik ekuivalen** : saat jumlah penitrasi tepat habis bereaksi. Biasanya ditandai dengan perubahan warna indikator.
    - **Indikator** : zat yang berubah warna di sekitar titik ekuivalen.
    
    ### Cara membaca hasil:
    - **pH** : menunjukkan tingkat keasaman. pH < 7 asam, > 7 basa, = 7 netral.
    - **pCa** : untuk titrasi kalsium; semakin tinggi pCa, semakin sedikit ion kalsium bebas.
    - **Potensial (V)** : untuk titrasi redoks; semakin tinggi potensial, larutan semakin oksidatif.
    
    ### Tips menggunakan simulator:
    - Pilih jenis titrasi sesuai percobaan yang ingin dipelajari.
    - Mulai dengan konsentrasi 0,1 M dan volume 50 mL.
    - Amati kurva: titik ekuivalen adalah bagian yang melonjak tajam.
    - Ganti indikator untuk melihat perubahan warna yang sesuai.
    - Jika kurva tidak menunjukkan lonjakan, naikkan volume maksimum buret.
    
    Selamat belajar! 🧪
    """)

st.caption("Simulator Titrasi Interaktif - Didesain untuk pembelajaran kimia yang mudah dipahami.")

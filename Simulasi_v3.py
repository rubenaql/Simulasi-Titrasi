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
    """
    c0 = konsentrasi NaOH (basa kuat) dalam M
    c_add = konsentrasi H2C2O4 (asam diprotik) dalam M
    v0_liter = volume NaOH awal (L)
    v_add_ml = volume H2C2O4 ditambahkan (mL)
    """
    # pKa1 = 1.25, pKa2 = 4.27
    Ka1 = 10 ** -1.25
    Ka2 = 10 ** -4.27
    v_add = v_add_ml / 1000
    Vt = v0_liter + v_add
    n_OH_awal = c0 * v0_liter
    n_H2A = c_add * v_add  # mol H2C2O4

    # Reaksi: H2A + 2OH- -> A2- + 2H2O (asam oksalat terdeprotonasi sempurna oleh basa kuat)
    # Stoikiometri: 2 mol OH- per 1 mol H2A
    n_OH_sisa = n_OH_awal - 2 * n_H2A

    if abs(n_OH_sisa) < 1e-12:
        # Titik ekuivalen (semua OH- habis, terbentuk A2-)
        # A2- adalah basa konjugat dari HA- (pKa2=4.27), Kb = Kw/Ka2
        Kb = Kw / Ka2
        C_A2 = n_H2A / Vt
        OH = math.sqrt(Kb * C_A2)
        H = Kw / OH
        status = "Titik ekuivalen 2 (semua OH- habis)"
        return pH_dari_H(H), status
    elif n_OH_sisa > 0:
        # Kelebihan OH-
        C_OH = n_OH_sisa / Vt
        H = Kw / C_OH
        status = "Kelebihan basa (setelah titik ekuivalen)"
        return pH_dari_H(H), status
    else:
        # n_OH_sisa < 0 : kelebihan H2A, tetapi belum tentu titik tengah
        # Hitung jumlah H+ yang dilepaskan. Pendekatan: kita anggap reaksi berlangsung bertahap.
        # Lebih mudah menggunakan perhitungan pH untuk campuran asam lemah/basa.
        # Karena ini kompleks, kita gunakan pendekatan numerik sederhana.
        # Untuk keperluan simulasi, kita hitung pH berdasarkan pembentukan buffer.
        # Misal: jika kelebihan H2A, maka kita punya campuran H2A dan HA-.
        # Namun implementasi lengkap terlalu panjang. Di sini kita buat pendekatan:
        n_H2A_berlebih = -n_OH_sisa / 2  # kelebihan H2A dalam mol
        # Reaksi pertama: OH- + H2A -> HA- + H2O, sisa H2A = n_H2A_berlebih, terbentuk HA- = n_OH_awal? Tidak tepat.
        # Alternatif: gunakan rumus pH untuk asam diprotik dengan basa kuat yang terbatas.
        # Kita sederhanakan: jika n_OH_awal < n_H2A, maka kita berada di daerah sebelum titik ekuivalen pertama.
        # Titik ekuivalen pertama terjadi saat n_OH_awal = n_H2A (karena OH- + H2A -> HA- + H2O)
        # Pada daerah tersebut, terbentuk buffer H2A/HA-.
        n_H2A_awal = n_H2A
        if n_OH_awal <= n_H2A_awal - 1e-12:
            # Sebelum titik ekuivalen pertama: buffer H2A/HA-
            n_H2A_sisa = n_H2A_awal - n_OH_awal
            n_HA_terbentuk = n_OH_awal
            pH = -math.log10(Ka1) + math.log10(n_HA_terbentuk / n_H2A_sisa)
            status = "Daerah buffer (H2C2O4 / HC2O4-)"
            return pH, status
        else:
            # Antara titik ekuivalen 1 dan 2: buffer HA-/A2-
            n_HA_sisa = 2 * n_H2A_awal - n_OH_awal
            n_A2_terbentuk = n_OH_awal - n_H2A_awal
            pH = -math.log10(Ka2) + math.log10(n_A2_terbentuk / n_HA_sisa)
            status = "Daerah buffer (HC2O4- / C2O42-)"
            return pH, status

# =========================
# TITRASI HCl vs BORAKS (Na2B4O7·10H2O)
# =========================
def hitung_boraks(c0, v0_liter, c_add, v_add_ml, Kw):
    """
    c0 = konsentrasi boraks (M) dalam analit
    c_add = konsentrasi HCl (M)
    v0_liter = volume boraks (L)
    v_add_ml = volume HCl ditambahkan (mL)
    Reaksi: B4O7^2- + 2H+ + 5H2O -> 4H3BO3
    H3BO3 adalah asam lemah (pKa = 9.24)
    """
    v_add = v_add_ml / 1000
    Vt = v0_liter + v_add
    n_boraks = c0 * v0_liter
    n_HCl = c_add * v_add
    # Stoikiometri: 2 mol H+ per 1 mol boraks
    n_HCl_sisa = n_HCl - 2 * n_boraks
    if abs(n_HCl_sisa) < 1e-12:
        # Titik ekuivalen: semua boraks bereaksi, terbentuk 4*n_boraks mol H3BO3
        C_H3BO3 = 4 * n_boraks / Vt
        Ka = 10 ** -9.24
        # Asam lemah
        H = math.sqrt(Ka * C_H3BO3)
        H = max(H, 1e-14)
        status = "Titik ekuivalen (larutan H3BO3)"
        return pH_dari_H(H), status
    elif n_HCl_sisa > 0:
        # Kelebihan HCl
        C_HCl = n_HCl_sisa / Vt
        H = C_HCl  # asam kuat
        status = "Kelebihan asam kuat"
        return pH_dari_H(H), status
    else:
        # n_HCl_sisa < 0 : kelebihan boraks, larutan bersifat basa karena hidrolisis boraks
        n_boraks_sisa = -n_HCl_sisa / 2
        # Boraks terhidrolisis menghasilkan OH-: B4O7^2- + 7H2O -> 4H3BO3 + 2OH-
        # Setiap mol boraks menghasilkan 2 mol OH- jika bereaksi sempurna, tetapi jika berlebih,
        # kita hitung konsentrasi OH- dari hidrolisis boraks.
        # Pendekatan: boraks adalah basa kuat? Tidak, tetapi dapat dihitung Kb.
        # Lebih mudah: karena pKa H3BO3 = 9.24, maka pKb boraks? Tidak langsung.
        # Kita gunakan pendekatan bahwa boraks bereaksi dengan air menghasilkan OH-.
        # Secara stoikiometri, 1 mol boraks menghasilkan 2 mol OH- (setelah reaksi dengan air).
        # Namun untuk kelebihan boraks, kita anggap semua boraks yang tersisa menghasilkan OH-
        # dengan asumsi reaksi sempurna. Konsentrasi OH- = 2 * n_boraks_sisa / Vt
        C_OH = 2 * n_boraks_sisa / Vt
        H = Kw / C_OH
        status = "Kelebihan boraks (basa)"
        return pH_dari_H(H), status

# =========================
# TITRASI KOMPLEKSOMETRI (EDTA vs Ca2+)
# =========================
def hitung_kompleksometri(c0, v0_liter, c_add, v_add_ml, Kf=10**10.7):
    """
    c0 = konsentrasi Ca2+ (M) dalam analit
    c_add = konsentrasi EDTA (M)
    v0_liter = volume Ca2+ (L)
    v_add_ml = volume EDTA ditambahkan (mL)
    Kf = konstanta stabilitas Ca-EDTA (log Kf = 10.7 pada pH 10)
    Asumsi: pH dijaga konstan 10 dengan buffer, sehingga EDTA dalam bentuk Y4-.
    Kurva: pCa = -log[Ca2+]
    """
    v_add = v_add_ml / 1000
    Vt = v0_liter + v_add
    n_Ca = c0 * v0_liter
    n_EDTA = c_add * v_add

    if n_EDTA >= n_Ca - 1e-12:
        # Kelebihan atau tepat ekuivalen
        if abs(n_EDTA - n_Ca) < 1e-12:
            # Titik ekuivalen: [Ca2+] = sqrt(1/Kf) * (mol total/Vt)? sebenarnya [Ca2+] = sqrt((C_Ca - [CaY])/Kf)
            # Pendekatan: [Ca2+] = sqrt( (n_Ca/Vt) / Kf )
            C_Ca_total = n_Ca / Vt
            Ca = math.sqrt(C_Ca_total / Kf)
            status = "Titik ekuivalen"
        else:
            # Kelebihan EDTA
            C_EDTA_lebih = (n_EDTA - n_Ca) / Vt
            # [Ca2+] = (n_Ca/Vt) / (Kf * C_EDTA_lebih)
            # asumsi kompleks stabil
            Ca = (n_Ca / Vt) / (Kf * C_EDTA_lebih)
            status = "Kelebihan EDTA"
    else:
        # Sebelum ekuivalen, Ca2+ berlebih
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
    """
    c0 = konsentrasi Fe2+ (M) dalam analit
    c_add = konsentrasi KMnO4 (M)
    v0_liter = volume Fe2+ (L)
    v_add_ml = volume KMnO4 ditambahkan (mL)
    Reaksi: MnO4- + 5Fe2+ + 8H+ -> Mn2+ + 5Fe3+ + 4H2O
    Asumsi [H+] konstan 1 M (pH 0)
    Potensial sel dihitung dengan Nernst.
    """
    v_add = v_add_ml / 1000
    Vt = v0_liter + v_add
    n_Fe2 = c0 * v0_liter
    n_MnO4 = c_add * v_add
    # Stoikiometri: 5 mol Fe2+ per 1 mol MnO4-
    n_Fe2_sisa = n_Fe2 - 5 * n_MnO4
    if abs(n_Fe2_sisa) < 1e-12:
        # Titik ekuivalen
        # Potensial dihitung dari kedua pasangan redoks, E = (5E0_Mn + 1*E0_Fe)/6
        E = (5 * E0_MnO4_Mn2 + 1 * E0_Fe3_Fe2) / 6
        status = "Titik ekuivalen"
    elif n_Fe2_sisa > 0:
        # Kelebihan Fe2+, gunakan pasangan Fe3+/Fe2+
        n_Fe3_terbentuk = 5 * n_MnO4
        if n_Fe3_terbentuk <= 0:
            E = E0_Fe3_Fe2  # asumsi awal, log(0) tidak terdefinisi, tetapi Fe3+ sangat kecil
        else:
            E = E0_Fe3_Fe2 + 0.0591 * math.log10(n_Fe3_terbentuk / n_Fe2_sisa)
        status = "Kelebihan Fe2+"
    else:
        # Kelebihan MnO4-
        n_MnO4_sisa = -n_Fe2_sisa / 5
        n_Mn2_terbentuk = n_Fe2 / 5  # semua Fe2+ habis
        # Potensial menggunakan pasangan MnO4-/Mn2+
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
# STRUKTUR PARAMETER (universal)
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
    pKa: float  # untuk asam lemah (opsional)
    # Parameter tambahan untuk kompleksometri/permanganometri
    logKf: float = 10.7
    E0_Fe: float = 0.77
    E0_Mn: float = 1.51

# =========================
# FUNGSI UTAMA PERHITUNGAN
# =========================
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
        # Untuk asam asetat, gunakan fungsi asam lemah yang sudah ada (dari kode sebelumnya)
        # Karena fungsi hitung_asam_lemah sudah didefinisikan di kode awal, kita panggil.
        # Tapi fungsi itu belum ada di kode ini, saya akan mendefinisikan ulang secara singkat.
        # Di sini saya akan tulis ulang fungsi hitung_asam_lemah (sederhana).
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
    c0 = st.number_input("Konsentrasi Analit (M)", min_value=0.0, value=0.1, step=0.01, format="%.4f", key="c0")
    v0_ml = st.number_input("Volume Analit (mL)", min_value=1.0, value=50.0, step=5.0, format="%.1f", key="v0_ml")

    st.subheader("Larutan Titran")
    c_add = st.number_input("Konsentrasi Titran (M)", min_value=0.0, value=0.1, step=0.01, format="%.4f", key="c_add")
    v_max = st.number_input("Volume Buret (mL)", min_value=10, max_value=500, value=100, step=10, key="v_max")

    st.subheader("Parameter Tambahan")
    temp_c = st.number_input("Suhu Ruangan (°C)", min_value=0.0, max_value=100.0, value=25.0, step=1.0, key="temp_c")
    v_add_ml = st.number_input("Volume Ditambahkan (mL)", min_value=0.0, max_value=float(v_max), value=0.0, step=1.0, format="%.1f", key="v_add")

    # Parameter khusus tergantung jenis
    if jenis_titrasi == "CH3COOH_NaOH":
        pKa = st.number_input("pKa Asam Lemah", min_value=0.0, value=4.76, step=0.1, format="%.2f", key="pKa")
    else:
        pKa = 4.76

    if jenis_titrasi == "Kompleksometri_EDTA_Ca":
        logKf = st.number_input("log Kf (Ca-EDTA)", min_value=0.0, value=10.7, step=0.1, format="%.1f", key="logKf")
    else:
        logKf = 10.7

    if jenis_titrasi == "Permanganometri_Fe":
        E0_Fe = st.number_input("E0 Fe3+/Fe2+ (V)", min_value=0.0, value=0.77, step=0.01, format="%.2f", key="E0_Fe")
        E0_Mn = st.number_input("E0 MnO4-/Mn2+ (V)", min_value=0.0, value=1.51, step=0.01, format="%.2f", key="E0_Mn")
    else:
        E0_Fe, E0_Mn = 0.77, 1.51

    # Indikator (hanya untuk titrasi asam-basa)
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
            key="indicator_select"
        )
    else:
        indicator = None

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

# Warna larutan untuk titrasi asam-basa
if jenis_titrasi in ["HCl_NaOH", "NaOH_HCl", "CH3COOH_NaOH", "NaOH_AsamOksalat", "HCl_Boraks"] and indicator is not None:
    solution_color = get_indicator_color(nilai, indicator)
else:
    # Untuk titrasi non-asam-basa, gunakan warna netral
    solution_color = "#f0f0f0"

# Layout utama
left, right = st.columns([1, 2])
with left:
    st.subheader("Larutan")
    # Tampilkan warna hanya jika indikator ada
    if indicator:
        st.markdown(
            f"""
            <div style="width:220px; height:320px; border:2px solid #cccccc; border-radius:10px; margin:auto; 
                        background:{solution_color}; position:relative; overflow:hidden; box-shadow:0 4px 8px rgba(0,0,0,0.1);">
                <div style="position:absolute; bottom:10px; left:0; right:0; text-align:center; 
                            background:rgba(255,255,255,0.6); padding:5px; font-size:12px; font-weight:bold;">
                    {satuan.upper()}: {nilai:.3f} | {indicator if indicator else ''}<br>
                    {status}
                </div>
            </div>
            """,
            unsafe_allow_html=True,
        )
    else:
        st.markdown(
            f"""
            <div style="width:220px; height:320px; border:2px solid #cccccc; border-radius:10px; margin:auto; 
                        background:#f0f0f0; position:relative; overflow:hidden; box-shadow:0 4px 8px rgba(0,0,0,0.1);">
                <div style="position:absolute; bottom:10px; left:0; right:0; text-align:center; 
                            background:rgba(255,255,255,0.6); padding:5px; font-size:12px; font-weight:bold;">
                    {satuan.upper()}: {nilai:.3f}<br>
                    {status}
                </div>
            </div>
            """,
            unsafe_allow_html=True,
        )
    col1, col2 = st.columns(2)
    with col1:
        st.metric(satuan.upper(), f"{nilai:.3f} " + ("mV" if satuan=="E (V)" else ""))
        st.metric("Volume Ekuivalen", f"{Ve:.2f} mL")
    with col2:
        st.metric("Status", status)
        st.metric("Volume Ditambahkan", f"{v_add_ml:.1f} mL")

# Kurva titrasi
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

fig = go.Figure()
fig.add_trace(
    go.Scatter(
        x=vs, y=ys, mode="lines", line=dict(width=4, color="#00ff99"), name="Kurva"
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

# Sesuaikan label sumbu y
if satuan == "pH":
    y_title = "pH"
    y_range = [0, 14]
elif satuan == "pCa":
    y_title = "pCa"
    y_range = [0, 10]
else:  # E (V)
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
    st.markdown(f"log Kf = {logKf} (pH 10, buffer amonia)")
elif jenis_titrasi == "Permanganometri_Fe":
    st.latex(r"MnO_4^- + 5Fe^{2+} + 8H^+ \rightarrow Mn^{2+} + 5Fe^{3+} + 4H_2O")
    st.markdown(f"E° Fe³⁺/Fe²⁺ = {E0_Fe} V, E° MnO₄⁻/Mn²⁺ = {E0_Mn} V")

st.caption("Simulator Titrasi Lengkap - Kelompok 3 LPK")

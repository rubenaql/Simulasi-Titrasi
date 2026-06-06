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
    page_icon="😜",
    initial_sidebar_state="expanded",
)

# =========================
# FUNGSI KIMIA (DIPERBAIKI)
# =========================

def hitung_kw(temp_c: float) -> float:
    return 10 ** (-14 + 0.031 * (temp_c - 25))

def pH_dari_H(H: float) -> float:
    H = max(H, 1e-14)
    return -math.log10(H)

def hitung_asam_kuat_basa_kuat(type_, c0, v0, c_add, v_add_ml, Kw):
    """Hitung pH titrasi asam kuat vs basa kuat dengan persamaan eksak."""
    v_add = v_add_ml / 1000
    Vt = v0 + v_add

    # --- Cegah ZeroDivisionError ---
    if Vt <= 0:
        return 7.0, "Volume total nol (periksa volume awal)"

    eps_mol = 1e-12

    if type_ == "strongA_strongB":
        n_asam = c0 * v0
        n_basa = c_add * v_add
        sisa_basa = n_basa - n_asam

        if abs(sisa_basa) < eps_mol:
            status = "Titik ekuivalen"
            # Pada titik ekuivalen, pH ditentukan oleh autoionisasi air
            H = math.sqrt(Kw)
        elif sisa_basa > 0:
            status = "Kelebihan basa"
            C_b = sisa_basa / Vt
            # Persamaan kuadrat untuk basa kuat dengan autoionisasi
            # [OH-] = C_b + [H+], dan [H+][OH-]=Kw
            # -> [H+]^2 + C_b[H+] - Kw = 0
            a = 1.0
            b = C_b
            c = -Kw
            H = (-b + math.sqrt(b*b - 4*a*c)) / (2*a)
            H = max(H, 1e-14)
        else:
            status = "Kelebihan asam"
            C_a = (-sisa_basa) / Vt
            # Persamaan kuadrat untuk asam kuat dengan autoionisasi
            a = 1.0
            b = C_a
            c = -Kw
            H = (-b + math.sqrt(b*b - 4*a*c)) / (2*a)
            H = max(H, 1e-14)

        return pH_dari_H(H), status

    if type_ == "strongB_strongA":
        n_basa = c0 * v0
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

def hitung_asam_lemah(c0, v0, c_add, v_add_ml, pKa, Kw):
    """Hitung pH titrasi asam lemah dengan basa kuat.
       Mengembalikan (pH, status)."""
    Ka = 10 ** (-pKa)
    v_add = v_add_ml / 1000
    Vt = v0 + v_add

    nHA = c0 * v0
    nOH = c_add * v_add
    eps = 1e-12

    # Kasus 1: belum ada basa atau sebelum ekuivalen
    if nOH < nHA - eps:
        sisa_HA = nHA - nOH
        terbentuk_A = nOH
        if terbentuk_A <= 0:
            # Hanya asam lemah, selesaikan persamaan kuadrat
            # [H+]^2 + Ka[H+] - Ka*Ca = 0
            Ca = nHA / Vt
            a = 1.0
            b = Ka
            c = -Ka * Ca
            H = (-b + math.sqrt(b*b - 4*a*c)) / (2*a)
            H = max(H, 1e-14)
            status = "Asam lemah (belum dititrasi)"
        else:
            # Daerah buffer, gunakan Henderson-Hasselbalch
            pH = pKa + math.log10(terbentuk_A / sisa_HA)
            status = "Daerah buffer"
            return pH, status
        return pH_dari_H(H), status

    # Kasus 2: titik ekuivalen
    if abs(nOH - nHA) < eps:
        C_garam = nHA / Vt
        Kb = Kw / Ka
        # Hidrolisis garam: [OH-] = sqrt(Kb * C_garam)
        OH = math.sqrt(Kb * C_garam)
        H = Kw / OH
        status = "Titik ekuivalen"
        return pH_dari_H(H), status

    # Kasus 3: kelebihan basa
    kelebihan = nOH - nHA
    C_b = kelebihan / Vt
    # Basa kuat + garam (garam tidak berpengaruh pada pH karena basa kuat dominan)
    # Tapi tetap perhitungkan autoionisasi air
    a = 1.0
    b = C_b
    c = -Kw
    H = (-b + math.sqrt(b*b - 4*a*c)) / (2*a)
    H = max(H, 1e-14)
    status = "Kelebihan basa"
    return pH_dari_H(H), status

@dataclass
class Parameter:
    type_: str
    temp_c: float
    c0: float
    v0: float
    c_add: float
    v_add_ml: float
    v_max: float
    pKa: float

def hitung_ph(params: Parameter):
    Kw = hitung_kw(params.temp_c)

    if params.type_ == "HCl_NaOH":
        pH, status = hitung_asam_kuat_basa_kuat(
            "strongA_strongB",
            params.c0,
            params.v0,
            params.c_add,
            params.v_add_ml,
            Kw,
        )
    elif params.type_ == "NaOH_HCl":
        pH, status = hitung_asam_kuat_basa_kuat(
            "strongB_strongA",
            params.c0,
            params.v0,
            params.c_add,
            params.v_add_ml,
            Kw,
        )
    else:  # CH3COOH_NaOH
        pH, status = hitung_asam_lemah(
            params.c0,
            params.v0,
            params.c_add,
            params.v_add_ml,
            params.pKa,
            Kw,
        )

    # Batasi pH antara 0 dan 14
    pH = max(0.0, min(14.0, pH))
    return pH, status, Kw

# =========================
# UI
# =========================

st.title(" Simulator Titrasi ")

st.markdown(
    """
Simulasi kimia interaktif menggunakan Streamlit.

Fitur:
- Asam kuat vs basa kuat
- Asam lemah vs basa kuat
- Titrasi otomatis
- Pembuatan kurva pH
- Simulasi warna indikator
- Ekspor CSV
"""
)

# =========================
# BILAH SAMPING
# =========================

with st.sidebar:
    st.header("⚙️ Pengaturan")

    type_ = st.selectbox(
        "Jenis Titrasi",
        [
            "HCl_NaOH",
            "NaOH_HCl",
            "CH3COOH_NaOH",
        ],
        format_func=lambda x: {
            "HCl_NaOH": "HCl + NaOH",
            "NaOH_HCl": "NaOH + HCl",
            "CH3COOH_NaOH": "CH3COOH + NaOH",
        }.get(x, x),
    )

    temp_c = st.slider("Suhu (°C)", 20.0, 30.0, 25.0)

    # --- Batasi nilai konsentrasi dan volume agar tidak negatif ---
    c0 = st.number_input("Konsentrasi Awal (M)", min_value=0.0, value=0.1, step=0.01, format="%.4f")
    v0 = st.number_input("Volume Awal (L)", min_value=0.001, value=0.05, step=0.01, format="%.4f")

    c_add = st.number_input("Konsentrasi Titran (M)", min_value=0.0, value=0.1, step=0.01, format="%.4f")

    v_max = st.slider("Volume Maksimum (mL)", 10, 100, 50)
    v_add_ml = st.slider("Volume Saat Ini (mL)", 0, int(v_max), 0)

    pKa = 4.76
    if type_ == "CH3COOH_NaOH":
        pKa = st.number_input("pKa", min_value=0.0, value=4.76, step=0.1, format="%.2f")

    indicator = st.selectbox(
        "Indikator",
        [
            "Phenolphthalein",
            "Methyl Orange",
            "Bromothymol Blue",
        ],
    )

params = Parameter(
    type_=type_,
    temp_c=temp_c,
    c0=c0,
    v0=v0,
    c_add=c_add,
    v_add_ml=v_add_ml,
    v_max=v_max,
    pKa=pKa,
)

# =========================
# HITUNG
# =========================

pH, status, Kw = hitung_ph(params)

# =========================
# VOLUME EKUIVALEN
# =========================
Ve = (c0 * v0 / c_add) * 1000 if c_add > 0 else 0

if type_ == "HCl_NaOH":
    indikator_rekomendasi = "Bromothymol Blue / Phenolphthalein"
elif type_ == "CH3COOH_NaOH":
    indikator_rekomendasi = "Phenolphthalein"
else:
    indikator_rekomendasi = "Methyl Orange"

# =========================
# WARNA INDIKATOR (DIPERBAIKI)
# =========================
# Gunakan border pada kotak agar indikator putih tetap terlihat
if indicator == "Phenolphthalein":
    solution_color = "#f0f0f0" if pH < 8.2 else "#ff69b4"  # abu-abu terang, bukan putih murni
elif indicator == "Methyl Orange":
    if pH < 3.1:
        solution_color = "#ff0000"
    elif pH > 4.4:
        solution_color = "#ffff00"
    else:
        solution_color = "#ff9900"
elif indicator == "Bromothymol Blue":
    if pH < 6:
        solution_color = "#ffff00"
    elif pH > 7.6:
        solution_color = "#0000ff"
    else:
        solution_color = "#00ff00"

# =========================
# TATA LETAK
# =========================

left, right = st.columns([1, 2])

with left:
    st.subheader("🧪 Larutan")

    # Tampilkan kotak warna dengan border abu-abu agar terlihat
    st.markdown(
        f"""
        <div style="
            width:220px;
            height:320px;
            border:4px solid #aaaaaa;
            border-radius:20px;
            background:{solution_color};
            margin:auto;
        "></div>
        """,
        unsafe_allow_html=True,
    )

    st.metric("pH", f"{pH:.2f}")
    st.metric("Status", status)
    st.metric("Volume Ekuivalen", f"{Ve:.2f} mL")
    st.metric("Volume Ditambahkan", f"{v_add_ml:.1f} mL")
    st.metric("Suhu", f"{temp_c:.1f} °C")

    st.markdown(
        f"""
        <div style="background:#111;color:#00ff00;font-size:30px;
        text-align:center;padding:10px;border-radius:10px;">
        pH {pH:.2f}
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.write(f"Kw ≈ {Kw:.2e}")

# =========================
# KURVA
# =========================

vs = np.linspace(0, v_max, 250)
phs = []

for v in vs:
    tmp_params = Parameter(
        type_=type_,
        temp_c=temp_c,
        c0=c0,
        v0=v0,
        c_add=c_add,
        v_add_ml=float(v),
        v_max=v_max,
        pKa=pKa,
    )
    curve_ph, _, _ = hitung_ph(tmp_params)
    phs.append(curve_ph)

fig = go.Figure()

fig.add_trace(
    go.Scatter(
        x=vs,
        y=phs,
        mode="lines",
        line=dict(width=4, color="#00ff99"),
        name="Kurva pH",
    )
)

if 0 <= Ve <= v_max:
    fig.add_vline(
        x=Ve,
        line_dash="dash",
        annotation_text="Titik Ekuivalen"
    )

fig.add_trace(
    go.Scatter(
        x=[v_add_ml],
        y=[pH],
        mode="markers",
        marker=dict(size=12, color="#4aa3ff"),
        name="Titik Saat Ini",
    )
)

fig.update_layout(
    template="plotly",
    height=500,
    xaxis_title="Volume Ditambahkan (mL)",
    yaxis_title="pH",
    yaxis=dict(range=[0, 14]),
)

with right:
    st.subheader("📈 Kurva Titrasi")
    st.plotly_chart(fig, use_container_width=True)

# =========================
# REAKSI
# =========================

st.subheader("⚗️ Reaksi Kimia")

if type_ == "HCl_NaOH":
    st.latex(r"HCl + NaOH \rightarrow NaCl + H_2O")
elif type_ == "NaOH_HCl":
    st.latex(r"NaOH + HCl \rightarrow NaCl + H_2O")
else:
    st.latex(r"CH_3COOH + NaOH \rightarrow CH_3COONa + H_2O")
    st.latex(r"pH = pK_a + \log\frac{[A^-]}{[HA]}")

st.subheader("🧮 Perhitungan")

with st.expander("Lihat Langkah Perhitungan"):
    mol_awal = c0 * v0
    mol_titran = c_add * (v_add_ml / 1000)

    st.write(f"Mol analit awal = {mol_awal:.5f} mol")
    st.write(f"Mol titran = {mol_titran:.5f} mol")

    if mol_awal > mol_titran + 1e-12:
        st.write(f"Sisa analit = {(mol_awal-mol_titran):.5f} mol")
    elif mol_titran > mol_awal + 1e-12:
        st.write(f"Kelebihan titran = {(mol_titran-mol_awal):.5f} mol")
    else:
        st.success("Titik ekuivalen")

st.sidebar.success(f"Indikator Disarankan: {indikator_rekomendasi}")

# =========================
# EKSPOR DATA
# =========================

st.subheader("📥 Ekspor Data")

export_df = pd.DataFrame({"Volume_mL": vs, "pH": phs})
with st.expander("📊 Tabel Data Kurva"):
    st.dataframe(export_df.round(3), use_container_width=True)

csv = export_df.to_csv(index=False)

st.download_button(
    label="Unduh CSV",
    data=csv,
    file_name="titration_curve.csv",
    mime="text/csv",
)

# =========================
# TITRASI OTOMATIS (DIPERBAIKI SEDIKIT)
# =========================

st.subheader("▶️ Titrasi Otomatis")

if st.button("Mulai Simulasi Otomatis"):
    progress = st.progress(0)
    chart_placeholder = st.empty()

    # Buat figure dasar sekali saja untuk efisiensi
    base_fig = go.Figure()
    base_fig.add_trace(go.Scatter(x=vs, y=phs, mode="lines", name="Kurva pH"))
    base_fig.update_layout(
        template="plotly",
        height=450,
        yaxis=dict(range=[0, 14]),
        xaxis_title="Volume Ditambahkan (mL)",
        yaxis_title="pH"
    )

    for i, vol in enumerate(np.linspace(0, v_max, 40)):
        tmp_params = Parameter(
            type_=type_,
            temp_c=temp_c,
            c0=c0,
            v0=v0,
            c_add=c_add,
            v_add_ml=float(vol),
            v_max=v_max,
            pKa=pKa,
        )
        temp_ph, _, _ = hitung_ph(tmp_params)

        # Salin figure agar tidak merusak yang asli
        temp_fig = go.Figure(base_fig)
        temp_fig.add_trace(
            go.Scatter(
                x=[vol],
                y=[temp_ph],
                mode="markers",
                marker=dict(size=12, color="red"),
                name="Titik saat ini"
            )
        )
        chart_placeholder.plotly_chart(temp_fig, use_container_width=True)
        progress.progress((i + 1) / 40)
        time.sleep(0.05)  # dikurangi sedikit agar lebih responsif

# =========================
# CATATAN AKHIR
# =========================

st.divider()

st.subheader("📋 Ringkasan")

ringkasan = pd.DataFrame(
    {
        "Parameter": ["pH", "Status", "Volume Ekuivalen", "Indikator"],
        "Nilai": [
            round(pH, 2),
            status,
            f"{Ve:.2f} mL",
            indikator_rekomendasi,
        ],
    }
)

st.table(ringkasan)

st.caption("Courtesy Of Kelompok 3 LPK")

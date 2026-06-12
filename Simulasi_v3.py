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
    page_icon="🧪",
    initial_sidebar_state="expanded",
)

# =========================
# CSS untuk menghapus border pada number input
# =========================
st.markdown(
    """
    <style>
    /* Hilangkan border dan background pada container utama Base Web input */
    div[data-testid="stNumberInput"] div[data-baseweb="input"],
    div[data-testid="stNumberInput"] div[data-baseweb="input"] > div {
        border: none !important;
        border-color: transparent !important;
        box-shadow: none !important;
        background-color: transparent !important;
    }
    
    /* Hilangkan ring/border saat input sedang fokus (diklik) */
    div[data-testid="stNumberInput"] div[data-baseweb="input"]:focus-within,
    div[data-testid="stNumberInput"] div[data-baseweb="input"] > div:focus-within {
        box-shadow: none !important;
        border: none !important;
    }

    /* Hilangkan border pada text input itu sendiri */
    div[data-testid="stNumberInput"] input {
        border: none !important;
        box-shadow: none !important;
        outline: none !important;
        background: transparent !important;
        padding: 0.25rem 0.5rem !important;
        width: auto !important;
        min-width: 80px !important;
    }
    
    /* Opsional: sesuaikan container agar tidak melebar */
    div[data-testid="stNumberInput"] {
        width: auto !important;
    }
    
    /* Style untuk wadah larutan */
    .solution-container {
        width: 220px;
        height: 320px;
        border: 2px solid #cccccc;
        border-radius: 10px;
        margin: auto;
        position: relative;
        overflow: hidden;
        transition: background-color 0.3s ease;
    }
    
    .solution-label {
        position: absolute;
        bottom: 10px;
        left: 0;
        right: 0;
        text-align: center;
        background: rgba(255,255,255,0.8);
        padding: 5px;
        font-size: 12px;
        font-weight: bold;
    }
    </style>
    """,
    unsafe_allow_html=True,
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

    if Vt <= 0:
        return 7.0, "Volume total nol (periksa volume awal)"

    eps_mol = 1e-12

    if type_ == "strongA_strongB":
        n_asam = c0 * v0
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
            H = (-b + math.sqrt(b*b - 4*a*c)) / (2*a)
            H = max(H, 1e-14)
        else:
            status = "Kelebihan asam"
            C_a = (-sisa_basa) / Vt
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

    if nOH < nHA - eps:
        sisa_HA = nHA - nOH
        terbentuk_A = nOH
        if terbentuk_A <= 0:
            Ca = nHA / Vt
            a = 1.0
            b = Ka
            c = -Ka * Ca
            H = (-b + math.sqrt(b*b - 4*a*c)) / (2*a)
            H = max(H, 1e-14)
            status = "Asam lemah (belum dititrasi)"
        else:
            pH = pKa + math.log10(terbentuk_A / sisa_HA)
            status = "Daerah buffer"
            return pH, status
        return pH_dari_H(H), status

    if abs(nOH - nHA) < eps:
        C_garam = nHA / Vt
        Kb = Kw / Ka
        OH = math.sqrt(Kb * C_garam)
        H = Kw / OH
        status = "Titik ekuivalen"
        return pH_dari_H(H), status

    kelebihan = nOH - nHA
    C_b = kelebihan / Vt
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
    else:
        pH, status = hitung_asam_lemah(
            params.c0,
            params.v0,
            params.c_add,
            params.v_add_ml,
            params.pKa,
            Kw,
        )

    pH = max(0.0, min(14.0, pH))
    return pH, status, Kw

# =========================
# FUNGSI WARNA INDIKATOR (DIPERBAIKI)
# =========================

def get_indicator_color(pH, indicator):
    """Mengembalikan kode warna hex berdasarkan pH dan indikator yang dipilih"""
    if indicator == "Phenolphthalein":
        # Phenolphthalein: tidak berwarna di pH < 8.2, pink di pH > 10, transisi 8.2-10
        if pH < 8.2:
            return "#f0f0f0"  # Tidak berwarna/bening
        elif pH < 10:
            # Gradasi pink
            ratio = (pH - 8.2) / 1.8
            r = int(255 - (255 - 255) * ratio)
            g = int(240 - (240 - 105) * ratio)
            b = int(240 - (240 - 180) * ratio)
            return f"#{r:02x}{g:02x}{b:02x}"
        else:
            return "#ff69b4"  # Pink terang
    elif indicator == "Methyl Orange":
        # Methyl Orange: merah di pH < 3.1, kuning di pH > 4.4, oranye di antaranya
        if pH < 3.1:
            return "#ff0000"  # Merah
        elif pH < 4.4:
            # Gradasi merah ke kuning
            ratio = (pH - 3.1) / 1.3
            r = int(255 - 255 * ratio)
            g = int(0 + 255 * ratio)
            b = 0
            return f"#{r:02x}{g:02x}{b:02x}"
        else:
            return "#ffff00"  # Kuning
    elif indicator == "Bromothymol Blue":
        # Bromothymol Blue: kuning di pH < 6.0, biru di pH > 7.6, hijau di antaranya
        if pH < 6.0:
            return "#ffff00"  # Kuning
        elif pH < 7.6:
            # Gradasi kuning ke biru
            ratio = (pH - 6.0) / 1.6
            r = int(255 - 255 * ratio)
            g = int(255 - (255 - 0) * ratio)
            b = int(0 + 255 * ratio)
            return f"#{r:02x}{g:02x}{b:02x}"
        else:
            return "#0000ff"  # Biru
    return "#f0f0f0"

def get_rgb_background(color_hex):
    """Mengubah warna hex ke format rgba untuk background gradient"""
    color_hex = color_hex.lstrip('#')
    r = int(color_hex[0:2], 16)
    g = int(color_hex[2:4], 16)
    b = int(color_hex[4:6], 16)
    return f"rgba({r}, {g}, {b}, 0.3)"

# =========================
# UI
# =========================

st.title("🧪 Simulator Titrasi Interaktif")

st.markdown(
    """
Simulasi kimia interaktif menggunakan Streamlit.  
Fitur: Asam kuat vs basa kuat | Asam lemah vs basa kuat | Kurva pH | Simulasi warna indikator | Ekspor CSV
"""
)

# =========================
# BILAH SAMPING (TATA LETAK BARU)
# =========================

with st.sidebar:
    st.header("⚙️ Pengaturan Titrasi")
    
    st.subheader("Jenis Titrasi")
    type_ = st.selectbox(
        "Pilih reaksi:",
        options=[
            "HCl_NaOH",
            "NaOH_HCl",
            "CH3COOH_NaOH",
        ],
        format_func=lambda x: {
            "HCl_NaOH": "Asam Kuat (HCl) + Basa Kuat (NaOH)",
            "NaOH_HCl": "Basa Kuat (NaOH) + Asam Kuat (HCl)",
            "CH3COOH_NaOH": "Asam Lemah (CH3COOH) + Basa Kuat (NaOH)",
        }.get(x, x),
        help="Pilih pasangan analit dan titran"
    )
    
    st.markdown("---")
    
    st.subheader("🧴 Larutan Analit")
    c0 = st.number_input(
        "Konsentrasi (M)", 
        min_value=0.0, 
        value=0.1, 
        step=0.01, 
        format="%.4f",
        help="Konsentrasi larutan yang akan dititrasi"
    )
    v0 = st.number_input(
        "Volume (L)", 
        min_value=0.001, 
        value=0.05, 
        step=0.01, 
        format="%.4f",
        help="Volume larutan analit (minimal 0.001 L)"
    )
    
    st.subheader("💧 Larutan Titran")
    c_add = st.number_input(
        "Konsentrasi (M)", 
        min_value=0.0, 
        value=0.1, 
        step=0.01, 
        format="%.4f",
        help="Konsentrasi larutan peniter"
    )
    v_max = st.slider(
        "Volume Maksimum (mL)", 
        min_value=10, 
        max_value=100, 
        value=50,
        help="Batas atas volume titran yang disimulasikan"
    )
    
    st.subheader("🌡️ Parameter Tambahan")
    temp_c = st.slider(
        "Suhu (°C)", 
        min_value=20.0, 
        max_value=30.0, 
        value=25.0,
        step=0.5,
        help="Mempengaruhi konstanta ionisasi air (Kw)"
    )
    v_add_ml = st.slider(
        "Volume Ditambahkan (mL)", 
        min_value=0, 
        max_value=int(v_max), 
        value=0,
        step=1,
        help="Volume titran yang sudah ditambahkan (simulasi interaktif)"
    )
    
    if type_ == "CH3COOH_NaOH":
        pKa = st.number_input(
            "pKa Asam Lemah", 
            min_value=0.0, 
            value=4.76, 
            step=0.1, 
            format="%.2f",
            help="Nilai pKa asam asetat adalah 4,76 pada 25°C"
        )
    else:
        pKa = 4.76
    
    st.markdown("---")
    
    st.subheader("🎨 Indikator pH")
    indicator = st.selectbox(
        "Pilih indikator untuk simulasi warna:",
        options=[
            "Phenolphthalein",
            "Methyl Orange",
            "Bromothymol Blue",
        ],
        format_func=lambda x: {
            "Phenolphthalein": "Phenolphthalein (trayek pH 8.2-10)",
            "Methyl Orange": "Methyl Orange (trayek pH 3.1-4.4)",
            "Bromothymol Blue": "Bromothymol Blue (trayek pH 6.0-7.6)",
        }.get(x, x),
        help="Warna larutan akan berubah sesuai pH dan indikator yang dipilih"
    )
    
    st.markdown("---")
    
    st.subheader("📊 Informasi Cepat")
    if c_add > 0:
        ve_calc = (c0 * v0 / c_add) * 1000
        st.metric("Volume Ekuivalen (teoritis)", f"{ve_calc:.2f} mL")
    else:
        st.warning("⚠️ Konsentrasi titran nol, tidak dapat menghitung volume ekuivalen.")
    
    if type_ == "HCl_NaOH":
        rec_ind = "Bromothymol Blue atau Phenolphthalein"
    elif type_ == "CH3COOH_NaOH":
        rec_ind = "Phenolphthalein"
    else:
        rec_ind = "Methyl Orange"
    st.info(f"💡 Indikator yang disarankan: {rec_ind}")
    
    st.caption("🔄 Geser slider dan amati perubahan pH serta warna larutan.")

# =========================
# PARAMETER OBJEK
# =========================

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
# HITUNG pH SAAT INI
# =========================

pH, status, Kw = hitung_ph(params)
Ve = (c0 * v0 / c_add) * 1000 if c_add > 0 else 0

# =========================
# WARNA INDIKATOR (MENGGUNAKAN FUNGSI YANG DIPERBAIKI)
# =========================
solution_color = get_indicator_color(pH, indicator)

# =========================
# TATA LETAK KOLOM UTAMA
# =========================

left, right = st.columns([1, 2])

with left:
    st.subheader("🧪 Larutan")
    
    # Menampilkan larutan dengan warna yang sesuai
    st.markdown(
        f"""
        <div class="solution-container" style="background: {solution_color};">
            <div class="solution-label">
                pH: {pH:.2f} | {indicator}
            </div>
        </div>
        """,
        unsafe_allow_html=True,
    )
    
    col1, col2 = st.columns(2)
    with col1:
        st.metric("📊 pH", f"{pH:.2f}")
        st.metric("⚖️ Volume Ekuivalen", f"{Ve:.2f} mL")
        st.metric("🌡️ Suhu", f"{temp_c:.1f} °C")
    with col2:
        st.metric("📝 Status", status)
        st.metric("💧 Volume Ditambahkan", f"{v_add_ml:.1f} mL")
        st.metric("🔬 Kw", f"{Kw:.2e}")
    
    # Menampilkan keterangan warna indikator
    with st.expander("🎨 Keterangan Warna Indikator"):
        if indicator == "Phenolphthalein":
            st.markdown("""
            - **pH < 8.2**: Tidak berwarna (bening)
            - **pH 8.2 - 10**: Pink (transisi)
            - **pH > 10**: Pink terang
            """)
        elif indicator == "Methyl Orange":
            st.markdown("""
            - **pH < 3.1**: Merah
            - **pH 3.1 - 4.4**: Oranye (transisi)
            - **pH > 4.4**: Kuning
            """)
        elif indicator == "Bromothymol Blue":
            st.markdown("""
            - **pH < 6.0**: Kuning
            - **pH 6.0 - 7.6**: Hijau (transisi)
            - **pH > 7.6**: Biru
            """)

# =========================
# KURVA TITRASI
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
fig.add_trace(go.Scatter(x=vs, y=phs, mode="lines", line=dict(width=4, color="#00ff99"), name="Kurva pH"))

if c_add > 0 and 0 <= Ve <= v_max:
    fig.add_vline(x=Ve, line_dash="dash", line_color="red", annotation_text="Titik Ekuivalen", annotation_position="top")

fig.add_trace(go.Scatter(x=[v_add_ml], y=[pH], mode="markers", marker=dict(size=15, color="red", symbol="circle"), name="Titik Saat Ini"))

# Menambahkan area warna indikator pada kurva
if indicator == "Phenolphthalein":
    fig.add_hrect(y0=8.2, y1=10, line_width=0, fillcolor="pink", opacity=0.2, annotation_text="Range Phenolphthalein")
elif indicator == "Methyl Orange":
    fig.add_hrect(y0=3.1, y1=4.4, line_width=0, fillcolor="orange", opacity=0.2, annotation_text="Range Methyl Orange")
elif indicator == "Bromothymol Blue":
    fig.add_hrect(y0=6.0, y1=7.6, line_width=0, fillcolor="green", opacity=0.2, annotation_text="Range Bromothymol Blue")

fig.update_layout(
    template="plotly_white",
    height=500,
    xaxis_title="Volume Ditambahkan (mL)",
    yaxis_title="pH",
    yaxis=dict(range=[0, 14], gridcolor='lightgray'),
    xaxis=dict(gridcolor='lightgray'),
    hovermode='closest',
)

with right:
    st.subheader("📈 Kurva Titrasi")
    st.plotly_chart(fig, use_container_width=True)

# =========================
# REAKSI DAN PERHITUNGAN
# =========================

st.subheader("⚗️ Reaksi Kimia")
if type_ == "HCl_NaOH":
    st.latex(r"HCl + NaOH \rightarrow NaCl + H_2O")
    st.markdown("**Asam Kuat + Basa Kuat → Garam + Air**")
elif type_ == "NaOH_HCl":
    st.latex(r"NaOH + HCl \rightarrow NaCl + H_2O")
    st.markdown("**Basa Kuat + Asam Kuat → Garam + Air**")
else:
    st.latex(r"CH_3COOH + NaOH \rightarrow CH_3COONa + H_2O")
    st.markdown("**Asam Lemah + Basa Kuat → Garam + Air**")
    st.latex(r"pH = pK_a + \log\frac{[A^-]}{[HA]} \quad \text{(Persamaan Henderson-Hasselbalch)}")

st.subheader("📝 Perhitungan")
with st.expander("🔍 Lihat Langkah Perhitungan"):
    mol_awal = c0 * v0
    mol_titran = c_add * (v_add_ml / 1000)
    st.write(f"**Mol analit awal** = {c0:.4f} M × {v0:.4f} L = {mol_awal:.6f} mol")
    st.write(f"**Mol titran** = {c_add:.4f} M × {v_add_ml/1000:.4f} L = {mol_titran:.6f} mol")
    
    if type_ in ["HCl_NaOH", "CH3COOH_NaOH"]:
        if mol_awal > mol_titran + 1e-12:
            sisa = mol_awal - mol_titran
            st.write(f"✅ **Sisa analit** = {sisa:.6f} mol")
        elif mol_titran > mol_awal + 1e-12:
            kelebihan = mol_titran - mol_awal
            st.write(f"⚠️ **Kelebihan titran** = {kelebihan:.6f} mol")
        else:
            st.success("🎯 **Titik ekuivalen tercapai!**")
    else:
        if mol_awal > mol_titran + 1e-12:
            sisa = mol_awal - mol_titran
            st.write(f"✅ **Sisa basa** = {sisa:.6f} mol")
        elif mol_titran > mol_awal + 1e-12:
            kelebihan = mol_titran - mol_awal
            st.write(f"⚠️ **Kelebihan asam** = {kelebihan:.6f} mol")
        else:
            st.success("🎯 **Titik ekuivalen tercapai!**")
    
    st.write(f"**pH terukur** = {pH:.4f}")
    st.write(f"**Status** = {status}")

# =========================
# EKSPOR DATA
# =========================

st.subheader("💾 Ekspor Data")
export_df = pd.DataFrame({"Volume_mL": vs, "pH": phs})
with st.expander("📊 Tabel Data Kurva"):
    st.dataframe(export_df.round(3), use_container_width=True, height=300)
csv = export_df.to_csv(index=False)
st.download_button(
    label="📥 Unduh CSV", 
    data=csv, 
    file_name=f"titration_curve_{type_}.csv", 
    mime="text/csv"
)

# =========================
# TITRASI OTOMATIS
# =========================

st.subheader("🤖 Titrasi Otomatis")
col1, col2 = st.columns([1, 3])
with col1:
    auto_speed = st.select_slider(
        "Kecepatan simulasi",
        options=["Lambat", "Normal", "Cepat"],
        value="Normal"
    )
    speed_map = {"Lambat": 0.15, "Normal": 0.08, "Cepat": 0.03}
    delay = speed_map[auto_speed]
    
if st.button("▶️ Mulai Simulasi Otomatis", type="primary"):
    progress = st.progress(0)
    status_text = st.empty()
    chart_placeholder = st.empty()
    
    base_fig = go.Figure()
    base_fig.add_trace(go.Scatter(x=vs, y=phs, mode="lines", name="Kurva pH", line=dict(color="#00ff99", width=3)))
    base_fig.update_layout(
        template="plotly_white",
        height=450,
        yaxis=dict(range=[0, 14]),
        xaxis_title="Volume Ditambahkan (mL)",
        yaxis_title="pH",
        showlegend=True
    )
    
    volumes = np.linspace(0, v_max, 50)
    for i, vol in enumerate(volumes):
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
        temp_ph, temp_status, _ = hitung_ph(tmp_params)
        
        temp_fig = go.Figure(base_fig)
        temp_fig.add_trace(
            go.Scatter(
                x=[vol], y=[temp_ph], 
                mode="markers", 
                marker=dict(size=12, color="red", symbol="circle", line=dict(color="darkred", width=2)), 
                name="Titik saat ini",
                showlegend=True if i == 0 else False
            )
        )
        chart_placeholder.plotly_chart(temp_fig, use_container_width=True)
        
        progress.progress((i + 1) / len(volumes))
        status_text.info(f"Volume: {vol:.1f} mL | pH: {temp_ph:.2f} | {temp_status}")
        time.sleep(delay)
    
    status_text.success("✅ Simulasi selesai!")
    time.sleep(1)
    status_text.empty()

# =========================
# RINGKASAN & FOOTER
# =========================

st.divider()
st.subheader("📋 Ringkasan")
ringkasan = pd.DataFrame(
    {
        "Parameter": ["pH", "Status", "Volume Ekuivalen (mL)", "Indikator", "Suhu (°C)", "Kw"],
        "Nilai": [round(pH, 2), status, f"{Ve:.2f}", rec_ind, f"{temp_c:.1f}", f"{Kw:.2e}"],
    }
)
st.table(ringkasan)

st.markdown("---")
st.caption("🧪 **Courtesy Of Kelompok 3 LPK** | Simulator Titrasi Interaktif v2.0")

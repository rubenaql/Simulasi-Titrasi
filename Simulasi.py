# Advanced Streamlit Titration Simulator (Ready Deploy)
import math
import time
from dataclasses import dataclass

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st


st.set_page_config(
    page_title="Simulasi Titrasi Dan Kadar pH",
    layout="wide",
    page_icon="🧪",
)


# =========================
# CHEMISTRY FUNCTIONS
# =========================


def calc_kw(temp_c: float) -> float:
    return 10 ** (-14 + 0.031 * (temp_c - 25))


def pH_from_H(H: float) -> float:
    return -math.log10(H)


def compute_strong_acid_strong_base(type_, c0, v0, c_add, v_add_ml, Kw):
    v_add = v_add_ml / 1000
    Vt = v0 + v_add

    if type_ == "strongA_strongB":
        nH = c0 * v0
        nOH = c_add * v_add
        excess = nOH - nH

        if abs(excess) < 1e-14:
            return 7.0, "Titik ekuivalen"

        if excess > 0:
            OH = excess / Vt
            H = Kw / OH
            return pH_from_H(H), "Kelebihan basa"
        else:
            H = (-excess) / Vt
            return pH_from_H(H), "Kelebihan asam"

    if type_ == "strongB_strongA":
        nOH = c0 * v0
        nH = c_add * v_add
        excess = nH - nOH

        if abs(excess) < 1e-14:
            return 7.0, "Titik ekuivalen"

        if excess > 0:
            H = excess / Vt
            return pH_from_H(H), "Kelebihan asam"
        else:
            OH = (-excess) / Vt
            H = Kw / OH
            return pH_from_H(H), "Kelebihan basa"

    raise ValueError(f"Unknown type_: {type_}")


def solve_weak_acid(c0, v0, c_add, v_add_ml, pKa, Kw):
    Ka = 10 ** (-pKa)

    v_add = v_add_ml / 1000
    Vt = v0 + v_add

    nHA = c0 * v0
    nOH = c_add * v_add

    if nOH < nHA:
        remaining_HA = nHA - nOH
        formed_A = nOH

        if formed_A <= 0:
            H = math.sqrt(Ka * c0)
            return pH_from_H(H)

        pH = pKa + math.log10(formed_A / remaining_HA)
        return pH

    elif abs(nOH - nHA) < 1e-12:
        Csalt = nHA / Vt
        Kb = Kw / Ka
        OH = math.sqrt(Kb * Csalt)
        H = Kw / OH
        return pH_from_H(H)

    else:
        excess = nOH - nHA
        OH = excess / Vt
        H = Kw / OH
        return pH_from_H(H)


@dataclass
class Params:
    type_: str
    temp_c: float
    c0: float
    v0: float
    c_add: float
    v_add_ml: float
    v_max: float
    pKa: float


def compute_ph(params):
    Kw = calc_kw(params.temp_c)

    if params.type_ in ["strongA_strongB", "strongB_strongA"]:
        pH, status = compute_strong_acid_strong_base(
            params.type_,
            params.c0,
            params.v0,
            params.c_add,
            params.v_add_ml,
            Kw,
        )
    else:
        pH = solve_weak_acid(
            params.c0,
            params.v0,
            params.c_add,
            params.v_add_ml,
            params.pKa,
            Kw,
        )

        nHA = params.c0 * params.v0
        nOH = params.c_add * (params.v_add_ml / 1000)

        if abs(nOH - nHA) < 1e-3:
            status = "Titik ekuivalen"
        elif nOH < nHA:
            status = "Daerah buffer"
        else:
            status = "Kelebihan basa"

    return max(0, min(14, pH)), status, Kw


# =========================
# UI
# =========================

st.title("🧪 Simulasi Titrasi Dan Grafik Kadar pH")

st.markdown(
    """
Alat Bantu Titrasi Interaktif.

Features:
- Asam Kuat vs Basa Kuat
- Asam Lemah vs Basa Kuat
- Auto Titrasi
- Kurva pH
- Simulasi Warna Indikator
- Data CSV Export
"""
)


# =========================
# SIDEBAR
# =========================

with st.sidebar:
    st.header("⚙️ Settings")

    type_ = st.selectbox(
        "Titration Type",
        [
            "strongA_strongB",
            "strongB_strongA",
            "weakA_strongB",
        ],
    )

    temp_c = st.slider("Temperature (°C)", 20.0, 30.0, 25.0)

    c0 = st.number_input("Initial Concentration (M)", value=0.1)
    v0 = st.number_input("Initial Volume (L)", value=0.05)

    c_add = st.number_input("Titrant Concentration (M)", value=0.1)

    v_max = st.slider("Maximum Volume (mL)", 10, 100, 50)

    v_add_ml = st.slider("Current Volume (mL)", 0, int(v_max), 0)

    pKa = 4.76
    if type_ == "weakA_strongB":
        pKa = st.number_input("pKa", value=4.76)

    indicator = st.selectbox(
        "Indicator",
        [
            "Phenolphthalein",
            "Methyl Orange",
            "Bromothymol Blue",
        ],
    )


params = Params(
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
# COMPUTE
# =========================

pH, status, Kw = compute_ph(params)


# =========================
# INDICATOR COLOR
# =========================

solution_color = "#00ff99"

if indicator == "Phenolphthalein":
    solution_color = "#ffffff" if pH < 8.2 else "#ff69b4"
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
# LAYOUT
# =========================

left, right = st.columns([1, 2])

with left:
    st.subheader("🧪 Solution")

    st.markdown(
        f"""
        <div style="
            width:220px;
            height:320px;
            border:4px solid white;
            border-radius:20px;
            background:{solution_color};
            margin:auto;
        "></div>
        """,
        unsafe_allow_html=True,
    )

    st.metric("pH", f"{pH:.2f}")
    st.metric("Status", status)
    st.write(f"Kw ≈ {Kw:.2e}")


# =========================
# CURVE
# =========================

vs = np.linspace(0, v_max, 250)
phs = []

for v in vs:
    temp_params = Params(
        type_=type_,
        temp_c=temp_c,
        c0=c0,
        v0=v0,
        c_add=c_add,
        v_add_ml=float(v),
        v_max=v_max,
        pKa=pKa,
    )
    curve_ph, _, _ = compute_ph(temp_params)
    phs.append(curve_ph)

fig = go.Figure()

fig.add_trace(
    go.Scatter(
        x=vs,
        y=phs,
        mode="lines",
        line=dict(width=4, color="#00ff99"),
        name="pH Curve",
    )
)

fig.add_trace(
    go.Scatter(
        x=[v_add_ml],
        y=[pH],
        mode="markers",
        marker=dict(size=12, color="#4aa3ff"),
        name="Current Point",
    )
)

fig.update_layout(
    template="plotly_dark",
    height=500,
    xaxis_title="Volume Added (mL)",
    yaxis_title="pH",
    yaxis=dict(range=[0, 14]),
)

with right:
    st.subheader("📈 Titration Curve")
    st.plotly_chart(fig, use_container_width=True)


# =========================
# REACTION
# =========================

st.subheader("⚗️ Chemical Reaction")

if type_ == "strongA_strongB":
    st.latex(r"HCl + NaOH \rightarrow NaCl + H_2O")
elif type_ == "strongB_strongA":
    st.latex(r"NaOH + HCl \rightarrow NaCl + H_2O")
else:
    st.latex(r"CH_3COOH + NaOH \rightarrow CH_3COONa + H_2O")
    st.latex(r"pH = pK_a + \log\frac{[A^-]}{[HA]}")


# =========================
# DATA EXPORT
# =========================

st.subheader("📥 Export Data")

export_df = pd.DataFrame({"Volume_mL": vs, "pH": phs})
csv = export_df.to_csv(index=False)

st.download_button(
    label="Download CSV",
    data=csv,
    file_name="titration_curve.csv",
    mime="text/csv",
)


# =========================
# AUTO TITRATION
# =========================

st.subheader("▶️ Auto Titration")

if st.button("Start Auto Simulation"):
    progress = st.progress(0)
    chart_placeholder = st.empty()

    for i, vol in enumerate(np.linspace(0, v_max, 40)):
        temp_params = Params(
            type_=type_,
            temp_c=temp_c,
            c0=c0,
            v0=v0,
            c_add=c_add,
            v_add_ml=float(vol),
            v_max=v_max,
            pKa=pKa,
        )
        temp_ph, _, _ = compute_ph(temp_params)

        temp_fig = go.Figure()
        temp_fig.add_trace(go.Scatter(x=vs, y=phs, mode="lines"))
        temp_fig.add_trace(
            go.Scatter(
                x=[vol],
                y=[temp_ph],
                mode="markers",
                marker=dict(size=12),
            )
        )

        temp_fig.update_layout(
            template="plotly_dark",
            height=450,
            yaxis=dict(range=[0, 14]),
        )

        chart_placeholder.plotly_chart(temp_fig, use_container_width=True)
        progress.progress((i + 1) / 40)
        time.sleep(0.08)


# =========================
# FOOTER
# =========================

st.divider()

st.caption("Credits To Kelompok 2,Streamlit, Plotly, NumPy, and Pandas.")

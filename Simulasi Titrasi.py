import math
from dataclasses import dataclass

import numpy as np
import plotly.graph_objects as go
import streamlit as st


def calc_kw(temp_c: float) -> float:
    """Pendekatan edukasi (20-30°C).

    log10(Kw) = -14 + 0.031*(tempC-25)
    sehingga Kw(25°C) ~ 1e-14.
    """
    return 10 ** (-14 + 0.031 * (temp_c - 25))


def pH_from_H(H: float) -> float:
    return -math.log10(H)


def compute_strong_acid_strong_base(
    type_: str,
    c0: float,
    v0: float,
    c_add: float,
    v_add_ml: float,
    Kw: float,
):
    """Titrasi asam kuat/basa kuat dengan koreksi Kw.

    type_:
      - strongA_strongB: asam kuat (awal) + basa kuat (ditambahkan)
      - strongB_strongA: basa kuat (awal) + asam kuat (ditambahkan)

    Catatan: diasumsikan ion kuat terdisosiasi sempurna.
    """
    v_add = v_add_ml / 1000.0
    Vt = v0 + v_add

    if type_ == "strongA_strongB":
        nH = c0 * v0
        nOH = c_add * v_add
        excess = nOH - nH  # + => OH berlebih

        if abs(excess) < 1e-14:
            return {"pH": 7.0, "status": "Titik ekuivalen"}

        if excess > 0:
            OH = excess / Vt
            H = Kw / OH
            return {"pH": pH_from_H(H), "status": "Kelebihan basa (OH⁻)"}
        else:
            H = (-excess) / Vt
            return {"pH": pH_from_H(H), "status": "Kelebihan asam (H⁺)"}

    if type_ == "strongB_strongA":
        nOH = c0 * v0
        nH = c_add * v_add
        excess = nH - nOH  # + => H berlebih

        if abs(excess) < 1e-14:
            return {"pH": 7.0, "status": "Titik ekuivalen"}

        if excess > 0:
            H = excess / Vt
            return {"pH": pH_from_H(H), "status": "Kelebihan asam (H⁺)"}
        else:
            OH = (-excess) / Vt
            H = Kw / OH
            return {"pH": pH_from_H(H), "status": "Kelebihan basa (OH⁻)"}

    raise ValueError(f"Unknown type_: {type_}")


def solve_ph_weak_acid_strong_base(
    c0: float,
    v0: float,
    c_add: float,
    v_add_ml: float,
    pKa: float,
    Kw: float,
):
    """Titrasi asam lemah (HA) + basa kuat (OH⁻) pada suhu dengan Kw.

    Kita pakai model distribusi asam lemah HA <-> H+ + A- dengan
    koefisien fraksi:
      HA = Ct * H / (H + Ka)
      A- = Ct * Ka / (H + Ka)

    Stoikiometri awal mempertimbangkan reaksi OH dengan HA:
      jika OH ditambahkan < OH yang dibutuhkan untuk menghabiskan HA,
      maka ada sisa HA.

    Persamaan keseimbangan muatan:
      H + Na = OH + A-
    dengan Na berasal dari basa kuat (diasumsikan lengkap terionisasi).
    """
    v_add = v_add_ml / 1000.0
    Vt = v0 + v_add
    Ka = 10 ** (-pKa)

    nHA0 = c0 * v0
    nOH_added = c_add * v_add

    # Berdasarkan stoikiometri netralisasi OH dengan HA:
    nA = min(nOH_added, nHA0)  # A- terbentuk
    nHA = max(0.0, nHA0 - nOH_added)  # sisa HA

    Ct = (nHA + nA) / Vt
    Na = nOH_added / Vt

    def alpha(H: float):
        denom = H + Ka
        # denom > 0 selalu
        return {"HA": Ct * H / denom, "A": Ct * Ka / denom}

    def f(pH: float) -> float:
        H = 10 ** (-pH)
        OH = Kw / H
        A_minus = alpha(H)["A"]
        return H + Na - OH - A_minus

    lo, hi = 0.01, 14.0
    flo, fhi = f(lo), f(hi)

    # Kalau langsung tepat
    if flo == 0:
        return lo
    if fhi == 0:
        return hi

    # Bracket scan jika perlu
    if flo * fhi > 0:
        prev_p, prev_f = lo, flo
        found = False
        for p in np.arange(0.05, 13.95 + 1e-12, 0.2):
            F = f(float(p))
            if prev_f * F <= 0:
                lo, hi = prev_p, float(p)
                flo, fhi = prev_f, F
                found = True
                break
            prev_p, prev_f = float(p), F

        if not found:
            # fallback: ambil pH yang meminimalkan |f|
            best_p = lo
            best = abs(flo)
            for p in np.arange(0.01, 14.01 + 1e-12, 0.05):
                F = abs(f(float(p)))
                if F < best:
                    best = F
                    best_p = float(p)
            return best_p

    # Bisection
    for _ in range(90):
        mid = (lo + hi) / 2
        fm = f(mid)
        if flo * fm <= 0:
            hi = mid
            fhi = fm
        else:
            lo = mid
            flo = fm

    return (lo + hi) / 2


@dataclass
class TitrationParams:
    type_: str
    temp_c: float
    Kw_override: float | None
    c0: float
    v0: float
    c_add: float
    v_add_max_ml: float
    v_add_ml: float
    pKa: float


def resolve_kw(temp_c: float, Kw_override: float | None) -> float:
    if Kw_override is not None and Kw_override > 0:
        return Kw_override
    return calc_kw(temp_c)


def compute_ph(params: TitrationParams):
    Kw = resolve_kw(params.temp_c, params.Kw_override)

    if params.type_ in ("strongA_strongB", "strongB_strongA"):
        res = compute_strong_acid_strong_base(
            params.type_,
            params.c0,
            params.v0,
            params.c_add,
            params.v_add_ml,
            Kw,
        )
        pH = float(res["pH"])
        status = str(res["status"])
        return pH, status, Kw

    # weakA_strongB
    pH = float(
        solve_ph_weak_acid_strong_base(
            params.c0,
            params.v0,
            params.c_add,
            params.v_add_ml,
            params.pKa,
            Kw,
        )
    )

    # Status (indikatif)
    nHA0 = params.c0 * params.v0
    nOH_added = params.c_add * (params.v_add_ml / 1000.0)
    eq_ratio = abs(nOH_added - nHA0) / max(1e-12, nHA0)

    if eq_ratio < 1e-3:
        status = "Mendekati titik ekuivalen"
    elif nOH_added < nHA0:
        status = "Daerah penyangga (HA/A⁻)"
    else:
        status = "Kelebihan basa (OH⁻)"

    return pH, status, Kw


def plot_curve(params: TitrationParams, steps: int = 80):
    vs = np.linspace(0, params.v_add_max_ml, steps + 1)
    phs = []
    for v in vs:
        p_tmp = TitrationParams(
            type_=params.type_,
            temp_c=params.temp_c,
            Kw_override=params.Kw_override,
            c0=params.c0,
            v0=params.v0,
            c_add=params.c_add,
            v_add_max_ml=params.v_add_max_ml,
            v_add_ml=float(v),
            pKa=params.pKa,
        )
        ph, _, _ = compute_ph(p_tmp)
        phs.append(float(np.clip(ph, 0, 14)))

    # Marker pada volume slider
    pH_current, _, _ = compute_ph(params)
    v_current = params.v_add_ml

    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=vs,
            y=phs,
            mode="lines",
            name="Kurva pH",
            line=dict(color="#2ee59d", width=3),
        )
    )
    fig.add_trace(
        go.Scatter(
            x=[v_current],
            y=[min(14, max(0, pH_current))],
            mode="markers",
            name="Titik saat ini",
            marker=dict(size=10, color="#4aa3ff"),
        )
    )

    fig.update_layout(
        template="plotly_dark",
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        margin=dict(l=50, r=20, t=30, b=40),
        xaxis_title="Volume penitrasi ditambahkan (mL)",
        yaxis_title="pH",
        yaxis=dict(range=[0, 14]),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="left", x=0),
        height=480,
    )
    return fig


st.set_page_config(page_title="Simulasi Titrasi (pH)", layout="wide")

st.title("Simulasi Titrasi dan Kurva pH")

st.caption(
    "Model edukasi (asam kuat/basa kuat) dan numerik (asam lemah + basa kuat) memakai pendekatan Kw terhadap suhu." \
    " Tidak mempertimbangkan aktivitas ion, efek ionik, atau faktor kompleks lainnya." \
)

with st.sidebar:
    st.header("Pengaturan")

    type_ = st.selectbox(
        "Jenis titrasi",
        options=[
            ("strongA_strongB", "Asam kuat + Basa kuat"),
            ("strongB_strongA", "Basa kuat + Asam kuat"),
            ("weakA_strongB", "Asam lemah (HA) + Basa kuat"),
        ],
        format_func=lambda x: x[1],
    )[0]

    temp_c = st.number_input("Suhu (°C)", min_value=20.0, max_value=30.0, value=25.0, step=0.1)

    Kw_override = st.number_input(
        "Konstanta air (Kw) override (opsional)",
        min_value=1e-18,
        value=float("nan"),
        step=1e-16,
        help="Kosongkan agar Kw dihitung otomatis dari suhu.",
    )
    Kw_override_val = None if (isinstance(Kw_override, float) and math.isnan(Kw_override)) else float(Kw_override)

    st.subheader("Konsentrasi & Volume")
    c0 = st.number_input("C awal (mol/L)", min_value=0.0, value=0.1, step=0.001)
    v0 = st.number_input("V awal (L)", min_value=0.0001, value=0.050, step=0.001)

    c_add = st.number_input("C penitrasi (mol/L)", min_value=0.0, value=0.1, step=0.001)
    v_add_max_ml = st.number_input("V maksimum (mL)", min_value=1.0, value=50.0, step=1.0)

    # Kurangi slider: pakai input angka untuk volume saat ini
    v_add_ml = st.number_input(
        "V saat ini (mL)",
        min_value=0.0,
        max_value=float(v_add_max_ml),
        value=0.0,
        step=1.0,
    )


    pKa = None
    if type_ == "weakA_strongB":
        pKa = st.number_input("pKa asam lemah (HA)", min_value=-10.0, max_value=20.0, value=4.76, step=0.01)

    st.divider()
    do_compute = st.button("Hitung pH", type="primary")
    do_plot = st.button("Plot kurva pH")


params = TitrationParams(
    type_=type_,
    temp_c=float(temp_c),
    Kw_override=Kw_override_val,
    c0=float(c0),
    v0=float(v0),
    c_add=float(c_add),
    v_add_max_ml=float(v_add_max_ml),
    v_add_ml=float(v_add_ml),
    pKa=float(pKa) if pKa is not None else 0.0,
)

# Validasi dasar
if params.v0 <= 0:
    st.error("Volume awal v0 harus > 0")
    st.stop()
if not (20.0 <= params.temp_c <= 30.0):
    st.error("Suhu harus di antara 20 sampai 30 °C")
    st.stop()
if params.v_add_ml < 0 or params.v_add_ml > params.v_add_max_ml:
    st.error("Volume penambahan vAdd di luar rentang")
    st.stop()

# Hitung pH
col1, col2 = st.columns([1, 1])

if do_compute or do_plot:
    with st.spinner("Menghitung pH..."):
        pH, status, Kw = compute_ph(params)

    col1.metric("pH (saat ini)", f"{pH:.2f}")
    col2.metric("Status", status)

    st.caption(f"Suhu: {params.temp_c:.1f} °C | Kw ≈ {Kw:.2e}")

# Plot kurva
if do_plot:
    with st.spinner("Mempersiapkan kurva..."):
        fig = plot_curve(params, steps=80)
    st.plotly_chart(fig, use_container_width=True)

# Default display (tanpa tombol plot)
if not do_plot:
    # Tampilkan kurva kecil saat belum plot (opsional: tampilkan tombol saja)
    st.info("Klik **Plot kurva pH** untuk melihat grafik pH vs volume.")


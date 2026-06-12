def get_indicator_color(pH, indicator):
    if indicator == "Phenolphthalein":
        if pH < 8.2:
            # tidak berwarna (kita gunakan putih transparan agar lebih natural)
            return "rgba(255, 255, 255, 0.1)"
        elif pH < 10.0:
            ratio = (pH - 8.2) / 1.8
            r = 255
            g = int(255 - 150 * ratio)   # 255 → 105
            b = int(255 - 75 * ratio)    # 255 → 180
            return f"#{r:02x}{g:02x}{b:02x}"
        else:
            return "#ff69b4"   # merah muda penuh

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
            return "#ffff00"            # kuning
        elif pH < 7.6:
            ratio = (pH - 6.0) / 1.6
            if ratio < 0.5:
                # kuning → hijau (R menurun, B tetap 0)
                r = int(255 * (1 - 2 * ratio))
                g = 255
                b = 0
            else:
                # hijau → biru (G menurun, B meningkat)
                r = 0
                g = int(255 * (2 - 2 * ratio))
                b = int(255 * (2 * ratio - 1))
            return f"#{r:02x}{g:02x}{b:02x}"
        else:
            return "#0000ff"            # biru

    return "#ffffff"

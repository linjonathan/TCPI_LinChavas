#!/usr/bin/env python3
"""
Local backend for the TC Potential Intensity tool.

Serves index.html and exposes /api/figure, which renders a publication-quality
SounderPy sounding with the potential-intensity overlays (SST parcel + markers +
PI annotation), computed with the real tcpyPI package.

Run:
    pip install -r requirements.txt
    python server.py
    # then open  http://localhost:5000  in your browser

The live, instant-update sounding in the page is pure JavaScript and needs no
server.  This backend is only used by the "Generate SounderPy figure" button,
which produces the high-resolution matplotlib figure for saving/publication.
"""
import io
import warnings
warnings.filterwarnings("ignore")

import matplotlib
matplotlib.use("Agg")  # headless
import matplotlib.pyplot as plt
import numpy as np

from flask import Flask, request, Response, send_from_directory
import os

# scientific stack
import sounderpy as spy
from tcpyPI import pi
from metpy.calc import (mixing_ratio_from_relative_humidity,
                        relative_humidity_from_dewpoint, moist_lapse,
                        saturation_mixing_ratio, virtual_temperature)
from metpy.units import units

HERE = os.path.dirname(os.path.abspath(__file__))
app = Flask(__name__, static_folder=None)

RD, G = 287.04, 9.80665


def hypsometric_heights(p_hpa, T_c, Td_c):
    """Geopotential heights (m) from the profile via the hypsometric equation."""
    p = np.asarray(p_hpa, float)
    Tk = np.asarray(T_c, float) + 273.15
    # mixing ratio (kg/kg) from dewpoint
    es = 6.112 * np.exp(17.67 * np.asarray(Td_c, float) / (243.5 + np.asarray(Td_c, float)))
    w = 0.622 * es / (p - es)
    Tv = Tk * (1 + 0.61 * w)
    z = np.zeros_like(p)
    z[0] = 10.0
    for i in range(1, len(p)):
        Tvm = 0.5 * (Tv[i] + Tv[i - 1])
        z[i] = z[i - 1] + RD * Tvm / G * np.log(p[i - 1] / p[i])
    return z


def synth_winds(z):
    """Gentle illustrative wind profile (winds are not part of PI physics)."""
    u = 2.5 + 0.0009 * z
    v = 1.0 + 0.0007 * z
    return u, v


def add_cors(resp):
    resp.headers["Access-Control-Allow-Origin"] = "*"
    resp.headers["Access-Control-Allow-Headers"] = "Content-Type"
    resp.headers["Access-Control-Allow-Methods"] = "POST, OPTIONS"
    return resp


@app.route("/")
def index():
    return send_from_directory(HERE, "index.html")


@app.route("/pi.js")
def pijs():
    return send_from_directory(HERE, "pi.js")


@app.route("/api/figure", methods=["POST", "OPTIONS"])
def figure():
    if request.method == "OPTIONS":
        return add_cors(Response(""))
    data = request.get_json(force=True)
    prof = data["prof"]                       # [{p,T,Td}, ...] surface-first
    SST = float(data.get("sst", 29.0))
    MSLP = float(data.get("mslp", 1010.0))
    diss = int(data.get("diss", 1))
    ckcd = float(data.get("ckcd", 0.9))
    ascent = int(data.get("ascent", 0))

    p = np.array([r["p"] for r in prof], float)
    T = np.array([r["T"] for r in prof], float)
    Td = np.array([r["Td"] for r in prof], float)
    z = hypsometric_heights(p, T, Td)
    u, v = synth_winds(z)

    # PI with the real package
    rh = relative_humidity_from_dewpoint(T * units.degC, Td * units.degC)
    mr = mixing_ratio_from_relative_humidity(p * units.hPa, T * units.degC, rh).magnitude * 1000
    VMAX, PMIN, IFL, TO, LNB = pi(SST, MSLP, p, T, mr,
                                  CKCD=ckcd, diss_flag=diss, ascent_flag=ascent)

    clean = {
        "p": p * units.hPa, "z": z * units.meter,
        "T": T * units.degC, "Td": Td * units.degC,
        "u": u * units("m/s"), "v": v * units("m/s"),
        "site_info": {
            "site-name": "TROPICAL", "site-lctn": "INTERACTIVE",
            "site-latlon": [15.0, -60.0], "site-elv": 10,
            "source": "USER INPUT DATA", "model": "no-model",
            "fcst-hour": "no-fcst-hour",
            "run-time": ["2026", "01", "01", "00"],
            "valid-time": ["2026", "01", "01", "00"],
        },
        "titles": {
            "top_title": "TROPICAL SOUNDING — POTENTIAL INTENSITY",
            "left_title": f"SST {SST:.1f}°C  |  MSLP {MSLP:.0f} hPa",
            "right_title": "Interactive PI tool",
        },
    }

    plt.close("all")
    spy.build_sounding(clean, style="full", color_blind=False,
                       dark_mode=False, special_parcels="simple")
    fig = plt.gcf()
    ax = fig.get_axes()[0]

    # ---- PI overlays (mirrors the notebook) ----
    valid = (not np.isnan(VMAX)) and IFL == 1
    plev = np.array([1000, 975, 950, 925, 900, 875, 850, 800, 750, 700, 650,
                     600, 550, 500, 450, 400, 350, 300, 250, 200, 150, 100]) * units.hPa
    sstT = moist_lapse(plev, SST * units.degC, reference_pressure=1000 * units.hPa)
    sstmr = saturation_mixing_ratio(plev, sstT)
    sstTv = virtual_temperature(sstT, sstmr)
    ax.plot(sstTv.magnitude, plev.magnitude, color="magenta", lw=3, ls="--",
            zorder=9, alpha=.9, label="SST parcel (Tv)")
    ax.plot(SST, 1000, marker="o", ms=9, color="red", mec="darkred", mew=1.5, zorder=10)
    ax.plot(sstTv[0].magnitude, 1000, marker="s", ms=9, color="magenta",
            mec="purple", mew=1.5, zorder=10)

    if valid:
        ax.plot(TO - 273.15, LNB, marker="o", ms=9, color="red",
                mec="darkred", mew=1.5, zorder=10)
        kt = float(VMAX) * 1.94384
        cat = ("Tropical Storm" if kt < 64 else "Category 1" if kt < 83 else
               "Category 2" if kt < 96 else "Category 3" if kt < 113 else
               "Category 4" if kt < 137 else "Category 5")
        ann = (f"TC Potential Intensity\nVMAX: {float(VMAX):.1f} m/s ({kt:.0f} kt)\n"
               f"PMIN: {float(PMIN):.1f} hPa\nCategory: {cat}\nSST: {SST:.1f}°C\n"
               f"(T0, LNB): ({TO-273.15:.0f}°C, {LNB:.0f} hPa)")
    else:
        ann = f"TC Potential Intensity\nVMAX = 0 (no storm supported)\nSST: {SST:.1f}°C"
    ax.text(0.97, 0.97, ann, transform=ax.transAxes, fontsize=16, va="top", ha="right",
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=.85, edgecolor="black", lw=2))

    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=int(data.get("dpi", 110)), bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return add_cors(Response(buf.read(), mimetype="image/png"))


if __name__ == "__main__":
    print("\nTC Potential Intensity — SounderPy backend")
    print("Open http://localhost:5000 in your browser.\n")
    # threaded=False keeps matplotlib's global state safe across requests
    app.run(host="127.0.0.1", port=5000, threaded=False, debug=False)

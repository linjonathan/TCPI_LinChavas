# Tropical Cyclone Potential Intensity — interactive sounding

An interactive, public-friendly tool for computing and explaining tropical-cyclone
**potential intensity (PI)** directly on an atmospheric sounding.

## Two ways to use it

### 1. Live tool (no install)
Just open **`index.html`** in any browser. The sounding and every PI quantity on it
update instantly as you move the sliders or edit the sounding table. All the physics
runs locally in JavaScript — a faithful port of Kerry Emanuel's *pcmin* algorithm
(the `tcpyPI` package), validated against the Python original to within 0.01 m/s / 0.01 hPa.

This file is fully self-contained and can be hosted anywhere (GitHub Pages, a shared
drive, etc.).

### 2. Publication-quality figure (optional backend)
The **"Generate SounderPy figure"** button renders the full, beautiful
[SounderPy](https://kylejgillett.github.io/sounderpy/) sounding with the PI overlays
(SST parcel, SST/outflow markers, PI annotation) using the *real* `tcpyPI` and
`SounderPy` packages. Because those are Python/matplotlib, this needs a small local server:

```bash
pip install -r requirements.txt
python server.py
# open http://localhost:5000 in your browser, then click "Generate SounderPy figure"
```

The button posts the current sounding/SST/options to the backend, which returns a
high-resolution PNG you can view and download. (You can also open `index.html` directly
as a file and point it at a separately-running `server.py`; the backend sends permissive
CORS headers so that works too.)

> Note: SounderPy needs a wind profile for its hodograph; the thermodynamics-only tool
> doesn't track winds, so the backend synthesizes a gentle illustrative wind profile.
> Winds play no role in the PI calculation.

## Files
| file | purpose |
|------|---------|
| `index.html` | the complete live tool (self-contained; embeds the PI engine) |
| `pi.js` | the standalone JavaScript PI engine (same code, reusable) |
| `server.py` | optional Flask backend for the SounderPy figure export |
| `requirements.txt` | Python deps for the backend only |

## References
Bister & Emanuel (2002, *JGR*); Emanuel (1986, *JAS*; 1994, *Atmospheric Convection*);
Gilford (2021, *GMD*) — pyPI/tcpyPI; Wing, Emanuel & Solomon (2015, *GRL*);
SounderPy — Kyle J. Gillett.

# TCPI Interactive Sounding — Development Notes

A summary of the bug fixes and visual/aesthetic changes made to the interactive
potential-intensity tool (`web/index.html`, `web/pi.js`) since the initial version.
The physics engine is a faithful JavaScript port of Emanuel's *pcmin* algorithm
(`tcpyPI`), validated to reproduce the Python reference to within 0.01 m/s / 0.01 hPa.

## Bugs fixed

- **Skew-T curves ran off the plot.** The original skew factor was too strong, so the
  temperature and parcel curves exited the right edge above ~550 hPa. Reduced the skew and
  widened the temperature window so the full profile stays on-chart.
- **Grid lines invisible on white.** Background isotherms/adiabats were too faint to read on
  the white panel; darkened them.
- **Canvas distortion.** Canvases used `width:100%` with a fixed height attribute, stretching
  the plot. Rewrote the setup to derive the pixel height from the rendered width, preserving
  aspect ratio and staying crisp on hi-DPI displays.
- **Info-box text overlap.** The Saffir–Simpson category line overlapped the T₀ line; rebuilt
  the box so each line gets its own row and the box auto-sizes.
- **Outflow marker cut off at the top.** With the chart top at 100 hPa the outflow level
  (~92 hPa) sat on the border; extended the pressure axis to 50 hPa so the marker and label
  are fully visible.
- **Shaded CAPE region disappeared.** The purple saturation-CAPE fill used a
  "contiguous-from-the-surface" method that broke whenever the ocean parcel had a little CIN
  near the surface (e.g., the Weisman–Klemp supercell at lower SST) — the entire fill above
  the first non-buoyant level vanished. Replaced it with a robust layer-by-layer
  positive-area fill.
- **Idealized sounding gave zero PI (should be nonzero).** The first "idealized moist adiabat"
  preset built the environment from the *saturated* 30 °C adiabat, which made the saturated
  SST parcel neutral and collapsed PI to ~0. Rebuilt it to lift the actual 80%-RH surface
  parcel (dry below its LCL, reversible-moist above), giving the intended result: ~zero CAPE
  but a real PI from the air–sea disequilibrium.
- **Pressure-dependence "off" baseline crash.** Evaluating the parcels at the true MSL
  (1010 hPa) for the toggle-off case launched them below the sounding's 1000-hPa surface and
  broke the reversible entropy inversion (CAPEMS → NaN at certain SSTs). Reverted to the safe
  surface-capped baseline.
- **Redundant condensate control.** A standalone "condensate loading" toggle turned out to be
  mathematically identical to the existing pseudo-adiabatic ascent option (in pcmin the ascent
  flag only toggles condensate loading, not the parcel temperature). Removed it and relabeled
  the ascent control instead.

## Aesthetic & visualization changes

- **Skew-T became the centerpiece** — enlarged, with the numeric results and the sensitivity
  plot moved below it.
- **Matched SounderPy styling** — reproduced MetPy's 47° skew transform, pressure range
  1050 → 50 hPa, surface-temperature-dependent temperature bounds, cornflower-blue dry/moist
  adiabats, black dashed mixing-ratio lines, blue 0 °C / −20 °C highlight lines, bold
  pressure/temperature tick labels, and km height labels up the left axis.
- **Curve colors matched SounderPy** — red temperature, green dewpoint, light-blue wet-bulb,
  dark-red dotted virtual/density temperature.
- **Inline curve labels replaced the legend** — each curve is labeled directly on the plot
  (T, Td, Tw, Tv, ocean parcel, near-surface parcel, and the "residual ∝ V²" energy region);
  the separate bottom legend was removed.
- **Precise parcel notation** — parcel curves use density-temperature subscript notation
  (T_ρ,ocean and T_ρ,parcel) drawn with proper lowered subscripts.
- **Two-region energy shading** — the purple wedge (ocean-parcel saturation CAPE, CAPEMS) with
  the subtracted near-surface-parcel CAPE (CAPEb) overlaid in red, so the residual purple is
  the fuel ∝ V².
- **Diagnostics strip and info box** — added under-plot stats (ocean-parcel CAPE, near-surface
  CAPE, LCL, PWAT, outflow T₀) and a ΔP = MSL − Pmin line in the on-plot PI box.
- **Text register adjusted** — the deep explainer was rewritten for a masters-level
  atmospheric-science audience, while the top intro was softened to an accessible
  "speed limit / maximum near-surface wind" framing (citing Rousseau-Rizzi & Emanuel 2019).
- **Marker & label polish** — repositioned the outflow, SST, and region labels to avoid
  overlaps; the SST marker and parcel-launch point now ride up to the R_max pressure when the
  pressure-dependence toggle is on.

## Note

Functional features added over the same period (not covered above) include the preset
soundings (Idealized moist adiabat, Jordan mean 1958, Weisman–Klemp supercell, plus the
tropical presets), the R_max pressure-dependence toggle, the virtual-temperature toggle, and
the optional SounderPy figure-export backend (`server.py`).

// pi.js — faithful JavaScript port of tcpyPI (Bister-Emanuel pcmin algorithm)
// Adapted from pcmin.m by Kerry Emanuel; pyPI/tcpyPI by Daniel Gilford.
// Ported to JS for the TC Potential Intensity web visualizer.
// Reference: Bister & Emanuel (2002); Emanuel (1994); Gilford (2021, JOSS).

const C = {
  CPD: 1005.7,      // J/kg/K specific heat dry air (const p)
  CPV: 1870.0,      // J/kg/K specific heat water vapor (const p)
  CL: 2500.0,       // J/kg/K modified specific heat liquid water
  RV: 461.5,        // J/kg/K gas constant water vapor
  RD: 287.04,       // J/kg/K gas constant dry air
  ALV0: 2.501e6,    // J/kg latent heat of vaporization at 0C
  A: 1669.0,        // pLCL empirical param
  B: 122.0,         // pLCL empirical param
  b: 2.0,           // azimuthal velocity exponent
};
C.CPVMCL = C.CPV - C.CL;     // -630
C.EPS = C.RD / C.RV;         // ~0.6220

// ---- thermodynamic helpers ----
const T_ktoC = (Tk) => Tk - 273.15;
const T_Ctok = (TC) => TC + 273.15;
const es_cc  = (TC) => 6.112 * Math.exp(17.67 * TC / (243.5 + TC));   // hPa
const Lv     = (TC) => C.ALV0 + C.CPVMCL * TC;                         // J/kg
const ev     = (R, P) => R * P / (C.EPS + R);                          // hPa
const rv     = (E, P) => C.EPS * E / (P - E);                          // g/g
// density temperature. TRHO_V toggles the vapor (virtual-temperature) buoyancy term; condensate
// loading (1+R)/(1+RT) is always applied, its magnitude set by RT (reversible vs pseudo-adiabatic).
let TRHO_V = true;
const Trho   = (T, RT, R) => T * (TRHO_V ? (1 + R / C.EPS) / (1 + R) : 1) * (1 + R) / (1 + RT); // K
const e_pLCL = (TP, RH, PP) => PP * Math.pow(RH, TP / (C.A - C.B * RH - TP));

function entropy_S(T, R, P) {
  const EV = ev(R, P);
  const ES = es_cc(T - 273.15);
  const RH = Math.min(EV / ES, 1.0);
  const ALV = Lv(T - 273.15);
  return (C.CPD + R * C.CL) * Math.log(T) - C.RD * Math.log(P - EV)
         + ALV * R / T - R * C.RV * Math.log(RH);
}

// Newton-Raphson invert saturated entropy for temperature
function solveTfromEntropy(S, P, RP, T_initial) {
  let TGNEW = T_initial;
  let TJC = T_ktoC(T_initial);
  let ES = es_cc(TJC);
  let RG = rv(ES, P);
  let NC = 0;
  let TG = 0;

  while (Math.abs(TGNEW - TG) > 0.001) {
    TG = TGNEW;
    const TC = T_ktoC(TG);
    const ENEW = es_cc(TC);
    RG = rv(ENEW, P);
    NC += 1;
    const ALV = Lv(TC);
    const SL = (C.CPD + RP * C.CL + ALV * ALV * RG / (C.RV * TG * TG)) / TG;
    const EM = ev(RG, P);
    const SG = (C.CPD + RP * C.CL) * Math.log(TG) - C.RD * Math.log(P - EM) + ALV * RG / TG;
    const AP = (NC < 3) ? 0.3 : 1.0;
    TGNEW = TG + AP * (S - SG) / SL;
    if (NC > 500 || ENEW > (P - 1)) {
      return { TG, RG, IFLAG: 2 };  // did not converge
    }
  }
  return { TG, RG, IFLAG: 1 };
}

function argminAbsToPtop(P, ptop) {
  let idx = 0, best = Infinity;
  for (let i = 0; i < P.length; i++) {
    const d = Math.abs(P[i] - ptop);
    if (d < best) { best = d; idx = i; }
  }
  return idx;
}

// CAPE of a parcel given parcel (TP K, RP g/g, PP hPa) and environment profiles
function cape(TP, RP, PP, T_in, R_in, P_in, ascent_flag = 0, ptop = 50, miss_handle = 1) {
  // missing-value handling
  const validIdx = T_in.map((t) => !Number.isNaN(t));
  let firstValid = validIdx.findIndex((v) => v);
  if (firstValid < 0) firstValid = 0;
  const numValid = validIdx.reduce((a, v) => a + (v ? 1 : 0), 0);
  let first_lvl;
  if (numValid !== P_in.length) {
    if (miss_handle !== 0) {
      return { CAPED: NaN, TOB: NaN, LNB: NaN, IFLAG: 3 };
    } else {
      let anyNanAbove = false;
      for (let i = firstValid; i < P_in.length; i++) if (Number.isNaN(T_in[i])) anyNanAbove = true;
      if (anyNanAbove) return { CAPED: NaN, TOB: NaN, LNB: NaN, IFLAG: 3 };
      first_lvl = firstValid;
    }
  } else {
    first_lvl = 0;
  }

  const N = argminAbsToPtop(P_in, ptop);
  const P = P_in.slice(first_lvl, N);
  const T = T_in.slice(first_lvl, N);
  const R = R_in.slice(first_lvl, N);
  const nlvl = P.length;
  const TVRDIF = new Array(nlvl).fill(0);

  // checks
  if (P[2] - P[1] > 0) return { CAPED: 0, TOB: NaN, LNB: NaN, IFLAG: 0 };
  if (RP < 1e-6 || TP < 200) return { CAPED: 0, TOB: NaN, LNB: NaN, IFLAG: 0 };

  const TPC = T_ktoC(TP);
  const ESP = es_cc(TPC);
  const EVP = ev(RP, PP);
  let RH = EVP / ESP;
  RH = Math.min(RH, 1.0);
  const S = entropy_S(TP, RP, PP);
  const PLCL = e_pLCL(TP, RH, PP);

  let CAPED = 0;
  let TOB = T[0];
  let IFLAG = 1;
  let jmin = 1e6;

  for (let j = 0; j < nlvl; j++) {
    jmin = Math.min(jmin, j);
    if (P[j] >= PLCL) {
      const TG = TP * Math.pow(P[j] / PP, C.RD / C.CPD);
      const RG = RP;
      const TLVR = Trho(TG, RG, RG);
      const TVENV = Trho(T[j], R[j], R[j]);
      TVRDIF[j] = TLVR - TVENV;
    } else {
      const res = solveTfromEntropy(S, P[j], RP, T[j]);
      if (res.IFLAG === 2) {
        return { CAPED: 0, TOB: T[0], LNB: P[0], IFLAG: 2 };
      }
      const TG = res.TG, RG = res.RG;
      const RMEAN = ascent_flag * RG + (1 - ascent_flag) * RP;
      const TLVR = Trho(TG, RMEAN, RG);
      const TENV = Trho(T[j], R[j], R[j]);
      TVRDIF[j] = TLVR - TENV;
    }
  }

  let NA = 0.0, PA = 0.0;
  let INB = 0;
  for (let j = nlvl - 1; j > jmin; j--) {
    if (TVRDIF[j] > 0) INB = Math.max(INB, j);
  }

  if (INB === 0) {
    return { CAPED: 0, TOB: T[0], LNB: 0, IFLAG };
  }

  for (let j = jmin + 1; j < INB + 1; j++) {
    const PFAC = C.RD * (TVRDIF[j] + TVRDIF[j - 1]) * (P[j - 1] - P[j]) / (P[j] + P[j - 1]);
    PA += Math.max(PFAC, 0.0);
    NA -= Math.min(PFAC, 0.0);
  }

  const PMA = PP + P[jmin];
  const PFAC = C.RD * (PP - P[jmin]) / PMA;
  PA += PFAC * Math.max(TVRDIF[jmin], 0.0);
  NA -= PFAC * Math.min(TVRDIF[jmin], 0.0);

  let PAT = 0.0;
  TOB = T[INB];
  let LNB = P[INB];
  if (INB < nlvl - 1) {
    const PINB = (P[INB + 1] * TVRDIF[INB] - P[INB] * TVRDIF[INB + 1]) / (TVRDIF[INB] - TVRDIF[INB + 1]);
    LNB = PINB;
    PAT = C.RD * TVRDIF[INB] * (P[INB] - PINB) / (P[INB] + PINB);
    TOB = (T[INB] * (PINB - P[INB + 1]) + T[INB + 1] * (P[INB] - PINB)) / (P[INB] - P[INB + 1]);
  }

  CAPED = PA + PAT - NA;
  CAPED = Math.max(CAPED, 0.0);
  return { CAPED, TOB, LNB, IFLAG: 1 };
}

// Maximum potential intensity. SSTC (C), MSL (hPa), P (hPa[]), TC (C[]), Rgkg (g/kg[])
function pi(SSTC, MSL, P, TC, Rgkg, opts = {}) {
  const CKCD = opts.CKCD ?? 0.9;
  const ascent_flag = opts.ascent_flag ?? 0;
  const diss_flag = opts.diss_flag ?? 1;
  const V_reduc = opts.V_reduc ?? 0.8;
  const ptop = opts.ptop ?? 50;
  const miss_handle = opts.miss_handle ?? 1;
  const pdep = opts.pressDep ?? true;   // parcels at R_max pressure (true) or background MSL (false)
  TRHO_V = opts.useVapor ?? true;       // toggle the vapor (virtual-temperature) buoyancy term

  const SSTK = T_Ctok(SSTC);
  const T = TC.map(T_Ctok);
  const R = Rgkg.map((r) => (Number.isNaN(r) ? 0 : r * 0.001));

  if (SSTC <= 5.0 || SSTC > 100) return { VMAX: NaN, PMIN: NaN, IFL: 0, TO: NaN, OTL: NaN };
  if (Number.isNaN(SSTC))         return { VMAX: NaN, PMIN: NaN, IFL: 0, TO: NaN, OTL: NaN };
  const minT = Math.min(...T), maxTC = Math.max(...TC);
  if (minT <= 100 || maxTC > 100) return { VMAX: NaN, PMIN: NaN, IFL: 0, TO: NaN, OTL: NaN };

  const ES0 = es_cc(SSTC);
  const NK = 0;

  // environmental CAPE
  let res = cape(T[NK], R[NK], P[NK], T, R, P, ascent_flag, ptop, miss_handle);
  const CAPEA = res.CAPED;
  let IFL = 1;
  if (res.IFLAG !== 1) IFL = res.IFLAG;

  let NP = 0;
  let PM = 970.0;
  let PMOLD = PM;
  let PNEW = 0.0;
  let TO = NaN, OTL = NaN, TVAV = NaN, RAT = NaN, CAPEM = NaN, CAPEMS = NaN;

  while (Math.abs(PNEW - PMOLD) > 0.5) {
    // parcels at the R_max pressure PM (pdep on) or the surface pressure (pdep off)
    const PPeval = pdep ? Math.min(PM, 1000.0) : Math.min(MSL, 1000.0);
    // CAPE at radius of max winds
    let TP = T[NK];
    let PP = PPeval;
    let RP = C.EPS * R[NK] * MSL / (PP * (C.EPS + R[NK]) - R[NK] * MSL);
    res = cape(TP, RP, PP, T, R, P, ascent_flag, ptop, miss_handle);
    CAPEM = res.CAPED;
    if (res.IFLAG !== 1) IFL = res.IFLAG;

    // saturation CAPE at radius of max winds
    TP = SSTK;
    PP = PPeval;
    RP = rv(ES0, PP);
    res = cape(TP, RP, PP, T, R, P, ascent_flag, ptop, miss_handle);
    CAPEMS = res.CAPED;
    if (res.IFLAG !== 1) IFL = res.IFLAG;
    TO = res.TOB;
    OTL = res.LNB;

    RAT = SSTK / TO;
    if (diss_flag === 0) RAT = 1.0;

    const RS0 = RP;
    const TV0 = Trho(T[NK], R[NK], R[NK]);
    const TVSST = Trho(SSTK, RS0, RS0);
    TVAV = 0.5 * (TV0 + TVSST);
    let CAT = (CAPEM - CAPEA) + 0.5 * CKCD * RAT * (CAPEMS - CAPEM);
    CAT = Math.max(CAT, 0.0);
    PNEW = MSL * Math.exp(-CAT / (C.RD * TVAV));

    PMOLD = PM;
    PM = PNEW;
    NP += 1;
    if (NP > 200 || PM < 400) return { VMAX: NaN, PMIN: NaN, IFL: 0, TO: NaN, OTL: NaN };
  }

  const CATFAC = 0.5 * (1 + 1 / C.b);
  let CAT = (CAPEM - CAPEA) + CKCD * RAT * CATFAC * (CAPEMS - CAPEM);
  CAT = Math.max(CAT, 0.0);
  const PMIN = MSL * Math.exp(-CAT / (C.RD * TVAV));
  const FAC = Math.max(0.0, CAPEMS - CAPEM);
  const VMAX = V_reduc * Math.sqrt(CKCD * RAT * FAC);

  return { VMAX, PMIN, IFL, TO, OTL, CAPEMS, CAPEA, PM };  // PM = converged R_max pressure
}

if (typeof module !== 'undefined' && module.exports) {
  module.exports = { pi, cape, es_cc, rv, ev, Trho, entropy_S, C };
}

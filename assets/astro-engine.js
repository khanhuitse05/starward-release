// Starward Clean-Room Ephemeris, Projection & Procedural Sky Model
// Pure JavaScript port of astronomical algorithms matching Starward's core conventions.

(function (root, factory) {
  if (typeof module === 'object' && module.exports) {
    module.exports = factory();
  } else {
    root.StarwardAstro = factory();
  }
}(typeof self !== 'undefined' ? self : this, function () {
  'use strict';

  const DEG = Math.PI / 180;
  const RAD = 180 / Math.PI;
  const TWO_PI = Math.PI * 2;
  const AU_KM = 149597870.7;
  const LIGHT_DAYS_PER_AU = 0.0057755183;
  const K_GAUSS = 0.01720209895;
  const EPS0_RAD = 23.43929111 * DEG;

  function clamp(v, lo, hi) {
    return v < lo ? lo : v > hi ? hi : v;
  }

  function normDeg(d) {
    let r = d % 360;
    return r < 0 ? r + 360 : r;
  }

  function normRad(r) {
    let x = r % TWO_PI;
    return x < 0 ? x + TWO_PI : x;
  }

  function smoothstep(e0, e1, x) {
    const t = clamp((x - e0) / (e1 - e0), 0, 1);
    return t * t * (3 - 2 * t);
  }

  function julianDate(date) {
    const ms = date instanceof Date ? date.getTime() : Number(date);
    const ut = ms / 86400000 + 2440587.5;
    const tt = ut + 69.2 / 86400;
    return {
      ut,
      tt,
      tUT: (ut - 2451545.0) / 36525.0,
      tTT: (tt - 2451545.0) / 36525.0
    };
  }

  function nutation(jd) {
    const T = jd.tTT;
    const om = (125.04452 - 1934.136261 * T + 0.0020708 * T * T) * DEG;
    const ls = (280.4665 + 36000.7698 * T) * DEG;
    const lm = (218.3165 + 481267.8813 * T) * DEG;
    const dpsiArcsec = -17.20 * Math.sin(om) - 1.32 * Math.sin(2 * ls) - 0.23 * Math.sin(2 * lm) + 0.21 * Math.sin(2 * om);
    const depsArcsec = 9.20 * Math.cos(om) + 0.57 * Math.cos(2 * ls) + 0.10 * Math.cos(2 * lm) - 0.09 * Math.cos(2 * om);
    const eps0Arcsec = 84381.448 - 46.8150 * T - 0.00059 * T * T + 0.001813 * T * T * T;
    return {
      dpsiRad: (dpsiArcsec / 3600) * DEG,
      depsRad: (depsArcsec / 3600) * DEG,
      meanEpsRad: (eps0Arcsec / 3600) * DEG,
      trueEpsRad: ((eps0Arcsec + depsArcsec) / 3600) * DEG
    };
  }

  function greenwichApparentSiderealRad(jd, nut) {
    const dt = jd.ut - 2451545.0;
    const T = jd.tUT;
    const gmstDeg = normDeg(
      280.46061837 +
      360.98564736629 * dt +
      0.000387933 * T * T -
      (T * T * T) / 38710000.0
    );
    const eqEqRad = nut.dpsiRad * Math.cos(nut.trueEpsRad);
    return normRad(gmstDeg * DEG + eqEqRad);
  }

  // Lieske (1979) 3x3 J2000 -> Equatorial of Date matrix
  function precessionMatrixFromJ2000(jd) {
    const T = jd.tTT;
    const T2 = T * T;
    const T3 = T2 * T;
    const s2r = DEG / 3600;
    const zeta = (2306.2181 * T + 0.30188 * T2 + 0.017998 * T3) * s2r;
    const z = (2306.2181 * T + 1.09468 * T2 + 0.018203 * T3) * s2r;
    const theta = (2004.3109 * T - 0.42665 * T2 - 0.041833 * T3) * s2r;

    const cz = Math.cos(z), sz = Math.sin(z);
    const ct = Math.cos(theta), st = Math.sin(theta);
    const cze = Math.cos(zeta), sze = Math.sin(zeta);

    return [
      cz * ct * cze - sz * sze, -cz * ct * sze - sz * cze, -cz * st,
      sz * ct * cze + cz * sze, -sz * ct * sze + cz * cze, -sz * st,
      st * cze,                 -st * sze,                  ct
    ];
  }

  function mat3MulVec(m, v) {
    return [
      m[0] * v[0] + m[1] * v[1] + m[2] * v[2],
      m[3] * v[0] + m[4] * v[1] + m[5] * v[2],
      m[6] * v[0] + m[7] * v[1] + m[8] * v[2]
    ];
  }

  function mat3Mul(a, b) {
    return [
      a[0]*b[0] + a[1]*b[3] + a[2]*b[6], a[0]*b[1] + a[1]*b[4] + a[2]*b[7], a[0]*b[2] + a[1]*b[5] + a[2]*b[8],
      a[3]*b[0] + a[4]*b[3] + a[5]*b[6], a[3]*b[1] + a[4]*b[4] + a[5]*b[7], a[3]*b[2] + a[4]*b[5] + a[5]*b[8],
      a[6]*b[0] + a[7]*b[3] + a[8]*b[6], a[6]*b[1] + a[7]*b[4] + a[8]*b[7], a[6]*b[2] + a[7]*b[5] + a[8]*b[8]
    ];
  }

  function eclipticVecToEquatorialJ2000(v) {
    const c = Math.cos(EPS0_RAD), s = Math.sin(EPS0_RAD);
    return [v[0], v[1] * c - v[2] * s, v[1] * s + v[2] * c];
  }

  function eqVectorToRaDec(v) {
    const r = Math.hypot(v[0], v[1], v[2]) || 1;
    const decRad = Math.asin(clamp(v[2] / r, -1, 1));
    const raRad = normRad(Math.atan2(v[1], v[0]));
    return { raRad, decRad, raDeg: raRad * RAD, decDeg: decRad * RAD };
  }

  function raDecDegToVec(raDeg, decDeg) {
    const ra = raDeg * DEG, dec = decDeg * DEG;
    const cd = Math.cos(dec);
    return [cd * Math.cos(ra), cd * Math.sin(ra), Math.sin(dec)];
  }

  // Transform J2000 ecliptic vector to apparent Equatorial (RA, Dec) of date
  function apparentFromEclipticJ2000(geoEcl, jd, precM, nut) {
    const eqJ2000 = eclipticVecToEquatorialJ2000(geoEcl);
    const eqDate = mat3MulVec(precM, eqJ2000);
    const base = eqVectorToRaDec(eqDate);
    const eps = nut.meanEpsRad;
    const ra = base.raRad, dec = base.decRad;
    const cosD = Math.cos(dec);
    const tanD = Math.abs(cosD) > 1e-4 ? Math.sin(dec) / cosD : 0;
    const dra = (Math.cos(eps) + Math.sin(eps) * Math.sin(ra) * tanD) * nut.dpsiRad - Math.cos(ra) * tanD * nut.depsRad;
    const ddec = Math.sin(eps) * Math.cos(ra) * nut.dpsiRad + Math.sin(ra) * nut.depsRad;
    const rRa = normRad(ra + dra);
    const rDec = clamp(dec + ddec, -Math.PI / 2, Math.PI / 2);
    return { raRad: rRa, decRad: rDec, raDeg: rRa * RAD, decDeg: rDec * RAD };
  }

  function createFrameContext(date, latDeg, lonDeg) {
    const jd = julianDate(date);
    const nut = nutation(jd);
    const gast = greenwichApparentSiderealRad(jd, nut);
    const last = normRad(gast + lonDeg * DEG);
    const latRad = latDeg * DEG;
    const precM = precessionMatrixFromJ2000(jd);

    const sLat = Math.sin(latRad), cLat = Math.cos(latRad);
    const sL = Math.sin(last), cL = Math.cos(last);

    // Date Equatorial unit vector -> Horizontal [East, North, Up]
    // H = LAST - RA =>
    // E = -sin(H)*cosD = -(sin(last)*cos(ra) - cos(last)*sin(ra))*cosD = -sin(last)*X + cos(last)*Y
    // N = -sin(lat)*cos(H)*cosD + cos(lat)*sinD
    // Z =  cos(lat)*cos(H)*cosD + sin(lat)*sinD
    const eqDateToHoriz = [
      -sL,         cL,          0,
      -sLat * cL, -sLat * sL,   cLat,
       cLat * cL,  cLat * sL,   sLat
    ];
    const j2000ToHoriz = mat3Mul(eqDateToHoriz, precM);

    return {
      date,
      jd,
      nut,
      gastRad: gast,
      lastRad: last,
      latDeg,
      lonDeg,
      latRad,
      precM,
      eqDateToHoriz,
      j2000ToHoriz
    };
  }

  function refractDeg(trueAltDeg) {
    if (trueAltDeg < -1.5) return 0;
    const h = Math.max(trueAltDeg, -0.5);
    const arg = (h + 10.3 / (h + 5.11)) * DEG;
    const rArcmin = 1.02 / Math.tan(arg);
    const atten = trueAltDeg < 0 ? smoothstep(-1.5, 0, trueAltDeg) : 1;
    return clamp((rArcmin / 60) * atten, 0, 0.7);
  }

  function apparentRaDecToHoriz(raDeg, decDeg, frame, applyRefraction) {
    const vEq = raDecDegToVec(raDeg, decDeg);
    const vH = mat3MulVec(frame.eqDateToHoriz, vEq);
    return horizVecToCoords(vH, applyRefraction);
  }

  function j2000VecToHoriz(vJ2000, frame, applyRefraction) {
    const vH = mat3MulVec(frame.j2000ToHoriz, vJ2000);
    return horizVecToCoords(vH, applyRefraction);
  }

  function horizVecToCoords(vH, applyRefraction) {
    const E = vH[0], N = vH[1], U = vH[2];
    const trueAltDeg = Math.asin(clamp(U, -1, 1)) * RAD;
    const altDeg = applyRefraction !== false ? trueAltDeg + refractDeg(trueAltDeg) : trueAltDeg;
    const azDeg = normDeg(Math.atan2(E, N) * RAD);
    if (applyRefraction !== false && altDeg !== trueAltDeg) {
      const aRad = altDeg * DEG;
      const ca = Math.cos(aRad);
      const azR = azDeg * DEG;
      return {
        vec: [ca * Math.sin(azR), ca * Math.cos(azR), Math.sin(aRad)],
        altDeg,
        trueAltDeg,
        azDeg
      };
    }
    return { vec: vH, altDeg, trueAltDeg, azDeg };
  }

  function horizToApparentRaDec(azDeg, altDeg, frame) {
    const az = azDeg * DEG, alt = altDeg * DEG;
    const E = Math.cos(alt) * Math.sin(az);
    const N = Math.cos(alt) * Math.cos(az);
    const U = Math.sin(alt);
    // Transpose of eqDateToHoriz
    const m = frame.eqDateToHoriz;
    const vx = m[0] * E + m[3] * N + m[6] * U;
    const vy = m[1] * E + m[4] * N + m[7] * U;
    const vz = m[2] * E + m[5] * N + m[8] * U;
    return eqVectorToRaDec([vx, vy, vz]);
  }

  // ---- Stereographic SkyProjection ----
  function createProjection(azDeg, altDeg, fovDeg, width, height) {
    const az = azDeg * DEG;
    const alt = altDeg * DEG;
    const fov = clamp(fovDeg, 0.35, 130) * DEG;

    const ca = Math.cos(alt), sa = Math.sin(alt);
    const caz = Math.cos(az), saz = Math.sin(az);

    const forward = [ca * saz, ca * caz, sa];
    const right = [caz, -saz, 0];
    // up = cross(right, forward)
    const up = [
      right[1] * forward[2] - right[2] * forward[1],
      right[2] * forward[0] - right[0] * forward[2],
      right[0] * forward[1] - right[1] * forward[0]
    ];

    const halfMin = Math.min(width, height) * 0.5;
    const halfMax = Math.max(width, height) * 0.5;
    const ppR = halfMin / (2 * Math.tan(fov * 0.25));
    const cx = width * 0.5;
    const cy = height * 0.5;
    const maxRadius = Math.hypot(width, height) * 0.65;

    function project(v) {
      const x = v[0], y = v[1], z = v[2];
      const c = x * forward[0] + y * forward[1] + z * forward[2];
      if (c <= -0.58) return null;
      const k = 2 / (1 + c);
      const rx = x * right[0] + y * right[1] + z * right[2];
      const uy = x * up[0] + y * up[1] + z * up[2];
      const px = cx + k * rx * ppR;
      const py = cy - k * uy * ppR;
      if (Math.abs(px - cx) > maxRadius && Math.abs(py - cy) > maxRadius) return null;
      return { x: px, y: py, scale: k, c };
    }

    function unproject(px, py) {
      const X = (px - cx) / ppR;
      const Y = (cy - py) / ppR;
      const rho2 = X * X + Y * Y;
      const denom = 4 + rho2;
      const kF = (4 - rho2) / denom;
      const kR = (4 * X) / denom;
      const kU = (4 * Y) / denom;
      const vx = kF * forward[0] + kR * right[0] + kU * up[0];
      const vy = kF * forward[1] + kR * right[1] + kU * up[1];
      const vz = kF * forward[2] + kR * right[2] + kU * up[2];
      const len = Math.hypot(vx, vy, vz) || 1;
      const nx = vx / len, ny = vy / len, nz = vz / len;
      return {
        vec: [nx, ny, nz],
        altDeg: Math.asin(clamp(nz, -1, 1)) * RAD,
        azDeg: normDeg(Math.atan2(nx, ny) * RAD)
      };
    }

    return {
      azDeg,
      altDeg,
      fovDeg,
      width,
      height,
      cx,
      cy,
      ppR,
      halfMin,
      halfMax,
      forward,
      right,
      up,
      project,
      unproject,
      pxPerRad: ppR * 2
    };
  }

  // ---- Solar Position (Meeus Ch. 25) ----
  function computeSun(jd, nut) {
    const T = jd.tTT;
    const L0 = normDeg(280.46646 + 36000.76983 * T + 0.0003032 * T * T);
    const M = normDeg(357.52911 + 35999.05029 * T - 0.0001537 * T * T) * DEG;
    const e = 0.016708634 - 0.000042037 * T - 0.0000001267 * T * T;
    const C =
      (1.914602 - 0.004817 * T - 0.000014 * T * T) * Math.sin(M) +
      (0.019993 - 0.000101 * T) * Math.sin(2 * M) +
      0.000289 * Math.sin(3 * M);
    const trueLonDeg = L0 + C;
    const v = M + C * DEG;
    const distanceAU = (1.000001018 * (1 - e * e)) / (1 + e * Math.cos(v));
    const omega = (125.04 - 1934.136 * T) * DEG;
    const lambdaDeg = trueLonDeg - 0.00569 - 0.00478 * Math.sin(omega);
    const lambda = lambdaDeg * DEG;
    const eps = nut.meanEpsRad + 0.00256 * DEG * Math.cos(omega);

    const sinLam = Math.sin(lambda), cosLam = Math.cos(lambda);
    const raRad = normRad(Math.atan2(Math.cos(eps) * sinLam, cosLam));
    const decRad = Math.asin(clamp(Math.sin(eps) * sinLam, -1, 1));
    return {
      id: 'sun',
      kind: 'sun',
      name: 'Sun',
      raRad,
      decRad,
      raDeg: raRad * RAD,
      decDeg: decRad * RAD,
      distanceAU,
      heliocentricAU: 0,
      phaseAngleDeg: 0,
      illuminatedFraction: 1,
      magnitude: -26.74,
      angularDiameterDeg: 2 * Math.atan(696340 / (distanceAU * AU_KM)) * RAD
    };
  }

  // ---- Lunar Position & Phase (Meeus Ch. 47) ----
  const MOON_LR = [
    [0,0,1,0,6288774,-20905355],[2,0,-1,0,1274027,-3699111],[2,0,0,0,658314,-2955968],
    [0,0,2,0,213618,-569925],[0,1,0,0,-185116,48888],[0,0,0,2,-114332,-3149],
    [2,0,-2,0,58793,246158],[2,-1,-1,0,57066,-152138],[2,0,1,0,53322,-170733],
    [2,-1,0,0,45758,-204586],[0,1,-1,0,-40923,-129620],[1,0,0,0,-34720,108743],
    [0,1,1,0,-30383,104755],[2,0,0,-2,15327,10321],[0,0,1,2,-12528,0],
    [0,0,1,-2,10980,79661],[4,0,-1,0,10675,-34782],[0,0,3,0,10034,-23210],
    [4,0,-2,0,8548,-21636],[2,1,-1,0,-7888,24208],[2,1,0,0,-6766,30824],
    [1,0,-1,0,-5163,-8379],[1,1,0,0,4987,-16675],[2,-1,1,0,4036,-12831],
    [2,0,2,0,3994,-10445],[4,0,0,0,3861,-11650],[2,0,-3,0,3665,14403],
    [0,1,-2,0,-2689,-7003],[2,0,-1,2,-2602,0],[2,-1,-2,0,2390,10056],
    [1,0,1,0,-2348,6322],[2,-2,0,0,2236,-9884],[0,1,2,0,-2120,5751],
    [0,2,0,0,-2069,0],[2,-2,-1,0,2048,-4950],[2,0,1,-2,-1773,4130],
    [2,0,0,2,-1595,0],[4,-1,-1,0,1215,-3958],[0,0,2,2,-1110,0],[3,0,-1,0,-892,3258]
  ];
  const MOON_B = [
    [0,0,0,1,5128122],[0,0,1,1,280602],[0,0,1,-1,277693],[2,0,0,-1,173237],
    [2,0,-1,1,55413],[2,0,-1,-1,46271],[2,0,0,1,32573],[0,0,2,1,17198],
    [2,0,1,-1,9266],[0,0,2,-1,8822],[2,-1,0,-1,8216],[2,0,-2,-1,4324],
    [2,0,1,1,4200],[2,1,0,-1,-3359],[2,-1,-1,1,2463],[2,-1,0,1,2211],
    [2,-1,-1,-1,2065],[0,1,-1,-1,-1870],[4,0,-1,-1,1828],[0,1,0,1,-1794]
  ];

  function computeMoon(jd, nut, sun) {
    const T = jd.tTT;
    const T2 = T * T, T3 = T2 * T, T4 = T3 * T;
    const Lp = normDeg(218.3164477 + 481267.88123421 * T - 0.0015786 * T2 + T3 / 538841 - T4 / 65194000);
    const D = normDeg(297.8501921 + 445267.1114034 * T - 0.0018819 * T2 + T3 / 545868 - T4 / 113065000);
    const M = normDeg(357.5291092 + 35999.0502909 * T - 0.0001536 * T2 + T3 / 24490000);
    const Mp = normDeg(134.9633964 + 477198.8675055 * T + 0.0087414 * T2 + T3 / 69699 - T4 / 14712000);
    const F = normDeg(93.2720950 + 483202.0175233 * T - 0.0036539 * T2 - T3 / 3526000 + T4 / 863310000);

    const A1 = normDeg(119.75 + 131.849 * T) * DEG;
    const A2 = normDeg(53.09 + 479264.290 * T) * DEG;
    const A3 = normDeg(313.45 + 481266.484 * T) * DEG;
    const E = 1 - 0.002516 * T - 0.0000074 * T2;
    const E2 = E * E;

    const Dr = D * DEG, Mr = M * DEG, Mpr = Mp * DEG, Fr = F * DEG;
    let sumL = 0, sumR = 0, sumB = 0;

    for (let i = 0; i < MOON_LR.length; i++) {
      const [cD, cM, cMp, cF, sL, sR] = MOON_LR[i];
      const arg = cD * Dr + cM * Mr + cMp * Mpr + cF * Fr;
      const eScale = Math.abs(cM) === 1 ? E : Math.abs(cM) === 2 ? E2 : 1;
      sumL += sL * eScale * Math.sin(arg);
      sumR += sR * eScale * Math.cos(arg);
    }
    for (let i = 0; i < MOON_B.length; i++) {
      const [cD, cM, cMp, cF, sB] = MOON_B[i];
      const arg = cD * Dr + cM * Mr + cMp * Mpr + cF * Fr;
      const eScale = Math.abs(cM) === 1 ? E : Math.abs(cM) === 2 ? E2 : 1;
      sumB += sB * eScale * Math.sin(arg);
    }

    sumL += 3958 * Math.sin(A1) + 1962 * Math.sin(Lp * DEG - Fr) + 318 * Math.sin(A2);
    sumB += -2235 * Math.sin(Lp * DEG) + 382 * Math.sin(A3) + 175 * Math.sin(A1 - Fr) + 175 * Math.sin(A1 + Fr);

    const lambda = (Lp + sumL / 1e6) * DEG + nut.dpsiRad;
    const beta = (sumB / 1e6) * DEG;
    const distanceKm = 385000.56 + sumR / 1000;
    const distanceAU = distanceKm / AU_KM;

    const eps = nut.trueEpsRad;
    const cb = Math.cos(beta), sb = Math.sin(beta);
    const cl = Math.cos(lambda), sl = Math.sin(lambda);
    const ce = Math.cos(eps), se = Math.sin(eps);

    const raRad = normRad(Math.atan2(sl * ce - Math.tan(beta) * se, cl));
    const decRad = Math.asin(clamp(sb * ce + cb * se * sl, -1, 1));

    // Phase angle (Meeus 48.2 / 48.4)
    const cosPsi =
      Math.sin(sun.decRad) * Math.sin(decRad) +
      Math.cos(sun.decRad) * Math.cos(decRad) * Math.cos(sun.raRad - raRad);
    const psi = Math.acos(clamp(cosPsi, -1, 1));
    const phaseAngleRad = Math.atan2(sun.distanceAU * Math.sin(psi), distanceAU - sun.distanceAU * Math.cos(psi));
    const phaseAngleDeg = phaseAngleRad * RAD;
    const illuminatedFraction = (1 + Math.cos(phaseAngleRad)) * 0.5;
    const isWaxing = Math.sin(Dr) >= 0;

    // Position angle of Moon's bright limb (Meeus 48.5)
    const chiRad = Math.atan2(
      Math.cos(sun.decRad) * Math.sin(sun.raRad - raRad),
      Math.sin(sun.decRad) * Math.cos(decRad) - Math.cos(sun.decRad) * Math.sin(decRad) * Math.cos(sun.raRad - raRad)
    );

    const mag = -12.73 + 0.026 * Math.abs(phaseAngleDeg) + 4e-9 * Math.pow(phaseAngleDeg, 4);
    const angularDiameterDeg = 2 * Math.atan(1737.4 / distanceKm) * RAD;

    let phaseName = 'Crescent';
    if (illuminatedFraction < 0.03) phaseName = 'New Moon';
    else if (illuminatedFraction < 0.45) phaseName = isWaxing ? 'Waxing Crescent' : 'Waning Crescent';
    else if (illuminatedFraction < 0.55) phaseName = isWaxing ? 'First Quarter' : 'Last Quarter';
    else if (illuminatedFraction < 0.96) phaseName = isWaxing ? 'Waxing Gibbous' : 'Waning Gibbous';
    else phaseName = 'Full Moon';

    return {
      id: 'moon',
      kind: 'moon',
      name: 'Moon',
      raRad,
      decRad,
      raDeg: raRad * RAD,
      decDeg: decRad * RAD,
      distanceKm,
      distanceAU,
      heliocentricAU: 1,
      phaseAngleDeg,
      illuminatedFraction,
      isWaxing,
      brightLimbChiRad: chiRad,
      phaseName,
      magnitude: mag,
      angularDiameterDeg
    };
  }

  // ---- Standish J2000 (1800–2050 AD) Keplerian Planet Elements ----
  const PLANET_ROWS = {
    mercury: { name: 'Mercury', radiusKm: 2439.7, a: 0.38709927, da: 0.00000037, e: 0.20563593, de: 0.00001906, I: 7.00497902, dI: -0.00594749, L: 252.25032350, dL: 149472.67411175, wb: 77.45779628, dwb: 0.16047689, Om: 48.33076593, dOm: -0.12534081 },
    venus:   { name: 'Venus',   radiusKm: 6051.8, a: 0.72333566, da: 0.00000390, e: 0.00677672, de: -0.00004107, I: 3.39467605, dI: -0.00078890, L: 181.97909950, dL: 58517.81538729, wb: 131.60246718, dwb: 0.00268329, Om: 76.67984255, dOm: -0.27769418 },
    emb:     { name: 'Earth',   radiusKm: 6371.0, a: 1.00000261, da: 0.00000562, e: 0.01671123, de: -0.00004392, I: -0.00001531, dI: -0.01294668, L: 100.46457166, dL: 35999.37244981, wb: 102.93768193, dwb: 0.32327364, Om: 0.0, dOm: 0.0 },
    mars:    { name: 'Mars',    radiusKm: 3389.5, a: 1.52371034, da: 0.00001847, e: 0.09339410, de: 0.00007882, I: 1.84969142, dI: -0.00813131, L: -4.55343205, dL: 19140.30268499, wb: -23.94362959, dwb: 0.44441088, Om: 49.55953891, dOm: -0.29257343 },
    jupiter: { name: 'Jupiter', radiusKm: 69911,  a: 5.20288700, da: -0.00011607, e: 0.04838624, de: -0.00013253, I: 1.30439695, dI: -0.00183714, L: 34.39644051, dL: 3034.74612775, wb: 14.72847983, dwb: 0.21252668, Om: 100.47390909, dOm: 0.20469106, b: -0.00012452, c: 0.06064060, s: -0.35635438, f: 38.35125000 },
    saturn:  { name: 'Saturn',  radiusKm: 58232,  a: 9.53667594, da: -0.00125060, e: 0.05386179, de: -0.00050991, I: 2.48599187, dI: 0.00193609, L: 49.95424423, dL: 1222.49362201, wb: 92.59887831, dwb: -0.41897216, Om: 113.66242448, dOm: -0.28867794, b: 0.00025899, c: -0.13434469, s: 0.87320147, f: 38.35125000 },
    uranus:  { name: 'Uranus',  radiusKm: 25362,  a: 19.18916464, da: -0.00196176, e: 0.04725744, de: -0.00004397, I: 0.77263783, dI: -0.00242939, L: 313.23810451, dL: 428.48202785, wb: 170.95427630, dwb: 0.40805281, Om: 74.01692503, dOm: 0.04240589, b: 0.00058331, c: -0.97731848, s: 0.17689245, f: 7.67025000 },
    neptune: { name: 'Neptune', radiusKm: 24622,  a: 30.06992276, da: 0.00026291, e: 0.00859048, de: 0.00005105, I: 1.77004347, dI: 0.00035372, L: -55.12002969, dL: 218.45945325, wb: 44.96476227, dwb: -0.32241464, Om: 131.78422574, dOm: -0.00508664, b: -0.00041348, c: 0.68346318, s: -0.10162547, f: 7.67025000 }
  };

  function solveKeplerElliptical(Mrad, e) {
    let E = Mrad + e * Math.sin(Mrad);
    for (let i = 0; i < 10; i++) {
      const f = E - e * Math.sin(E) - Mrad;
      const fp = 1 - e * Math.cos(E);
      const dE = f / fp;
      E -= dE;
      if (Math.abs(dE) < 1e-11) break;
    }
    return E;
  }

  function planetHeliocentricEclJ2000(row, jd) {
    const T = jd.tTT;
    const a = row.a + row.da * T;
    const e = row.e + row.de * T;
    const I = (row.I + row.dI * T) * DEG;
    const L = row.L + row.dL * T;
    const wb = row.wb + row.dwb * T;
    const Om = (row.Om + row.dOm * T) * DEG;
    const w = (wb - (row.Om + row.dOm * T)) * DEG;

    let Mdeg = normDeg(L - wb);
    if (row.b) {
      const fT = row.f * T * DEG;
      Mdeg += row.b * T * T + row.c * Math.cos(fT) + row.s * Math.sin(fT);
    }
    const Mrad = normRad(Mdeg * DEG);
    const E = solveKeplerElliptical(Mrad, e);

    const xOrb = a * (Math.cos(E) - e);
    const yOrb = a * Math.sqrt(1 - e * e) * Math.sin(E);

    const cw = Math.cos(w), sw = Math.sin(w);
    const cO = Math.cos(Om), sO = Math.sin(Om);
    const cI = Math.cos(I), sI = Math.sin(I);

    return [
      (cw * cO - sw * sO * cI) * xOrb + (-sw * cO - cw * sO * cI) * yOrb,
      (cw * sO + sw * cO * cI) * xOrb + (-sw * sO + cw * cO * cI) * yOrb,
      (sw * sI) * xOrb + (cw * sI) * yOrb
    ];
  }

  function planetMagnitude(id, r, delta, iDeg) {
    const base = 5 * Math.log10(r * delta);
    const i = iDeg, i2 = i * i, i3 = i2 * i;
    switch (id) {
      case 'mercury': return -0.42 + base + 0.0380 * i - 0.000273 * i2 + 0.000002 * i3;
      case 'venus': return -4.40 + base + 0.0009 * i + 0.000239 * i2 - 0.00000065 * i3;
      case 'mars': return -1.52 + base + 0.016 * i;
      case 'jupiter': return -9.40 + base + 0.005 * i;
      case 'saturn': return -8.88 + base + 0.044 * i;
      case 'uranus': return -7.19 + base;
      case 'neptune': return -6.87 + base;
      default: return 0;
    }
  }

  // Approximate Galilean satellites in Jovian equatorial plane for zoomed-in Jupiter view
  const GALILEAN = [
    { name: 'Io', aRJ: 5.905, periodD: 1.7691378, phase0Deg: 106.1, mag: 5.02 },
    { name: 'Europa', aRJ: 9.397, periodD: 3.5511810, phase0Deg: 175.7, mag: 5.29 },
    { name: 'Ganymede', aRJ: 14.992, periodD: 7.1545530, phase0Deg: 120.6, mag: 4.61 },
    { name: 'Callisto', aRJ: 26.330, periodD: 16.6890184, phase0Deg: 84.4, mag: 5.65 }
  ];

  function computeGalileanMoons(jd, jupApparent, jupDistAU) {
    const tDays = jd.tt - 2451545.0;
    const jupRadiusRad = Math.atan(69911 / (jupDistAU * AU_KM));
    const ra0 = jupApparent.raRad, dec0 = jupApparent.decRad;
    // Jupiter pole ~ RA 268.05 deg, Dec +64.49 deg -> P ~ -25..+25 deg tilt along ecliptic
    const pa = -22 * DEG;
    const cp = Math.cos(pa), sp = Math.sin(pa);
    const invCosDec = 1 / Math.max(Math.cos(dec0), 0.1);
    return GALILEAN.map(g => {
      const u = (g.phase0Deg + (360 / g.periodD) * tDays) * DEG;
      const xR = g.aRJ * Math.sin(u);
      const yR = -g.aRJ * Math.cos(u) * 0.06;
      const dx = (xR * cp - yR * sp) * jupRadiusRad;
      const dy = (xR * sp + yR * cp) * jupRadiusRad;
      return {
        name: g.name,
        raDeg: normDeg((ra0 + dx * invCosDec) * RAD),
        decDeg: (dec0 + dy) * RAD,
        inFront: Math.cos(u) > 0,
        mag: g.mag
      };
    });
  }

  function computePlanets(jd, precM, nut) {
    const earthHelio = planetHeliocentricEclJ2000(PLANET_ROWS.emb, jd);
    const earthDist = Math.hypot(earthHelio[0], earthHelio[1], earthHelio[2]);
    const list = ['mercury', 'venus', 'mars', 'jupiter', 'saturn', 'uranus', 'neptune'];

    return list.map(id => {
      const row = PLANET_ROWS[id];
      let helio = planetHeliocentricEclJ2000(row, jd);
      let geo = [helio[0] - earthHelio[0], helio[1] - earthHelio[1], helio[2] - earthHelio[2]];
      let delta = Math.hypot(geo[0], geo[1], geo[2]);
      const tau = delta * LIGHT_DAYS_PER_AU;
      const retarded = { ut: jd.ut - tau, tt: jd.tt - tau, tUT: (jd.ut - tau - 2451545) / 36525, tTT: (jd.tt - tau - 2451545) / 36525 };
      helio = planetHeliocentricEclJ2000(row, retarded);
      geo = [helio[0] - earthHelio[0], helio[1] - earthHelio[1], helio[2] - earthHelio[2]];

      const r = Math.hypot(helio[0], helio[1], helio[2]);
      delta = Math.hypot(geo[0], geo[1], geo[2]);
      const cosI = clamp((r * r + delta * delta - earthDist * earthDist) / (2 * r * delta), -1, 1);
      const iDeg = Math.acos(cosI) * RAD;
      const place = apparentFromEclipticJ2000(geo, jd, precM, nut);
      const angDiamDeg = 2 * Math.atan(row.radiusKm / (delta * AU_KM)) * RAD;

      let ringTiltDeg = 0;
      if (id === 'saturn') {
        // Approximate Saturn pole J2000: RA 40.589 deg, Dec 83.537 deg
        const pr = 40.589 * DEG, pd = 83.537 * DEG;
        const pole = [Math.cos(pd) * Math.cos(pr), Math.cos(pd) * Math.sin(pr), Math.sin(pd)];
        const geoEq = eclipticVecToEquatorialJ2000(geo);
        const uGeo = [geoEq[0] / delta, geoEq[1] / delta, geoEq[2] / delta];
        ringTiltDeg = Math.asin(clamp(-(pole[0]*uGeo[0] + pole[1]*uGeo[1] + pole[2]*uGeo[2]), -1, 1)) * RAD;
      }

      return {
        id,
        kind: 'planet',
        name: row.name,
        raRad: place.raRad,
        decRad: place.decRad,
        raDeg: place.raDeg,
        decDeg: place.decDeg,
        distanceAU: delta,
        heliocentricAU: r,
        phaseAngleDeg: iDeg,
        illuminatedFraction: (1 + Math.cos(iDeg * DEG)) * 0.5,
        magnitude: planetMagnitude(id, r, delta, iDeg),
        angularDiameterDeg: angDiamDeg,
        ringTiltDeg,
        galileanMoons: id === 'jupiter' ? computeGalileanMoons(jd, place, delta) : null
      };
    });
  }

  // ---- Stumpff Universal Variable Small Body Orbit Solver ----
  function stumpff(z) {
    if (z > 1e-6) {
      const sz = Math.sqrt(z);
      return { c2: (1 - Math.cos(sz)) / z, c3: (sz - Math.sin(sz)) / (z * sz) };
    }
    if (z < -1e-6) {
      const sz = Math.sqrt(-z);
      return { c2: (Math.cosh(sz) - 1) / (-z), c3: (Math.sinh(sz) - sz) / ((-z) * sz) };
    }
    return { c2: 0.5 - z / 24 + (z * z) / 720, c3: 1 / 6 - z / 120 + (z * z) / 5040 };
  }

  // Solve perifocal (p, q, dp, dq) from perihelion qAU, eccentricity e, and dtDays = jd - periJD
  function perifocalAt(qAU, e, dtDays) {
    const mu = K_GAUSS * K_GAUSS;
    const sqrtMu = K_GAUSS;
    const alpha = (1 - e) / qAU;
    const h = Math.sqrt(mu * qAU * (1 + e));

    // Initial chi guess
    let chi;
    if (Math.abs(alpha) < 1e-5) {
      const s = 3 * sqrtMu * dtDays / (2 * Math.pow(qAU, 1.5));
      chi = Math.sqrt(qAU) * (Math.cbrt(s + Math.sqrt(s * s + 1)) + Math.cbrt(s - Math.sqrt(s * s + 1))) * 0.8;
    } else if (alpha > 0) {
      chi = sqrtMu * dtDays * alpha;
    } else {
      const sgn = dtDays >= 0 ? 1 : -1;
      const a = 1 / alpha;
      const val = Math.max(1.1, Math.abs(-2 * mu * alpha * dtDays / ((e - 1) * sqrtMu * Math.sqrt(-a) + 1e-6)));
      chi = sgn * Math.sqrt(-a) * Math.log(val);
    }
    if (!Number.isFinite(chi)) chi = 0;

    for (let iter = 0; iter < 30; iter++) {
      const z = alpha * chi * chi;
      const { c2, c3 } = stumpff(z);
      const chi2 = chi * chi;
      const chi3 = chi2 * chi;
      const r = qAU * (1 - z * c2) + chi2 * c2;
      const F = qAU * chi * (1 - z * c3) + chi3 * c3 - sqrtMu * dtDays;
      const dF = r;
      const d2F = chi * (1 - z * c3) * (1 - alpha * qAU);
      const dFsafe = Math.abs(dF) > 1e-12 ? dF : 1e-12;
      let step = F / dFsafe;
      if (Math.abs(step) > 0.5) {
        const disc = Math.max(0, 16 * dFsafe * dFsafe - 20 * F * d2F);
        step = (5 * F) / (dFsafe + (dFsafe >= 0 ? 1 : -1) * Math.sqrt(disc) || 1);
      }
      chi -= step;
      if (Math.abs(step) < 1e-10) break;
    }

    const z = alpha * chi * chi;
    const { c2, c3 } = stumpff(z);
    const f = 1 - (chi * chi / qAU) * c2;
    const g = dtDays - (chi * chi * chi / sqrtMu) * c3;
    const p = f * qAU;
    const q = g * Math.sqrt(mu * (1 + e) / qAU);
    return [p, q];
  }

  function computeOrbitalElementsHelio(qAU, e, iDeg, nodeDeg, wDeg, dtDays) {
    const [p, q] = perifocalAt(qAU, e, dtDays);
    const w = wDeg * DEG, Om = nodeDeg * DEG, inc = iDeg * DEG;
    const cw = Math.cos(w), sw = Math.sin(w);
    const cO = Math.cos(Om), sO = Math.sin(Om);
    const cI = Math.cos(inc), sI = Math.sin(inc);

    return [
      (cw * cO - sw * sO * cI) * p + (-sw * cO - cw * sO * cI) * q,
      (cw * sO + sw * cO * cI) * p + (-sw * sO + cw * cO * cI) * q,
      (sw * sI) * p + (cw * sI) * q
    ];
  }

  function computeAsteroid(ast, jd, earthHelio, earthDist, precM, nut) {
    // [number, name, kind, a_au, e, i_deg, node_deg, peri_deg, m_deg, epoch_jd, h, g]
    const [, name, kind, a, e, iDeg, nodeDeg, wDeg, mDeg, epochJD, hMag, gSlope] = ast;
    const nDegDay = (K_GAUSS * RAD) / Math.pow(a, 1.5);
    // Shift to perihelion time nearest to jd
    const mAtJd = normDeg(mDeg + nDegDay * (jd.tt - epochJD));
    const mSigned = mAtJd > 180 ? mAtJd - 360 : mAtJd;
    const dtDays = mSigned / nDegDay;
    const qAU = a * (1 - e);

    let helio = computeOrbitalElementsHelio(qAU, e, iDeg, nodeDeg, wDeg, dtDays);
    let geo = [helio[0] - earthHelio[0], helio[1] - earthHelio[1], helio[2] - earthHelio[2]];
    const delta0 = Math.hypot(geo[0], geo[1], geo[2]);
    const tau = delta0 * LIGHT_DAYS_PER_AU;
    helio = computeOrbitalElementsHelio(qAU, e, iDeg, nodeDeg, wDeg, dtDays - tau);
    geo = [helio[0] - earthHelio[0], helio[1] - earthHelio[1], helio[2] - earthHelio[2]];

    const r = Math.hypot(helio[0], helio[1], helio[2]);
    const delta = Math.hypot(geo[0], geo[1], geo[2]);
    const cosPhi = clamp((r * r + delta * delta - earthDist * earthDist) / (2 * r * delta), -1, 1);
    const alpha = Math.acos(cosPhi);
    const phi1 = Math.exp(-3.33 * Math.pow(Math.tan(alpha * 0.5), 0.63));
    const phi2 = Math.exp(-1.87 * Math.pow(Math.tan(alpha * 0.5), 1.22));
    const vmag = hMag + 5 * Math.log10(Math.max(r * delta, 1e-4)) - 2.5 * Math.log10((1 - gSlope) * phi1 + gSlope * phi2);
    const place = apparentFromEclipticJ2000(geo, jd, precM, nut);

    return {
      id: 'ast:' + ast[0],
      kind: 'asteroid',
      subkind: kind,
      name,
      raRad: place.raRad,
      decRad: place.decRad,
      raDeg: place.raDeg,
      decDeg: place.decDeg,
      distanceAU: delta,
      heliocentricAU: r,
      magnitude: vmag
    };
  }

  function computeComet(comet, jd, earthHelio, earthDist, precM, nut) {
    // [id, name, desig, q_au, e, i_deg, node_deg, w_deg, peri_jd, h0, slope_n, peak_iso, note]
    const [id, name, desig, qAU, e, iDeg, nodeDeg, wDeg, periJD, h0, slopeN, peakISO, note] = comet;
    const dt = jd.tt - periJD;
    let helio = computeOrbitalElementsHelio(qAU, e, iDeg, nodeDeg, wDeg, dt);
    let geo = [helio[0] - earthHelio[0], helio[1] - earthHelio[1], helio[2] - earthHelio[2]];
    let delta = Math.hypot(geo[0], geo[1], geo[2]);
    const tau = delta * LIGHT_DAYS_PER_AU;
    helio = computeOrbitalElementsHelio(qAU, e, iDeg, nodeDeg, wDeg, dt - tau);
    geo = [helio[0] - earthHelio[0], helio[1] - earthHelio[1], helio[2] - earthHelio[2]];
    delta = Math.hypot(geo[0], geo[1], geo[2]);
    const r = Math.hypot(helio[0], helio[1], helio[2]);

    const place = apparentFromEclipticJ2000(geo, jd, precM, nut);

    // 3D anti-solar ion tail endpoint & lagged dust tail endpoint projected on sky
    const uHelio = [helio[0] / (r || 1), helio[1] / (r || 1), helio[2] / (r || 1)];
    const tailLenAU = clamp(0.14 / Math.pow(Math.max(r, 0.25), 0.75), 0.03, 0.38);
    const ionTipGeo = [geo[0] + uHelio[0] * tailLenAU, geo[1] + uHelio[1] * tailLenAU, geo[2] + uHelio[2] * tailLenAU];
    const ionTipPlace = apparentFromEclipticJ2000(ionTipGeo, jd, precM, nut);

    // Dust tail lags behind anti-solar vector along negative orbital velocity
    const helioPrev = computeOrbitalElementsHelio(qAU, e, iDeg, nodeDeg, wDeg, dt - tau - 0.5);
    const vel = [helio[0] - helioPrev[0], helio[1] - helioPrev[1], helio[2] - helioPrev[2]];
    const vNorm = Math.hypot(vel[0], vel[1], vel[2]) || 1;
    const dustDir = [
      uHelio[0] * 0.78 - (vel[0] / vNorm) * 0.55,
      uHelio[1] * 0.78 - (vel[1] / vNorm) * 0.55,
      uHelio[2] * 0.78 - (vel[2] / vNorm) * 0.55
    ];
    const dNorm = Math.hypot(dustDir[0], dustDir[1], dustDir[2]) || 1;
    const dustLenAU = tailLenAU * 0.72;
    const dustTipGeo = [
      geo[0] + (dustDir[0] / dNorm) * dustLenAU,
      geo[1] + (dustDir[1] / dNorm) * dustLenAU,
      geo[2] + (dustDir[2] / dNorm) * dustLenAU
    ];
    const dustTipPlace = apparentFromEclipticJ2000(dustTipGeo, jd, precM, nut);

    const cosElong = clamp(-(geo[0] * earthHelio[0] + geo[1] * earthHelio[1] + geo[2] * earthHelio[2]) / ((delta || 1) * (earthDist || 1)), -1, 1);
    const elongationDeg = Math.acos(cosElong) * RAD;
    const vmag = h0 + 5 * Math.log10(Math.max(delta, 0.02)) + 2.5 * slopeN * Math.log10(Math.max(r, 0.08));

    return {
      id: 'comet:' + id,
      rawId: id,
      kind: 'comet',
      name,
      desig,
      peakISO,
      note,
      raRad: place.raRad,
      decRad: place.decRad,
      raDeg: place.raDeg,
      decDeg: place.decDeg,
      ionTipRaDeg: ionTipPlace.raDeg,
      ionTipDecDeg: ionTipPlace.decDeg,
      dustTipRaDeg: dustTipPlace.raDeg,
      dustTipDecDeg: dustTipPlace.decDeg,
      distanceAU: delta,
      heliocentricAU: r,
      elongationDeg,
      magnitude: vmag
    };
  }

  // ---- Ballesteros (2012) Star B-V Blackbody Color ----
  function bvToRgb(bv) {
    const c = clamp(bv == null ? 0.4 : bv, -0.4, 2.0);
    const t = 4600 * (1 / (0.92 * c + 1.7) + 1 / (0.92 * c + 0.62));
    // Planck locus approximation -> sRGB mixed 45% toward pure starlight white
    const k = t / 100;
    let r, g, b;
    if (k <= 66) {
      r = 255;
      g = clamp(99.4708025861 * Math.log(k) - 161.1195681661, 0, 255);
      b = k <= 19 ? 0 : clamp(138.5177312231 * Math.log(k - 10) - 305.0447927307, 0, 255);
    } else {
      r = clamp(329.698727446 * Math.pow(k - 60, -0.1332047592), 0, 255);
      g = clamp(288.1221695283 * Math.pow(k - 60, -0.0755148492), 0, 255);
      b = 255;
    }
    // Mix 45% toward white (StarAppearance.swift)
    return [
      Math.round(r * 0.55 + 255 * 0.45),
      Math.round(g * 0.55 + 255 * 0.45),
      Math.round(b * 0.55 + 255 * 0.45)
    ];
  }

  // ---- Deterministic Milky Way Procedural Blobs (Galactic -> J2000) ----
  const GAL_TO_J2000 = [
    [-0.0548755604, -0.8734370902, -0.4838350155],
    [ 0.4941094279, -0.4448296300,  0.7469822445],
    [-0.8676661490, -0.1980763734,  0.4559837762]
  ];

  function mulberry32(seed) {
    let a = seed >>> 0;
    return function () {
      a |= 0;
      a = (a + 0x6D2B79F5) | 0;
      let t = Math.imul(a ^ (a >>> 15), 1 | a);
      t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }

  function randGauss(rng) {
    const u1 = Math.max(1e-7, rng()), u2 = rng();
    return Math.sqrt(-2 * Math.log(u1)) * Math.cos(TWO_PI * u2);
  }

  function galToJ2000(lDeg, bDeg) {
    const l = lDeg * DEG, b = bDeg * DEG;
    const cb = Math.cos(b);
    const g = [cb * Math.cos(l), cb * Math.sin(l), Math.sin(b)];
    // Transpose of GAL_TO_J2000 rows? Wait: rows above ARE J2000->Galactic rows, so transpose converts Gal -> J2000:
    return [
      GAL_TO_J2000[0][0]*g[0] + GAL_TO_J2000[1][0]*g[1] + GAL_TO_J2000[2][0]*g[2],
      GAL_TO_J2000[0][1]*g[0] + GAL_TO_J2000[1][1]*g[1] + GAL_TO_J2000[2][1]*g[2],
      GAL_TO_J2000[0][2]*g[0] + GAL_TO_J2000[1][2]*g[1] + GAL_TO_J2000[2][2]*g[2]
    ];
  }

  function generateMilkyWayBlobs() {
    const rng = mulberry32(0x5A7B9C1D);
    const blobs = [];

    // 1. Bright galactic bulge & spiral arms (Sagittarius, Scutum, Norma, Cygnus, Carina)
    for (let i = 0; i < 620; i++) {
      let lDeg;
      const roll = rng();
      if (roll < 0.34) {
        // Sagittarius–Scutum bulge concentration (-32..+32 deg)
        lDeg = randGauss(rng) * 18;
      } else if (roll < 0.50) {
        // Cygnus star cloud (~75 deg)
        lDeg = 76 + randGauss(rng) * 15;
      } else if (roll < 0.66) {
        // Carina-Vela bright stretch (~288 deg)
        lDeg = 288 + randGauss(rng) * 20;
      } else {
        lDeg = rng() * 360;
      }
      const lRad = lDeg * DEG;
      const bulge = Math.exp(-Math.pow(Math.sin(lRad * 0.5), 2) * 8);
      const bDeg = randGauss(rng) * (4.2 + bulge * 4.5);
      const radiusDeg = 3.2 + rng() * (5.2 + bulge * 2.8);
      const intensity = (0.042 + rng() * 0.068) * (0.58 + 0.85 * bulge);
      blobs.push({
        vec: galToJ2000(lDeg, bDeg),
        radiusDeg,
        intensity,
        warm: bulge * 0.85 + rng() * 0.15,
        dark: false
      });
    }

    // 2. Broader outer Milky Way haze
    for (let i = 0; i < 240; i++) {
      const lDeg = rng() * 360;
      const bDeg = randGauss(rng) * 11.5;
      blobs.push({
        vec: galToJ2000(lDeg, bDeg),
        radiusDeg: 7.5 + rng() * 7.0,
        intensity: 0.025 + rng() * 0.032,
        warm: 0.25,
        dark: false
      });
    }

    // 3. Great Rift & Coalsack dark dust lanes
    for (let i = 0; i < 120; i++) {
      const inRift = rng() < 0.85;
      const lDeg = inRift ? (-8 + rng() * 92) : (300 + rng() * 14);
      const bDeg = inRift ? (0.8 + randGauss(rng) * 1.9) : (-0.6 + randGauss(rng) * 1.4);
      blobs.push({
        vec: galToJ2000(lDeg, bDeg),
        radiusDeg: inRift ? 2.6 + rng() * 3.4 : 2.2 + rng() * 2.0,
        intensity: 0.16 + rng() * 0.18,
        warm: 0,
        dark: true
      });
    }
    return blobs;
  }

  // ---- Jagged Horizon Terrain Profile matching Starward ----
  function terrainAltitudeDeg(azDeg) {
    const t = azDeg * DEG;
    const ridge1 = 1.45 * (1 - Math.abs(Math.sin(t * 2.5 + 0.4)));
    const ridge2 = 0.88 * (1 - Math.abs(Math.sin(t * 5.0 - 1.15)));
    const ridge3 = 0.56 * (1 - Math.abs(Math.sin(t * 9.0 + 0.7)));
    const ridge4 = 0.28 * (1 - Math.abs(Math.cos(t * 17.0 - 0.3)));
    const crags = 0.14 * Math.sin(t * 31.0 + 1.2) + 0.08 * (1 - Math.abs(Math.sin(t * 47.0)));
    return Math.max(0.45, 0.7 + ridge1 + ridge2 + ridge3 + ridge4 + crags);
  }

  // ---- Darkness State & Windows (DarknessPill) ----
  function darknessStateForSunAlt(altDeg) {
    if (altDeg >= -0.833) return 'day';
    if (altDeg >= -6) return 'civil';
    if (altDeg >= -12) return 'nautical';
    if (altDeg >= -18) return 'astronomical';
    return 'dark';
  }

  function sunAltAtDate(date, latDeg, lonDeg) {
    const frame = createFrameContext(date, latDeg, lonDeg);
    const sun = computeSun(frame.jd, frame.nut);
    return apparentRaDecToHoriz(sun.raDeg, sun.decDeg, frame, true).altDeg;
  }

  function formatTimeHHMM(date, tzOffsetMin) {
    const ms = date.getTime() + tzOffsetMin * 60000;
    const d = new Date(ms);
    const hh = String(d.getUTCHours()).padStart(2, '0');
    const mm = String(d.getUTCMinutes()).padStart(2, '0');
    return `${hh}:${mm}`;
  }

  function computeDarknessStatus(date, latDeg, lonDeg, tzOffsetMin) {
    const altNow = sunAltAtDate(date, latDeg, lonDeg);
    const state = darknessStateForSunAlt(altNow);
    const baseMs = date.getTime();
    const stepMs = 5 * 60000;
    const maxSteps = 288; // 24 hours

    function findTransition(targetPredicate) {
      let prevMs = baseMs;
      let prevVal = sunAltAtDate(new Date(prevMs), latDeg, lonDeg);
      if (!targetPredicate(prevVal)) return null;
      for (let i = 1; i <= maxSteps; i++) {
        const curMs = baseMs + i * stepMs;
        const curVal = sunAltAtDate(new Date(curMs), latDeg, lonDeg);
        if (!targetPredicate(curVal)) {
          // Binary bisect to <30s
          let lo = prevMs, hi = curMs;
          for (let b = 0; b < 8; b++) {
            const mid = (lo + hi) * 0.5;
            if (targetPredicate(sunAltAtDate(new Date(mid), latDeg, lonDeg))) lo = mid;
            else hi = mid;
          }
          return new Date((lo + hi) * 0.5);
        }
        prevMs = curMs;
      }
      return null;
    }

    function findFirstEntry(targetPredicate) {
      let prevMs = baseMs;
      for (let i = 1; i <= maxSteps; i++) {
        const curMs = baseMs + i * stepMs;
        const curVal = sunAltAtDate(new Date(curMs), latDeg, lonDeg);
        if (targetPredicate(curVal)) {
          let lo = prevMs, hi = curMs;
          for (let b = 0; b < 8; b++) {
            const mid = (lo + hi) * 0.5;
            if (targetPredicate(sunAltAtDate(new Date(mid), latDeg, lonDeg))) hi = mid;
            else lo = mid;
          }
          return new Date((lo + hi) * 0.5);
        }
        prevMs = curMs;
      }
      return null;
    }

    let label = 'Dark sky';
    let detail = '';

    if (state === 'dark') {
      label = 'Dark sky';
      const endsAt = findTransition(a => a < -18);
      detail = endsAt ? `until ${formatTimeHHMM(endsAt, tzOffsetMin)}` : 'all night';
    } else if (state === 'astronomical' || state === 'nautical' || state === 'civil') {
      label = state === 'astronomical' ? 'Astronomical twilight' : state === 'nautical' ? 'Nautical twilight' : 'Civil twilight';
      const darkStarts = findFirstEntry(a => a < -18);
      const sunRises = findFirstEntry(a => a >= -0.833);
      if (darkStarts && (!sunRises || darkStarts.getTime() < sunRises.getTime())) {
        detail = `dark at ${formatTimeHHMM(darkStarts, tzOffsetMin)}`;
      } else if (sunRises) {
        detail = `sunrise ${formatTimeHHMM(sunRises, tzOffsetMin)}`;
      }
    } else {
      label = 'Daylight';
      const sunset = findFirstEntry(a => a < -0.833);
      detail = sunset ? `sunset ${formatTimeHHMM(sunset, tzOffsetMin)}` : 'sun up';
    }

    return { state, sunAltDeg: altNow, label, detail };
  }

  return {
    DEG,
    RAD,
    TWO_PI,
    clamp,
    normDeg,
    normRad,
    smoothstep,
    julianDate,
    nutation,
    greenwichApparentSiderealRad,
    createFrameContext,
    refractDeg,
    raDecDegToVec,
    mat3MulVec,
    apparentRaDecToHoriz,
    j2000VecToHoriz,
    horizToApparentRaDec,
    createProjection,
    computeSun,
    computeMoon,
    computePlanets,
    planetHeliocentricEclJ2000,
    PLANET_ROWS,
    computeAsteroid,
    computeComet,
    bvToRgb,
    generateMilkyWayBlobs,
    terrainAltitudeDeg,
    darknessStateForSunAlt,
    computeDarknessStatus,
    formatTimeHHMM
  };
}));

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Ledicia Díaz Lago
MD one_run analysis

Lee:
    outputs/one_run/instantaneous_energies.dat
    outputs/one_run/corr_*.dat
    outputs/one_run/corrmean_*.dat
    outputs/one_run/md_final_results.txt
    outputs/one_run/rva.dat   (Fortran unformatted sequential)  [opcional]

Produce:
    outputs/one_run/analysis/*.png
    outputs/one_run/analysis/summary.json
    outputs/one_run/analysis/summary.txt

Plots:
  1) Series temporales: epot/ekin/etot/T/P
  2) Medias móviles
  3) Histogramas: T, P, Etot
  4) Dispersión: P vs T, Epot vs Ekin, Etot vs T
  5) ACF desde corr_*.dat + corrmean_*.dat (normalizada)
  6) Estimación de tiempos de correlación (tau_int, 1/e)
  7) Desde rva.dat (si existe):
        - MSD(τ) promediada en orígenes temporales (time-origin averaged)
        - VACF(τ) promediada en orígenes temporales
        - RDF g(r)
        - Difusión (coeficiente D) por dos vías físicas estándar:
            * Einstein:  MSD(t) ≈ 6 D t  (régimen difusivo a tiempos largos)
            * Green–Kubo: D = (1/3) ∫_0^∞ VACF(t) dt

Por qué hay "fit_frac" y "tmax_frac" (y por qué NO usar siempre todo):
  - MSD:
      A tiempos muy cortos el movimiento es balístico (MSD ~ t^2).
      La relación MSD ≈ 6Dt solo vale en el régimen difusivo (tiempos largos, MSD ~ t).
      Si ajustas con todo, mezclas t^2 y t y sesgas D.
      -> Por eso se ajusta solo la “cola” (p.ej. último 50%).
  - VACF:
      En teoría la integral es hasta ∞, pero en datos finitos:
        * la VACF decae a ~0
        * después queda dominada por ruido estadístico alrededor de 0
      Integrar muy lejos puede acumular ruido y hacer D inestable.
      -> Por eso puede interesar cortar la integral a una fracción del tiempo.

Flags CLI (opcionales; si NO las pones, se comporta como antes):
  --msd-fit-frac     fracción final del MSD usada para el ajuste lineal (default: 0.5)
  --vacf-tmax-frac   fracción del tiempo total usada para integrar VACF (default: 1.0 = todo, como antes)
  --rva-max-lag      max lag en snapshots para MSD/VACF (default: auto como antes)
  --origin-stride    salto entre orígenes temporales para promediar (default: auto como antes)

Uso:
  python3 scripts/md_one_run_analysis.py

Opcional:
  python3 md_one_run_analysis.py --root ..
  python3 md_one_run_analysis.py --skip-rva
  python3 md_one_run_analysis.py --msd-fit-frac 0.4 --vacf-tmax-frac 0.8
"""

from __future__ import annotations

import argparse
import json
import math
import re
import struct
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import matplotlib.pyplot as plt


# -----------------------------
# Helpers: safe IO + parsing
# -----------------------------
def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def read_text_table(path: Path, comment_prefix: str = "#") -> np.ndarray:
    """
    Reads whitespace-separated numeric table with optional comment lines.
    Returns a 2D float array (nrows, ncols). Raises if empty.
    """
    if not path.exists():
        raise FileNotFoundError(f"Missing file: {path}")
    rows: List[List[float]] = []
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith(comment_prefix):
                continue
            parts = s.split()
            try:
                rows.append([float(x) for x in parts])
            except ValueError:
                continue
    if not rows:
        raise ValueError(f"No numeric data found in: {path}")
    return np.array(rows, dtype=float)


def parse_md_final_results(path: Path) -> Dict[str, float]:
    """
    Very tolerant parser for md_final_results.txt (append-mode file).
    Extracts last block of results and returns numeric key->value for common fields.
    """
    out: Dict[str, float] = {}
    if not path.exists():
        return out

    text = path.read_text(encoding="utf-8", errors="replace")
    blocks = text.split("************** MD PRODUCTION RESULTS **************")
    if len(blocks) < 2:
        return out
    last = blocks[-1]

    for line in last.splitlines():
        line = line.strip()
        if not line:
            continue

        # key: value
        m = re.match(r"^([A-Za-z0-9_<>\-]+)\s*:\s*([+\-0-9.eEdD]+)\s*$", line)
        if m:
            k = m.group(1)
            v = m.group(2).replace("D", "e").replace("d", "e")
            try:
                out[k] = float(v)
            except ValueError:
                pass
            continue

        # <Epot>: val std: val
        m2 = re.match(
            r"^(<[^>]+>|[A-Za-z_0-9]+)\s*:\s*([+\-0-9.eEdD]+)\s+std\s*:\s*([+\-0-9.eEdD]+)",
            line,
        )
        if m2:
            k = m2.group(1)
            v = m2.group(2).replace("D", "e").replace("d", "e")
            s = m2.group(3).replace("D", "e").replace("d", "e")
            try:
                out[k] = float(v)
                out[k + "_std"] = float(s)
            except ValueError:
                pass
            continue

        # Temperature: X Pressure: Y
        m3 = re.match(
            r"^Temperature\s*:\s*([+\-0-9.eEdD]+)\s+Pressure\s*:\s*([+\-0-9.eEdD]+)",
            line,
        )
        if m3:
            try:
                out["Temperature"] = float(m3.group(1).replace("D", "e").replace("d", "e"))
                out["Pressure"] = float(m3.group(2).replace("D", "e").replace("d", "e"))
            except ValueError:
                pass
            continue

    return out


# -----------------------------
# ACF utilities
# -----------------------------
@dataclass
class CorrInfo:
    name: str
    lag: np.ndarray
    c: np.ndarray
    cn: np.ndarray


def load_corr_files(one_run_dir: Path) -> Tuple[Dict[str, CorrInfo], Dict[str, CorrInfo]]:
    """
    Loads corr_*.dat and corrmean_*.dat if present.
    Returns: (corr, corrmean) dicts by observable name.
    """
    corr: Dict[str, CorrInfo] = {}
    corrmean: Dict[str, CorrInfo] = {}

    for obs in ["epot", "ekin", "etot", "temp", "press"]:
        p = one_run_dir / f"corr_{obs}.dat"
        if p.exists():
            arr = read_text_table(p)
            lag = arr[:, 0].astype(int)
            c = arr[:, 1]
            cn = arr[:, 2]
            corr[obs] = CorrInfo(obs, lag, c, cn)

        pm = one_run_dir / f"corrmean_{obs}.dat"
        if pm.exists():
            arrm = read_text_table(pm)
            lagm = arrm[:, 0].astype(int)
            cm = arrm[:, 1]
            cnm = arrm[:, 2]
            corrmean[obs] = CorrInfo(obs, lagm, cm, cnm)

    return corr, corrmean


def estimate_tau_int(cn: np.ndarray, dt_sample: float, max_lag: int) -> float:
    """
    Integrated autocorrelation time (rough):
      tau_int ≈ dt * (1 + 2 * sum_{l=1..L*} cn[l])
    We stop at first non-positive cn or at max_lag (whichever first).
    """
    if len(cn) < 2:
        return float("nan")
    L = min(max_lag, len(cn) - 1)
    s = 0.0
    for l in range(1, L + 1):
        if cn[l] <= 0:
            break
        s += cn[l]
    return dt_sample * (1.0 + 2.0 * s)


def estimate_tau_1e(cn: np.ndarray, dt_sample: float) -> float:
    """Time when cn drops below 1/e. Linear interpolation between nearest lags."""
    target = 1.0 / math.e
    if len(cn) < 2:
        return float("nan")
    for i in range(1, len(cn)):
        if cn[i] <= target:
            x0, x1 = (i - 1) * dt_sample, i * dt_sample
            y0, y1 = cn[i - 1], cn[i]
            if y1 == y0:
                return x1
            t = (target - y0) / (y1 - y0)
            return x0 + t * (x1 - x0)
    return float("nan")


# -----------------------------
# Fortran unformatted sequential reader (rva.dat)
# -----------------------------
class FortranSequentialReader:
    """
    Reads Fortran unformatted sequential records with record markers.
    Supports 4-byte or 8-byte record markers (auto-detect per file).
    """

    def __init__(self, path: Path):
        self.path = path
        self.f = path.open("rb")
        self.marker_size = self._detect_marker_size()

    def close(self) -> None:
        try:
            self.f.close()
        except Exception:
            pass

    def _detect_marker_size(self) -> int:
        """
        Heuristic:
          - read first 32 bytes, try interpret first 4 or 8 as record length
          - record length should be positive and "reasonable"
          - trailing marker should match
        """
        pos0 = self.f.tell()
        head = self.f.read(64)
        self.f.seek(pos0)

        # try 4-byte
        if len(head) >= 8:
            n4 = struct.unpack("<i", head[0:4])[0]
            if 0 < n4 < 1_000_000_000 and len(head) >= 4 + n4 + 4:
                tail = struct.unpack("<i", head[4 + n4 : 4 + n4 + 4])[0]
                if tail == n4:
                    return 4

        # try 8-byte
        if len(head) >= 16:
            n8 = struct.unpack("<q", head[0:8])[0]
            if 0 < n8 < 1_000_000_000 and len(head) >= 8 + n8 + 8:
                tail = struct.unpack("<q", head[8 + n8 : 8 + n8 + 8])[0]
                if tail == n8:
                    return 8

        return 4  # common with gfortran

    def read_record_bytes(self) -> bytes:
        ms = self.marker_size
        hdr = self.f.read(ms)
        if not hdr:
            raise EOFError
        if len(hdr) != ms:
            raise IOError("Truncated record header")
        n = struct.unpack("<i" if ms == 4 else "<q", hdr)[0]
        if n < 0:
            raise IOError(f"Negative record length: {n}")
        data = self.f.read(n)
        if len(data) != n:
            raise IOError("Truncated record payload")
        ftr = self.f.read(ms)
        if len(ftr) != ms:
            raise IOError("Truncated record footer")
        n2 = struct.unpack("<i" if ms == 4 else "<q", ftr)[0]
        if n2 != n:
            raise IOError(f"Record marker mismatch: {n} vs {n2}")
        return data

    def read_record(self, fmt: str) -> Tuple:
        data = self.read_record_bytes()
        size = struct.calcsize(fmt)
        if len(data) != size:
            raise IOError(f"Record size mismatch: got {len(data)} bytes, expected {size}")
        return struct.unpack(fmt, data)

    def read_record_array(self, dtype: np.dtype, count: int) -> np.ndarray:
        data = self.read_record_bytes()
        arr = np.frombuffer(data, dtype=dtype, count=count)
        if arr.size != count:
            raise IOError(f"Array record size mismatch: got {arr.size}, expected {count}")
        return arr.copy()


@dataclass
class RVAData:
    n: int
    L: float
    dt: float
    output_interval: int
    n_snapshots: int
    rx: np.ndarray
    ry: np.ndarray
    rz: np.ndarray
    rux: np.ndarray
    ruy: np.ndarray
    ruz: np.ndarray
    vx: np.ndarray
    vy: np.ndarray
    vz: np.ndarray


def read_rva(path: Path) -> RVAData:
    """
    Expects file format written by your md_simulation:
      header record: n, box_length, dt, output_interval, n_snapshots_expected
      per snapshot:
        record1: rx, ry, rz      (3*n doubles)
        record2: rux, ruy, ruz   (3*n doubles)
        record3: vx, vy, vz      (3*n doubles)
        record4: ax, ay, az      (3*n doubles)  -> ignored
    """
    r = FortranSequentialReader(path)
    try:
        n, L, dt, out_int, n_snap = r.read_record("<i d d i i")
        n = int(n)
        out_int = int(out_int)
        n_snap = int(n_snap)

        rx = np.empty((n_snap, n), dtype=np.float64)
        ry = np.empty((n_snap, n), dtype=np.float64)
        rz = np.empty((n_snap, n), dtype=np.float64)
        rux = np.empty((n_snap, n), dtype=np.float64)
        ruy = np.empty((n_snap, n), dtype=np.float64)
        ruz = np.empty((n_snap, n), dtype=np.float64)
        vx = np.empty((n_snap, n), dtype=np.float64)
        vy = np.empty((n_snap, n), dtype=np.float64)
        vz = np.empty((n_snap, n), dtype=np.float64)

        count3 = 3 * n

        for s in range(n_snap):
            a = r.read_record_array(np.float64, count3)
            rx[s, :] = a[0:n]
            ry[s, :] = a[n:2 * n]
            rz[s, :] = a[2 * n:3 * n]

            a = r.read_record_array(np.float64, count3)
            rux[s, :] = a[0:n]
            ruy[s, :] = a[n:2 * n]
            ruz[s, :] = a[2 * n:3 * n]

            a = r.read_record_array(np.float64, count3)
            vx[s, :] = a[0:n]
            vy[s, :] = a[n:2 * n]
            vz[s, :] = a[2 * n:3 * n]

            _ = r.read_record_array(np.float64, count3)  # accel ignored

        return RVAData(
            n=n, L=float(L), dt=float(dt), output_interval=out_int, n_snapshots=n_snap,
            rx=rx, ry=ry, rz=rz, rux=rux, ruy=ruy, ruz=ruz, vx=vx, vy=vy, vz=vz
        )
    finally:
        r.close()


# -----------------------------
# MD observables from rva
#   - time-origin averaged MSD(τ) and VACF(τ)
# -----------------------------
def compute_msd_tau_timeorig(
    rux: np.ndarray, ruy: np.ndarray, ruz: np.ndarray,
    max_lag: Optional[int] = None,
    origin_stride: int = 1,
) -> np.ndarray:
    """
    MSD(τ) time-origin averaged:
      MSD(τ) = < |ru(t+τ) - ru(t)|^2 >_{particles, time origins}

    Arrays shape: (n_snap, n)
    Returns msd_tau shape: (max_lag+1,)
    """
    n_snap, n = rux.shape
    if n_snap < 2:
        return np.array([0.0], dtype=np.float64)

    if max_lag is None:
        max_lag = n_snap - 1
    max_lag = int(min(max_lag, n_snap - 1))
    origin_stride = max(1, int(origin_stride))

    msd = np.zeros(max_lag + 1, dtype=np.float64)
    counts = np.zeros(max_lag + 1, dtype=np.int64)

    for t0 in range(0, n_snap - 1, origin_stride):
        L = min(max_lag, (n_snap - 1) - t0)
        if L <= 0:
            continue
        dx = rux[t0:t0 + L + 1, :] - rux[t0, :][None, :]
        dy = ruy[t0:t0 + L + 1, :] - ruy[t0, :][None, :]
        dz = ruz[t0:t0 + L + 1, :] - ruz[t0, :][None, :]
        d2 = np.mean(dx*dx + dy*dy + dz*dz, axis=1)
        msd[:L + 1] += d2
        counts[:L + 1] += 1

    mask = counts > 0
    msd[mask] /= counts[mask]
    return msd


def compute_vacf_tau_timeorig(
    vx: np.ndarray, vy: np.ndarray, vz: np.ndarray,
    max_lag: Optional[int] = None,
    origin_stride: int = 1,
) -> np.ndarray:
    """
    VACF(τ) time-origin averaged:
      VACF(τ) = < v(t) · v(t+τ) >_{particles, time origins}

    Arrays shape: (n_snap, n)
    Returns vacf_tau shape: (max_lag+1,)
    """
    n_snap, n = vx.shape
    if n_snap < 2:
        v0 = np.mean(vx[0]*vx[0] + vy[0]*vy[0] + vz[0]*vz[0])
        return np.array([v0], dtype=np.float64)

    if max_lag is None:
        max_lag = n_snap - 1
    max_lag = int(min(max_lag, n_snap - 1))
    origin_stride = max(1, int(origin_stride))

    vacf = np.zeros(max_lag + 1, dtype=np.float64)
    counts = np.zeros(max_lag + 1, dtype=np.int64)

    for t0 in range(0, n_snap - 1, origin_stride):
        L = min(max_lag, (n_snap - 1) - t0)
        if L <= 0:
            continue
        v0x = vx[t0, :][None, :]
        v0y = vy[t0, :][None, :]
        v0z = vz[t0, :][None, :]

        dot = np.mean(
            vx[t0:t0 + L + 1, :]*v0x +
            vy[t0:t0 + L + 1, :]*v0y +
            vz[t0:t0 + L + 1, :]*v0z,
            axis=1
        )

        vacf[:L + 1] += dot
        counts[:L + 1] += 1

    mask = counts > 0
    vacf[mask] /= counts[mask]
    return vacf


def diffusion_from_msd_linear_fit(t: np.ndarray, msd: np.ndarray, fit_frac: float = 0.5) -> float:
    """
    Einstein relation (3D):
      MSD(t) ≈ 6 D t   (solo en régimen difusivo a tiempos largos)

    Estimamos D ajustando linealmente la cola de MSD:
      - usamos la fracción final fit_frac del intervalo temporal
      - slope ≈ 6D  => D = slope/6

    fit_frac=0.5 significa: usar el último 50% de puntos para el ajuste.
    """
    if len(t) < 10:
        return float("nan")
    fit_frac = float(fit_frac)
    fit_frac = max(0.05, min(fit_frac, 1.0))
    i0 = int((1.0 - fit_frac) * len(t))
    i0 = max(0, min(i0, len(t) - 5))
    x = t[i0:]
    y = msd[i0:]
    A = np.vstack([x, np.ones_like(x)]).T
    slope, _ = np.linalg.lstsq(A, y, rcond=None)[0]
    return float(slope / 6.0)


def diffusion_from_vacf_integral(t: np.ndarray, vacf: np.ndarray, tmax_frac: float = 1.0) -> float:
    """
    Green–Kubo (3D):
      D = (1/3) ∫_0^∞ < v(0)·v(t) > dt

    En datos finitos integramos hasta un tiempo máximo disponible.
    tmax_frac permite recortar la cola ruidosa:
      - tmax_frac=1.0: integra todo (como el script original)
      - tmax_frac=0.8: integra hasta el 80% del tiempo total

    Integramos con trapecios.
    """
    if len(t) < 2:
        return float("nan")
    tmax_frac = float(tmax_frac)
    tmax_frac = max(0.05, min(tmax_frac, 1.0))
    imax = int(math.floor(tmax_frac * (len(t) - 1)))
    imax = max(1, min(imax, len(t) - 1))
    return float((1.0 / 3.0) * np.trapz(vacf[: imax + 1], t[: imax + 1]))


def compute_rdf(
    rx: np.ndarray, ry: np.ndarray, rz: np.ndarray, L: float,
    nbins: int = 200, rmax: Optional[float] = None
) -> Tuple[np.ndarray, np.ndarray]:
    """
    RDF g(r) computed from wrapped positions using MIC distances.
    Uses snapshots and all unique pairs i<j (O(n^2) per snapshot).
    Includes automatic subsampling for performance.

    Returns: r_centers, g_r
    """
    n_snap, n = rx.shape
    if rmax is None:
        rmax = 0.5 * L

    dr = rmax / nbins
    hist = np.zeros(nbins, dtype=np.float64)

    snap_idx = np.arange(n_snap)
    if n_snap > 200:
        snap_idx = np.linspace(0, n_snap - 1, 200, dtype=int)

    part_idx = np.arange(n)
    if n > 800:
        part_idx = np.linspace(0, n - 1, 800, dtype=int)

    n_eff = len(part_idx)
    if n_eff < 2:
        raise ValueError("Not enough particles for RDF after subsampling.")

    vol = L ** 3
    rho = n_eff / vol

    for s in snap_idx:
        x = rx[s, part_idx]
        y = ry[s, part_idx]
        z = rz[s, part_idx]
        for i in range(n_eff - 1):
            dx = x[i + 1:] - x[i]
            dy = y[i + 1:] - y[i]
            dz = z[i + 1:] - z[i]
            dx -= L * np.rint(dx / L)
            dy -= L * np.rint(dy / L)
            dz -= L * np.rint(dz / L)
            r = np.sqrt(dx*dx + dy*dy + dz*dz)
            valid = r < rmax
            bins = (r[valid] / dr).astype(int)
            np.add.at(hist, bins, 2.0)  # i-j and j-i

    r_edges = np.linspace(0.0, rmax, nbins + 1)
    r_centers = 0.5 * (r_edges[:-1] + r_edges[1:])
    shell_vol = (4.0 / 3.0) * math.pi * (r_edges[1:] ** 3 - r_edges[:-1] ** 3)

    n_samples = len(snap_idx)
    norm = n_samples * n_eff * rho * shell_vol
    g = np.zeros_like(r_centers)
    mask = norm > 0
    g[mask] = hist[mask] / norm[mask]
    return r_centers, g


# -----------------------------
# Plotting (matplotlib)
# -----------------------------
def savefig(path: Path) -> None:
    plt.tight_layout()
    plt.savefig(path, dpi=180)
    plt.close()


def plot_timeseries(t: np.ndarray, data: Dict[str, np.ndarray], outdir: Path) -> None:
    for k, y in data.items():
        plt.figure()
        plt.plot(t, y)
        plt.xlabel("time")
        plt.ylabel(k)
        plt.title(f"Time series: {k}")
        savefig(outdir / f"timeseries_{k}.png")

    plt.figure()
    if "epot" in data:
        plt.plot(t, data["epot"], label="epot")
    if "ekin" in data:
        plt.plot(t, data["ekin"], label="ekin")
    if "etot" in data:
        plt.plot(t, data["etot"], label="etot")
    plt.xlabel("time")
    plt.ylabel("energy")
    plt.title("Energies vs time")
    plt.legend()
    savefig(outdir / "timeseries_energies.png")


def rolling_mean(x: np.ndarray, w: int) -> np.ndarray:
    if w <= 1:
        return x.copy()
    w = min(w, len(x))
    kernel = np.ones(w) / w
    return np.convolve(x, kernel, mode="same")


def plot_rolling(t: np.ndarray, data: Dict[str, np.ndarray], outdir: Path, window_frac: float = 0.02) -> None:
    n = len(t)
    w = max(5, int(window_frac * n))
    for k, y in data.items():
        plt.figure()
        plt.plot(t, y, alpha=0.35, label="raw")
        plt.plot(t, rolling_mean(y, w), label=f"rolling mean (w={w})")
        plt.xlabel("time")
        plt.ylabel(k)
        plt.title(f"Rolling mean: {k}")
        plt.legend()
        savefig(outdir / f"rolling_{k}.png")


def plot_histograms(data: Dict[str, np.ndarray], outdir: Path) -> None:
    for k in ["T", "P", "etot"]:
        if k not in data:
            continue
        plt.figure()
        plt.hist(data[k], bins=60)
        plt.xlabel(k)
        plt.ylabel("count")
        plt.title(f"Histogram: {k}")
        savefig(outdir / f"hist_{k}.png")


def plot_scatter(data: Dict[str, np.ndarray], outdir: Path) -> None:
    def scatter(a: str, b: str) -> None:
        if a not in data or b not in data:
            return
        plt.figure()
        plt.scatter(data[a], data[b], s=8, alpha=0.5)
        plt.xlabel(a)
        plt.ylabel(b)
        plt.title(f"Scatter: {b} vs {a}")
        savefig(outdir / f"scatter_{a}_vs_{b}.png")

    scatter("T", "P")
    scatter("epot", "ekin")
    scatter("T", "etot")


def plot_corr(
    corr: Dict[str, CorrInfo],
    corrmean: Dict[str, CorrInfo],
    dt_sample: float,
    outdir: Path
) -> Dict[str, Dict[str, float]]:
    """
    Plots normalized ACF and returns tau estimates per observable.
    """
    stats: Dict[str, Dict[str, float]] = {}

    if corr:
        plt.figure()
        for obs, ci in corr.items():
            plt.plot(ci.lag * dt_sample, ci.cn, label=obs)
        plt.xlabel("tau")
        plt.ylabel("C_norm(tau)")
        plt.title("Normalized autocorrelation (full series)")
        plt.legend()
        savefig(outdir / "acf_fullseries_overlay.png")

    if corrmean:
        plt.figure()
        for obs, ci in corrmean.items():
            plt.plot(ci.lag * dt_sample, ci.cn, label=obs)
        plt.xlabel("tau")
        plt.ylabel("<C_norm(tau)>_blocks")
        plt.title("Normalized autocorrelation (block-averaged)")
        plt.legend()
        savefig(outdir / "acf_blockmean_overlay.png")

    for obs in sorted(set(list(corr.keys()) + list(corrmean.keys()))):
        plt.figure()
        if obs in corr:
            plt.plot(corr[obs].lag * dt_sample, corr[obs].cn, label="full series")
        if obs in corrmean:
            plt.plot(corrmean[obs].lag * dt_sample, corrmean[obs].cn, label="block mean")
        plt.xlabel("tau")
        plt.ylabel("C_norm(tau)")
        plt.title(f"Normalized autocorrelation: {obs}")
        plt.legend()
        savefig(outdir / f"acf_{obs}.png")

        source = corr.get(obs) or corrmean.get(obs)
        if source is not None:
            tau_int = estimate_tau_int(source.cn, dt_sample, max_lag=int(source.lag[-1]))
            tau_1e = estimate_tau_1e(source.cn, dt_sample)
            stats[obs] = {"tau_int": float(tau_int), "tau_1e": float(tau_1e)}

    return stats


def plot_msd_vacf_rdf(
    rva: RVAData,
    outdir: Path,
    msd_fit_frac: float = 0.5,
    vacf_tmax_frac: float = 1.0,
    rva_max_lag: Optional[int] = None,
    origin_stride: Optional[int] = None,
) -> Dict[str, float]:
    """
    Produce:
      - msd_tau.png  (time-origin averaged)
      - vacf_tau.png (time-origin averaged)
      - rdf.png
    and return diffusion estimates.

    Parametrización (con defaults “como antes”):
      - msd_fit_frac=0.5: ajusta la cola final del MSD (Einstein) usando el último 50%
      - vacf_tmax_frac=1.0: integra toda la VACF disponible (Green–Kubo) (como antes)
      - rva_max_lag=None: usa el mismo heurístico que antes (min(n-1,3000))
      - origin_stride=None: usa el mismo heurístico que antes (según longitud)
    """
    dt_sample = rva.dt * rva.output_interval

    # --- max_lag: por defecto como antes
    if rva_max_lag is None:
        max_lag = min(rva.n_snapshots - 1, 3000)
    else:
        max_lag = int(max(1, min(rva_max_lag, rva.n_snapshots - 1)))

    # --- origin_stride: por defecto como antes
    if origin_stride is None:
        origin_stride_used = 1
        if rva.n_snapshots > 5000:
            origin_stride_used = 5
        if rva.n_snapshots > 20000:
            origin_stride_used = 20
    else:
        origin_stride_used = int(max(1, origin_stride))

    msd_tau = compute_msd_tau_timeorig(
        rva.rux, rva.ruy, rva.ruz,
        max_lag=max_lag,
        origin_stride=origin_stride_used
    )
    vacf_tau = compute_vacf_tau_timeorig(
        rva.vx, rva.vy, rva.vz,
        max_lag=max_lag,
        origin_stride=origin_stride_used
    )

    tau = np.arange(len(msd_tau)) * dt_sample

    # MSD(τ)
    plt.figure()
    plt.plot(tau, msd_tau)
    plt.xlabel("Tau [reduced time]")
    plt.ylabel("MSD(τ)")
    plt.title("Mean squared displacement (time-origin averaged)")
    savefig(outdir / "msd_tau.png")

    # VACF(τ)
    plt.figure()
    plt.plot(tau, vacf_tau)
    plt.xlabel("Tau [reduced time]")
    plt.ylabel("VACF(τ)")
    plt.title("Velocity autocorrelation (time-origin averaged)")
    savefig(outdir / "vacf_tau.png")

    # --- Difusión ---
    # Einstein: MSD ~ 6 D t (régimen difusivo)
    D_msd = diffusion_from_msd_linear_fit(tau, msd_tau, fit_frac=msd_fit_frac)

    # Green–Kubo: D = (1/3) ∫ VACF dt (con posible recorte de cola)
    D_vacf = diffusion_from_vacf_integral(tau, vacf_tau, tmax_frac=vacf_tmax_frac)

    # RDF
    try:
        r, g = compute_rdf(rva.rx, rva.ry, rva.rz, rva.L, nbins=200, rmax=0.5 * rva.L)
        plt.figure()
        plt.plot(r, g)
        plt.xlabel("r [reduced length]")
        plt.ylabel("g(r)")
        plt.title("Radial distribution function (RDF)")
        savefig(outdir / "rdf.png")
    except Exception as e:
        (outdir / "rdf_error.txt").write_text(str(e), encoding="utf-8")

    return {
        "dt_sample_rva": float(dt_sample),
        "max_lag_used": int(len(msd_tau) - 1),
        "origin_stride_used": int(origin_stride_used),
        "msd_fit_frac_used": float(msd_fit_frac),
        "vacf_tmax_frac_used": float(vacf_tmax_frac),
        "D_from_MSD_tau_fit": float(D_msd),
        "D_from_VACF_tau_int": float(D_vacf),
    }


# -----------------------------
# Root autodetect
# -----------------------------
def default_project_root() -> Path:
    """
    If this script is in <root>/scripts/md_one_run_analysis.py,
    return <root>. Otherwise fallback to cwd.
    """
    here = Path(__file__).resolve()
    parent = here.parent  # scripts/
    cand = parent.parent  # root/
    if (cand / "outputs").exists() or (cand / "inputs").exists():
        return cand
    return Path.cwd().resolve()


# -----------------------------
# Main analysis pipeline
# -----------------------------
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--root",
        type=str,
        default=None,
        help="Project root (contains outputs/ and inputs/). Default: auto-detect from script location.",
    )
    ap.add_argument(
        "--one-run-subdir",
        type=str,
        default="outputs/one_run",
        help="Path from root to one_run outputs",
    )
    ap.add_argument("--skip-rva", action="store_true", help="Skip reading rva.dat and RDF/VACF/MSD")

    # --- NEW (opcionales): si no los pones, se usa el comportamiento original ---
    ap.add_argument(
        "--msd-fit-frac",
        type=float,
        default=0.5,
        help="Fraction of the MSD tail used for linear fit (Einstein). Default: 0.5 (same as before).",
    )
    ap.add_argument(
        "--vacf-tmax-frac",
        type=float,
        default=1.0,
        help="Fraction of total time used to integrate VACF (Green–Kubo). Default: 1.0 (use all, same as before).",
    )
    ap.add_argument(
        "--rva-max-lag",
        type=int,
        default=None,
        help="Max lag (in snapshots) for MSD/VACF. Default: auto (same heuristic as before).",
    )
    ap.add_argument(
        "--origin-stride",
        type=int,
        default=None,
        help="Stride between time origins for MSD/VACF averaging. Default: auto (same heuristic as before).",
    )

    args = ap.parse_args()

    root = Path(args.root).resolve() if args.root is not None else default_project_root()
    one_run = (root / args.one_run_subdir).resolve()
    outdir = one_run / "analysis"
    ensure_dir(outdir)

    # -----------------------------
    # 1) instantaneous_energies.dat
    # -----------------------------
    inst_path = one_run / "instantaneous_energies.dat"
    inst = read_text_table(inst_path)

    if inst.shape[1] < 6:
        raise ValueError(f"Expected >= 6 columns in {inst_path}, got {inst.shape[1]}")

    t = inst[:, 0]
    epot = inst[:, 1]
    ekin = inst[:, 2]
    etot = inst[:, 3]
    T = inst[:, 4]
    P = inst[:, 5]

    data = {"epot": epot, "ekin": ekin, "etot": etot, "T": T, "P": P}

    plot_timeseries(t, data, outdir)
    plot_rolling(t, data, outdir, window_frac=0.02)
    plot_histograms(data, outdir)
    plot_scatter(data, outdir)

    dt_sample = float(np.median(np.diff(t))) if len(t) > 1 else float("nan")

    # -----------------------------
    # 2) correlations
    # -----------------------------
    corr, corrmean = load_corr_files(one_run)
    acf_stats: Dict[str, Dict[str, float]] = {}
    if (corr or corrmean) and np.isfinite(dt_sample):
        acf_stats = plot_corr(corr, corrmean, dt_sample, outdir)

    # -----------------------------
    # 3) md_final_results.txt summary
    # -----------------------------
    final_path = one_run / "md_final_results.txt"
    final_stats = parse_md_final_results(final_path)

    # -----------------------------
    # 4) rva.dat: MSD/VACF/RDF (+ diffusion)
    # -----------------------------
    rva_stats: Dict[str, float] = {}
    rva_path = one_run / "rva.dat"
    if (not args.skip_rva) and rva_path.exists():
        rva = read_rva(rva_path)
        rva_stats = plot_msd_vacf_rdf(
            rva,
            outdir,
            msd_fit_frac=args.msd_fit_frac,
            vacf_tmax_frac=args.vacf_tmax_frac,
            rva_max_lag=args.rva_max_lag,
            origin_stride=args.origin_stride,
        )
        if not np.isfinite(dt_sample):
            dt_sample = rva.dt * rva.output_interval

    # -----------------------------
    # 5) write summary
    # -----------------------------
    summary = {
        "root": str(root),
        "one_run_dir": str(one_run),
        "n_samples_instantaneous": int(len(t)),
        "dt_sample_from_instantaneous": float(dt_sample),
        "final_results_parsed": final_stats,
        "acf_tau_estimates": acf_stats,
        "rva_stats": rva_stats,
        "files_used": {
            "instantaneous_energies": str(inst_path),
            "corr_files_present": sorted([str((one_run / f"corr_{o}.dat")) for o in corr.keys()]),
            "corrmean_files_present": sorted([str((one_run / f"corrmean_{o}.dat")) for o in corrmean.keys()]),
            "md_final_results": str(final_path) if final_path.exists() else None,
            "rva": str(rva_path) if rva_path.exists() else None,
        },
        "plots_dir": str(outdir),
        "cli_diffusion_params": {
            "msd_fit_frac": float(args.msd_fit_frac),
            "vacf_tmax_frac": float(args.vacf_tmax_frac),
            "rva_max_lag": None if args.rva_max_lag is None else int(args.rva_max_lag),
            "origin_stride": None if args.origin_stride is None else int(args.origin_stride),
        },
    }

    (outdir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    lines: List[str] = []
    lines.append("MD one_run analysis summary")
    lines.append("=" * 28)
    lines.append(f"one_run dir: {one_run}")
    lines.append(f"samples (instantaneous): {len(t)}")
    lines.append(f"dt_sample (from time series): {dt_sample}")
    lines.append("")
    if final_stats:
        lines.append("Parsed md_final_results (last block):")
        for k in sorted(final_stats.keys()):
            lines.append(f"  {k}: {final_stats[k]}")
        lines.append("")
    if acf_stats:
        lines.append("ACF tau estimates (from normalized ACF):")
        for obs, d in acf_stats.items():
            lines.append(f"  {obs}: tau_int={d.get('tau_int')}, tau_1e={d.get('tau_1e')}")
        lines.append("")
    if rva_stats:
        lines.append("Trajectory-derived stats (time-origin MSD/VACF + diffusion):")
        # D por dos vías físicas:
        # - Einstein (MSD): robusto si hay tramo lineal claro a tiempos largos
        # - Green–Kubo (VACF): sensible al ruido de la cola si integras demasiado
        for k, v in rva_stats.items():
            lines.append(f"  {k}: {v}")
        lines.append("")
    lines.append(f"Plots saved to: {outdir}")

    (outdir / "summary.txt").write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
import os
import re
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ---------- helpers ----------
def read_params(hpp_path="src/parameters.hpp"):
    if not os.path.exists(hpp_path):
        return {}
    txt = open(hpp_path, "r").read()

    def grab(name, cast=float):
        m = re.search(rf"{name}\s*=\s*([0-9.eE+-]+)", txt)
        return cast(m.group(1)) if m else None

    return {
        "N": grab("N", int),
        "a": grab("a", float),
        "eta": grab("eta", float),  # may fail if eta is now 'inline double eta'
    }


def semiclassical_density_1loop(eta):
    """
    One-loop (Gaussian) instanton density used as reference line.
    This uses S0 = 4 eta^3 / 3 and the prefactor from Schäfer's lecture.
    """
    S0 = 4.0 * eta**3 / 3.0
    return 8.0 * eta**2.5 * np.sqrt(2.0 / np.pi) * np.exp(-S0)


def semiclassical_density_2loop(eta):
    """
    Two-loop correction: multiply 1-loop result by exp(-71/(72 S0)).
    """
    S0 = 4.0 * eta**3 / 3.0
    return semiclassical_density_1loop(eta) * np.exp(-(71.0 / 72.0) / S0)


def S0_classical(eta):
    return 4.0 * eta**3 / 3.0


def load_density_csv(path):
    """
    Expect one of:
      - new format: columns [n_cool, n_inst, density]
      - legacy format: columns [n_cool, n_inst, action] (no density)
    Returns dict with keys:
      n_cool, density (or None), n_inst, action (or None)
    """
    df = pd.read_csv(path)
    cols = {c.lower(): c for c in df.columns}  # case-insensitive map
    need = ["n_cool", "n_inst"]
    for k in need:
        if k not in cols:
            raise ValueError(f"{path} missing column '{k}'")
    n_cool = df[cols["n_cool"]].to_numpy()
    n_inst = df[cols["n_inst"]].to_numpy()
    density = df[cols["density"]].to_numpy() if "density" in cols else None
    action = df[cols["action"]].to_numpy() if "action" in cols else None
    return {"n_cool": n_cool, "n_inst": n_inst, "density": density, "action": action}


def label_from_file(path, params):
    eta = params.get("eta", None)
    base = os.path.basename(path)
    stem = os.path.splitext(base)[0]
    if eta is not None:
        return rf"{stem}  ($\eta={eta}$)"
    return stem


# ---------- NEW: ensemble Fig.7 mode ----------
if os.path.exists("data/fig7_long.csv"):
    df = pd.read_csv("data/fig7_long.csv")

    need_cols = {"eta", "conf", "n_cool", "density", "s_per_inst"}
    missing = need_cols.difference(set(df.columns))
    if missing:
        raise ValueError(f"data/fig7_long.csv missing columns: {sorted(missing)}")

    # group statistics
    g = df.groupby(["eta", "n_cool"], sort=True)

    dens_mean = g["density"].mean()
    dens_err = g["density"].std(ddof=1) / np.sqrt(g.size())

    s_mean = g["s_per_inst"].mean()
    s_err = g["s_per_inst"].std(ddof=1) / np.sqrt(g["s_per_inst"].count())

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(12, 5), sharex=False)

    for eta in sorted(df["eta"].unique()):
        # (a) density with error bars
        nc = dens_mean.loc[eta].index.to_numpy()
        y = dens_mean.loc[eta].to_numpy()
        yerr = dens_err.loc[eta].to_numpy()

        axA.errorbar(
            nc, y, yerr=yerr,
            marker="o", ms=3, linestyle="-", capsize=2,
            label=rf"$\eta={eta}$"
        )

        # reference lines (one- and two-loop)
        axA.axhline(semiclassical_density_1loop(eta), linestyle="-", color="#000000", linewidth=1)
        axA.axhline(semiclassical_density_2loop(eta), linestyle="--", color="#000000", linewidth=1)

        # (b) S/Ninst with error bars (NaNs already ignored by mean/count)
        nc2 = s_mean.loc[eta].index.to_numpy()
        y2 = s_mean.loc[eta].to_numpy()
        y2err = s_err.loc[eta].to_numpy()

        axB.errorbar(
            nc2, y2, yerr=y2err,
            marker="s", ms=3, linestyle="-", capsize=2,
            label=rf"$\eta={eta}$"
        )

        # classical line 
        axB.axhline(S0_classical(eta), linestyle="-", color="#000000", linewidth=1)

    axA.set_xscale("log")
    axB.set_xscale("log")

    axA.set_xlabel(r"$n_{\mathrm{cool}}$")
    axA.set_ylabel(r"instanton density  $n_{I{+}A} = N_{I{+}A}/\beta$")
    axA.set_title("Fig.7(a): Instanton density vs cooling sweeps")
    axA.grid(True, alpha=0.3)
    axA.legend()

    axB.set_xlabel(r"$n_{\mathrm{cool}}$")
    axB.set_ylabel(r"$S / N_{\mathrm{inst}}$")
    axB.set_title("Fig.7(b): Action per instanton vs cooling sweeps")
    axB.grid(True, alpha=0.3)
    axB.legend()

    plt.tight_layout()
    plt.show()
    raise SystemExit


# ---------- LEGACY: single-run files mode ----------
# ---------- gather files ----------
files = sorted(glob.glob("data/instanton_density_vs_ncool*.csv"))
files += sorted(glob.glob("data/instanton_density_*.csv"))

if not files:
    raise FileNotFoundError(
        "No instanton density CSVs found. Expected e.g. "
        "'data/fig7_long.csv' (new) or 'data/instanton_density_vs_ncool.csv' (legacy)."
    )

params = read_params()
beta = None
if params.get("N") is not None and params.get("a") is not None:
    beta = params["N"] * params["a"]

series = []
for f in files:
    try:
        s = load_density_csv(f)
        s["file"] = f
        series.append(s)
    except Exception as e:
        print(f"[skip] {f}: {e}")

if not series:
    raise RuntimeError("No valid CSV could be parsed for figure 7.")

fig, (axA, axB) = plt.subplots(1, 2, figsize=(12, 5), sharex=False)

# (a) density vs n_cool
any_density = False
for s in series:
    label = label_from_file(s["file"], params)
    n_cool = s["n_cool"]
    dens = s["density"]
    if dens is None:
        if beta is not None and np.all(np.isfinite(s["n_inst"])):
            dens = s["n_inst"] / beta
        else:
            print(f"[warn] {s['file']} has no 'density' and beta unknown; skipping.")
            continue
    axA.plot(n_cool, dens, marker="o", ms=3, label=label)
    any_density = True

axA.set_xlabel(r"$n_{\mathrm{cool}}$")
axA.set_ylabel(r"instanton density  $n_{I{+}A} = N_{I{+}A}/\beta$")
axA.set_title("Instanton density vs cooling sweeps (legacy)")
axA.grid(True, alpha=0.3)
if any_density:
    axA.legend()

# (b) action per instanton vs n_cool (legacy: needs 'action')
any_action = False
for s in series:
    if s["action"] is None:
        continue
    n_cool = s["n_cool"]
    n_inst = s["n_inst"].astype(float)
    action = s["action"].astype(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        s_per_inst = np.where(n_inst > 0, action / n_inst, np.nan)
    label = label_from_file(s["file"], params)
    axB.plot(n_cool, s_per_inst, marker="s", ms=3, label=label)
    any_action = True

axB.set_xlabel(r"$n_{\mathrm{cool}}$")
axB.set_ylabel(r"$S / N_{\mathrm{inst}}$")
axB.set_title("Action per instanton vs cooling sweeps (legacy)")
axB.grid(True, alpha=0.3)
if any_action:
    axB.legend()

plt.tight_layout()
plt.show()
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
        "N":   grab("N", int),
        "a":   grab("a", float),
        "eta": grab("eta", float),
    }

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
    need = ["n_cool","n_inst"]
    for k in need:
        if k not in cols:
            raise ValueError(f"{path} missing column '{k}'")
    n_cool = df[cols["n_cool"]].to_numpy()
    n_inst = df[cols["n_inst"]].to_numpy()
    density = None
    action = None
    if "density" in cols:
        density = df[cols["density"]].to_numpy()
    if "action" in cols:
        action = df[cols["action"]].to_numpy()
    return {"n_cool": n_cool, "n_inst": n_inst, "density": density, "action": action}

def label_from_file(path, params):
    # Try to include eta in label if available, else use filename stem
    eta = params.get("eta", None)
    base = os.path.basename(path)
    stem = os.path.splitext(base)[0]
    if eta is not None:
        return rf"{stem}  ($\eta={eta}$)"
    return stem

# ---------- gather files ----------
# Preferred (new) export
files = sorted(glob.glob("data/instanton_density_vs_ncool*.csv"))
# Also accept legacy exports that include action: "instanton_density_*.csv"
files += sorted(glob.glob("data/instanton_density_*.csv"))

if not files:
    raise FileNotFoundError(
        "No instanton density CSVs found. Expected e.g. "
        "'data/instanton_density_vs_ncool.csv' or 'data/instanton_density_*.csv'.")

params = read_params()
beta = None
if params.get("N") is not None and params.get("a") is not None:
    beta = params["N"] * params["a"]

# ---------- load all series ----------
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

# ---------- plot ----------
fig, (axA, axB) = plt.subplots(1, 2, figsize=(12,5), sharex=False)

# (a) density vs n_cool
any_density = False
for s in series:
    label = label_from_file(s["file"], params)
    n_cool = s["n_cool"]
    dens = s["density"]
    if dens is None:
        # Try to compute density = n_inst / beta if beta is known
        if beta is not None and np.all(np.isfinite(s["n_inst"])):
            dens = s["n_inst"] / beta
        else:
            print(f"[warn] {s['file']} has no 'density' and beta unknown; skipping density plot for this file.")
            continue
    axA.plot(n_cool, dens, marker="o", ms=3, label=label)
    any_density = True

axA.set_xlabel(r"$n_{\mathrm{cool}}$")
axA.set_ylabel(r"instanton density  $n_{I{+}A} = N_{I{+}A}/\beta$")
titleA = "Instanton density vs cooling sweeps"
if params.get("eta") is not None:
    titleA += rf"  ($\eta={params['eta']}$)"
axA.set_title(titleA)
axA.grid(True, alpha=0.3)
if any_density:
    axA.legend()

# (b) action per instanton vs n_cool (requires 'action' column NOT IMPLEMENTED YET)
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
axB.set_title("Action per instanton vs cooling sweeps")
axB.grid(True, alpha=0.3)
if any_action:
    axB.legend()

plt.tight_layout()
plt.show()
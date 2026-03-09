import numpy as np
import matplotlib.pyplot as plt

d = np.genfromtxt("data/fig15_paths.csv", delimiter=",", names=True)

if d.size == 0:
    raise RuntimeError("fig15_paths.csv is empty")

if d.shape == ():
    d = np.array([d], dtype=d.dtype)

ratios = np.unique(d["ratio"])[::-1]

def zero_crossings_interp(tau, x):
    z = []
    n = len(x)
    for i in range(n - 1):
        xi, xj = x[i], x[i + 1]
        if xi == 0.0:
            z.append(tau[i])
        elif xi * xj < 0.0:
            t = xi / (xi - xj)
            z0 = tau[i] + t * (tau[i + 1] - tau[i])
            z.append(z0)
    return np.array(z)

beta = np.max(d["tau"]) + (d["tau"][1] - d["tau"][0])

fig, ax = plt.subplots(1, 2, figsize=(11, 4))

for r in ratios:
    m = np.isclose(d["ratio"], r)
    tau = d["tau"][m].copy()
    x = d["x"][m].copy()
    s = d["s"][m].copy()

    order0 = np.argsort(tau)
    tau = tau[order0]
    x = x[order0]
    s = s[order0]

    z = zero_crossings_interp(tau, x)

    if len(z) >= 2:
        tau0 = 0.5 * (z[0] + z[1])
    else:
        tau0 = tau[np.argmin(np.abs(x))]

    tau_shift = tau - tau0
    tau_shift = (tau_shift + 0.5 * beta) % beta - 0.5 * beta

    order = np.argsort(tau_shift)

    ax[0].plot(tau_shift[order], x[order], label=f"{r:.1f}")
    ax[1].plot(tau_shift[order], s[order], label=f"{r:.1f}")

ax[0].set_xlim(-2.0, 2.0)
ax[1].set_xlim(-2.0, 2.0)

ax[0].set_xlabel(r"$\tau$")
ax[0].set_ylabel(r"$x(\tau)$")
ax[0].set_title("Relaxed IA path")
ax[0].grid(True, alpha=0.3)
ax[0].legend(title=r"$S/S_0$")

ax[1].set_xlabel(r"$\tau$")
ax[1].set_ylabel(r"$s(\tau)$")
ax[1].set_title("Action density")
ax[1].grid(True, alpha=0.3)
ax[1].legend(title=r"$S/S_0$")

plt.tight_layout()
plt.savefig("figures/Figure_15.png", dpi=200)
plt.show()
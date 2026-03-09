import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

sum_tauIA_df = pd.read_csv("data/fig14_sum_tauIA.csv")
sum_tauZ_df = pd.read_csv("data/fig14_sum_tauZ.csv")
stream_df = pd.read_csv("data/fig14_streamline_tauZ.csv")
mc_df = pd.read_csv("data/fig14_mc_tauZ.csv")

print(mc_df[["tauZ","Sint_over_S0"]].describe())

# keep only the physically interesting window
sum_tauIA_df = sum_tauIA_df[
    (sum_tauIA_df["tauIA"] >= 0.0) &
    (sum_tauIA_df["tauIA"] <= 4.0) &
    (sum_tauIA_df["Sint_over_S0"] >= -2.05) &
    (sum_tauIA_df["Sint_over_S0"] <= 0.05)
]

sum_tauZ_df = sum_tauZ_df[
    (sum_tauZ_df["tauZ"] >= 0.0) &
    (sum_tauZ_df["tauZ"] <= 4.0) &
    (sum_tauZ_df["Sint_over_S0"] >= -2.05) &
    (sum_tauZ_df["Sint_over_S0"] <= 0.05)
]

stream_df = stream_df[
    (stream_df["tauZ"] >= 0.0) &
    (stream_df["tauZ"] <= 4.0) &
    (stream_df["Sint_over_S0"] >= -2.05) &
    (stream_df["Sint_over_S0"] <= 0.05)
]

#mc_df = mc_df[
#    (mc_df["tauZ"] >= 0.0) &
#    (mc_df["tauZ"] <= 4.0) &
#    (mc_df["Sint_over_S0"] >= -2.05) &
#    (mc_df["Sint_over_S0"] <= 0.05)
#]
print("raw MC rows =", len(mc_df))
print(mc_df.head())
print(mc_df.describe())

fig, ax = plt.subplots(figsize=(7, 5))

# sum ansatz vs tau_IA
ax.plot(
    sum_tauIA_df["tauIA"],
    sum_tauIA_df["Sint_over_S0"],
    linewidth=2,
    label="sum ansatz"
)

# sum ansatz (zcr): same data reparametrized by zero-crossing distance
ax.scatter(
    sum_tauZ_df["tauZ"],
    sum_tauZ_df["Sint_over_S0"],
    marker="^",
    s=28,
    label="sum ansatz (zcr)"
)

# streamline points
if len(stream_df) > 0:
    ax.scatter(
        stream_df["tauZ"],
        stream_df["Sint_over_S0"],
        color="tab:orange",
        s=18,
        label="streamline"
    )

# Monte Carlo points
#if len(mc_df) > 0:
    ax.scatter(
        mc_df["tauZ"],
        mc_df["Sint_over_S0"],
        color="#bc36da",
        s=14,
        label="Monte Carlo"
    )

# analytic large-separation interaction
eta = 1.4
S0 = 4.0 * eta**3 / 3.0

tau = np.linspace(0.0, 4.0, 300)
Sint = 2.0 * S0 * (1.0 - 6.0 * np.exp(-eta * tau))
Sint_over_S0 = Sint / S0 - 2.0

ax.plot(
    tau,
    Sint_over_S0,
    "-",
    color="black",
    linewidth=2,
    label=r"analytic tail $2S_0(1-6e^{-\eta\tau})$"
)

ax.set_xlim(0.0, 4.0)
ax.set_ylim(-2.5, 0.5)

ax.set_xlabel(r"$\tau$")
ax.set_ylabel(r"$(S_{IA}-2S_0)/S_0$")
ax.set_title("Fig. 14")

ax.grid(True, alpha=0.3)
ax.legend()

plt.tight_layout()
plt.savefig("figures/Figure_14.png", dpi=200)
plt.show()

# to be remove
plt.scatter(mc_df["tauZ"], mc_df["Sint_over_S0"], s=10)
plt.xlabel("tauZ")
plt.ylabel("Sint_over_S0")
plt.show()
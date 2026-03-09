import numpy as np
import matplotlib.pyplot as plt

d = np.genfromtxt("data/fig17_positions.csv", delimiter=",", names=True)

mI = d["kind"] == 1
mA = d["kind"] == -1

plt.figure(figsize=(10,6))

plt.scatter(d["conf"][mI], d["tau"][mI], s=1, label="I", alpha=0.9)
plt.scatter(d["conf"][mA], d["tau"][mA], s=1, label="AI", alpha=0.9)

plt.xlabel("configuration index")
plt.ylabel(r"$\tau$")
plt.title("Fig.17: Typical IILM instanton configuration (positions vs MC time)")

plt.grid(True, linestyle="--", alpha=0.35)

plt.legend()
plt.tight_layout()
plt.show()
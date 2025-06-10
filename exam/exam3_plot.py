import numpy as np
import matplotlib.pyplot as plt

# Load data, skipping the first column (row numbers) if present
data = np.loadtxt("exam3_results.txt")

# Extract columns
x = data[:, 0]
y1 = data[:, 1]
y2 = data[:, 2]
y3 = data[:, 3]

fig, axs = plt.subplots(3, 1, figsize=(3, 4), sharex=True)

axs[0].plot(x, y1, label="x(t)")
axs[0].set_ylabel("x(t)")
axs[0].legend()
axs[0].grid(True)

axs[1].plot(x, y2, label="x'(t)", color="orange")
axs[1].set_ylabel("x'(t)")
axs[1].legend()
axs[1].grid(True)

axs[2].plot(x, y3, label="X_F(t)-x(t)", color="green")
axs[2].set_xlabel("t")
axs[2].set_ylabel("X_F(t)-x(t)")
axs[2].legend()
axs[2].grid(True)

fig.suptitle("Exam exercise 3 Results")
# plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.show()

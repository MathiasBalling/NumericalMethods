import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("2022jun3_results.txt")
t = np.arange(data.shape[0])
y1 = data[:, 0]
y2 = data[:, 1]

plt.figure(figsize=(8, 5))
plt.plot(t, y1, label="y1(t)")
plt.plot(t, y2, label="y2(t)")
plt.xlabel("Step")
plt.ylabel("Value")
plt.title("Solutions y1(t) and y2(t) from 2022jun3_results.txt")
plt.legend()
plt.tight_layout()
plt.show()

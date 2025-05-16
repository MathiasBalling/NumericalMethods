import numpy as np
import matplotlib.pyplot as plt

# Results from mandatory_5/main.cpp
data = np.loadtxt("results.txt")
u_vals = data[:, 0]
v_vals = data[:, 1]


# Compute the 3D coordinates of the geodesic
def surface_r(u, v):
    x = u
    y = v
    z = v**4 - 2 * u**4
    return x, y, z


x_geo, y_geo, z_geo = surface_r(u_vals, v_vals)

# Create the correct surface
U, V = np.meshgrid(np.linspace(-1, 1, 100), np.linspace(-1, 1, 100))
X = U
Y = V
Z = V**4 - 2 * U**4

# Find midpoint
u_mid = (u_vals[0] + u_vals[-1]) / 2
# Find the closest point in u_vals to u_mid
mid_idx = np.argmin(np.abs(u_vals - u_mid))
x_mid, y_mid, z_mid = surface_r(u_vals[mid_idx], v_vals[mid_idx])

# Plotting
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection="3d")

# Plot surface
ax.plot_surface(X, Y, Z, alpha=0.6, cmap="viridis", edgecolor="none")

# Overlay geodesic curve
ax.plot3D(
    x_geo, y_geo, z_geo, color="red", linewidth=2, label="From mandatory_5/main.cpp"
)

# Mark midpoint
ax.scatter(x_mid, y_mid, z_mid, color="blue", s=50, label="Midpoint V((uA + uB)/2)")

# Labels and legend
ax.set_xlabel("u")
ax.set_ylabel("v")
ax.set_zlabel("z")
ax.set_title("Geodesic Curve on the Parametric Surface")
ax.legend()

plt.tight_layout()
plt.show()

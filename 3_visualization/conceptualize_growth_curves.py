import numpy as np
import matplotlib.pyplot as plt

def growth_model_1(t, B0, nearest_mature, k, theta):
    return B0 + (nearest_mature - B0) * (1 - np.exp(-k * t))**theta

def growth_model_2(t, B0, nearest_mature, k, k0, theta):
    return B0 + (nearest_mature - B0) * (1 - np.exp(-(k0 + k * t)))**theta

def growth_model_3(t, B0, nearest_mature, k, k0, theta):
    return B0 + (nearest_mature - B0) * (1 - np.exp(-t * (k0 + k)))**theta

def growth_model_4(t, B0, nearest_mature, k, theta):
    return B0 + nearest_mature * (1 - np.exp(-k * t))**theta

# Set up parameters
B0 = 40
nearest_mature = 100
k = 0.1
theta = 3

# Create a figure with three subplots side by side
fig, (ax1, ax4, ax2, ax3) = plt.subplots(1, 4, figsize = (24, 6))

# Plot 1: Original Growth Model
t = np.linspace(0, 50, 500)
y_pred1 = growth_model_1(t, B0, nearest_mature, k, theta)

ax1.plot(t, y_pred1, 'b-', linewidth = 2)
ax1.set_title("Original Growth Model")
ax1.set_xlabel("Time")
ax1.set_ylabel("Size")
ax1.grid(True)

# Plot 2: Growth Model with k0 added to k*t
k_values = [0, 0.5, 1]
colors = ['r', 'g', 'b']

for k0, color in zip(k_values, colors):
    y2 = growth_model_2(t, B0, nearest_mature, k, k0, theta)
    ax2.plot(t, y2, color = color, label = f'k0 = {k0}')

ax2.set_xlabel('Time')
ax2.set_ylabel('Size')
ax2.set_title('Growth Model with k0 added to k*t')
ax2.legend()
ax2.grid(True)

# Plot 3: New Growth Model with t*(k0 + k)

for k0, color in zip(k_values, colors):
    y3 = growth_model_3(t, B0, nearest_mature, k, k0, theta)
    ax3.plot(t, y3, color = color, label = f'k0 = {k0}')

ax3.set_xlabel('Time')
ax3.set_ylabel('Size')
ax3.set_title('New Growth Model with t*(k0 + k)')
ax3.legend()
ax3.grid(True)

# Plot 4: Growth Model with B0 removed
y_pred4 = growth_model_4(t, B0, nearest_mature, k, theta)

ax4.plot(t, y_pred4, 'b-', linewidth = 2)
ax4.set_xlabel('Time')
ax4.set_ylabel('Size')
ax4.set_title('Growth Model with B0 removed from asymptote')
ax4.legend()
ax4.grid(True)

# Adjust layout and display
plt.tight_layout()
plt.show()
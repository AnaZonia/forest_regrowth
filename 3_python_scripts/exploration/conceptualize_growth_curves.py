import numpy as np
import matplotlib.pyplot as plt

def growth_model_1(t, B0, nearest_mature, k, theta):
    return B0 + (nearest_mature - B0) * (1 - np.exp(-k * t))**theta

def growth_model_2(t, B0, nearest_mature, k, k0, theta):
    return B0 + (nearest_mature - B0) * (1 - np.exp(-(k0 + k * t)))**theta

def growth_model_3(t, B0, nearest_mature, k, k0, theta):
    return B0 + (nearest_mature - B0) * (1 - np.exp(-t * (k0 + k)))**theta

# Set up parameters
B0 = 40
nearest_mature = 100
k = 0.1
theta = 3

# Create a figure with three subplots side by side
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize = (24, 6))

# Plot 1: Original Growth Model
t1 = np.linspace(0, 50, 500)
y_pred1 = growth_model_1(t1, B0, nearest_mature, k, theta)

ax1.plot(t1, y_pred1, 'b-', linewidth = 2)
ax1.set_title("Original Growth Model")
ax1.set_xlabel("Time")
ax1.set_ylabel("Size")
ax1.grid(True)
ax1.axhline(y = nearest_mature, color = 'r', linestyle = '--', label = 'Maximum size')
ax1.axhline(y = B0, color = 'g', linestyle = '--', label = 'Initial size')
ax1.legend()

# Plot 2: Growth Model with k0 added to k*t
t2 = np.linspace(0, 50, 500)
k_values = [0, 0.5, 1]
colors = ['r', 'g', 'b']

for k0, color in zip(k_values, colors):
    y2 = growth_model_2(t2, B0, nearest_mature, k, k0, theta)
    ax2.plot(t2, y2, color = color, label = f'k0 = {k0}')

ax2.set_xlabel('Time')
ax2.set_ylabel('Size')
ax2.set_title('Growth Model with k0 added to k*t')
ax2.legend()
ax2.grid(True)

# Plot 3: New Growth Model with t*(k0 + k)
t3 = np.linspace(0, 50, 500)

for k0, color in zip(k_values, colors):
    y3 = growth_model_3(t3, B0, nearest_mature, k, k0, theta)
    ax3.plot(t3, y3, color = color, label = f'k0 = {k0}')

ax3.set_xlabel('Time')
ax3.set_ylabel('Size')
ax3.set_title('New Growth Model with t*(k0 + k)')
ax3.legend()
ax3.grid(True)

# Adjust layout and display
plt.tight_layout()
plt.show()

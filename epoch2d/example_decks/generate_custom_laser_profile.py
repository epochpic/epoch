# Generates amplitude.dat and phase.dat for custom_laser_profile.deck.
# Run this once, with both output files placed in the same data directory
# as the deck, before starting epoch2d.
import numpy as np

N_T = 50
N_Y = 30
Y_MIN, Y_MAX = -10e-6, 10e-6
T_START, T_END = 0.0, 30e-15

y = np.linspace(Y_MIN, Y_MAX, N_Y)
t = np.linspace(T_START, T_END, N_T)

# meshgrid with indexing="ij" gives arr[t_index, y_index]
Tg, Yg = np.meshgrid(t, y, indexing="ij")

# --- amplitude envelope (values in [0, 1]) ---
t0, tau = 15e-15, 5e-15
w0 = 3e-6
amplitude = np.exp(-((Tg - t0) / tau) ** 2) * np.exp(-(Yg / w0) ** 2)

assert amplitude.shape == (N_T, N_Y)
amplitude.astype(np.float64).tofile("amplitude.dat")

# --- phase envelope (radians) ---
phase = np.zeros((N_T, N_Y))  # flat phase for a plane wave
phase.astype(np.float64).tofile("phase.dat")

print("Wrote amplitude.dat and phase.dat "
      f"({N_T} x {N_Y} points, {amplitude.nbytes} bytes each)")

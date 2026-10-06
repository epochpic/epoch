# Generates amplitude.dat and phase.dat for custom_laser_profile.deck.
# Run this once, with both output files placed in the same data directory
# as the deck, before starting epoch3d.
import numpy as np

N_Y, N_Z, N_T = 48, 40, 16  # n_tr1, n_tr2, n_t in the deck
Y_MIN, Y_MAX = -8e-6, 8e-6
Z_MIN, Z_MAX = -8e-6, 8e-6
T_START, T_END = 0.0, 60e-15

y = np.linspace(Y_MIN, Y_MAX, N_Y)
z = np.linspace(Z_MIN, Z_MAX, N_Z)
t = np.linspace(T_START, T_END, N_T)

# indexing="ij" gives arr[t_index, z_index, y_index] = (n_t, n_tr2, n_tr1)
Tg, Zg, Yg = np.meshgrid(t, z, y, indexing="ij")

# --- amplitude envelope (values in [0, 1]) ---
t0, tau = 30e-15, 10e-15
w0 = 3e-6
amplitude = (np.exp(-((Tg - t0) / tau) ** 2)
             * np.exp(-(Yg ** 2 + Zg ** 2) / w0 ** 2))

assert amplitude.shape == (N_T, N_Z, N_Y)
amplitude.astype(np.float64).tofile("amplitude.dat")

# --- phase envelope (radians) ---
phase = np.zeros((N_T, N_Z, N_Y))  # flat phase for a plane wave
phase.astype(np.float64).tofile("phase.dat")

print("Wrote amplitude.dat and phase.dat "
      f"({N_T} x {N_Z} x {N_Y} points, {amplitude.nbytes} bytes each)")

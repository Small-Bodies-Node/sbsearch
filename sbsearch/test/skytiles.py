# Licensed with the 3-clause BSD license.  See LICENSE for details.
import numpy as np

# tile 1/2 the sky
N_tiles = 10
ra_steps = np.linspace(0, 2 * np.pi, N_tiles + 1)
dec_steps = np.linspace(-np.pi / 2, np.pi / 2, N_tiles + 1)
sky_tiles = np.zeros((N_tiles**2, 10))
for i in range(N_tiles):
    for j in range(N_tiles):
        sky_tiles[i * N_tiles + j] = (
            np.mean(ra_steps[i:i+2]),
            np.mean(dec_steps[j:j+2]),
            ra_steps[i], dec_steps[j],
            ra_steps[i], dec_steps[j+1],
            ra_steps[i+1], dec_steps[j+1],
            ra_steps[i+1], dec_steps[j])
del ra_steps, dec_steps, i, j

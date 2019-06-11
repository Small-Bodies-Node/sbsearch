# Licensed with the 3-clause BSD license.  See LICENSE for details.
import numpy as np

N_tiles = 10
ra_steps = np.linspace(0, 360, N_tiles + 1)
dec_steps = np.linspace(-90, 90, N_tiles + 1)
sky_tiles = []
for i in range(N_tiles):
    for j in range(N_tiles):
        tile = 'SRID=40001;POLYGON(({} {},{} {},{} {},{} {},{} {}))'.format(
            ra_steps[i], dec_steps[j],
            ra_steps[i], dec_steps[j+1],
            ra_steps[i+1], dec_steps[j+1],
            ra_steps[i+1], dec_steps[j],
            ra_steps[i], dec_steps[j])
        sky_tiles.append(tile)
del ra_steps, dec_steps, i, j, tile

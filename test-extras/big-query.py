# Licensed with the 3-clause BSD license.  See LICENSE for details.

# generate a 1,000,000 observation database and query it

import os
import numpy as np
from numpy.random import rand
from sbsearch import SBSearch
from sbsearch.model import Observation
from sbsearch.spatial import offset_by
from astropy.coordinates.angle_utilities import angular_separation


class Config:
    mjd_start = 59253
    survey_length = 365 * 5
    exp = 30 / 86400
    fov_size = np.radians(0.2)
    points_per_day = 1  # ephemeris points


def uniform_latitude(u):
    """picked between 0 and pi/2"""
    return np.pi / 2 - np.arccos(np.sqrt(1 - u))


def random_point(th0=0, th1=np.pi):
    u = 2 * rand() - 1
    th = np.sign(u) * uniform_latitude(np.abs(u))
    return rand() * 2 * np.pi, th


def random_obs(N):
    n = 0
    rho = Config.fov_size / np.sqrt(2)
    observations = []
    while(n < N):
        ph, th = random_point()
        fov = np.array(
            [offset_by(ph, th, np.pi * i / 2, rho)
             for i in range(4)]
        )
        mjd = rand() * Config.survey_length + Config.mjd_start
        obs = Observation(
            mjd_start=mjd,
            mjd_stop=mjd + Config.exp,
        )
        obs.set_fov(np.degrees(fov[:, 0]), np.degrees(fov[:, 1]))
        observations.append(obs)
        n += 1
        if n % 10000 == 0:
            print(n)

    return observations


def random_eph():
    # proper motion
    mu = rand()**2 * 0.040 * 24 / 206245  # rad/day
    mu *= np.sign(rand() - 0.5)
    # pa = rand() * 2 * np.pi  # radians E of N

    # starting point
    ph, th = random_point()
    mjd = (Config.mjd_start +
           np.linspace(0, Config.survey_length,
                       Config.survey_length * Config.points_per_day + 1))
    dt = mjd[1] - mjd[0]
    ra = [ph]
    dec = [th]
    for i in range(len(mjd) - 1):
        ph, th = offset_by(ra[-1], dec[-1], np.pi / 2, mu * dt)
        ra.append(ph)
        dec.append(th)

    return ra, dec, mjd


db_existed = os.path.exists('test.db')
sbs = SBSearch(1e-3, 'sqlite:///./temp/test.db')
if not db_existed:
    for i in range(10):
        sbs.add_observations(random_obs(100000))

# query N random "ephemerides"
speed_test = True
if speed_test:
    for i in range(10):
        ra, dec, mjd = random_eph()
        obs = sbs.find_observations_intersecting_line(ra, dec, mjd)
        print(len(obs), 'found')
else:
    # accuracy test
    obs = []
    # get best matches
    while len(obs) == 0:
        ra, dec, mjd = random_eph()
        N = len(ra)
        obs = sbs.find_observations_intersecting_line(
            ra[:N], dec[:N], mjd[:N])
    print(len(obs), 'found using time')
    # get all matches to line
    obs1 = sorted(
        sbs.find_observations_intersecting_line(ra[:N], dec[:N]),
        key=lambda o: o.mjd_start
    )
    print(len(obs1), 'found without')

    obs2 = [o for o in obs1
            if (o.mjd_start <= mjd[N - 1]) and (o.mjd_stop >= mjd[0])
            ]
    print(len(obs2), 'of those within full time limits:', mjd[0], mjd[N - 1])

    dt = []
    dr = []
    for i in range(N):
        t = mjd[i:i+2].mean()
        _dt = np.abs([(o.mjd_start + o.mjd_stop) / 2 - t for o in obs2])
        # pick the closest in time, and all within exp time
        j = _dt.argmin()
        dt.append(_dt[j])
        _ra = np.mean(ra[i:i+2])
        _dec = np.mean(dec[i:i+2])
        cen = np.mean(np.radians(np.array(
            [c.split(':') for c in obs2[j].fov.split(',')],
            float
        )), 0)
        dr.append(angular_separation(_ra, _dec, cen[0], cen[1]))

    import matplotlib.pyplot as plt
    plt.clf()
    ax = plt.gca()
    ax.scatter(np.abs(dt), dr, marker='.')
    plt.setp(ax, yscale='log', xscale='log', xlabel='dt', ylabel='dr')
    plt.axhline(Config.fov_size / np.sqrt(2))
    plt.axvline(Config.exp)

# Licensed with the 3-clause BSD license.  See LICENSE for details.

# generate a 1,000,000 observation database and query it

import argparse
import numpy as np
from numpy.random import rand
from astropy.time import Time
import astropy.units as u
from sbsearch import SBSearch
from sbsearch.model import ExampleSurvey
from sbsearch.spatial import offset_by
from astropy.coordinates.angle_utilities import angular_separation


parser = argparse.ArgumentParser()
parser.add_argument('action', choices=['add', 'reindex', 'speed', 'accuracy',
                                       'encke', 'ceres'],
                    help=('add: 1,000,000 observations to the database; '
                          'reindex: after changing indexing parameters '
                          '(must edit this file or the sbsearch source code); '
                          'speed: test 10 random 5-year searches; '
                          'accuracy: verifies observation matches; '
                          'encke: search for this comet'
                          )
                    )
args = parser.parse_args()


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
        obs = ExampleSurvey(
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
    mu = rand()**2 * 200
    print(mu, 'arcsec/hr')
    mu *= (u.arcsec / u.hr).to(u.rad / u.day)
    mu *= np.sign(rand() - 0.5)

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


try:
    sbs = SBSearch('postgresql://@/big_query_test', min_edge_length=0.005)
except Exception as e:
    raise ValueError(
        "Failed db connection, is it created and online?"
    ) from e

sbs.source = 'example_survey'

if args.action == 'add':
    sbs.db.drop_spatial_index()
    for i in range(10):
        sbs.add_observations(random_obs(100000))
    sbs.db.create_spatial_index()
elif args.action == 'reindex':
    sbs.re_index()
    sbs.__exit__()
elif args.action == 'speed':
    # query N random "ephemerides"
    for i in range(10):
        ra, dec, mjd = random_eph()
        t = Time.now()
        obs = sbs.find_observations_intersecting_line_at_time(ra, dec, mjd)
        print(len(obs), 'found in {:.1f}'.format((Time.now() - t).to('s')))
elif args.action == 'accuracy':
    # accuracy test
    obs = []
    # generate random ephemerides until at least one match is found
    # get best matches
    while len(obs) == 0:
        ra, dec, mjd = random_eph()
        N = len(ra)
        obs = sbs.find_observations_intersecting_line_at_time(
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
    plt.savefig('big-query-accuracy.pdf')
elif args.action == 'encke':
    encke = sbs.get_designation('2P', add=True)
    eph = encke.ephemeris_over_date_range(
        Time(Config.mjd_start, format='mjd'),
        Time(Config.mjd_start + Config.survey_length, format='mjd'),
    )

    t = Time.now()
    obs = sbs.find_observations_by_ephemeris(eph)
    print('{} found in {:.1f}'.format(len(obs), (Time.now() - t).to('s')))

    # sbs.padding = 1

    # t = Time.now()
    # obs = sbs.find_observations_by_ephemeris(eph)
    # print('{} found in {:1f} with {} arcmin padding'.format(
    #     len(obs), (Time.now() - t).to('s'), sbs.padding))
elif args.action == 'ceres':
    ceres = sbs.get_designation('1', add=True)
    eph = ceres.ephemeris_over_date_range(
        Time(Config.mjd_start, format='mjd'),
        Time(Config.mjd_start + Config.survey_length, format='mjd'),
    )

    t = Time.now()
    obs = sbs.find_observations_by_ephemeris(eph, approximate=True)
    print('{} found in {:.1f}'.format(len(obs), (Time.now() - t).to('s')))

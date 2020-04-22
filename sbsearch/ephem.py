# Licensed with the 3-clause BSD license.  See LICENSE for details.
from itertools import groupby
import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.table import vstack
from sbpy.data import Ephem, Orbit, Names, QueryError
from . import util

STEP_LIMIT = 300


def generate(desg, location, epochs, source='jpl', orbit=None, cache=False):
    """Generate ephemeris at specific epochs from external source.

    Parameters
    ----------
    desg : string
        Object designation.

    location : string
        Observer location.

    epochs : array-like, Time, or dict
        Compute ephemeris at these epochs.  Arrays must be
        floats (for Julian date) or else parsable by
        `~astropy.time.Time`.  Dictionaries are passed to the
        ephemeris source as is.  Set ``step=None`` for an
        adaptable time step that mitigates interpolation errors.

    source : string, optional
        Source to use: 'mpc', 'jpl', or 'oorb'.  'oorb' requires
        ``orbit`` parameter.

    orbit : `~sbpy.data.Orbit`, optional
        Orbital elements for ``source=oorb``.

    cache : bool, optional
        Use cached ephemerides; primarily for testing.

    Returns
    -------
    eph : `~sbpy.data.ephem.Ephem`
        Ephemeris.

    """

    if source not in ['mpc', 'jpl', 'oorb']:
        raise ValueError('Source must be "mpc", "jpl", or "oorb".')

    _epochs = _format_epochs(epochs)
    if isinstance(_epochs, dict):
        if epochs.get('step') is None:
            eph = _get_adaptable_steps(desg, location, _epochs,
                                       source=source, orbit=orbit,
                                       cache=cache)
    else:
        eph = _get_fixed_steps(desg, location, _epochs, source=source,
                               orbit=orbit, cache=cache)
    return eph


def generate_orbit(desg, epochs, cache=False):
    """Generate orbital parameters from JPL Horizons at specific epochs.

    Parameters
    ----------
    desg: string
        Object designation.

    epochs: array-like or dict
        Compute orbital elements at these epochs.  For arrays,
        must be floats(for Julian date) or else parsable by
        `~astropy.time.Time`.

    cache: bool, optional
        Use cached ephemerides; primarily for testing.

    Returns
    -------
    orb: `~sbpy.data.Orbit`
        Orbital elements.

    """

    _epochs = _format_epochs(epochs)
    if not isinstance(_epochs, dict):
        if len(_epochs) > STEP_LIMIT:
            orb = None
            N = np.ceil(len(_epochs) / STEP_LIMIT)
            for e in np.array_split(_epochs, N):
                _orb = generate_orbit(desg, e, cache=cache)
                if orb:
                    orb.add_rows(_orb)
                else:
                    orb = _orb
            return orb
        else:
            pass
            # ... and continue

    kwargs = dict(epochs=_epochs, cache=cache)
    if Names.asteroid_or_comet(desg) == 'comet':
        kwargs['id_type'] = 'designation'
        if desg.strip()[0] != 'A':
            cap_limit = closest_apparition_limit(_epochs)
            kwargs.update(closest_apparition=cap_limit,
                          no_fragments=True)

    orb = Orbit.from_horizons(desg, **kwargs)

    return orb


def closest_apparition_limit(epochs):
    if isinstance(epochs, dict):
        start = epochs['start'].jd
    else:
        start = epochs[0].jd
    return '<{}'.format(start)


def _format_epochs(epochs):
    if isinstance(epochs, dict):
        start, stop = util.epochs_to_time((epochs['start'], epochs['stop']))
        step = (None if epochs.get('step') is None
                else u.Quantity(epochs.get('step')))
        e = {
            'start': start,
            'stop': stop,
            'step': step
        }
    else:
        e = util.epochs_to_time(epochs)
        if len(epochs) > 1:
            d = np.diff(e.jd)
            if any(d <= 0):
                raise ValueError('Epoch dates must be increasing and unique.')

    return e


def _get_fixed_steps(desg, location, epochs, source='jpl', orbit=None,
                     cache=False):
    if not isinstance(epochs, dict):
        # list of specific dates, divide into chunks, as needed
        if len(epochs) > STEP_LIMIT:
            eph = None
            N = np.ceil(len(epochs) / STEP_LIMIT)
            for e in np.array_split(epochs, N):
                _eph = _get_fixed_steps(desg, location, list(e),
                                        source=source, cache=cache)
                if eph:
                    eph.add_rows(_eph)
                else:
                    eph = _eph
            return eph
        else:
            pass
            # and proceed...

    if source == 'mpc':
        eph = Ephem.from_mpc(desg, epochs=epochs,
                             location=location,
                             proper_motion='sky',
                             proper_motion_unit='rad/s',
                             cache=cache)

        z = np.zeros(len(eph))
        if 'Uncertainty 3sig' not in eph.table.colnames:
            eph.table.add_column(u.Quantity(z, 'arcsec'),
                                 name='SMAA_3sigma')
            eph.table.add_column(u.Quantity(z, 'arcsec'),
                                 name='SMIA_3sigma')
            eph.table.add_column(u.Quantity(z, 'rad'),
                                 name='Theta_3sigma')
        else:
            # MPC's ephemeris uncertainty is a line, rather than
            # an ellipse
            eph.table.add_column(u.Quantity(z, 'arcsec'),
                                 name='SMIA_3sigma')
            eph.table.add_column(eph['Uncertainty 3sig'],
                                 name='SMAA_3sigma')
            eph.table.add_column(eph['Unc. P.A.'],
                                 name='Theta_3sigma')

        # convert date from Time to Julian date, this helps with table
        # manipulation since Time columns cannot be appended to.
        eph['jd'] = eph['date'].jd
    elif source == 'jpl':
        # column 7 (sidereal time) isn't needed, but adding to avoid sbpy-0.2.1 crash
        kwargs = dict(
            epochs=epochs,
            location=location,
            quantities='1,3,8,9,7,19,20,23,24,27,36,37',
            cache=cache
        )
        if Names.asteroid_or_comet(desg) == 'comet':
            kwargs['id_type'] = 'designation'
            if desg.strip()[0] != 'A':
                cap_limit = closest_apparition_limit(epochs)
                kwargs.update(
                    closest_apparition=cap_limit,
                    no_fragments=True
                )

        try:
            eph = Ephem.from_horizons(desg, **kwargs)
        except QueryError:
            # Dual-listed objects should be queried without CAP/NOFRAG. If
            # this was a comet query, try again.
            if Names.asteroid_or_comet(desg) == 'comet':
                del kwargs['closest_apparition'], kwargs['no_fragments']
                eph = Ephem.from_horizons(desg, **kwargs)
            else:
                raise

        # Horizon's theta is measured N of E.  We want E of N.
        eph['Theta_3sigma'] = (np.pi / 2) * u.rad - eph['Theta_3sigma']

        # create a plain Julian date column, helps with table
        # manipulation in _get_adaptable_steps
        eph['jd'] = eph['epoch'].jd
    elif source == 'oorb':
        eph = Ephem.from_oo(orbit, epochs=epochs, location=location)
        # no uncertainties from oorb
        z = np.zeros(len(eph))
        eph.table.add_column(u.Quantity(z, 'arcsec'),
                             name='SMAA_3sigma')
        eph.table.add_column(u.Quantity(z, 'arcsec'),
                             name='SMIA_3sigma')
        eph.table.add_column(u.Quantity(z, 'rad'),
                             name='Theta_3sigma')
        eph['jd'] = eph['epoch'].jd  # same note as for jpl

    return eph


def _get_adaptable_steps(desg, location, epochs, source='jpl', orbit=None,
                         cache=False):
    """Ephemerides with adaptable time step.

    Based on ZChecker analysis, Oct 2018, the interpolation error is
    ~2" / delta for 6h time step.

    """

    _epochs = epochs.copy()

    # daily ephemeris for delta > 1
    _epochs['step'] = 1 * u.day
    eph = _get_fixed_steps(desg, location, _epochs, source=source,
                           cache=cache)

    for limit, substep in ((1 * u.au, 4 * u.hr), (0.25 * u.au, 1 * u.hr)):
        updated = []
        groups = groupby(zip(eph['delta'], eph['jd']),
                         lambda e: e[0] < limit)
        for inside, epochs in groups:
            if not inside:
                continue

            jd = [e[1] for e in epochs]
            if len(jd) > 1:
                sub_epochs = _format_epochs({
                    'start': jd[0],
                    'stop': jd[-1],
                    'step': substep
                })
                update = _get_fixed_steps(desg, location, sub_epochs,
                                          source=source, cache=cache)
                updated.append(update.table)

        eph = Ephem.from_table(vstack([eph.table] + updated))

        # remove duplicates
        eph.table.sort('jd')
        d = np.diff(eph['jd'])
        duplicates = np.flatnonzero(np.isclose(d, 0)) + 1
        eph.table.remove_rows(duplicates)

    return eph

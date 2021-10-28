# Licensed with the 3-clause BSD license.  See LICENSE for details.

from typing import Dict, Iterable, List, Optional, Union
from abc import ABC, abstractmethod
from itertools import groupby

import numpy as np
from astropy.table import Table, Row
from astropy.time import Time
import astropy.units as u
from astroquery import jplhorizons
from sbpy.data import Names
from sbpy.data.names import TargetNameParseError

from .target import MovingTarget
from .model import Ephemeris

__all__: List[str] = [
    'get_ephemeris_generator',
    'set_ephemeris_generator',
    'EphemerisGenerator',
    'Horizons'
]

Epochs = Union[Dict[str, Union[Time, u.Quantity]], Time]


class EphemerisGenerator(ABC):
    _CHUNK_SIZE = 300  # maximum number of dates to request

    @classmethod
    def target_at_dates(cls, observer: str, target: MovingTarget, dates: Time,
                        cache: bool = True) -> List[Ephemeris]:
        """Generate ephemeris for target at specific dates.


        Parameters
        ----------
        observer : str
            The observatory code.

        target : MovingTarget
            The target to observe.

        dates : Time
            Observe at these times.

        cache : bool, optional
            Cache the results.


        Returns
        -------
        eph : list of Ephemeris
            Ephemerides.

        """

        if dates.isscalar:
            return cls._ephemeris(observer, target, dates, cache=cache)

        n_dates: int = len(dates)
        n_chunks: int = int(np.ceil(n_dates / cls._CHUNK_SIZE))

        eph: List[Ephemeris] = []
        i: np.ndarray
        for i in np.array_split(np.arange(n_dates), n_chunks):
            eph.extend(
                cls._ephemeris(observer, target, dates[i], cache=cache)
            )

        # set the retrieved date
        retrieved: str = Time.now().iso
        for row in eph:
            row.retrieved = retrieved

        return eph

    @classmethod
    def target_over_date_range(cls, observer: str, target: MovingTarget,
                               start: Time, stop: Time,
                               step: Optional[u.Quantity] = None,
                               cache: bool = False) -> List[Ephemeris]:
        """Generate ephemeris for target at specific dates.

        The adaptable step sizes are based on ZChecker analysis, Oct 2018, 
        which showed an interpolation error of ~2" / delta for 6h time steps.


        Parameters
        ----------
        observer : str
            The observatory code.

        target : MovingTarget
            The target to observe.

        start, stop : Time
            Start and stop dates.

        step : Quantity, optional
            Step size.  If `None`, then use an adaptable step size.

        cache : bool, optional
            Use cached results, if possible, otherwise cache results.  For
            ephemerides generated via astroquery.

        Returns
        -------
        eph : list of Ephemeris
            Ephemerides.

        """

        epochs: Dict[str, Union[Time, u.Quantity]] = {
            'start': start,
            'stop': stop,
            'step': step
        }
        eph: List[Ephemeris]
        if step is None:
            # adaptable time steps
            # start with daily ephemeris, which is good for delta > 1
            epochs['step'] = 1 * u.day
            eph = cls._ephemeris(observer, target, epochs, cache=cache)

            def grouper(eph):
                return eph.delta < limit

            # progressively check interior to 1 au
            for limit, substep in ((1, 4 * u.hr), (0.25, 1 * u.hr)):
                updated: List[Ephemeris] = []
                groups: Iterable = groupby(eph, grouper)
                inside: bool
                epochs: Iterable
                for inside, epochs in groups:
                    if not inside:
                        continue

                    dates: List[Time] = Time(
                        [e.mjd for e in epochs], format='mjd')
                    if len(dates) > 1:
                        sub_epochs: Dict[str, Union[Time, u.Quantity]] = {
                            'start': dates[0],
                            'stop': dates[-1],
                            'step': substep
                        }
                        updated.extend(
                            cls._ephemeris(observer, target, sub_epochs,
                                           cache=cache)
                        )

                eph.extend(updated)

                # sort and remove duplicates
                if len(eph) > 1:
                    eph = list(np.unique(sorted(eph)))

                    # mjd: np.ndarray = np.array([e.mjd for e in eph])
                    # duplicate: np.ndarray = np.r_[
                    #     False, np.isclose(np.diff(mjd), 0)]
                    # eph = [e for i, e in enumerate(eph) if ~duplicate[i]]

            # finally, set the retrieved date
            retrieved: str = Time.now().iso
            for row in eph:
                row.retrieved = retrieved

            return eph

        # otherwise, return fixed steps
        eph = cls._ephemeris(observer, target, epochs, cache=cache)

        # set the retrieved date
        retrieved: str = Time.now().iso
        for row in eph:
            row.retrieved = retrieved

        return eph

    @abstractmethod
    def _ephemeris(self, *args, **kwargs):
        """Actual ephemeris code goes here."""


class Horizons(EphemerisGenerator):
    """Generate ephemerides with JPL Horizons."""

    generator_name = 'jpl'
    _QUANTITIES: str = '1,3,8,9,7,19,20,23,24,27,36,37,41'
    # _KEEP_COLUMNS: List[str] = [
    #     'datetime_jd', 'RA', 'DEC', 'RA_rate', 'DEC_rate', 'r', 'r_rate',
    #     'delta', 'V', 'alpha', 'sunTargetPA', 'velocityPA', 'SMAA_3sigma',
    #     'SMIA_3sigma', 'Theta_3sigma', 'true_anom'
    # ]  # all columns except magnitude, which varies by object type

    @classmethod
    def _none_if_masked(cls, value):
        if getattr(value, 'mask', False):
            return None
        else:
            return value

    @classmethod
    def _ephemeris(cls, observer: str, target: MovingTarget,
                   epochs: Epochs, cache: bool = True) -> Table:
        """Get ephemeris from JPL.


        Parameters
        ----------
        observer : str
            Observer location code or name.

        target : MovingTarget
            Target to observe.

        epochs: Time or dict
            Specific dates or a dictionary with 'start', 'stop', 'step'.

        cache : bool, optional
            Use cached results, if possible, otherwise cache results.


        Returns
        -------
        eph : Table

        """

        id_type: str = 'smallbody'
        closest_apparition: Union[float, bool] = False
        no_fragments: bool = False

        # if this is a comet, use the closeset apparition and do not match
        # fragments
        comet: bool = False
        try:
            comet = Names.parse_comet(target.primary_designation)
        except TargetNameParseError:
            pass

        if comet:
            # but A/ objects are asteroids
            if target.primary_designation.strip()[0] != 'A':
                id_type = 'designation'
                closest_apparition = True
                no_fragments = True

        _epochs: Union[Dict[str, str], np.ndarray]
        sort_index: Union[np.ndarray, None] = None
        unsort_index: Union[np.ndarray, None] = None
        if isinstance(epochs, Time):
            if epochs.isscalar:
                _epochs = np.array([epochs.utc.jd])
            else:
                _epochs = epochs.utc.jd

            sort_index = np.argsort(_epochs)
            unsort_index = np.empty_like(sort_index)
            unsort_index[sort_index] = np.arange(len(sort_index))

            # search for apparition closest to but no later than start date
            if closest_apparition:
                closest_apparition = '<' + str(_epochs[0])
        else:
            _epochs = {
                'start': epochs['start'].utc.iso,
                'stop': epochs['stop'].utc.iso,
                'step': '{:.0f}{}'.format(
                    epochs['step'].value,
                    str(epochs['step'].unit)[0]
                )
            }
            if closest_apparition:
                closest_apparition = '<' + _epochs['start'][:4]

        query: jplhorizons.Horizons = jplhorizons.Horizons(
            id=target.primary_designation,
            id_type=id_type,
            location=observer,
            epochs=_epochs
        )

        eph: Table
        try:
            eph = query.ephemerides(closest_apparition=closest_apparition,
                                    no_fragments=no_fragments,
                                    quantities=cls._QUANTITIES,
                                    cache=cache)
        except ValueError:
            # Dual-listed objects should be queried without CAP/NOFRAG.  If
            # this was a comet query, try again without them.
            eph = query.ephemerides(cache=cache, quantities=cls._QUANTITIES)

        # returned columns
        # eph: targetname, datetime_str, datetime_jd, M1, solar_presence, k1,
        #      flags, RA, DEC, RA_rate, DEC_rate, airmass, magextinct, Tmag,
        #      Nmag, siderealtime, r, r_rate, delta, delta_rate, elong,
        #      elongFlag, alpha, sunTargetPA, velocityPA, RA_3sigma,
        #      DEC_3sigma, SMAA_3sigma, SMIA_3sigma, Theta_3sigma, Area_3sigma
        # ele: targetname, datetime_jd, datetime_str, M1, e, k1, q, incl,
        #      Omega, w, Tp_jd, n, M, nu, a, Q, P

        if sort_index is not None:
            eph = [eph[i] for i in unsort_index]

        # convert to list of Ephemeris objects
        result: List[Ephemeris] = []
        for row in eph:
            result.append(
                Ephemeris(
                    object_id=target.object_id,
                    mjd=row['datetime_jd'] - 2400000.5,
                    rh=row['r'],
                    delta=row['delta'],
                    phase=row['alpha'],
                    drh=row['r_rate'],
                    true_anomaly=row['true_anom'],
                    ra=row['RA'],
                    dec=row['DEC'],
                    dra=row['RA_rate'],
                    ddec=row['DEC_rate'],
                    unc_a=cls._none_if_masked(row['SMAA_3sigma']),
                    unc_b=cls._none_if_masked(row['SMIA_3sigma']),
                    unc_theta=cls._none_if_masked(row['Theta_3sigma']),
                    elong=row['elong'],
                    sangle=(row['sunTargetPA'] - 180) % 360,
                    vangle=(row['velocityPA'] - 180) % 360,
                    vmag=cls._vmag(row)
                )
            )

        return result

        # eph['V'] = cls._vmag(eph)

        # eph.keep_columns(cls._KEEP_COLUMNS)

        # eph['datetime_jd'].name = 'date'
        # eph['date'] = Time(eph['date'], format='jd')
        # eph['RA'].name = 'ra'
        # eph['DEC'].name = 'dec'
        # eph['RA_rate'].name = 'dra/dt'
        # eph['DEC_rate'].name = 'ddec/dt'
        # eph['r'].name = 'rh'
        # eph['r_rate'].name = 'drh/dt'
        # eph['alpha'].name = 'phase'
        # eph['sunTargetPA'] = (eph['sunTargetPA'] - 180) % 360
        # eph['sunTargetPA'].name = 'sangle'
        # eph['velocityPA'] = (eph['velocityPA'] - 180) % 360
        # eph['velocityPA'].name = 'vangle'
        # eph['SMAA_3sigma'].name = 'unc_a'
        # eph['SMIA_3sigma'].name = 'unc_b'
        # eph['Theta_3sigma'].name = 'unc_theta'

        # return eph

    @staticmethod
    def _vmag(row: Row, missing=99):
        """Get most relevant magnitude estimate from ephemeris.

        ``Ephem`` does not support masking, so ephemerides are populated with
        zeros.  They are replaced in order of preference: Tmag, Nmag, APmag, V.


        Parameters
        ----------
        row : astropy Table Row
            Ephemeris (meta)data from Horizons.

        missing : float, optional
            Use this value for missing magnitudes.

        """
        vmag = missing
        for k in ['V', 'APmag', 'Nmag', 'Tmag']:
            if k in row.colnames:
                if row[k] != 0:
                    vmag = row[k]
        return vmag


_generator: str = 'jpl'
_generators: Dict[str, EphemerisGenerator] = {
    g.generator_name: g for g in EphemerisGenerator.__subclasses__()
}


def get_ephemeris_generator() -> EphemerisGenerator:
    """Return the currently set ephemeris generator."""
    global _generator, _generators
    return _generators[_generator]


def set_ephemeris_generator(name: str):
    """Set the ephemeris generator."""
    global _generator, _generators
    if name in _generators:
        _generator = name
    else:
        raise ValueError(f"Invalid ephemeris generator: {name}.")

# Licensed with the 3-clause BSD license.  See LICENSE for details.

import pytest

from typing import List
import numpy as np
import sqlalchemy as sa
from astropy.table import Table
from astropy.time import Time

from ..sbsearch import SBSearch
from ..model import (Ephemeris, Observation, ObservationSpatialTerm,
                     UnspecifiedSurvey)
from ..ephemeris import get_ephemeris_generator
from ..target import MovingTarget
from ..exceptions import UnknownSource


@pytest.fixture
def sbs() -> SBSearch:
    engine: sa.engine.Engine = sa.create_engine('sqlite://')
    sessionmaker: sa.orm.sessionmaker = sa.orm.sessionmaker(bind=engine)
    with SBSearch(3e-4, sessionmaker()) as sbs:
        yield sbs


@pytest.fixture
def observations() -> List[Observation]:
    observations: List[Observation] = [
        Observation(
            mjd_start=59252.1,
            mjd_stop=59252.2,
        ),
        Observation(
            mjd_start=59252.21,
            mjd_stop=59252.31,
        )
    ]
    observations[0].set_fov([1, 2, 2, 1], [3, 3, 4, 4])
    observations[1].set_fov([2, 3, 3, 2], [3, 3, 4, 4])
    return observations


class TestSBSearch:
    def test_source(self, sbs: SBSearch) -> None:
        assert sbs.source == Observation

        sbs.source = UnspecifiedSurvey
        assert sbs.source == UnspecifiedSurvey

        sbs.source = 'observation'
        assert sbs.source == Observation

        sbs.source = 'unspecified_survey'
        assert sbs.source == UnspecifiedSurvey

    def test_source_error(self, sbs: SBSearch) -> None:
        with pytest.raises(UnknownSource):
            sbs.source = 'NEAT'

        with pytest.raises(UnknownSource):
            sbs.source = ObservationSpatialTerm

    def test_sources(self, sbs: SBSearch) -> None:
        assert sbs.sources == {
            'observation': Observation,
            'unspecified_survey': UnspecifiedSurvey
        }

    def test_uncertainty_ellipse(self, sbs: SBSearch) -> None:
        sbs.padding = 1
        sbs.uncertainty_ellipse = True
        assert sbs.uncertainty_ellipse

    def test_padding(self, sbs: SBSearch) -> None:
        sbs.uncertainty_ellipse = True
        sbs.padding = 1
        assert sbs.padding == 1

    def test_get_designation(self, sbs: SBSearch) -> None:
        object_id: int = sbs.add_designation('2P').object_id
        assert sbs.get_designation('2P').object_id == object_id
        assert sbs.get_object_id(object_id).designations[0] == '2P'

    def test_add_get_ephemeris(self, sbs: SBSearch) -> None:
        target: MovingTarget = sbs.add_designation('2P')
        sbs.add_ephemeris('500', target, '2021-01-01', '2021-02-01',
                          cache=True)
        eph: Table = sbs.get_ephemeris(target, '2021-01-01', '2021-02-01')
        assert len(eph) == 32
        assert np.isclose(eph[0].mjd, 59215.0)
        assert np.isclose(eph[-1].mjd, 59246.0)

    def test_add_get_observation(self, sbs):
        obs: Observation = Observation(
            mjd_start=59252.1,
            mjd_stop=59252.2,
        )
        obs.set_fov([1, 1, 2, 2], [1, 2, 2, 1])
        sbs.add_observations([obs])

        observations = sbs.get_observations()

        observations: List[Observation] = sbs.get_observations(
            source='observation', mjd=[59252, 59253])
        assert obs.observation_id == observations[0].observation_id

    def test_add_observations_terms(self, sbs):
        obs: Observation = Observation(
            mjd_start=59252.1,
            mjd_stop=59252.2,
        )
        obs.set_fov([1, 2, 2, 1], [3, 3, 4, 4])
        sbs.add_observations([obs])

        terms = (
            sbs.db.session.query(ObservationSpatialTerm)
            .filter(ObservationSpatialTerm.observation_id
                    == obs.observation_id)
            .all()
        )
        assert all(
            [term.observation_id == obs.observation_id for term in terms]
        )
        assert (len(set([term.term for term in terms])
                    - set([b'$10195', b'10195', b'10194', b'1019', b'101c',
                           b'101', b'$10197', b'10197', b'$10199', b'10199',
                           b'1019c', b'$101b', b'101b', b'$101c1', b'101c1',
                           b'101c4', b'101d', b'$101c7', b'101c7', b'$101c9',
                           b'101c9', b'101cc', b'$101eb', b'101eb', b'101ec',
                           b'101f'])
                    ) == 0)

    def test_find_observations_intersecting_polygon(self, sbs):
        observations: List[Observation] = [
            Observation(
                mjd_start=59252.1,
                mjd_stop=59252.2,
            ),
            Observation(
                mjd_start=59252.21,
                mjd_stop=59252.31,
            )
        ]
        observations[0].set_fov([1, 2, 2, 1], [3, 3, 4, 4])
        observations[1].set_fov([2, 3, 3, 2], [3, 3, 4, 4])
        sbs.add_observations(observations)

        found: List[Observation] = sbs.find_observations_intersecting_polygon(
            np.radians([0.5, 1.5, 1.5, 0.5]),
            np.radians([2.5, 2.5, 3.5, 3.5])
        )

        assert len(found) == 1
        assert found[0].observation_id == observations[0].observation_id

    @pytest.mark.parametrize(
        "ra, dec, n, indices",
        [
            ([0, 1], [2.5, 2.5], 0, []),
            ([0.5, 1.5], [3.5, 3.5], 1, [0]),
            ([0, 3], [3, 4], 2, [0, 1])
        ]
    )
    def test_find_observations_intersecting_line(
            self, ra, dec, n, indices, sbs, observations):
        sbs.add_observations(observations)

        found: List[Observation] = sbs.find_observations_intersecting_line(
            np.radians(ra), np.radians(dec))
        assert len(found) == n

        # verify that the expected observations are being matched (if any)
        obs_ids: List[int] = [
            observations[i].observation_id for i in range(len(observations))
            if i in indices
        ]
        assert all([f.observation_id in obs_ids for f in found])

        # verify that no unexpected observations are being matched (if any)
        obs_ids: List[int] = [
            observations[i].observation_id for i in range(len(observations))
            if i not in indices
        ]
        assert not any([f.observation_id in obs_ids for f in found])

    def test_find_observations_intersecting_line_with_padding(
            self, sbs, observations):
        sbs.add_observations(observations)
        ra: np.ndarray = np.radians((0.95, 0.95))
        dec: np.ndarray = np.radians((0.0, 5.0))

        # line missed the field
        found: List[Observation] = sbs.find_observations_intersecting_line(
            ra, dec)
        assert len(found) == 0

        # add a small amount of padding, it should still miss
        padding: np.ndarray = np.radians((0.03, 0.03))

        # get approximate results
        found: List[Observation] = sbs.find_observations_intersecting_line(
            ra, dec, approximate=True)
        assert len(found) == 1

        # get detailed results
        found = sbs.find_observations_intersecting_line(
            ra, dec, a=padding, b=padding)
        assert len(found) == 0

        # add enough padding to intersect 1 field
        padding *= 2
        found = sbs.find_observations_intersecting_line(
            ra, dec, a=padding, b=padding)
        assert len(found) == 1
        assert found[0].observation_id == observations[0].observation_id

        # add time to the query, should miss
        found = sbs.find_observations_intersecting_line_at_time(
            ra, dec, mjd=[59252.3, 59253.0], a=padding, b=padding)
        assert len(found) == 0

        # adjust time, should hit
        found = sbs.find_observations_intersecting_line_at_time(
            ra, dec, mjd=[59252.1, 59253.2], a=padding, b=padding)
        assert len(found) == 1
        assert found[0].observation_id == observations[0].observation_id

        # add a ridiculous amount of padding and cover both fields
        padding *= 20
        found = sbs.find_observations_intersecting_line(
            ra, dec, a=padding, b=padding)
        assert len(found) == 2

        # parameter validation testing
        with pytest.raises(ValueError):
            sbs.find_observations_intersecting_line(ra, dec, a=padding, b=[1])

    def test_find_observations_intersecting_line_at_time(
            self, sbs, observations):
        sbs.add_observations(observations)

        # approximate check
        args = (np.radians([1, 1.9]), np.radians(
            [3.5, 3.5]), [59252.1, 59252.3])
        found = sbs.find_observations_intersecting_line_at_time(
            *args, approximate=True
        )
        assert len(found) == 2

        # detailed check
        found = sbs.find_observations_intersecting_line_at_time(*args)
        assert len(found) == 1
        assert found[0] == observations[0]

    def test_find_observations_intersecting_line_at_time_errors(self, sbs):
        with pytest.raises(ValueError):
            sbs.find_observations_intersecting_line_at_time(
                [0, 1], [1, 0, -1], [59400, 59401])

        with pytest.raises(ValueError):
            sbs.find_observations_intersecting_line_at_time(
                [0, 1], [1, 0], [59400, 59401], a=[1, 1], b=[2])

    def test_find_observations_by_ephemeris(self, sbs, observations) -> None:
        target: MovingTarget = MovingTarget('2P')

        start: Time = Time(59252, format='mjd')
        stop: Time = Time(59253, format='mjd')
        encke: List[Ephemeris] = (
            get_ephemeris_generator()
            .target_over_date_range('500@', target, start, stop, cache=True)
        )
        sbs.add_observations(observations)

        encke_image: Observation = Observation(
            mjd_start=encke[1].mjd - 10 / 86400,
            mjd_stop=encke[1].mjd + 20 / 86400
        )
        encke_image.set_fov(
            np.array([-1, 1, 1, -1]) * 0.2 + encke[1].ra,
            np.array([-1, -1, 1, 1]) * 0.2 + encke[1].dec,
        )
        sbs.add_observations([encke_image])

        encke_offset: Observation = Observation(
            mjd_start=encke[1].mjd - 10 / 86400,
            mjd_stop=encke[1].mjd + 20 / 86400
        )
        encke_offset.set_fov(
            np.array([1, 2, 2, 1]) * 0.2 + encke[1].ra,
            np.array([-1, -1, 1, 1]) * 0.2 + encke[1].dec,
        )
        sbs.add_observations([encke_offset])

        found: List[Observation] = sbs.find_observations_by_ephemeris(encke)
        assert len(found) == 1
        assert found[0] == encke_image

        # add some padding to get the offset image
        sbs.padding = 0.2 * 60
        found = sbs.find_observations_by_ephemeris(encke)
        assert len(found) == 2
        assert encke_image in found
        assert encke_offset in found

        # remove padding, test uncertainty ellipse
        sbs.padding = 0
        sbs.uncertainty_ellipse = True
        found = sbs.find_observations_by_ephemeris(encke)
        assert len(found) == 1
        assert encke_image in found

        # embiggen uncertainty ellipse
        sbs.uncertainty_ellipse = True
        for i in range(len(encke)):
            encke[i].unc_a = 0.2 * 3600
            encke[i].unc_b = 0.2 * 3600

        found = sbs.find_observations_by_ephemeris(encke)
        assert len(found) == 2
        assert encke_image in found
        assert encke_offset in found

        # that was too much
        sbs.uncertainty_ellipse = True
        for i in range(len(encke)):
            encke[i].unc_a = 0.1 * 3600
            encke[i].unc_b = 0.1 * 3600

        found = sbs.find_observations_by_ephemeris(encke)
        assert len(found) == 1
        assert encke_image in found

        # pad it out
        sbs.padding = 6
        found = sbs.find_observations_by_ephemeris(encke)
        assert len(found) == 2
        assert encke_image in found
        assert encke_offset in found

    def test_ephemeris_uncertainty_offsets(self):
        eph: List[Ephemeris] = [
            Ephemeris(ra=0, dec=0, dra=0, ddec=1e-3,
                      unc_a=30, unc_b=10, unc_theta=0)
        ]

        with pytest.raises(ValueError):
            SBSearch._ephemeris_uncertainty_offsets(eph)

        eph.append(Ephemeris(ra=0, dec=1, dra=0, ddec=1e-3,
                             unc_a=30, unc_b=10, unc_theta=0))
        a: np.ndarray  # parallel to target motion
        b: np.ndarray  # perpendicular to target motion
        a, b = SBSearch._ephemeris_uncertainty_offsets(eph)

        # unc_a is 30 along PA=0, unc_a --> parallel --> a
        assert np.allclose(a, np.radians(30 / 3600))
        assert np.allclose(b, np.radians(10 / 3600))

        # rotate unc ellipse
        for e in eph:
            e.unc_theta = 45

        # now unc_a dominates:
        a, b = SBSearch._ephemeris_uncertainty_offsets(eph)
        assert np.allclose(a, np.radians(30 * np.cos(np.radians(45)) / 3600))
        assert np.allclose(b, np.radians(30 * np.cos(np.radians(45)) / 3600))

        # rotate unc ellipse
        for e in eph:
            e.unc_theta = 90

        # unc_a is 30 along PA=90, unc_a --> perpendicular --> b
        a, b = SBSearch._ephemeris_uncertainty_offsets(eph)
        assert np.allclose(a, np.radians(10 / 3600))
        assert np.allclose(b, np.radians(30 / 3600))

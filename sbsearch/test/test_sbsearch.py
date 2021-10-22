# Licensed with the 3-clause BSD license.  See LICENSE for details.

import pytest

from typing import List
import numpy as np
import sqlalchemy as sa
from astropy.table import Table
from astropy.time import Time
from astropy.tests.helper import remote_data

from ..sbsearch import SBSearch
from ..model import Ephemeris, Observation
from ..model.example_survey import ExampleSurveySpatialTerm, ExampleSurvey
from ..ephemeris import get_ephemeris_generator
from ..target import MovingTarget
from ..exceptions import UnknownSource, DesignationError
from ..config import Config


@pytest.fixture(name='sbs')
def fixture_sbs() -> SBSearch:
    engine: sa.engine.Engine = sa.create_engine('sqlite://')
    sessionmaker: sa.orm.sessionmaker = sa.orm.sessionmaker(bind=engine)
    with SBSearch(sessionmaker()) as sbs:
        yield sbs


@pytest.fixture(name='observations')
def fixture_observations() -> List[Observation]:
    observations: List[ExampleSurvey] = [
        ExampleSurvey(
            mjd_start=59252.1,
            mjd_stop=59252.2,
        ),
        ExampleSurvey(
            mjd_start=59252.21,
            mjd_stop=59252.31,
        )
    ]
    observations[0].set_fov([1, 2, 2, 1], [3, 3, 4, 4])
    observations[1].set_fov([2, 3, 3, 2], [3, 3, 4, 4])
    return observations


class TestSBSearch:
    def test_with_config(self) -> None:
        config: Config = Config(database='sqlite://',
                                uncertainty_ellipse=True)
        sbs: SBSearch = SBSearch.with_config(config)
        assert sbs.uncertainty_ellipse == True

    def test_source(self, sbs: SBSearch) -> None:
        with pytest.raises(ValueError):
            sbs.source

        sbs.source = ExampleSurvey
        assert sbs.source == ExampleSurvey

        sbs.source = 'example_survey'
        assert sbs.source == ExampleSurvey

    def test_source_error(self, sbs: SBSearch) -> None:
        with pytest.raises(UnknownSource):
            sbs.source = 'NEAT'

        with pytest.raises(UnknownSource):
            sbs.source = ExampleSurveySpatialTerm

    def test_sources(self, sbs: SBSearch) -> None:
        assert sbs.sources == {
            'example_survey': ExampleSurvey
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

        with pytest.raises(DesignationError):
            sbs.get_designation('1P')

        target: MovingTarget = sbs.get_designation('1P', add=True)
        assert sbs.get_designation('1P').object_id == target.object_id

    @remote_data
    def test_add_get_ephemeris(self, sbs: SBSearch) -> None:
        target: MovingTarget = sbs.add_designation('2P')
        sbs.add_ephemeris('500', target, '2021-01-01', '2021-02-01',
                          cache=True)
        eph: Table = sbs.get_ephemeris(target, '2021-01-01', '2021-02-01')
        assert len(eph) == 32
        assert np.isclose(eph[0].mjd, 59215.0)
        assert np.isclose(eph[-1].mjd, 59246.0)

    def test_add_get_observations(self, sbs, observations):
        sbs.add_observations(observations[:1])

        sbs.source = 'example_survey'
        observations = sbs.get_observations()

        obs: List[Observation] = sbs.get_observations(
            source='example_survey', mjd=[59252, 59253])
        assert obs[0].observation_id == observations[0].observation_id

    def test_get_observations_by_source(self, sbs, observations):
        sbs.add_observations(observations)
        sbs.source = 'example_survey'
        assert len(sbs.get_observations()) == 2

        new_obs = [ExampleSurvey(
            mjd_start=59215.0,
            mjd_stop=59215.1,
            fov='1:2,1:3,2:3,2:2'
        )]
        sbs.add_observations(new_obs)
        assert len(sbs.get_observations()) == 3

        # we wouldn't normally add an observation like this, but for the
        # purposes of testing, it is OK
        sbs.db.session.add(Observation(
            mjd_start=59215.0,
            mjd_stop=59215.1,
            fov='1:2,1:3,2:3,2:2'
        ))
        sbs.db.session.commit()

        assert len(sbs.get_observations()) == 3

    @remote_data
    def test_add_observations_terms(self, sbs, observations):
        sbs.add_observations(observations[:1])

        terms = (
            sbs.db.session.query(ExampleSurveySpatialTerm)
            .filter(ExampleSurveySpatialTerm.source_id
                    == observations[0].id)
            .all()
        )
        assert all(
            [term.source_id == observations[0].id
             for term in terms]
        )
        assert (len(set([term.term for term in terms])
                    - set(['$10195', '10195', '10194', '1019', '101c',
                           '101', '$10197', '10197', '$10199', '10199',
                           '1019c', '$101b', '101b', '$101c1', '101c1',
                           '101c4', '101d', '$101c7', '101c7', '$101c9',
                           '101c9', '101cc', '$101eb', '101eb', '101ec',
                           '101f'])
                    ) == 0)

    @remote_data
    def test_add_get_found(self, sbs, observations):
        # targets not really found, but we can still exercise the code
        sbs.add_observations(observations)

        found: list = sbs.get_found()
        assert len(found) == 0

        halley: MovingTarget = sbs.add_designation('1P')
        sbs.add_found(halley, observations[:1], cache=True)

        encke: MovingTarget = sbs.add_designation('2P')
        sbs.add_found(encke, observations, cache=True)

        found = sbs.get_found(target=encke)
        assert len(found) == 2
        obs_ids = [obs.observation_id for obs in observations]
        assert all([f.observation_id in obs_ids for f in found])

        found = sbs.get_found()
        assert len(found) == 3
        assert (set([f.object_id for f in found])
                == set((halley.object_id, encke.object_id)))

        assert len(sbs.get_found(mjd=[50000, 51000])) == 0
        assert len(sbs.get_found(mjd=[59252, 59252.2])) == 2

        # we wouldn't normally use an Observation like this, but it is
        # sufficient for this test
        another: Observation = Observation(
            mjd_start=59252.1,
            mjd_stop=59252.2,
        )
        another.set_fov([1, 2, 2, 1], [3, 3, 4, 4])
        sbs.db.session.add(another)
        sbs.db.session.commit()
        with pytest.raises(ValueError):
            sbs.add_found(encke, observations + [another])

    def test_re_index(self, sbs, observations):
        sbs.source = 'example_survey'
        sbs.add_observations(observations)
        count = sbs.db.session.query(ExampleSurveySpatialTerm).count()
        sbs.re_index()
        assert count == sbs.db.session.query(ExampleSurveySpatialTerm).count()

        obs = ExampleSurvey(
            mjd_start=59252.32,
            mjd_stop=59252.42,
        )
        obs.set_fov([2, 3, 3, 2], [3, 3, 4, 4])
        sbs.add_observations([obs])

        sbs.source = ExampleSurvey
        count2 = (
            sbs.db.session.query(ExampleSurveySpatialTerm)
            .join(sbs.source)
            .count()
        )
        assert count != count2
        sbs.re_index(drop_index=True)
        assert count2 == (
            sbs.db.session.query(ExampleSurveySpatialTerm)
            .join(sbs.source)
            .count()
        )

    def test_find_observations_intersecting_polygon(self, sbs, observations):
        sbs.source = 'example_survey'
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

        sbs.source = 'example_survey'
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
        sbs.source = 'example_survey'

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
        sbs.source = 'example_survey'

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
        sbs.source = 'example_survey'
        with pytest.raises(ValueError):
            sbs.find_observations_intersecting_line_at_time(
                [0, 1], [1, 0, -1], [59400, 59401])

        with pytest.raises(ValueError):
            sbs.find_observations_intersecting_line_at_time(
                [0, 1], [1, 0], [59400, 59401], a=[1, 1], b=[2])

    @remote_data
    def test_find_observations_by_ephemeris(self, sbs, observations) -> None:
        target: MovingTarget = MovingTarget('2P')

        start: Time = Time(59252, format='mjd')
        stop: Time = Time(59253, format='mjd')
        encke: List[Ephemeris] = (
            get_ephemeris_generator()
            .target_over_date_range('500@', target, start, stop, cache=True)
        )
        sbs.add_observations(observations)
        sbs.source = 'example_survey'

        encke_image: ExampleSurvey = ExampleSurvey(
            mjd_start=encke[1].mjd - 10 / 86400,
            mjd_stop=encke[1].mjd + 20 / 86400
        )
        encke_image.set_fov(
            np.array([-1, 1, 1, -1]) * 0.2 + encke[1].ra,
            np.array([-1, -1, 1, 1]) * 0.2 + encke[1].dec,
        )
        sbs.add_observations([encke_image])

        encke_offset: ExampleSurvey = ExampleSurvey(
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

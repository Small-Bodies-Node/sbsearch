# Licensed with the 3-clause BSD license.  See LICENSE for details.

import pytest

from ..model import Designation, Observation, ObservationSpatialTerm, UnspecifiedSurvey


def test_designation_repr() -> None:
    d: Designation = Designation(object_id=1, desg='1P', primary=True)
    assert repr(d) == "<Designation(desg='1P', object_id=1, primary=True)>"


def test_observation_set_fov() -> None:
    obs: Observation = Observation()
    obs.set_fov([1, 2, 2, 1], [1, 1, 2, 2])
    assert obs.fov == '1.000000:1.000000,2.000000:1.000000,2.000000:2.000000,1.000000:2.000000'


def test_observation_set_fov_error() -> None:
    obs: Observation = Observation()
    with pytest.raises(ValueError):
        obs.set_fov([1, 2, 2, 1, 1], [1, 1, 2, 2, 1])


def test_observation_repr() -> None:
    obs: Observation = Observation(
        observation_id=1,
        mjd_start=59400.0,
        mjd_stop=59400.1
    )
    obs.set_fov([1, 2, 2, 1], [1, 1, 2, 2])
    assert repr(obs) == "<Observation observation_id=1, source='observation', fov='1.000000:1.000000,2.000000:1.000000,2.000000:2.000000,1.000000:2.000000' mjd_start=59400.0 mjd_stop=59400.1>"

    obs: UnspecifiedSurvey = UnspecifiedSurvey(
        observation_id=2,
        mjd_start=59400.0,
        mjd_stop=59400.1
    )
    obs.set_fov([1, 2, 2, 1], [1, 1, 2, 2])
    assert repr(obs) == "<UnspecifiedSurvey observation_id=2, source='unspecified_survey', fov='1.000000:1.000000,2.000000:1.000000,2.000000:2.000000,1.000000:2.000000' mjd_start=59400.0 mjd_stop=59400.1>"


def test_observationspatialterm_repr() -> None:
    term: ObservationSpatialTerm = ObservationSpatialTerm(
        term_id=1,
        observation_id=1,
        term='e'
    )
    return repr(term) == '<ObservationSpatialTerm term_id=1 observation_id=1 term="e">'

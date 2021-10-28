# Licensed with the 3-clause BSD license.  See LICENSE for details.

import pytest

from ..model import Designation, Observation
from ..model.example_survey import ExampleSurvey


def test_designation_repr() -> None:
    d: Designation = Designation(object_id=1, desg='1P', primary=True)
    assert repr(d) == "<Designation(desg='1P', object_id=1, primary=True)>"


def test_observation_set_fov() -> None:
    obs: ExampleSurvey = ExampleSurvey()
    obs.set_fov([1, 2, 2, 1], [1, 1, 2, 2])
    assert obs.fov == '1.000000:1.000000,2.000000:1.000000,2.000000:2.000000,1.000000:2.000000'


def test_observation_set_fov_error() -> None:
    obs: ExampleSurvey = ExampleSurvey()
    with pytest.raises(ValueError):
        obs.set_fov([1, 2, 2, 1, 1], [1, 1, 2, 2, 1])


def test_observation_repr() -> None:
    # one wouldn't normally create an observation like this, but for the purposes
    # of testing this is OK:
    obs: Observation = Observation(
        observation_id=1,
        mjd_start=59400.0,
        mjd_stop=59400.1
    )
    obs.set_fov([1, 2, 2, 1], [1, 1, 2, 2])
    assert repr(obs) == "<Observation observation_id=1, source='observation', fov='1.000000:1.000000,2.000000:1.000000,2.000000:2.000000,1.000000:2.000000' mjd_start=59400.0 mjd_stop=59400.1>"

    obs: ExampleSurvey = ExampleSurvey(
        observation_id=2,
        mjd_start=59400.0,
        mjd_stop=59400.1
    )
    obs.set_fov([1, 2, 2, 1], [1, 1, 2, 2])
    assert repr(obs) == "<ExampleSurvey observation_id=2, source='example_survey', fov='1.000000:1.000000,2.000000:1.000000,2.000000:2.000000,1.000000:2.000000' mjd_start=59400.0 mjd_stop=59400.1>"

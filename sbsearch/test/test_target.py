# Licensed with the 3-clause BSD license.  See LICENSE for details.

from typing import List
from astropy.coordinates.sky_coordinate import SkyCoord
import pytest
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time

from ..model import Ephemeris
from ..target import *
from ..sbsdb import SBSDatabase
from ..ephemeris import set_ephemeris_generator
from ..exceptions import (
    DatabaseNotConnected,
    PrimaryDesignationError,
    SecondaryDesignationError,
    ObjectError,
    DesignationError,
)


from . import fixture_db


class TestFixedTarget:
    def test_init(self):
        target: FixedTarget = FixedTarget(SkyCoord(123 * u.deg, 45.6 * u.deg))
        assert target.ra.deg == 123
        assert target.dec.deg == 45.6

        target = FixedTarget.from_radec(123, 45.6, u.deg)
        assert target.ra.deg == 123
        assert target.dec.deg == 45.6

    def test_str(self):
        target: FixedTarget = FixedTarget(SkyCoord(123 * u.deg, 45.6 * u.deg))
        assert str(target) == "fixed(123 45.6)"

    def test_init_error(self):
        with pytest.raises(ValueError):
            FixedTarget(SkyCoord([1, 2] * u.deg, [3, 4] * u.deg))

    def test_ephemeris_at_dates(self):
        dates: Time = Time(("2020-06-01", "2020-07-01"))
        eph: List[Ephemeris] = FixedTarget(
            SkyCoord(123 * u.deg, 45.6 * u.deg)
        ).ephemeris(dates=dates)
        assert len(eph) == 2
        assert np.allclose([e.ra for e in eph], 123)
        assert np.allclose([e.dec for e in eph], 45.6)

    def test_ephemeris_over_date_range(self):
        start: Time = Time("2020-06-01")
        stop: Time = Time("2020-07-01")
        step: u.Quantity = 10 * u.day
        eph: List[Ephemeris] = FixedTarget(
            SkyCoord(123 * u.deg, 45.6 * u.deg)
        ).ephemeris(start=start, stop=stop, step=step)
        assert len(eph) == 4
        assert np.allclose([e.ra for e in eph], 123)
        assert np.allclose([e.dec for e in eph], 45.6)


class TestMovingTarget:
    def test_init(self):
        target = MovingTarget("2P")
        assert target.primary_designation == "2P"

        with pytest.raises(DatabaseNotConnected):
            target.db

    def test_str(self):
        assert str(MovingTarget("2P")) == "2P"

    @pytest.mark.parametrize("desg,object_id", (("1P", 1), ("2P", None)))
    def test_init_with_db(self, desg, object_id, db):
        target = MovingTarget(desg, db)
        assert target.primary_designation == desg
        assert target.db == db
        assert target.object_id is None

    def test_db(self):
        target = MovingTarget("2P")
        with pytest.raises(DatabaseNotConnected):
            target.db

    def test_designations(self):
        target = MovingTarget("2P")
        target.designations = ["2P", "Encke"]
        assert target.primary_designation == "2P"
        assert target.secondary_designations == ["Encke"]

    def test_designations_error(self):
        target = MovingTarget("2P")
        with pytest.raises(PrimaryDesignationError):
            target.designations = ["2P", "Encke", "2P"]

    def test_secondary_designations_error(self):
        target = MovingTarget("2P")
        with pytest.raises(SecondaryDesignationError):
            target.secondary_designations = ["Encke", "2P"]

    def test_from_id(self, db):
        target = MovingTarget.from_id(2, db)
        assert target.primary_designation == "C/1995 O1"
        assert target.secondary_designations == ["Hale-Bopp"]

    def test_from_designation(self, db):
        target = MovingTarget.from_designation("C/1995 O1", db)
        assert target.object_id == 2

    @pytest.mark.parametrize(
        "desg,object_id", (("1P", 1), ("C/1995 O1", 2), ("Hale-Bopp", 2))
    )
    def test_resolve_designation(self, desg, object_id, db):
        assert object_id == MovingTarget.resolve_designation(desg, db)

    def test_resolve_designation_error(self, db):
        with pytest.raises(DesignationError):
            MovingTarget.resolve_designation("2P", db)

    @pytest.mark.parametrize(
        "object_id,desgs", ((1, ["1P"]), (2, ["C/1995 O1", "Hale-Bopp"]))
    )
    def test_resolve_id(self, object_id, desgs, db):
        resolved = MovingTarget.resolve_id(object_id, db)
        assert len(resolved) == len(desgs)
        assert all([d in resolved for d in desgs])

    def test_resolve_id_object_error(self, db):
        with pytest.raises(ObjectError):
            target = MovingTarget.resolve_id(100, db)

    def test_resolve_id_primarydesignationerror(self, db):
        db.session.execute('UPDATE designation SET "primary"=TRUE WHERE object_id=2')
        with pytest.raises(PrimaryDesignationError):
            target = MovingTarget.resolve_id(2, db)

    def test_add(self, db):
        target = MovingTarget("2P", db)
        assert target.object_id is None
        target.add()
        assert target.object_id is not None

    def test_add_secondary_designation(self, db):
        target = MovingTarget("2P", db)
        target.add()
        target.add_secondary_designation("Encke")
        target.update()
        target_id = target.object_id
        del target
        target = MovingTarget.from_id(target_id, db)
        assert target.primary_designation == "2P"
        assert target.secondary_designations == ["Encke"]

    def test_add_designationerror(self, db):
        target = MovingTarget("C/1995 O1", db)
        with pytest.raises(DesignationError):
            target.add()

    def test_remove(self, db):
        target = MovingTarget.from_id(1, db)
        target.remove()

        assert target.object_id is None

        with pytest.raises(DesignationError):
            MovingTarget.from_designation(target.primary_designation, db)

        with pytest.raises(ObjectError):
            MovingTarget.from_id(target.object_id, db)

    def test_remove_object_id(self, db):
        MovingTarget.remove_object_id(1, db)

        with pytest.raises(DesignationError):
            MovingTarget.from_designation("1P", db)

        with pytest.raises(ObjectError):
            MovingTarget.from_id(1, db)

    def test_remove_object_id_error(self, db):
        with pytest.raises(ObjectError):
            MovingTarget.remove_object_id(100, db)

    def test_remove_designation(self, db):
        MovingTarget.remove_designation("1P", db)

        with pytest.raises(DesignationError):
            MovingTarget.from_designation("1P", db)

    def test_remove_designation_error(self, db):
        with pytest.raises(DesignationError):
            MovingTarget.remove_designation("103P", db)

    def test_update(self, db):
        target = MovingTarget.from_id(1, db)
        target.primary_designation = "1P/1982 U1"
        target.secondary_designations += ["1P", "Halley"]
        target.update()
        del target
        target = MovingTarget.from_id(1, db)
        assert target.primary_designation == "1P/1982 U1"
        assert "1P" in target.secondary_designations
        assert "Halley" in target.secondary_designations

    def test_update_error(self, db):
        target = MovingTarget("2P", db)
        with pytest.raises(ObjectError):
            target.update()

    def test_update_primary_designation(self, db):
        target = MovingTarget.from_id(1, db)
        target.update_primary_designation("1P/1982 U1")
        del target
        target = MovingTarget.from_id(1, db)
        assert target.object_id == 1
        assert target.primary_designation == "1P/1982 U1"
        assert target.secondary_designations == ["1P"]

        target.update_primary_designation("1P")
        del target
        target = MovingTarget.from_id(1, db)
        assert target.object_id == 1
        assert target.primary_designation == "1P"
        assert target.secondary_designations == ["1P/1982 U1"]

    @pytest.mark.remote_data
    def test_ephemeris_at_dates(self):
        set_ephemeris_generator("jpl")
        target: MovingTarget = MovingTarget("2P")
        dates: Time = Time(("2020-06-01", "2020-07-01"))
        eph: List[Ephemeris] = target.ephemeris(dates=dates)
        assert np.allclose([e.ra for e in eph], [65.4021, 120.97483], rtol=1e-3)
        assert np.allclose([e.dec for e in eph], [26.36761, 17.72957], rtol=1e-3)
        assert np.allclose([e.mjd for e in eph], dates.mjd)

    @pytest.mark.remote_data
    def test_ephemeris_over_date_range(self):
        set_ephemeris_generator("jpl")
        target: MovingTarget = MovingTarget("2P")
        start: Time = Time("2020-06-01")
        stop: Time = Time("2020-07-01")
        step: u.Quantity = 10 * u.day
        eph: List[Ephemeris] = target.ephemeris(start=start, stop=stop, step=step)
        assert np.allclose(
            [e.ra for e in eph], [65.4021, 81.36382, 100.97247, 120.97483], rtol=1e-3
        )
        assert np.allclose(
            [e.dec for e in eph], [26.36761, 26.87021, 24.42504, 17.72957], rtol=1e-3
        )
        assert np.allclose([e.mjd for e in eph], [59001.0, 59011.0, 59021.0, 59031.0])

    def test_ephemeris_dates_errors(self):
        target: MovingTarget = MovingTarget("2P")
        with pytest.raises(ValueError):
            target.ephemeris()

        with pytest.raises(ValueError):
            target.ephemeris(dates=1, start=1)

# Licensed with the 3-clause BSD license.  See LICENSE for details.

from abc import ABC, abstractmethod
from typing import List, Optional, Set, Tuple, Type, TypeVar, Union

from sqlalchemy import desc

import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time

from .exceptions import (
    ObjectError, DesignationError, PrimaryDesignationError,
    SecondaryDesignationError, DatabaseNotConnected
)
from .model import Ephemeris, Obj, Designation
from .sbsdb import SBSDatabase

__all__: str = [
    'Target',
    'MovingTarget',
    'FixedTarget'
]


class Target(ABC):
    """Abstract base class for all astronomical targets."""

    def ephemeris(self, observer: str = '500@', dates: Optional[Time] = None,
                  start: Optional[Time] = None, stop: Optional[Time] = None,
                  step: Optional[u.Quantity] = None) -> List[Ephemeris]:
        """Generate ephemeris for requested dates.


        Parameters
        ----------
        observer : string, optional
            Observatory code.

        dates : Time, optional
            Request these specific dates.

        start, stop : Time, optional
            Request this date range.

        step : Quantity, optional
            Step size for date range.


        Returns
        -------
        eph : list of Ephemeris
            Ephemeris at requested dates.

        """

        if dates is not None and (start is not None or stop is not None):
            raise ValueError(
                'dates and start/stop cannot be simultaneously requested.')

        if dates is not None:
            return self.ephemeris_at_dates(dates, observer=observer)
        else:
            if None in [start, stop]:
                raise ValueError(
                    'Both stop and start must be specified for date ranges.')
            return self.ephemeris_over_date_range(
                start, stop, step=step, observer=observer)

    @abstractmethod
    def ephemeris_at_dates(self, dates: Time, observer: str = '500@'
                           ) -> List[Ephemeris]:
        """Generate ephemeris at requested dates.


        Parameters
        ----------
        dates : Time
            Request these specific dates.

        observer : string, optional
            Observatory code.


        Returns
        -------
        eph : list of Ephemeris
            Ephemeris at requested dates.

        """

    @abstractmethod
    def ephemeris_over_date_range(
        self, start: Time, stop: Time, step: Optional[u.Quantity] = None,
        observer: str = '500@'
    ) -> List[Ephemeris]:
        """Generate ephemeris over requested date range.


        Parameters
        ----------
        start, stop : Time
            Request this date range.

        observer : string, optional
            Observatory code.


        Returns
        -------
        eph : list of Ephemeris
            Ephemeris at requested dates.

        """


class FixedTarget(Target):
    """Fixed (Celestial) target.


    Parameters
    ----------
    coords : SkyCoord

    """

    def __init__(self, coords: SkyCoord):
        if not coords.isscalar and len(coords) != 1:
            raise ValueError('Multiple coordinates supplied.')

        self._coords: SkyCoord = coords.icrs

    def coordinates(self, *args) -> List[Ephemeris]:
        """This target's coordinates.


        Returns
        -------
        coords : SkyCoord

        """

        return self._coords

    def ephemeris_at_dates(self, dates: Time, observer: str = '500@'
                           ) -> List[Ephemeris]:
        return [
            Ephemeris(mjd=mjd, ra=self._coords.ra.deg,
                      dec=self._coords.dec.deg)
            for mjd in dates.mjd
        ]
    ephemeris_at_dates.__doc__ = Target.ephemeris_at_dates.__doc__

    def ephemeris_over_date_range(self, start: Time, stop: Time,
                                  step: Optional[u.Quantity] = None,
                                  observer: str = '500@') -> List[Ephemeris]:
        days: float = (stop - start).jd
        _step: float = 1 if step is None else step.to_value('day')
        n: int = int(days / _step)
        dates: Time = start + (
            u.Quantity(np.arange(n + 1) * _step, 'day')
        )
        return [
            Ephemeris(mjd=mjd, ra=self._coords.ra.deg,
                      dec=self._coords.dec.deg)
            for mjd in dates.mjd
        ]
    ephemeris_over_date_range.__doc__ = (
        Target.ephemeris_over_date_range.__doc__
    )


MT = TypeVar('MT', bound='MovingTarget')


class MovingTarget(Target):
    """Moving target.


    Parameters
    ----------
    primary_designation: string
        Object primary designation.

    db: SBSDatabase, optional
        Database for storing persistent objects.

    secondary_designations: list of strings, optional
        Secondary(alternate) designations.


    Examples
    --------
    >>> target = MovingTarget('2P')

    >>> db = SBSDatabase('sqlite://')
    >>> target = MovingTarget('2P', db=db)

    # use ``from_id`` or ``from_desg`` when you want an object copied
    # from the database:
    >>> target = MovingTarget.from_id(1, db)
    >>> target = MovingTarget.from_desg('2P', db)

    """

    def __init__(self, primary_designation: str,
                 db: Optional[SBSDatabase] = None,
                 secondary_designations: List[str] = []):
        self.primary_designation: str = primary_designation
        self.secondary_designations: List[str] = secondary_designations
        self._db: Union[SBSDatabase, None] = db
        self._object_id: Union[int, None] = None

    @property
    def db(self) -> SBSDatabase:
        """SBSDatabase object for object persistence."""
        if self._db is None:
            raise DatabaseNotConnected
        return self._db

    @property
    def object_id(self) -> Union[int, None]:
        """Object ID.

        Generally, if this is set, the object has been stored in the database.

        """
        return self._object_id

    @property
    def designations(self) -> List[str]:
        """All object designations, the primary is first."""
        return [self.primary_designation] + list(self._secondary_designations)

    @designations.setter
    def designations(self, desigs: List[str]):
        if desigs[0] in desigs[1:]:
            raise PrimaryDesignationError(f'{desigs[0]} repeated')
        self.primary_designation = desigs[0]
        self.secondary_designations = desigs[1:]

    @property
    def secondary_designations(self) -> List[str]:
        """Secondary / alternate designations."""
        return list(self._secondary_designations)

    @secondary_designations.setter
    def secondary_designations(self, desgs: List[str]):
        # check that none match the primary
        if self.primary_designation in desgs:
            raise SecondaryDesignationError(
                f'{self.primary_designation} is already defined as the '
                'primary designation')
        self._secondary_designations: Set[str] = set(desgs)

    @classmethod
    def from_id(cls: Type[MT], object_id: int, db: SBSDatabase) -> MT:
        """Initialize from database based on object ID.


        Parameters
        ----------
        object_id: int
            Object ID.

        db: SBSDatabase
            Database interface.


        Returns
        -------
        target: MovingTarget


        Raises
        ------
        ObjectError

        """

        desgs: List[str] = MovingTarget.resolve_id(object_id, db)
        target: MovingTarget = cls(
            desgs[0], db=db, secondary_designations=desgs[1:])
        target._object_id = object_id
        return target

    @classmethod
    def from_designation(cls: Type[MT], desg: str, db: SBSDatabase) -> MT:
        """Initialize from database based on object designation.


        Parameters
        ----------
        desg: string
            Object designation.

        db: SBSDatabase
            Database interface.


        Returns
        -------
        target: MovingTarget


        Raises
        ------
        DesignationError

        """

        object_id: int = cls.resolve_designation(desg, db)
        return MovingTarget.from_id(object_id, db)

    @staticmethod
    def resolve_designation(desg: str, db: SBSDatabase) -> int:
        """Resolve designation to database object ID.


        Parameters
        ----------
        desg: string
            Object designation.

        db: SBSDatabase
            Database interface.


        Returns
        -------
        object_id: int
            Database object ID.


        Raises
        ------
        DesignationError

        """

        object_id: Union[None, Tuple[int]] = (
            db.session.query(Designation.object_id)
            .filter(Designation.desg == desg)
            .first())
        if object_id is None:
            raise DesignationError('')

        return object_id[0]

    @staticmethod
    def resolve_id(object_id: int, db: SBSDatabase) -> List[str]:
        """Resolve object ID to designations.

        Parameters
        ----------
        object_id: int
            Database object ID.

        db: SBSDatabase
            Database interface.

        Returns
        -------
        desgs: list of strings
            Object designation.  The first is the "primary" designation.

        """

        designations: Tuple[Designation] = (
            db.session.query(Designation)
            .filter(Designation.object_id == object_id)
            .order_by(desc(Designation.primary))
            .all()
        )

        if len(designations) == 0:
            raise ObjectError(f'{object_id} not found.')

        # verify that there is only one primary designation
        test: int = sum([d.primary for d in designations])
        if test != 1:
            raise PrimaryDesignationError(
                f'{object_id} has {test} primary designation')

        # unpack into a single list
        desgs: List[str] = [d.desg for d in designations]

        return desgs

    def add(self) -> None:
        """Add to database and set object ID.

        Also adds secondary designations, if any are defined.

        Raises
        ------
        DesignationError if any of the designations are already in use.

        """

        # verify that the object's designations are not already in the database
        n: int = (self.db.session.query(Designation)
                  .filter(Designation.desg.in_(self.designations))
                  .count())
        if n > 0:
            raise DesignationError(
                f'{n} of the target designations are already used')

        obj: Obj = Obj()
        self.db.session.add(obj)
        self.db.session.commit()

        i: int
        for i in range(len(self.designations)):
            desg: Designation = Designation(
                object_id=obj.object_id, desg=self.designations[i],
                primary=i == 0)
            self.db.session.add(desg)

        self.db.session.commit()
        self._object_id = obj.object_id

    def add_secondary_designation(self, desg: str) -> None:
        """Add secondary designation to the object.

        Parameters
        ----------
        desg: string

        """

        self.secondary_designations += [desg]

    def remove(self) -> None:
        """Remove this object from the database.

        Raises
        ------
        ObjectError

        """

        MovingTarget.remove_object_id(self.object_id, self.db)
        self._object_id = None

    @staticmethod
    def remove_object_id(object_id: int, db: SBSDatabase) -> None:
        """Remove given object ID from the database.

        Parameters
        ----------
        object_id: int
            Database object ID.

        db: SBSDatabase
            SBSearch database object.

        Raises
        ------
        ObjectError if it does not exist.

        """

        obj: Union[Obj, None] = (
            db.session.query(Obj)
            .filter(Obj.object_id == object_id)
            .first()
        )
        if obj is None:
            raise ObjectError(f'{object_id} not in database.')

        db.session.delete(obj)
        db.session.commit()

    @staticmethod
    def remove_designation(desg: str, db: SBSDatabase) -> None:
        """Remove given designation from the database.

        The object ID associated with the designation will always persist
        after designation removal.

        Parameters
        ----------
        desg: string
            Designation in database.

        db: SBSDatabase
            SBSearch database object.

        Raises
        ------
        DesignationError if it does not exist.

        """

        d: Union[Designation, None] = (
            db.session.query(Designation)
            .filter(Designation.desg == desg)
            .first()
        )
        if d is None:
            raise DesignationError(f'{desg} not in database.')

        db.session.delete(d)
        db.session.commit()

    def update(self) -> None:
        """Sync object designations with database.

        Raises
        ------
        ObjectError

        """

        if self.object_id is None:
            raise ObjectError('Object ID not defined; consider add() instead.')

        # get current designations
        current: List[str] = MovingTarget.resolve_id(self.object_id, self.db)

        # delete them
        desg: str
        for desg in current:
            self.remove_designation(desg, self.db)
        self.db.session.commit()

        i: int
        for i in range(len(self.designations)):
            desg: Designation = Designation(
                object_id=self.object_id, desg=self.designations[i],
                primary=i == 0)
            self.db.session.add(desg)

        self.db.session.commit()

    def update_primary_designation(self, desg: str) -> None:
        """Update primary designation in the database.

        Moves the old primary designation to a secondary designation.

        Parameters
        ----------
        desg: string
            The new primary designation.  Does not need to be a current
            designation of the object.

        """

        # move current primary to secondary
        d: Designation = (
            self.db.session.query(Designation)
            .filter(Designation.object_id == self.object_id)
            .filter(Designation.primary)
            .one()
        )
        d.primary = False
        self.db.session.merge(d)

        # if desg does not exist, add it, else update
        d: Union[Designation, None] = (
            self.db.session.query(Designation)
            .filter(Designation.desg == desg)
            .first()
        )

        if d is None:
            d = Designation(
                object_id=self.object_id,
                desg=desg,
                primary=True
            )
        else:
            d.primary = True

        self.db.session.merge(d)
        self.db.session.commit()

        # sync designations from the database
        self.designations = MovingTarget.resolve_id(self.object_id, self.db)

    def ephemeris_at_dates(
        self, dates: Time, observer: str = '500@', cache=True
    ) -> List[Ephemeris]:
        from .ephemeris import EphemerisGenerator, get_ephemeris_generator
        g: EphemerisGenerator = get_ephemeris_generator()
        return g.target_at_dates(observer, self, dates, cache=cache)

    ephemeris_at_dates.__doc__ = Target.ephemeris_at_dates.__doc__

    def ephemeris_over_date_range(
        self, start: Time, stop: Time, step: Optional[u.Quantity] = None,
        observer: str = '500@', cache=True
    ) -> List[Ephemeris]:
        from .ephemeris import EphemerisGenerator, get_ephemeris_generator
        g: EphemerisGenerator = get_ephemeris_generator()
        return g.target_over_date_range(
            observer, self, start, stop, step, cache=cache)

    ephemeris_over_date_range.__doc__ = (
        Target.ephemeris_over_date_range.__doc__
    )

# Licensed with the 3-clause BSD license.  See LICENSE for details.

from typing import Union, List
import numpy as np
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import (Column, BigInteger, Integer, Float, String, ForeignKey,
                        Boolean, Text)
from sqlalchemy.dialects.postgresql import ARRAY
from ..spatial import SpatialIndexer  # pylint: disable=E0611

__all__: List[str] = [
    'Base',
    'Designation',
    'Obj',
    'Observation',
    'Ephemeris',
    'Found',
]

Base = declarative_base()


class Obj(Base):
    __tablename__: str = 'obj'
    object_id: int = Column(Integer, primary_key=True)


class Orbit(Base):
    __tablename__: str = 'orbit'
    orbit_id: int = Column(Integer, primary_key=True)
    # remaining columns are TBD


class Designation(Base):
    __tablename__: str = 'designation'
    desg_id: int = Column(Integer, primary_key=True)
    desg: str = Column(String(64), unique=True, index=True, nullable=False)
    object_id: int = Column(
        Integer, ForeignKey('obj.object_id', onupdate='CASCADE',
                            ondelete='CASCADE'),
        index=True, nullable=False)
    primary: bool = Column(Boolean, nullable=False)

    def __repr__(self) -> str:
        return (f"<Designation(desg='{self.desg}', "
                f"object_id={self.object_id}, primary={self.primary})>")


class Ephemeris(Base):
    __tablename__: str = 'ephemeris'
    ephemeris_id: int = Column(BigInteger, primary_key=True)
    object_id: int = Column(
        Integer, ForeignKey('obj.object_id', onupdate='CASCADE',
                            ondelete='CASCADE'),
        index=True, nullable=False)
    # keep synced with Found
    mjd: float = Column(Float(32), index=True, nullable=False,
                        doc='Modified Julian date, UT')
    rh: float = Column(Float(32), doc='heliocentric distance, au')
    delta: float = Column(Float(32), doc='observer-target distance, au')
    phase: float = Column(Float(32), doc='Sun-comet-observer angle, deg')
    drh: float = Column(Float(32), doc='heliocentric velocity, km/s')
    true_anomaly: float = Column(Float(32), doc='true anomaly angle, deg')
    ra: float = Column(Float(32), nullable=False,
                       doc='Right Ascension ICRF, deg')
    dec: float = Column(Float(32), nullable=False,
                        doc='Declination ICRF, deg')
    dra: float = Column(Float(32), doc=('sky motion, includes cos(Dec) factor,'
                                        ' arcsec/hr'))
    ddec: float = Column(Float(32), doc='arcsec/hr')
    unc_a: float = Column(
        Float(32), doc='error ellipse semi-major axis, arcsec')
    unc_b: float = Column(
        Float(32), doc='error ellipse semi-minor axis, arcsec')
    unc_theta: float = Column(
        Float(32), doc='error ellipse position angle (E of N), deg')
    elong: float = Column(Float(32), doc='solar elongation')
    sangle: float = Column(
        Float(32), doc='projected comet-Sun position angle, deg')
    vangle: float = Column(
        Float(32), doc='projected comet velocity vector position angle, deg')
    vmag: float = Column(Float(32), doc='predicted visual brightness, mag')
    retrieved: str = Column(String(64))

    # Ephemeris comparisons
    def __lt__(self, other) -> bool:
        """Compare dates."""
        return self.mjd < other.mjd

    def __le__(self, other) -> bool:
        """Compare dates."""
        return self.mjd <= other.mjd

    def __eq__(self, other) -> bool:
        """Compare dates using numpy.isclose."""
        return self.mjd == other.mjd

    def __ne__(self, other) -> bool:
        """Compare dates using numpy.isclose, see __eq__."""
        return not self.__eq__(other)

    def __gt__(self, other) -> bool:
        """Compare dates."""
        return self.mjd > other.mjd

    def __ge__(self, other) -> bool:
        """Compare dates."""
        return self.mjd >= other.mjd


class Observation(Base):
    """Observation data class.

    Derive all surveys from this class.  See ExampleSurvey for an example.

    """
    # Override the following in derived classes.
    __tablename__: str = 'observation'
    __data_source_name__ = 'All data sources'
    __obscode__ = '500'  # MPC observatory code

    # Required attributes for basic sbsearch functionality.
    observation_id: int = Column(BigInteger, primary_key=True)
    source: str = Column(String(64), default='observation',
                         doc='source survey')
    mjd_start: float = Column(Float(32), nullable=False, index=True,
                              doc='shutter open, modified Julian date, UTC')
    mjd_stop: float = Column(Float(32), nullable=False, index=True,
                             doc='shutter close, modified Julian date, UTC')
    fov: str = Column(String(128), nullable=False, doc=(
        'field of view as set of comma-separated RA:Dec pairs in degrees,'
        'e.g., "1:1, 1:2, 2:1" (tip: see set_fov)'))
    spatial_terms: List[str] = Column(ARRAY(Text), nullable=False,
                                      doc='spatial index terms')

    # Common attributes.  Additional attributes and data-set-specific
    # attributes are defined in each data set object.
    filter: str = Column(String(16), doc='filter/bandpass')
    exposure: float = Column(Float(32), doc='exposure time, s')
    seeing: float = Column(Float(32), doc='point source FWHM, arcsec')
    airmass: float = Column(Float(32), doc='observation airmass')
    maglimit: float = Column(Float(32), doc='detection limit, mag')

    # Class methods.
    def set_fov(self, ra: Union[List[float], np.ndarray],
                dec: Union[List[float], np.ndarray]) -> None:
        """Set ``fov``  with these vertices.

        Parameters
        ----------
        ra, dec : arrays of float
            Vertices expressed in degrees.

        """

        if len(ra) > 4:
            raise ValueError('No more than 4 vertices are allowed.')
        values = []
        for i in range(len(ra)):
            values.append(f'{ra[i]:.6f}:{dec[i]:.6f}')
        self.fov = ','.join(values)

    def __repr__(self) -> str:
        return (f'<{self.__class__.__name__}'
                f' observation_id={self.observation_id},'
                f' source={repr(self.source)}, fov={repr(self.fov)}'
                f' mjd_start={self.mjd_start} mjd_stop={self.mjd_stop}>')

    __mapper_args__ = {
        "polymorphic_identity": "observation",
        "polymorphic_on": source
    }


class Found(Base):
    __tablename__ = 'found'
    found_id = Column(BigInteger, primary_key=True)
    object_id = Column(Integer,
                       ForeignKey('obj.object_id',
                                  onupdate='CASCADE',
                                  ondelete='CASCADE'),
                       index=True)
    orbit_id = Column(Integer,
                      ForeignKey('orbit.orbit_id',
                                 onupdate='CASCADE',
                                 ondelete='CASCADE'),
                      index=True)
    observation_id = Column(BigInteger,
                            ForeignKey('observation.observation_id',
                                       onupdate='CASCADE',
                                       ondelete='CASCADE'))
    # keep synced with Ephemeris
    mjd: float = Column(Float(32), index=True, nullable=False,
                        doc='Modified Julian date, UTC')
    rh: float = Column(Float(32), doc='heliocentric distance, au')
    delta: float = Column(Float(32), doc='observer-target distance, au')
    phase: float = Column(Float(32), doc='Sun-target-observer angle, deg')
    drh: float = Column(Float(32), doc='heliocentric velocity, km/s')
    true_anomaly: float = Column(Float(32), doc='true anomaly angle, deg')
    ra: float = Column(Float(32), nullable=False,
                       doc='Right Ascension ICRF, deg')
    dec: float = Column(Float(32), nullable=False,
                        doc='Declination ICRF, deg')
    dra: float = Column(Float(32), doc=('sky motion, includes cos(Dec) factor,'
                                        ' arcsec/hr'))
    ddec: float = Column(Float(32), doc='arcsec/hr')
    unc_a: float = Column(Float(32),
                          doc='error ellipse semi-major axis, arcsec')
    unc_b: float = Column(Float(32),
                          doc='error ellipse semi-minor axis, arcsec')
    unc_theta: float = Column(Float(32),
                              doc='error ellipse position angle (E of N), deg')
    elong: float = Column(Float(32), doc='solar elongation')
    sangle: float = Column(Float(32),
                           doc='projected comet-Sun position angle, deg')
    vangle: float = Column(Float(32),
                           doc=('projected comet velocity vector position '
                                'angle, deg'))
    vmag: float = Column(Float(32), doc='predicted visual brightness, mag')
    retrieved: str = Column(String(64))

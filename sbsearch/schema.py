# Licensed with the 3-clause BSD license.  See LICENSE for details.

__all__ = [
    'Obj',
    'Eph',
    'Obs',
    'GenericObs',
    'Found'
]

import sqlalchemy as sa
from sqlalchemy import (MetaData, Table, Column, Integer, Float, String,
                        BigInteger, LargeBinary, ForeignKey, Index)
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.ext.compiler import compiles
from geoalchemy2 import Geography
from .exceptions import MissingSpatialReferenceSystemError

Base = declarative_base()


class Obj(Base):
    __tablename__ = 'obj'
    objid = Column(Integer, primary_key=True)
    desg = Column(String(64), unique=True)


class Eph(Base):
    __tablename__ = 'eph'
    ephid = Column(BigInteger, primary_key=True)
    objid = Column(Integer,
                   ForeignKey(
                       'obj.objid',
                       onupdate='CASCADE',
                       ondelete='CASCADE'),
                   index=True)
    jd = Column(Float(32), index=True, doc='Julian date, UT')
    rh = Column(Float(32), doc='heliocentric distance, au')
    delta = Column(Float(32), doc='observer-target distance, au')
    ra = Column(Float(32), doc='Right Ascension ICRF, deg')
    dec = Column(Float(32), doc='Declination ICRF, deg')
    dra = Column(Float(32), doc=('sky motion, includes cos(Dec) factor,'
                                 ' arcsec/hr'))
    ddec = Column(Float(32), doc='arcsec/hr')
    unc_a = Column(Float(32), doc='error ellipse semi-major axis, arcsec')
    unc_b = Column(Float(32), doc='error ellipse semi-minor axis, arcsec')
    unc_theta = Column(Float(32), doc='error ellipse position angle, deg')
    vmag = Column(Float(32), doc='predicted visual brightness, mag')
    segment = Column(Geography('LINESTRING', 40001))
    retrieved = Column(String(64))


class Obs(Base):
    __tablename__ = 'obs'
    obsid = Column(BigInteger, primary_key=True)
    source = Column(String(64), index=True, default='unspecified',
                    doc='source survey')
    jd_start = Column(Float(32), index=True, doc=(
        'shutter open, Julian date, UT'))
    jd_stop = Column(Float(32), index=True, doc=(
        'shutter close, Julian date, UT'))
    fov = Column(Geography('POLYGON', 40001), doc='image field of view')
    filter = Column(String(16), doc='filter/bandpass')
    exposure = Column(Float(32), doc='exposure time, s')
    seeing = Column(Float(32), doc='point source FWHM, arcsec')
    airmass = Column(Float(32))
    maglimit = Column(
        Float(32), doc='detection limit for point sources (mag)')

    __mapper_args__ = {
        "polymorphic_identity": "obs",
        "polymorphic_on": source
    }


class GenericObs(Obs):
    __tablename__ = 'generic_obs'
    id = Column(BigInteger, primary_key=True)
    obsid = Column(ForeignKey('obs.obsid', onupdate='CASCADE',
                              ondelete='CASCADE'))
    __mapper_args__ = {
        'polymorphic_identity': 'unspecified'
    }


class Found(Base):
    __tablename__ = 'found'
    foundid = Column(BigInteger, primary_key=True)
    objid = Column(BigInteger,
                   ForeignKey(
                       'obj.objid',
                       onupdate='CASCADE',
                       ondelete='CASCADE'),
                   index=True)
    obsid = Column(BigInteger,
                   ForeignKey(
                       'obs.obsid',
                       onupdate='CASCADE',
                       ondelete='CASCADE'),
                   index=True)
    jd = Column(Float(32), doc='observation mid-time, Julian date, UT')
    ra = Column(Float(32), doc='Right Ascension ICRF, deg')
    dec = Column(Float(32), doc='Declination ICRF, deg')
    dra = Column(Float(32), doc='arcsec/hr')
    dra = Column(Float(32), doc=('sky motion, includes cos(Dec) factor,'
                                 ' arcsec/hr'))
    ddec = Column(Float(32), doc='arcsec/hr')
    unc_a = Column(Float(32), doc='error ellipse semi-major axis, arcsec')
    unc_b = Column(Float(32), doc='error ellipse semi-minor axis, arcsec')
    unc_theta = Column(Float(32), doc='error ellipse position angle, deg')
    vmag = Column(Float(32), doc='predicted visual brightness, mag')
    rh = Column(Float(32), doc='heliocentric distance, au')
    rdot = Column(Float(32), doc='heliocentric radial velocity, km/s')
    delta = Column(Float(32), doc='observer-target distance, au')
    phase = Column(Float(32), doc='Sun-observer-target angle, deg')
    selong = Column(Float(32), doc='solar elongation, deg')
    sangle = Column(Float(32), doc=('projected target-Sun vector position'
                                    ' angle, deg'))
    vangle = Column(Float(32), doc=('projected target velocity vector'
                                    ' position angle, deg'))
    trueanomaly = Column(Float(32), doc='deg')
    tmtp = Column(Float(32), doc='time from perihelion, T-T_P, days')

    Index('found_obj_obs', 'objid', 'obsid', unique=True)


def _exists_sbsearch_spatial_ref_sys(con):
    # verify that a plain spherical coordinate system exists
    try:
        proj4text = con.execute('''
        SELECT proj4text FROM public.spatial_ref_sys
        WHERE srid=40001
        ''').fetchone()[0].strip()
        return proj4text == '+proj=longlat +ellps=sphere +no_defs'
    except TypeError:
        return False


def _add_sbsearch_spatial_ref_sys(con):
    # based on answer from PolyGeo at:
    # https://gis.stackexchange.com/questions/2459/what-coordinate-system-should-be-used-to-store-geography-data-for-celestial-coor
    con.execute('''
    INSERT INTO public.spatial_ref_sys
    VALUES(
      40001,
      'SBSearch',
      1,
      'GEOGCS["Normal Sphere (r=6370997)",DATUM["unknown",SPHEROID["sphere",6370997,0]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]]',
      '+proj=longlat +ellps=sphere +no_defs'
    );''')
    return _exists_sbsearch_spatial_ref_sys(con)


def create(engine):
    con = engine.connect()

    if not _exists_sbsearch_spatial_ref_sys(con):
        verified = _add_sbsearch_spatial_ref_sys(con)
        if not verified:
            raise MissingSpatialReferenceSystemError(
                'Missing spatial reference system.  Are the PostGIS '
                'extensions loaded?')

    Base.metadata.create_all(con)
    con.close()

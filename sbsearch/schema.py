# Licensed with the 3-clause BSD license.  See LICENSE for details.

__all__ = [
    'Obj',
    'Eph',
    'eph_tree',
    'Obs',
    'obs_tree',
    'Found'
]

import sqlalchemy as sa
from sqlalchemy import (MetaData, Table, Column, Integer, Float, String,
                        LargeBinary, ForeignKey)
from sqlalchemy.ext.declarative import declarative_base

Base = declarative_base()


class Obj(Base):
    __tablename__ = 'obj'
    objid = Column(Integer, primary_key=True)
    name = Column(String(64), unique=True)


class Eph(Base):
    __tablename__ = 'eph'
    ephid = Column(Integer, primary_key=True)
    objid = Column(Integer,
                   ForeignKey('obj.objid', onupdate='CASCADE',
                              ondelete='CASCADE'),
                   index=True)
    jd = Column(Float(64))
    rh = Column(Float(32))
    delta = Column(Float(32))
    ra = Column(Float(64))
    dec = Column(Float(64))
    dra = Column(Float(32))
    ddec = Column(Float(32))
    unc_a = Column(Float(32))
    unc_b = Column(Float(32))
    unc_theta = Column(Float(32))
    vmag = Column(Float(32))
    retrieved = Column(String(64))


class Obs(Base):
    __tablename__ = 'obs'
    obsid = Column(Integer, primary_key=True)
    source = Column(String(64), index=True)
    jd_start = Column(Float(64))
    jd_stop = Column(Float(64))
    fov = Column(LargeBinary(),
                 comment="RA, Dec center and corners in radians")


class Found(Base):
    __tablename__ = 'found'
    foundid = Column(Integer, primary_key=True)
    objid = Column(Integer,
                   ForeignKey('obj.objid', onupdate='CASCADE',
                              ondelete='CASCADE'),
                   index=True)
    obsid = Column(Integer, ForeignKey('obs.obsid', onupdate='CASCADE',
                                       ondelete='CASCADE'))
    obsjd = Column(Float(64))
    ra = Column(Float(64))
    dec = Column(Float(64))
    dra = Column(Float(32))
    ddec = Column(Float(32))
    unc_a = Column(Float(32))
    unc_b = Column(Float(32))
    unc_theta = Column(Float(32))
    vmag = Column(Float(32))
    rh = Column(Float(32))
    rdot = Column(Float(32))
    delta = Column(Float(32))
    phase = Column(Float(32))
    selong = Column(Float(32))
    sangle = Column(Float(32))
    vangle = Column(Float(32))
    trueanomaly = Column(Float(32))
    tmtp = Column(Float(32))


sqlite_eph_tree = Table('eph_tree', base.metadata
'''
CREATE VIRTUAL TABLE IF NOT EXISTS eph_tree USING RTREE(
    ephid INTEGER PRIMARY KEY,
    mjd0 FLOAT,
    mjd1 FLOAT,
    x0 FLOAT,
    x1 FLOAT,
    y0 FLOAT,
    y1 FLOAT,
    z0 FLOAT,
    z1 FLOAT
);
/* observation rtree */
CREATE VIRTUAL TABLE IF NOT EXISTS obs_tree USING RTREE(
  obsid INTEGER PRIMARY KEY,
  mjd0 FLOAT,
  mjd1 FLOAT,
  x0 FLOAT,
  x1 FLOAT,
  y0 FLOAT,
  y1 FLOAT,
  z0 FLOAT,
  z1 FLOAT
);

tree triggers:

    CREATE TRIGGER IF NOT EXISTS delete_eph_from_tree BEFORE DELETE ON eph
    BEGIN
      DELETE FROM eph_tree WHERE ephid=old.ephid;
    END;

    CREATE TRIGGER IF NOT EXISTS delete_object_from_eph BEFORE DELETE ON obj
    BEGIN
      DELETE FROM eph WHERE objid=old.objid;
    END;

    CREATE TRIGGER IF NOT EXISTS delete_obs_from_obs_tree
    BEFORE DELETE ON obs
    BEGIN
      DELETE FROM obs_tree WHERE obsid=old.obsid;
    END;
'''

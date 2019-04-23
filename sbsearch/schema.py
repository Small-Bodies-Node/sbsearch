# Licensed with the 3-clause BSD license.  See LICENSE for details.

__all__ = [
    'obj',
    'eph', 
    'eph_tree',
    'obs',
    'obs_tree',
    'found'
]

import sqlalchemy as sa
from sqlalchemy import (MetaData, Table, Column, Integer, Float, String,
                        ForeignKey)

metadata = MetaData()

obj = Table('obj', metadata,
            Column('objid', Integer, primary_key=True),
            Column('name', String(64), unique=True))

eph = Table('eph', metadata,
            Column('ephid', Integer, primary_key=True),
            Column('objid', Integer,
                   ForeignKey('obj.objid', onupdate='CASCADE',
                              ondelete='CASCADE'),
                   index=True),
            Column('jd', Float(32)),
            Column('rh', Float(32)),
            Column('delta', Float(32)),
            Column('ra', Float(32)),
            Column('dec', Float(32)),
            Column('dra', Float(32)),
            Column('ddec', Float(32)),
            Column('unc_a', Float(16)),
            Column('unc_b', Float(16)),
            Column('unc_theta', Float(16)),
            Column('vmag', Float(16)),
            Column('retrieved', String(64)),
)

sa.event.listen(mysql_eph_tree, 'after_create',
                sa.DDL("CREATE SPATIAL INDEX mjd ON eph_tree(mjd)"))
            

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

/* delete from eph tree before delete from eph table */
CREATE TRIGGER IF NOT EXISTS delete_eph_from_tree BEFORE DELETE ON eph
BEGIN
    DELETE FROM eph_tree WHERE ephid=old.ephid;
END;

/* observation table, RA, Dec in radians */
CREATE TABLE IF NOT EXISTS obs(
  obsid INTEGER PRIMARY KEY,
  source TEXT,
  jd_start FLOAT,
  jd_stop FLOAT,
  fov BLOB
);

CREATE INDEX IF NOT EXISTS obs_sources ON obs (source);

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

/* delete from obs tree before delete from obs table */
CREATE TRIGGER IF NOT EXISTS delete_obs_from_obs_tree BEFORE DELETE ON obs
BEGIN
  DELETE FROM obs_tree WHERE obsid=old.obsid;
END;

/* objects found in obs table */
CREATE TABLE IF NOT EXISTS found(
  foundid INTEGER PRIMARY KEY,
  objid INTEGER,
  obsid INTEGER,
  obsjd FLOAT,
  ra FLOAT,
  dec FLOAT,
  dra FLOAT,
  ddec FLOAT,
  ra3sig FLOAT,
  dec3sig FLOAT,
  vmag FLOAT,
  rh FLOAT,
  rdot FLOAT,
  delta FLOAT,
  phase FLOAT,
  selong FLOAT,
  sangle FLOAT,
  vangle FLOAT,
  trueanomaly FLOAT,
  tmtp FLOAT,
  FOREIGN KEY(objid) REFERENCES obj(objid) ON DELETE CASCADE,
  FOREIGN KEY(obsid) REFERENCES obs(obsid) ON DELETE CASCADE
);

'''

    db.executescript(schema)


def create(engine):
    metadata = sa.MetaData()
    conn = engine.connect()
    if not engine.dialect.has_table(conn, 'obj'):
        pass

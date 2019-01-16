# Licensed with the 3-clause BSD license.  See LICENSE for details.


def create(db):
    schema = '''
CREATE TABLE IF NOT EXISTS obj(
    objid INTEGER PRIMARY KEY,
    desg TEXT UNIQUE
);

CREATE TABLE IF NOT EXISTS eph(
    ephid INTEGER PRIMARY KEY,
    objid INTEGER,
    jd FLOAT,
    rh FLOAT,
    delta FLOAT,
    ra FLOAT,
    dec FLOAT,
    dra FLOAT,
    ddec FLOAT,
    unc_a FLOAT,
    unc_b FLOAT,
    unc_theta FLOAT,
    vmag FLOAT,
    retrieved TEXT,
    FOREIGN KEY(objid) REFERENCES obj(objid)
);

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

CREATE TRIGGER IF NOT EXISTS delete_eph_from_tree BEFORE DELETE ON eph
BEGIN
    DELETE FROM eph_tree WHERE ephid=old.ephid;
END;

CREATE TRIGGER IF NOT EXISTS delete_object_from_eph BEFORE DELETE ON obj
BEGIN
    DELETE FROM eph WHERE objid=old.objid;
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
  FOREIGN KEY(objid) REFERENCES obj(objid)
);

CREATE TRIGGER IF NOT EXISTS delete_obs_from_found BEFORE DELETE ON obs
BEGIN
  DELETE FROM found WHERE obsid=old.obsid;
END;

CREATE TRIGGER IF NOT EXISTS delete_object_from_found BEFORE DELETE ON obj
BEGIN
  DELETE FROM found WHERE objid=old.objid;
END;
'''

    db.executescript(schema)

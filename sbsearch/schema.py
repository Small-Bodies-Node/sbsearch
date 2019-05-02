# Licensed with the 3-clause BSD license.  See LICENSE for details.

__all__ = [
    'Obj',
    'Eph',
    'EphTree',
    'Obs',
    'ObsTree',
    'Found'
]

import sqlalchemy as sa
from sqlalchemy import (MetaData, Table, Column, Integer, Float, String,
                        LargeBinary, ForeignKey)
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.ext.compiler import compiles

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
    ra = Column(Float(64), doc='rad')
    dec = Column(Float(64), doc='rad')
    dra = Column(Float(32), doc='arcsec/hr')
    ddec = Column(Float(32), doc='arcsec/hr')
    unc_a = Column(Float(32), doc='rad')
    unc_b = Column(Float(32), doc='rad')
    unc_theta = Column(Float(32), doc='rad')
    vmag = Column(Float(32))
    retrieved = Column(String(64))


class Obs(Base):
    __tablename__ = 'obs'
    obsid = Column(Integer, primary_key=True)
    source = Column(String(64), index=True)
    jd_start = Column(Float(64))
    jd_stop = Column(Float(64))
    fov = Column(
        LargeBinary(), doc="RA, Dec center and corners in radians")
    filter = Column(String(16))
    seeing = Column(Float(32))
    airmass = Column(Float(32))
    maglimit = Column(
        Float(32), doc='detection limit for point sources (mag)')

    def coords_to_fov(self, *coords):
        """Save coordinate list as fov array.


        Parameters
        ----------
        *coords : lists or tuples
            RA, Dec pairs for the field of view center and
            corners in radians.


        Notes
        -----
        The order of the coordinates is not important

        """

        radec = sum(coords, type(coords[0])())
        self.fov = struct.pack('{}d'.format(len(coords) * 2), *radec)

    def fov_to_coords(self):
        """Convert FOV to array of coodinates.


        Returns
        -------
        coords : ndarray
            Nx2 array of RA and Dec (radians).

        """

        N = len(self.fov) // 80
        c = struct.unpack('{}d'.format(10 * N), fov)
        return np.array(c).reshape((N, 2))

    def fov_to_xyz(self):
        """Convert FOV to array of Cartesian coordinates.

        Returns
        -------
        xyz : ndarray
            Nx3 array of x, y, z.

        """
        ra, dec = self.fov_to_coords().T
        return np.array((np.cos(dec) * np.cos(ra),
                         np.cos(dec) * np.sin(ra),
                         np.sin(dec))).T

    __mapper_args__ = {
        "polymorphic_identity": "obs",
        "polymorphic_on": source
    }


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


class RTree(Base):
    __abstract__ = True


class EphTree(RTree):
    __tablename__ = 'eph_tree'
    ephid = Column(Integer, primary_key=True)
    mjd0 = Column(Float(32))
    mjd1 = Column(Float(32))
    x0 = Column(Float(32))
    x1 = Column(Float(32))
    y0 = Column(Float(32))
    y1 = Column(Float(32))
    z0 = Column(Float(32))
    z1 = Column(Float(32))

    __mapper_args__ = {
        'polymorphic_identity': 'eph_tree',
    }


class ObsTree(RTree):
    __tablename__ = 'obs_tree'
    obsid = Column(Integer, primary_key=True)
    mjd0 = Column(Float(32))
    mjd1 = Column(Float(32))
    x0 = Column(Float(32))
    x1 = Column(Float(32))
    y0 = Column(Float(32))
    y1 = Column(Float(32))
    z0 = Column(Float(32))
    z1 = Column(Float(32))

    __mapper_args__ = {
        'polymorphic_identity': 'obs_tree',
    }


@compiles(sa.schema.CreateTable)
def _compile_tables(element, compiler, **kwargs):
    cmd = compiler.visit_create_table(element, **kwargs)
    rtrees = [table.__tablename__ for table in RTree.__subclasses__()]
    if element.element.name in rtrees:
        # sqlite R-trees:
        table = element.element
        name = table.name
        pk = table.primary_key.columns.values()[0].name
        cmd = ' '.join(cmd.split())
        cmd = (cmd.replace('CREATE TABLE', 'CREATE VIRTUAL TABLE')
               .replace(name, name + ' USING RTREE')
               .replace(pk + ' INTEGER NOT NULL', pk + ' INTEGER PRIMARY KEY')
               .replace(', PRIMARY KEY ({})'.format(pk), ''))

    return cmd


# R-tree foreign keys not supported, use triggers
triggers = {
    'delete_eph_from_tree': '''
    CREATE TRIGGER IF NOT EXISTS delete_eph_from_tree BEFORE DELETE ON eph
    BEGIN
      DELETE FROM eph_tree WHERE ephid=old.ephid;
    END;
    ''',

    'delete_object_from_eph': '''
    CREATE TRIGGER IF NOT EXISTS delete_object_from_eph BEFORE DELETE ON obj
    BEGIN
      DELETE FROM eph WHERE objid=old.objid;
    END;
    ''',

    'delete_obs_from_obs_tree': '''
    CREATE TRIGGER IF NOT EXISTS delete_obs_from_obs_tree
    BEFORE DELETE ON obs
    BEGIN
      DELETE FROM obs_tree WHERE obsid=old.obsid;
    END;
    '''
}


def create(engine):
    con = engine.connect()
    Base.metadata.create_all(con)
    for trigger in triggers.values():
        con.execute(trigger)
    con.close()

# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""config
=========

The configuration is a JSON-formatted file.  See ``_config_example``.
Applications may add whatever variables to the configuration file.

"""

import os
import json

__all__ = ['Config']

_config_example = '''
{
  "database": "/path/to/sbsearch.db",
  "log": "/path/to/sbsearch.log",
  "location": "IAU observatory code",
  "obstable": "obs_table_name"
}
'''


class Config:
    """SBSearch configuration.

    Controls database location, log file location, object cutout and
    stack locations.  Parameters are stored as object keys:
    ``Config['user']``, ``Config['log']``, etc..


    Parameters
    ----------
    filename : string
        The file to load, or `None` to load the default location
        '~/.config/sbsearch.config'.

    **kwargs
      Additional or updated configuration parameters and values.


    Attributes
    ----------
    obs_table : Observation table name.
    obs_tree : Observation search tree name.
    found_table : Found object table name.

    """

    def __init__(self, filename=None, **kwargs):
        if filename is None:
            filename = os.path.expanduser('~/.config/sbsearch.config')

        with open(filename) as f:
            self.config = json.load(f)

        self.config.update(kwargs)

    def __getitem__(self, k):
        return self.config[k]

    @property
    def obs_table(self):
        return self.config['obstable']

    @property
    def obs_tree(self):
        return self.config['obstable'] + '_tree'

    @property
    def found_table(self):
        return self.config['obstable'] + '_found'

    @classmethod
    def from_args(cls, args):
        """Initialize from command-line arguments.

        Parameters
        ----------
        args : result from argparse.ArgumentParser.parse_args()
          Options checked: --config, --database, --log, --location,
          --obstable.

        Returns
        -------
        config : Config

        """

        updates = {}

        config_file = getattr(args, 'config', None)

        db = getattr(args, 'database', None)
        if db is not None:
            updates['database'] = db

        log = getattr(args, 'log', None)
        if log is not None:
            updates['log'] = log

        path = getattr(args, 'location', None)
        if path is not None:
            updates['location'] = path

        path = getattr(args, 'obstable', None)
        if path is not None:
            updates['obstable'] = path

        return cls(config_file, **updates)

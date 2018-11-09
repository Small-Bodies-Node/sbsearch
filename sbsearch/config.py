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
  "location": "IAU observatory code"
}
'''


class Config:
    """SBSearch configuration.

    Controls database location, log file location, object cutout and
    stack locations.  Parameters are stored as object keys:
    ``Config['user']``, ``Config['log']``, etc..


    Parameters
    ----------
    **kwargs
        Configuration parameters and values.

    """

    def __init__(self, **kwargs):
        self.config = {
            "database": ":memory:",
            "log": "/dev/null",
            "location": "500"
        }
        self.config.update(kwargs)

    def __getitem__(self, k):
        return self.config[k]

    @classmethod
    def from_args(cls, args):
        """Initialize from command-line arguments.

        Parameters
        ----------
        args : result from argparse.ArgumentParser.parse_args()
          Options checked: --config, --database, --log, --location

        read_defaults : bool, optional
            Set to ``True`` to read the default configuration file.
            If --config is define, this option is ignored.

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

        return cls.from_file(config_file, **updates)

    @classmethod
    def from_file(cls, filename=None, **kwargs):
        """Initialize from JSON-formatted file.

        Parameters
        ----------
        filename : string, optional
            Name of the file to read, or ``None`` for the default
            file: ~/.config/sbsearch.config .

        **kwargs
            Override saved parameters with these values.

        Returns
        -------
        config : Config

        """
        if filename is None:
            filename = os.path.expanduser('~/.config/sbsearch.config')

        with open(filename) as f:
            config = json.load(f)

        config.update(**kwargs)

        return cls(**config)

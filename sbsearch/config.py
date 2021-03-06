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
  "database": "postgresql://host/database",
  "log": "/path/to/sbsearch.log"
}
'''


class Config:
    """SBSearch configuration.

    Controls database location and log file location.  Parameters are stored
    as object keys: ``Config['database']``, ``Config['log']``.

    Parameters
    ----------
    **kwargs
        Configuration parameters and values.

    """

    # list of default files in order of precedence
    DEFAULT_FILES = ['sbsearch.cfg', '.sbsearch.cfg',
                     os.path.expanduser('~/.config/sbsearch.config')]

    def __init__(self, **kwargs):
        self.config = {
            "database": "postgresql:///sbsearch_test",
            "log": "/dev/null"
        }
        self.config.update(kwargs)

    def __getitem__(self, k):
        return self.config[k]

    @classmethod
    def from_args(cls, args, **updates):
        """Initialize from command-line arguments.

        The configuration file specified by --config (or the default,
        if ``None``) is read first.

        Parameters
        ----------
        args : result from argparse.ArgumentParser.parse_args()
          Options checked: --config for a configuration file,
          --option, where option is a configuration item, replacing
          spaces with underscores.

        **updates
            Any other configuration items.  However, `args` will take
            precedence.

        Returns
        -------
        config : Config

        """

        config_file = getattr(args, 'config', None)

        for k in _config_example:
            v = getattr(args, k.replace(' ', '_'), None)
            if v is not None:
                updates[k] = v

        return cls.from_file(config_file, **updates)

    @classmethod
    def from_file(cls, filename=None, **kwargs):
        """Initialize from JSON-formatted file.

        Parameters
        ----------
        filename : string, optional
            Name of the file to read, or ``None`` for the default
            file (in order of precedence):
                {}

        **kwargs
            Override saved parameters with these values.

        Returns
        -------
        config : Config

        """.format('\n                '.join(cls.DEFAULT_FILES))

        if filename is None:
            for filename in cls.DEFAULT_FILES:
                if os.path.exists(filename):
                    break

        try:
            with open(filename) as f:
                config = json.load(f)
        except IOError:
            print(_config_example)
            raise

        config.update(**kwargs)

        return cls(**config)

    def update(self, **kwargs):
        self.config.update(**kwargs)

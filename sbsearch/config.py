# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""config
=========

The configuration is a JSON-formatted file.  See ``_config_example``.
Applications may add whatever variables to the configuration file.

"""

import os
import json
import argparse
from typing import Dict, List, Optional, TypeVar, Union

__all__: List[str] = ['Config']

ConfigObject = TypeVar('ConfigObject', bound='Config')

_config_example: str = '''
{
  "database": "sqlite://test.db",
  "log": "/path/to/sbsearch.log",
  "min_edge_length": 3e-4,
  "uncertainty_ellipse": false,
  "padding": 0,
  "debug": false
}
'''


class Config:
    """SBSearch configuration.

    Controls database location, log file location, and ``SBSearch`` options.
    Parameters are stored as object keys: ``Config['user']``,
    ``Config['log']``, etc.


    Parameters
    ----------
    **kwargs
        Configuration parameters and values.

    """

    # list of default files in order of precedence
    DEFAULT_FILES: List[str] = [
        'sbsearch.cfg',
        '.sbsearch.cfg',
        os.path.expanduser('~/.config/sbsearch.config')
    ]

    DEFAULT_PARAMETERS: Dict[str, Union[str, float, int, bool]] = {
        "database": "sqlite://",
        "log": "/dev/null",
        "min_edge_length": 3e-4,
        "uncertainty_ellipse": False,
        "padding": 0,
        "debug": False
    }

    def __init__(self, **kwargs) -> None:
        self.config: Dict[str, Union[str, float, int, bool]] = (
            self.DEFAULT_PARAMETERS.copy()
        )
        self.update(kwargs)

    def __getitem__(self, k: str) -> Union[str, float, int, bool]:
        return self.config[k]

    @classmethod
    def find_default_file(cls) -> Union[str, None]:
        """Find the default configuration file, if any exists.

        Returns
        -------
        filename : string
            The found file or ``None``, if one could not be found.
            Searches the following locations in order of precedence:
                {}

        """.format('\n                '.join(cls.DEFAULT_FILES))

        filename: Union[str, None] = None
        for filename in cls.DEFAULT_FILES:
            if os.path.exists(filename):
                break
        else:
            filename = None
        return filename

    @classmethod
    def from_args(cls, args: argparse.Namespace, find_default: bool = False,
                  **updates) -> ConfigObject:
        """Initialize from command-line arguments.

        Parameters
        ----------
        args: argparse.Namespace
            For example, a result from argparse.ArgumentParser.parse_args().
            Options checked:
                - --config for a configuration file,
                - --option, where option is a configuration key (e.g., log,
                    min_edge_length), replacing spaces with underscores.

        find_default: bool, optional
            Attempt to find and read a default configuration file when no
            config file name is provided in the options.  See
            ``find_default_file``.

        **updates
            Any other configuration items.  However, `args` will take
            precedence.

        Returns
        -------
        config: Config

        """

        # config file treated separately
        config_file: Union[str, None] = getattr(args, 'config')
        if config_file is None:
            config_file = cls.find_default_file()

        k: str
        for k in cls.DEFAULT_PARAMETERS:
            v: Union[str, float, int, bool, None] = getattr(
                args, k.replace(' ', '_'), None)
            if v is not None:
                updates[k] = v

        if find_default:
            return cls.from_file(config_file, **updates)
        elif config_file is not None:
            return cls.from_file(config_file, **updates)
        else:
            return cls(**updates)

    @ classmethod
    def from_file(cls, filename: Optional[str] = None, **kwargs) -> ConfigObject:
        """Initialize from JSON-formatted file.

        Parameters
        ----------
        filename: string, optional
            Name of the file to read, or ``None`` for the default
            file ( in order of precedence):
                {}

        **kwargs
            Override saved parameters with these values.

        Returns
        -------
        config: Config

        """.format('\n                '.join(cls.DEFAULT_FILES))

        fn: str = filename
        if fn is None:
            for fn in cls.DEFAULT_FILES:
                if os.path.exists(fn):
                    break

        try:
            with open(fn) as f:
                config = json.load(f)
        except IOError:
            print(_config_example)
            raise

        config.update(**kwargs)

        return cls(**config)

    def update(self, kwargs) -> None:
        self.config.update(kwargs)

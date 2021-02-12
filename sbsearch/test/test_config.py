# Licensed with the 3-clause BSD license.  See LICENSE for details.
import argparse
import tempfile
import pytest

from ..config import Config


@pytest.fixture(name="config_file")
def fixture_config_file():
    with tempfile.NamedTemporaryFile('w+') as f:
        f.write('''{
  "database": "sqlite:///path/to/sbsearch.db",
  "log": "/path/to/sbsearch.log",
  "min_edge_length": 1e-4
}''')
        f.seek(0)
        yield f.name


class TestConfig:
    def test_init(self):
        config = Config(database='sqlite:///tmp/test.db')
        assert config['database'] == 'sqlite:///tmp/test.db'

    def test_from_args(self, config_file, monkeypatch):
        parser = argparse.ArgumentParser()
        parser.add_argument('--config')
        parser.add_argument('--min-edge-length', type=float)
        parser.add_argument('--database')
        parser.add_argument('--log')
        args = parser.parse_args(
            ['--database=sqlite:///tmp/other.db', '--log=alternate.log'])
        config = Config.from_args(args, find_default=False)
        assert config['database'] == args.database
        assert config['log'] == args.log
        assert config['min_edge_length'] == Config.DEFAULT_PARAMETERS['min_edge_length']

        # this time, read the config file
        args = parser.parse_args(
            ['--database=sqlite:///tmp/other.db', '--log=alternate.log',
             '--config=' + config_file])
        config = Config.from_args(args, find_default=False)
        assert config['database'] == args.database
        assert config['log'] == args.log
        assert config['min_edge_length'] == 1e-4

        # this time, no config file, but let Config discover a default
        monkeypatch.setattr(Config, "DEFAULT_FILES", [config_file])
        args = parser.parse_args(
            ['--database=sqlite:///tmp/other.db', '--log=alternate.log'])
        config = Config.from_args(args, find_default=True)
        assert config['database'] == args.database
        assert config['log'] == args.log
        assert config['min_edge_length'] == 1e-4

    def test_from_file(self, config_file):
        config = Config.from_file(config_file)
        assert config['database'] == 'sqlite:///path/to/sbsearch.db'
        assert config['min_edge_length'] == 1e-4
        assert config['log'] == '/path/to/sbsearch.log'

    def test_from_file_not_found(self):
        with pytest.raises(IOError):
            Config.from_file('/this_file_does_not_exist.config')

    def test_find_default_file(self, config_file, monkeypatch):
        monkeypatch.setattr(Config, "DEFAULT_FILES", [config_file])
        assert Config.find_default_file() == config_file
        monkeypatch.setattr(Config, "DEFAULT_FILES", [
                            '/this_file_does_not_exist.config'])
        assert Config.find_default_file() is None

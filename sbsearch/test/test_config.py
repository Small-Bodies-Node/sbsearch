# Licensed with the 3-clause BSD license.  See LICENSE for details.
import argparse
import tempfile
import pytest

from ..config import Config


@pytest.fixture
def config_file():
    with tempfile.NamedTemporaryFile('w+') as f:
        f.write('''{
  "database": "/path/to/sbsearch.db",
  "log": "/path/to/sbsearch.log",
  "location": "IAU observatory code"
}''')
        f.seek(0)
        yield f.name


class TestConfig:
    def test_init(self):
        config = Config(database=':memory:')
        assert config['database'] == ':memory:'

    def test_from_args(self, config_file):
        parser = argparse.ArgumentParser()
        parser.add_argument('--config')
        parser.add_argument('--location')
        parser.add_argument('--database')
        parser.add_argument('--log')
        args = parser.parse_args(
            ['--database=:memory:', '--location=I41', '--log=a.log',
             '--config=' + config_file])
        config = Config.from_args(args)
        assert config['database'] == ':memory:'
        assert config['location'] == 'I41'
        assert config['log'] == 'a.log'

    def test_from_file(self, config_file):
        config = Config.from_file(config_file)
        assert config['database'] == '/path/to/sbsearch.db'
        assert config['location'] == 'IAU observatory code'
        assert config['log'] == '/path/to/sbsearch.log'

    def test_default_file(self, config_file):
        DEFAULT_FILE = Config.DEFAULT_FILE
        Config.DEFAULT_FILE = config_file
        config = Config.from_file()
        assert config['database'] == '/path/to/sbsearch.db'
        assert config['location'] == 'IAU observatory code'
        assert config['log'] == '/path/to/sbsearch.log'

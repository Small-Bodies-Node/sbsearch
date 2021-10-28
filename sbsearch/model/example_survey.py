# Licensed with the 3-clause BSD license.  See LICENSE for details.

from typing import List
from sqlalchemy import Column, ForeignKey

# replace "from . import" with "from sbsearch.model import"
from .core import BigIntegerType, Base, Observation

__all__: List[str] = [
    'ExampleSurvey'
]


class ExampleSurvey(Observation):
    """Example survey data object.

    This class serves as an example data source.  It is used in sbsearch
    testing, but should not be used elsewhere.

    Copy this example to your code base and edit as needed.  See the catch tool
    for examples: github.com/Small-Bodies-Node/catch

    """
    __tablename__ = 'example_survey'
    __data_source_name__ = 'Example survey'
    __obscode__ = '500'  # MPC observatory code

    source_id = Column(BigIntegerType, primary_key=True)
    observation_id = Column(BigIntegerType,
                            ForeignKey('observation.observation_id',
                                       onupdate='CASCADE',
                                       ondelete='CASCADE'),
                            nullable=False,
                            index=True)

    # Add any additional attributes here, e.g.:
    # instrument = Column(String(64), doc='Instrument / detector name')

    __mapper_args__ = {
        'polymorphic_identity': 'example_survey'
    }

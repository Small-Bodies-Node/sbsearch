# Licensed with the 3-clause BSD license.  See LICENSE for details.

from typing import List
import sqlalchemy as sa
from sqlalchemy import Column, String, ForeignKey, Index

# replace "from . import" with "from sbsearch.model import"
from .core import BigIntegerType, Base, Observation

__all__: List[str] = [
    'ExampleSurvey'
]


class ExampleSurvey(Observation):
    """Example survey data object.

    This class serves as an example data source.  It is used in sbsearch
    testing, but should not be used elsewhere.

    Three objects must be defined to add a new data source:

        1. The observation data object (this class).
        2. The spatial term data object (see ExampleSurveySpatialTerm).
        3. The spatial term index (see ExampleSurveySpatialTermIndex).

    Copy each of the three examples to your code base and edit as needed.
    See the catch tool for examples: github.com/Small-Bodies-Node/catch

    """
    __tablename__ = 'example_survey'
    __data_source_name__ = 'Example survey'
    __obscode__ = '500'  # MPC observatory code

    id = Column(BigIntegerType, primary_key=True)
    observation_id = Column(BigIntegerType,
                            ForeignKey('observation.observation_id',
                                       onupdate='CASCADE',
                                       ondelete='CASCADE'),
                            nullable=False,
                            index=True)

    # Replace ExampleSurveySpatialTerm with an appropriate and unique name for
    # this data source, and define a new table based on the example below.
    terms = sa.orm.relationship("ExampleSurveySpatialTerm",
                                back_populates='source')

    # Add any additional attributes here, e.g.:
    # instrument = Column(String(64), doc='Instrument / detector name')

    __mapper_args__ = {
        'polymorphic_identity': 'example_survey'
    }


class ExampleSurveySpatialTerm(Base):
    """Example survey data source spatial term.

    The spatial term encodes the footprint of the observation on the sky.

    Replace instances of "example_survey" with the string used your data
    source's ``__tablename__``.

    """
    __tablename__ = 'example_survey_spatial_terms'
    term_id = Column(BigIntegerType, primary_key=True)
    source_id = Column(BigIntegerType,
                       ForeignKey('example_survey.id', onupdate='CASCADE',
                                  ondelete='CASCADE'),
                       nullable=False, index=True)
    term = Column(String(32), nullable=False)

    # Replace ExampleSurvey with the name of the data object for this source
    source = sa.orm.relationship("ExampleSurvey", back_populates="terms")

    def __repr__(self) -> str:
        return (f'<{self.__class__.__name__} term_id={self.term_id}'
                f' observation_id={self.source_id},'
                f' term={repr(self.term)}>')


# Define this index outside of the table to make it easier to drop/create.
# The naming scheme must be strictly followed: ix_(source table name)_spatial_terms
# Manually create with:
# CREATE INDEX ix_example_survey_spatial_terms ON example_survey_spatial_terms (term);
ExampleSurveySpatialTermIndex = Index(
    "ix_example_survey_spatial_terms",
    ExampleSurveySpatialTerm.term
)

# Licensed with the 3-clause BSD license.  See LICENSE for details.
class SBSException(Exception):
    pass


class BadObjectID(SBSException):
    pass


class NoEphemerisError(SBSException):
    pass


class SourceNotFoundError(SBSException):
    pass

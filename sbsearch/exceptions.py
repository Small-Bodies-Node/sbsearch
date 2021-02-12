# Licensed with the 3-clause BSD license.  See LICENSE for details.
class SBSException(Exception):
    pass


class ObjectError(SBSException):
    pass


class DesignationError(SBSException):
    pass


class PrimaryDesignationError(DesignationError):
    pass


class SecondaryDesignationError(DesignationError):
    pass


class DatabaseNotConnected(SBSException):
    pass


class UnknownSource(SBSException):
    pass

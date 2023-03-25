#ifndef SBS_EXCEPTIONS_H_
#define SBS_EXCEPTIONS_H_

#include <stdexcept>

namespace sbsearch
{
    class SBSException : public std::runtime_error
    {
    public:
        SBSException(const std::string &what_arg) : std::runtime_error(what_arg) {}
    };

    class MovingTargetError : public SBSException
    {
    public:
        MovingTargetError(const std::string &what_arg) : SBSException("Moving target error (" + what_arg + ")") {}
    };

    class EphemerisError : public SBSException
    {
    public:
        EphemerisError(const std::string &what_arg) : SBSException("Ephemeris error (" + what_arg + ")") {}
    };

    class ObservatoryError : public SBSException
    {
    public:
        ObservatoryError(const std::string &what_arg) : SBSException("Observatory error (" + what_arg + ")") {}
    };

}

#endif // SBS_EXCEPTIONS_H_
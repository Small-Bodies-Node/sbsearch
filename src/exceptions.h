#ifndef SBS_EXCEPTIONS_H_
#define SBS_EXCEPTIONS_H_

#include <exception>

namespace sbsearch
{
    class SBSException : public std::exception
    {
    };

    class EphemerisPaddingNotEnabled : public SBSException
    {
    };
}

#endif // SBS_EXCEPTIONS_H_
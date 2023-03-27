#ifndef SBS_FOUND_H_
#define SBS_FOUND_H_

#include "config.h"
#include "ephemeris.h"
#include "observation.h"

namespace sbsearch
{
    struct Found
    {
        Observation observation;
        Ephemeris ephemeris;

        Found(Observation o, Ephemeris e) : observation(o), ephemeris(e){};

        const bool operator==(const Found &other) const
        {
            return ((observation == other.observation) & (ephemeris == other.ephemeris));
        }

        const bool operator!=(const Found &other) const
        {
            return !(*this == other);
        }
    };
}

#endif // SBS_FOUND_H_

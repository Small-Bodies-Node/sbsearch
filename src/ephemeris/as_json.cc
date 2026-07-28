#include "config.h"

#include "ephemeris.h"

namespace sbsearch
{
    json::array Ephemeris::as_json()
    {
        json::array data_array;
        for (Datum datum : data())
            data_array.emplace_back(datum.as_json());
        return data_array;
    }
}
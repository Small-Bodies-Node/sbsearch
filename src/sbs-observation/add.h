#ifndef SBS_OBSERVATION_ADD_H_
#define SBS_OBSERVATION_ADD_H_

#include <iostream>
#include <utility>
#include <boost/json.hpp>

#include "arguments.h"
#include "../csv.h"
#include "../observation.h"
#include "../sbsearch.h"

namespace sbsearch::sbs_observation
{
    // Convert JSON object to Observation
    Observation tag_invoke(json::value_to_tag<Observation>, json::value const &jv);

    // Get Observation object from input stream in CSV format.  The first time
    // this function is used, the CSV header line will be read in and saved.
    // Sets failbit on is if a line could not be converted to an observation.
    std::istream &operator>>(std::istream &is, Observation &observation);

    // add observations from a JSON stream
    template <typename DB>
    void add_from_json(const Arguments &args, SBSearch<DB> &sbs, std::istream *input);

    // add observations from a CSV stream
    template <typename DB>
    void add_from_csv(const Arguments &args, SBSearch<DB> &sbs, std::istream *input);

    // add observations to the database
    template <typename DB>
    void add(const Arguments &args, SBSearch<DB> &sbs);
}

#endif

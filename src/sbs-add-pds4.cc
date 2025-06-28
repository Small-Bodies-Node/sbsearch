#include <cstdlib>
#include <iostream>
#include <set>
#include <string>
#include <libxml++/libxml++.h>

#include "config.h"
#include "env.h"
#include "logging.h"
#include "observation.h"
#include "sbsearch/sbsearch.h"
#include "sbsdb/postgresql.h"

#define TESTING false

using namespace sbsearch;
using std::cerr;
using std::cout;
using std::set;
using std::string;

class PDS4Observation : public Observation
{
public:
    friend std::istream &operator>>(std::istream &is, PDS4Observation &observation)
    {
        string line;
        char band;
        string visit, lldec, llra, lrdec, lrra, uldec, ulra, urdec, urra;
        double mjd_mid, exptime, mjd_start;

        if (std::getline(is, line))
        {
            std::istringstream iss(line);
            iss >> band >> visit >> mjd_mid >> exptime >> lldec >> llra >> lrdec >> lrra >> mjd_start >> uldec >> ulra >> urdec >> urra;
        }
        observation.source("lsst-dp0.2");
        observation.observatory("X05");
        observation.product_id(visit);
        observation.mjd_start(mjd_start);
        observation.mjd_stop(mjd_start + exptime / 86400);
        observation.fov(llra + ":" + lldec + "," +
                        lrra + ":" + lrdec + "," +
                        urra + ":" + urdec + "," +
                        ulra + ":" + uldec);

        return is;
    }
};

template <typename DB>
void add(SBSearch<DB> &sbs, const std::set<string> &inventory, std::istream &input)
{
    string line, fov;
    ProgressTriangle progress;

    Observations observations;
    observations.data.reserve(10000);

    sbs.drop_observations_indices();

    // skip the header
    for (int i = 0; i < 7; i++)
        std::getline(input, line);

    std::istream_iterator<PDS4Observation> start(input), end;
    while (start != end)
    {
        size_t count = 0;
        while (start != end && count < 10000)
        {
            observations.append(*start);
            start++;
            count++;
        }

        if (!TESTING)
        {
            sbs.add_observations(observations);
        }
        progress.update(count);
        observations.data.clear();
    }

    if (TESTING)
        Logger::info() << "Processed " << progress.count() << " observations." << std::endl;
    else
        Logger::info() << "Added " << progress.count() << " observations." << std::endl;

    sbs.create_observations_indices();
}

const set<string> &get_inventory(std::string_view collection)
{
}

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        cerr << R"(
Usage: sbs-add-pds4 <collection> <path>

Add observations from PDS4 labels.

<collection>    A PDS4 collection label.  Only files in this collection
                will be added.

<path>          Directory within which to search for label files
                (*xml). Nested directories are not recursively
                searched.

)" << std::endl;
        return 1;
    }

    if (!ENV.database | !ENV.log_file)
    {
        cerr << "Missing SBS_DATABASE and/or SBS_LOG_FILE environment variables.\n";
        return 1;
    }

    string collection(argv[1]);
    string path(argv[2]);

    try
    {
        Logger::info() << "sbs-add-pds4" << std::endl;

        SBSearch<sbsdb::Postgresql> sbs(ENV.database.value(), {.log_file = ENV.log_file.value()});

        Logger::info() << "Reading PDS4 labels from " << collection << " and " << path << std::endl;

        set<string> inventory = get_inventory(collection);
        add(sbs, inventory, path);
    }
    catch (std::exception &e)
    {
        cerr << "\nError: " << e.what() << "\n";
        return 1;
    }
}
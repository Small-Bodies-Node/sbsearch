#include <iostream>
#include <string>

#include "config.h"
#include "logging.h"
#include "observation.h"
#include "sbsearch.h"
#include "sbsdb/postgresql.h"

#define TESTING false

using namespace sbsearch;
using std::cerr;
using std::cout;
using std::string;

class LSSTObservation : public Observation
{
public:
    friend std::istream &operator>>(std::istream &is, LSSTObservation &observation)
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

/* add observations from an LSST DP0.2 table output

\col.band.arraySize = *
\col.band.type = char
\
|band|ccdVisitId|expMidptMJD  |expTime|llcdec     |llcra     |lrcdec     |lrcra     |obsStartMJD  |ulcdec     |ulcra     |urcdec     |urcra     |
|char|long      |double       |double |double     |double    |double     |double    |double       |double     |double    |double     |double    |
|    |          |d            |s      |deg        |deg       |deg        |deg       |             |deg        |deg       |deg        |deg       |
|    |          |             |       |           |          |           |          |             |           |          |           |          |
 g     178142072 59817.3294062    30.0 -44.6348735 49.8294117 -44.5055577 50.0828991 59817.3292326 -44.4511016 49.6455085 -44.3221413 49.8986892
 g     178142075 59817.3294062    30.0  -44.498233 50.0970101 -44.3682984 50.3492978 59817.3292326 -44.3148262  49.912809 -44.1852584 50.1647986
 g     178142076 59817.3294062    30.0 -44.3079991 49.9058018 -44.1784619  50.157799 59817.3292326 -44.1241311 49.7226748 -43.9949586 49.9743502

*/
template <typename DB>
void add(SBSearch<DB> &sbs, std::istream &input)
{
    string line;
    ProgressTriangle progress;

    Observations observations;
    observations.data.reserve(10000);

    sbs.drop_observations_indices();

    // skip the header
    for (int i = 0; i < 7; i++)
        std::getline(input, line);

    std::istream_iterator<LSSTObservation> start(input), end;
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

int main(int argc, char *argv[])
{
    string filename;
    if (argc > 1)
        filename = string(argv[1]);
    else
    {
        cerr << "Usage: sbs-add-lsst <filename>" << std::endl;
        return 1;
    }

    try
    {
        // SBSearch sbs("sqlite3://lsst.db");
        SBSearch<sbsdb::Postgresql> sbs("postgres:///lsst");
        Logger::info() << "sbs-add-lsst" << std::endl;

        Logger::info() << "Reading observations from " << filename << std::endl;
        std::ifstream input(filename);
        if (!input)
            throw std::runtime_error("Error opening file: " + filename);
        add(sbs, input);
    }
    catch (std::exception &e)
    {
        cerr << "\nError: " << e.what() << "\n";
        return 1;
    }
}
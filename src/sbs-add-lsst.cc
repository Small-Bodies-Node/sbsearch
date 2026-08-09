/*

Adds LSST Data Preview 1 images to an sbsearch database.

The input is a CSV formatted table from the Rubin Science Platform.  Query the
CCDVisit table, and get all rows (about 16,000).  Expected column names:

ccdVisitId,expTime,llcdec,llcra,lrcdec,lrcra,ulcdec,ulcra,urcdec,urcra,obsStartMJD

*/

#include <iostream>
#include <string>
#include "sofa/sofa.h"

#include "csv.h"
#include "logging.h"
#include "observation.h"
#include "progress_widgets.h"
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
    friend CsvStream &operator>>(CsvStream &csv, LSSTObservation &observation)
    {
        observation = LSSTObservation();

        // nothing to do
        if (!csv.good())
            return csv;

        vector<string> cells;
        csv >> cells;

        // skip lines without data
        if (cells.size() == 0)
            return csv >> observation;

        auto columns = csv.columns();
        try
        {
            // ccdVisitId,expTime,llcdec,llcra,lrcdec,lrcra,ulcdec,ulcra,urcdec,urcra,obsStartMJD
            observation.source("lsst dp1");
            observation.observatory("X05");
            observation.product_id(cells[columns["ccdVisitId"]]);

            // convert TAI to UTC
            double mjd_start_tai = std::stod(cells[columns["obsStartMJD"]]);
            double mjd_start, utc0;
            iauTaiutc(2400000.5, mjd_start_tai, &utc0, &mjd_start);

            observation.mjd_start(mjd_start);

            double exptime = std::stod(cells[columns["expTime"]]);
            observation.mjd_stop(mjd_start + exptime / 86400.);

            std::stringstream fov;
            fov << cells[columns["llcra"]] << ":" << cells[columns["llcdec"]] << ","
                << cells[columns["lrcra"]] << ":" << cells[columns["lrcdec"]] << ","
                << cells[columns["urcra"]] << ":" << cells[columns["urcdec"]] << ","
                << cells[columns["ulcra"]] << ":" << cells[columns["ulcdec"]];
            observation.fov(fov.str());

            observation.observation_id({});
        }
        catch (std::invalid_argument &exc)
        {
            throw SBSException("Invalid Observation data on CSV line " + std::to_string(csv.line()) + ": " + exc.what());
        }

        return csv;
    }
};

template <typename DB>
void add(SBSearch<DB> &sbs, std::istream &input)
{
    string line;
    ProgressTriangle progress;

    sbs.drop_observations_indices();

    Observations observations;
    const int batch_size = 10000;
    observations.data.reserve(batch_size);

    CsvStream csv(input);

    // peek first to set eof as needed
    while (csv.peek() && csv.good())
    {
        observations.data.clear();

        int count = 0;
        while (csv.peek() && csv.good() && count < batch_size)
        {
            LSSTObservation obs;
            csv >> obs;
            observations.append(obs);
            count++;
        }

        if (!TESTING)
        {
            sbs.add_observations(observations);
        }

        progress.update(observations.size());
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
        SBSearch<sbsdb::Postgresql> sbs("postgres:///lsst", {"lsst.log", LogLevel::INFO, true});
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

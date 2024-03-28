#include <iostream>
#include <string>

#include "config.h"
#include "logging.h"
#include "observation.h"
#include "sbsearch.h"
#include "util.h"

#define TESTING false

using namespace sbsearch;
using std::cerr;
using std::cout;
using std::string;

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
void add(SBSearch &sbs, std::istream &input)
{
    string line, fov;
    ProgressTriangle progress;

    Observations observations;
    observations.reserve(10000);

    Logger::info() << "Dropping observations indices." << std::endl;
    sbs.drop_observations_indices();

    // skip the header
    for (int i = 0; i < 7; i++)
        std::getline(input, line);

    while (input)
    {
        std::getline(input, line);
        if (line.size() < 145)
            continue;

        double mjd_start = std::stod(line.substr(85, 13));

        fov = line.substr(51, 10) + ":" + line.substr(39, 11) + "," +
              line.substr(74, 10) + ":" + line.substr(62, 11) + "," +
              line.substr(134, 10) + ":" + line.substr(122, 11) + "," +
              line.substr(111, 10) + ":" + line.substr(99, 11);

        observations.push_back(
            Observation(
                "lsst-dp0.2",
                "X05",
                line.substr(6, 10),
                mjd_start,
                mjd_start + std::stod(line.substr(31, 7)) / 86400,
                fov));

        if (observations.size() == 10000)
        {
            if (!TESTING)
            {
                sbs.add_observations(observations);
            }
            progress.update(10000);
            observations.clear();
        }
    }

    if (!TESTING)
    {
        sbs.add_observations(observations);
    }
    else
    {
        observations[0].format.show_fov = true;
        cout << observations << std::endl;
    }
    progress.update(observations.size());
    observations.clear();

    if (TESTING)
        cout << "Processed " << progress.count() << " observations.\n\n";
    else
        cout << "Added " << progress.count() << " observations.\n\n";

    Logger::info() << "Building observations indices." << std::endl;
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
        SBSearch sbs(SBSearch::sqlite3, "lsst.db");
        Logger::info() << "sbs-add-lsst" << std::endl;

        Logger::info() << "Reading observations from " << filename << std::endl;
        std::ifstream input(filename);
        if (!input)
            throw std::runtime_error("Error opening file: " + filename);
        add(sbs, input);
        input.close();
    }
    catch (std::exception &e)
    {
        cerr << "\nError: " << e.what() << "\n";
        return 1;
    }
}
// Configure sbsearch database
#include "config.h"

#include <iostream>
#include <string>
#include <getopt.h>

#include "indexer.h"
#include "sbsearch.h"

using sbsearch::Indexer;
using sbsearch::SBSearch;
using std::cout;
using std::string;

Indexer::Options parse_arguments(int argc, char **argv, Indexer::Options options)
{
    int c;

    while (1)
    {
        static struct option long_options[] =
            {
                {"max-spatial-cells", required_argument, 0, 'c'},
                {"max-spatial-resolution", required_argument, 0, 'M'},
                {"min-spatial-resolution", required_argument, 0, 'm'},
                {"temporal-resolution", required_argument, 0, 't'},
                {0, 0, 0, 0}};
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long(argc, argv, "c:x:n:t:",
                        long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {
        case 'c':
            options.max_spatial_cells(std::stoi(optarg));
            break;

        case 'M':
            options.max_spatial_resolution(std::stod(optarg));
            break;

        case 'm':
            options.min_spatial_resolution(std::stod(optarg));
            break;

        case 't':
            options.temporal_resolution(std::stoi(optarg));
            break;

        case '?':
            /* getopt_long already printed an error message. */
            break;

        default:
            abort();
        }
    }
}

int main(int argc, char **argv)
{
    SBSearch sbs(SBSearch::sqlite3, "sbsearch_test.db");
    const Indexer::Options previous_options = sbs.indexer_options();
    Indexer::Options updated_options = parse_arguments(argc, argv, previous_options);

    cout << "\nPrevious index setup:"
         << "\n  Maximum spatial cells: " << previous_options.max_spatial_cells()
         << "\n  Maximum spatial resolution (deg) / level: "
         << previous_options.max_spatial_resolution() / DEG
         << " / " << previous_options.min_spatial_level()
         << "\n  Minimum spatial resolution (deg) / level: "
         << previous_options.min_spatial_resolution() / DEG
         << " / " << previous_options.max_spatial_level()
         << "\n  Temporal resolution (1/day): " << previous_options.temporal_resolution()
         << "\n\n"
         << "\nNew index setup:"
         << "\n  Maximum spatial cells: " << updated_options.max_spatial_cells()
         << "\n  Maximum spatial resolution (deg) / level: "
         << updated_options.max_spatial_resolution() / DEG
         << " / " << updated_options.min_spatial_level()
         << "\n  Minimum spatial resolution (deg) / level: "
         << updated_options.min_spatial_resolution() / DEG
         << " / " << updated_options.max_spatial_level()
         << "\n  Temporal resolution (1/day): " << updated_options.temporal_resolution()
         << "\n\n";

    sbs.reindex(updated_options);

    return 0;
}
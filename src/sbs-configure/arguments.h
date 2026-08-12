#ifndef SBS_CONFIGURE_ARGUMENTS_H_
#define SBS_CONFIGURE_ARGUMENTS_H_

#include "config.h"
#include "indexer.h"
#include "cli.h"

using namespace sbsearch;
using namespace sbsearch::cli;

namespace sbsearch::sbs_configure
{
    struct Arguments : CommonArguments
    {
        bool create;
        Indexer::Options indexer_options;
        bool reconfigured;
        bool drop_indices;
        bool reindex;
    };

    Arguments get_arguments(int argc,
                            char *argv[],
                            Indexer::Options current_options = Indexer::Options());
}

#endif
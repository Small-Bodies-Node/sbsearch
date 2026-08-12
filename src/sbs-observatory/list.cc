#include <iostream>
#include <string>
#include <string_view>

#include "config.h"
#include "logging.h"
#include "observatory.h"
#include "sbsearch.h"
#include "table.h"
#include "sbsdb/get.h"
#include "sbsdb/postgresql.h"
#include "sbs-observatory/list.h"

using namespace sbsearch;
using namespace sbsearch::table;
using std::cerr;
using std::cout;
using std::string;
using std::string_view;

namespace sbsearch::sbs_observatory
{
    template <typename DB>
    void list(string_view output_filename, SBSearch<DB> &sbs)
    {
        Observatories observatories = sbsdb::get::all_observatories(sbs.db());

        std::ostream *os;
        std::ofstream outf;

        if (output_filename.empty())
            os = &cout;
        else
        {
            outf.open(string(output_filename));
            os = &outf;
        }

        if (observatories.size() == 0)
            cout << "# No observatories in the database.\n";
        else
        {
            int N = observatories.size();
            vector<string> names;
            vector<double> lon, rho_cos_phi, rho_sin_phi;

            for (const auto &item : observatories)
            {
                names.push_back(item.first);
                lon.push_back(item.second.longitude);
                rho_cos_phi.push_back(item.second.rho_cos_phi);
                rho_sin_phi.push_back(item.second.rho_sin_phi);
            }

            Table table;
            table.add(Column("name", "%s", names));
            table.add(Column("longitude", "%.3lf", lon));
            table.add(Column("rho cos(phi)", "%.6lf", rho_cos_phi));
            table.add(Column("rho sin(phi)", "%.6lf", rho_sin_phi));
            *os << table;
        }
    }

    template void list(string_view, SBSearch<sbsdb::Postgresql> &);
}

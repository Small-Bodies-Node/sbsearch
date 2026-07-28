#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <boost/json.hpp>

#include "./add.h"
#include "./arguments.h"
#include "./csv.h"
#include "logging.h"
#include "observation.h"
#include "sbsdb/postgresql.h"
#include "sbsearch.h"
#include "util/string.h"

using namespace sbsearch::cli;
using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::vector;

namespace sbsearch
{
    Observation tag_invoke(json::value_to_tag<Observation>, json::value const &jv)
    {
        json::object const &obj = jv.as_object();

        Observation obs(
            json::value_to<string>(obj.at("source")),
            json::value_to<string>(obj.at("observatory")),
            json::value_to<string>(obj.at("product_id")),
            json::value_to<double>(obj.at("mjd_start")),
            json::value_to<double>(obj.at("mjd_stop")),
            json::value_to<string>(obj.at("fov")));

        if (obj.contains("observation_id"))
        {
            auto value = obj.at("observation_id");
            if (!value.is_null())
                obs.observation_id(json::value_to<int64_t>(value));
        }

        return obs;
    }
}

namespace sbsearch::sbs_observation
{
    template <typename DB>
    void add_from_json(const Arguments &args, SBSearch<DB> &sbs, std::istream *input)
    {
        sbsearch::ProgressTriangle progress;

        Observations observations;
        observations.data.reserve(args.batch_size);

        json::stream_parser parser;
        boost::system::error_code error;
        parser.reset();

        size_t buffered = 0;
        while (!input->eof())
        {
            char buffer[4096];
            input->read(buffer, sizeof(buffer));

            size_t n = parser.write_some(buffer, input->gcount(), error);
            buffered += n;
            if (error)
                throw std::runtime_error("Error parsing JSON:" + error.message());

            int remainder = input->gcount() - n; // number of unconsumed characters

            // when JSON objects are ready, process them
            while (parser.done())
            {
                json::value data = parser.release();
                parser.reset();
                buffered = 0;

                Observations observations(json::value_to<vector<Observation>>(data));
                if (!args.noop)
                {
                    if (args.action == "update")
                        sbs.update_observations(observations);
                    else
                        sbs.add_observations(observations);
                }
                progress += observations.size();

                if (remainder)
                {
                    int m = parser.write_some(buffer + n, remainder);
                    remainder -= m; // update unconsumed characters
                    buffered += m;
                    n += m;
                }
            }
        }

        if (buffered > 0)
            // input is done, and parser was done parsing objects, but there is
            // still data in the buffer, this is an error:
            throw std::runtime_error("Processing complete, but " + std::to_string(buffered) + " bytes remain in the buffer.");

        string action = (args.action == "update") ? "Updated " : "Added ";
        cout << (args.noop ? "Read in " : action) << progress.count() << " observations.\n\n";
    }

    template <typename DB>
    void add_from_csv(const Arguments &args, SBSearch<DB> &sbs, std::istream *input)
    {
        string line;
        ProgressTriangle progress;

        Observations observations;
        observations.data.reserve(args.batch_size);

        CsvStream csv(*input);
        // peek first to set eof as needed
        while (csv.peek() && csv.good())
        {
            observations.data.clear();

            int count = 0;
            while (csv.peek() && csv.good() && count < args.batch_size)
            {
                Observation obs;
                csv >> obs;
                observations.append(obs);
                count++;
            }

            if (!args.noop)
            {
                if (args.action == "update")
                    sbs.update_observations(observations);
                else
                    sbs.add_observations(observations);
            }

            progress += observations.size();
        }

        string action = (args.action == "update") ? "Updated " : "Added ";
        cout << (args.noop ? "Read in " : action) << progress.count() << " observations.\n\n";
    }

    template <typename DB>
    void add(const Arguments &args, SBSearch<DB> &sbs)
    {
        if (args.drop_indices)
        {
            Logger::info() << "Dropping observations indices." << endl;
            sbs.db()->drop_observations_indices();
        }

        OutputFormat file_format = args.file_format;
        if (file_format == AUTO)
        {
            const int len = args.file.size();
            if ((len >= 5) && (args.file.compare(len - 5, 5, ".json") == 0))
                file_format = JSON;
            else if ((len >= 4) && (args.file.compare(len - 4, 4, ".csv") == 0))
                file_format = CSV;
        }

        std::istream *input;
        std::ifstream inf;
        if (args.file == "-")
        {
            cout << "\nReading observations from stdin.\n";
            input = &std::cin;
        }
        else
        {
            cout << "Reading observations from \"" + args.file + "\".\n";
            inf.open(args.file);
            if (!inf)
                throw std::runtime_error("Error opening file: " + args.file);
            input = &inf;
        }

        if (file_format == JSON)
            add_from_json(args, sbs, input);
        else
            add_from_csv(args, sbs, input);

        if (args.drop_indices)
        {
            Logger::info() << "Building observations indices." << endl;
            sbs.db()->create_observations_indices();
        }
    }

    template void add_from_json(const Arguments &, SBSearch<sbsearch::sbsdb::Postgresql> &, std::istream *);
    template void add_from_csv(const Arguments &, SBSearch<sbsearch::sbsdb::Postgresql> &, std::istream *);
    template void add(const Arguments &, SBSearch<sbsearch::sbsdb::Postgresql> &);
}
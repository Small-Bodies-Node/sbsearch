#include <iostream>
#include <iterator>
#include <map>
#include <string>
#include <vector>
#include <boost/json.hpp>

#include "add.h"
#include "arguments.h"
#include "../logging.h"
#include "../observation.h"
#include "../sbsdb/postgresql.h"
#include "../sbsearch.h"
#include "../util/string.h"

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

    std::istream &operator>>(std::istream &is, Observation &observation)
    {
        static std::map<string, int> columns;
        static int64_t line = 0;

        if (!is.good())
            return is;

        // get next set of cells
        vector<string> cells = sbsearch::util::get_csv_cells(is);
        line++;

        // skip lines without data
        if (cells.size() == 0)
            return is >> observation;

        // has the column header been read in yet?
        if (columns.size() == 0)
        {
            int index = 0;
            for (auto const &column : cells)
                columns[column] = index++;
            return is;
        }

        static bool has_observation_id = columns.find("observation_id") != columns.end();

        // it is an error to read partial or too much data
        if (cells.size() != columns.size())
        {
            throw SBSException("Invalid Observation data on CSV line " + std::to_string(line));
            is.setstate(std::ios_base::failbit);
            return is;
        }

        try
        {
            observation.source(cells[columns["source"]]);
            observation.observatory(cells[columns["observatory"]]);
            observation.product_id(cells[columns["product_id"]]);
            observation.mjd_start(std::stod(cells[columns["mjd_start"]]));
            observation.mjd_stop(std::stod(cells[columns["mjd_stop"]]));
            observation.fov(cells[columns["fov"]]);
            if (has_observation_id)
                observation.observation_id((int64_t)std::stoll(cells[columns["observation_id"]]));
        }
        catch (std::invalid_argument &exc)
        {
            for (auto const &[key, value] : columns)
                cerr << key << ": " << value << endl;
            cerr << util::join(cells, "/") << endl;
            throw SBSException("Invalid Observation data on CSV line " + std::to_string(line) + ": " + exc.what());
        }

        return is;
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
        while (input)
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
                    sbs.add_observations(observations);
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

        cout << "Added " << progress.count() << " observations.\n\n";
    }

    template <typename DB>
    void add_from_csv(const Arguments &args, SBSearch<DB> &sbs, std::istream *input)
    {
        string line;
        ProgressTriangle progress;

        Observations observations;
        observations.data.reserve(10000);

        std::istream_iterator<Observation> start(*input), end;
        start++; // the first line read was the column header
        while (start != end)
        {
            int count = 0;
            while (start != end && count < 10000)
            {
                observations.append(*start);
                count++;
                start++;
            }

            if (!args.noop)
                sbs.add_observations(observations);

            progress += observations.size();
        }

        cout << "Added " << progress.count() << " observations.\n\n";
    }

    template <typename DB>
    void add(const Arguments &args, SBSearch<DB> &sbs)
    {
        if (args.drop_indices)
        {
            Logger::info() << "Dropping observations indices." << endl;
            sbs.db()->drop_observations_indices();
        }

        FileFormat file_format = args.file_format;
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
#include <iostream>
#include <exception>
#include <map>
#include <string>
#include <vector>

#include "csv.h"
#include "exceptions.h"
#include "observation.h"
#include "util/string.h"

using sbsearch::util::strip;
using std::string;
using std::vector;

namespace sbsearch
{
    void CsvStream::init_columns()
    {
        // get next set of cells as the column header, skipping comments and
        // empty lines
        vector<string> header;
        do
        {
            header = get_cells();
        } while (good() && (header.size() == 0));

        int index = 0;
        for (auto const &column : header)
            columns_[column] = index++;
        has_observation_id = columns_.find("observation_id") != columns_.end();
        has_meta = columns_.find("meta") != columns_.end();
    }

    CsvStream &operator>>(CsvStream &csv, vector<string> &cells)
    {
        csv.cells_.clear();
        cells.clear();

        // if stream is not good, then nothing to do
        if (!csv.good())
            return csv;

        csv.cells_ = csv.get_cells();

        // skip lines without data
        if (csv.cells_.size() == 0)
            return csv >> cells;

        // it is an error to read partial or too much data
        if (csv.cells_.size() != csv.columns_.size())
        {
            throw SBSException("Invalid data on CSV line " + std::to_string(csv.line()));
            csv.setstate(std::ios_base::failbit);
        }

        cells = {csv.cells_.begin(), csv.cells_.end()};
        return csv;
    }

    CsvStream &operator>>(CsvStream &csv, Observation &observation)
    {
        observation = Observation();

        // nothing to do
        if (!csv.good())
            return csv;

        vector<string> cells;
        csv >> cells;

        // skip lines without data
        if (cells.size() == 0)
            return csv >> observation;

        try
        {
            observation.source(cells[csv.columns_["source"]]);
            observation.observatory(cells[csv.columns_["observatory"]]);
            observation.product_id(cells[csv.columns_["product_id"]]);
            observation.mjd_start(std::stod(cells[csv.columns_["mjd_start"]]));
            observation.mjd_stop(std::stod(cells[csv.columns_["mjd_stop"]]));
            observation.fov(cells[csv.columns_["fov"]]);

            if (csv.has_observation_id)
            {
                int64_t observation_id = (int64_t)std::stoll(cells[csv.columns_["observation_id"]]);
                if (observation_id > 0)
                    observation.observation_id(observation_id);
            }
            else
                observation.observation_id({});

            if (csv.has_meta)
            {
                string meta = cells[csv.columns_["meta"]];
                if (meta == "")
                    observation.meta({});
                else
                    observation.meta(std::move(meta));
            }
            else
                observation.meta({});
        }
        catch (std::invalid_argument &exc)
        {
            throw SBSException("Invalid Observation data on CSV line " + std::to_string(csv.line()) + ": " + exc.what());
        }

        return csv;
    }

    string CsvStream::get_cell()
    {
        string cell;
        cell.reserve(1025);

        // discard leading spaces
        while (peek() == ' ')
            get();

        // Begin cell processing: is this cell a quoted string?
        if (peek() == '"')
        {
            *this >> std::setw(1025) >> std::quoted(cell);

            // skip to , discarding extra characters but not \n
            while ((peek() != '\n') && good()) // peek first to set eof as needed
            {
                if (get() == ',')
                    break;
            }
        }
        else
        {
            // this cell is not a quoted string, get all characters up to ,
            // or \n
            while ((peek() != '\n') && good()) // peek first to set eof as needed
            {
                char next = get();
                if (next == ',')
                    break;
                cell.push_back(next);
            }

            // only strip spaces around an unquoted cell
            cell = string(strip(cell));
        }

        if (cell.size() > 1024)
            throw std::runtime_error("CSV cell length must be <=1024 characters.");

        return cell;
    }

    vector<string> CsvStream::get_cells()
    {
        vector<string> cells;

        // nothing to do
        if (!good())
            return cells;

        // skip blank lines
        if (peek() == '\n')
        {
            get();
            line_++;
            return cells;
        }

        // skip comment lines
        if (peek() == '#')
        {
            while ((peek() != '\n') && good()) // peek first to set eof as needed
                get();

            // consume trailing newline
            if (peek() == '\n')
                get();

            line_++;
            return cells;
        }

        while (good())
        {
            string cell = get_cell();
            cells.emplace_back(cell);

            // if the next character is a newline, then the line is done
            if (peek() == '\n')
                break;
        }

        line_++;

        // remove trailing newline so another line may be processed
        if (peek() == '\n')
            get();

        return cells;
    }
}
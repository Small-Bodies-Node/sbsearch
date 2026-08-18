#ifndef SBS_CSV_H_
#define SBS_CSV_H_

#include <iostream>
#include <map>
#include <string>
#include <vector>
#include "observation.h"

using std::string;
using std::vector;

namespace sbsearch
{
    /**
     * @brief Read CSV (character separated value) formatted stream into a
     * vector of strings or an SBSearch object.
     *
     * Supported objects:
     *   - Observation
     *
     * Blank lines and lines starting with # are skipped.  To read a line with a
     * single empty cell, use quotes, e.g., "".
     *
     * Strings may be "quoted."  Leading and trailing whitespace within a cell
     * is stripped unless within a quoted string.
     *
     * The maximum cell length is 1024 characters.
     *
     * Upon initialization, the first valid CSV line is read in as the column
     * header.  See the relevant `operator>>` functions for column definitions
     * and formats.  It is an error to have a mismatch between the number of
     * cells in a line and in the column header.
     *
     * std::fstream inf("observations.csv"); CsvStream csv(inf); Observation
     * obs; csv >> obs;
     */
    class CsvStream : public std::istream
    {
    public:
        // Constructor, automatically initializes the column definitions
        CsvStream(std::istream &is, const char delim = ',') : std::istream(is.rdbuf()), delimiter(delim)
        {
            init_columns();
        };

        // initialize the column definitions with the next valid CSV line
        void init_columns();

        // column definitions
        const std::map<string, int> &columns() { return columns_; };

        // number of lines read, including blank and commented lines
        size_t line() { return line_; };

        // last set of cells read and parsed
        const vector<string> &last_cells() { return cells_; };

        // Read a CSV line into a vector of strings.
        friend CsvStream &operator>>(CsvStream &csv, vector<string> &cells);

        // Read a CSV line and save the result to an Observation object.  Column
        // labels and example data:
        //
        // source,observatory,product_id,mjd_start,mjd_stop,fov,observation_id,meta
        // "Big Survey Project","I41","unique product ID",60000.00,60000.01,"0:0, 1:0, 1:1, 0:1",1,"{\"maglim\": 25.0}"
        //
        // * "fov" is a string of comma-separated RA:Dec pairs in units of degrees.
        // * "observation_id" is optional, values < 0 are treated as null.
        // * "meta" is optional.  Empty strings are treated as null.
        // * The column order is arbitrary.
        friend CsvStream &operator>>(CsvStream &csv, Observation &observation);

    private:
        // value delimiter
        char delimiter;

        // associate column definition with column index
        std::map<string, int> columns_;

        // last cells read
        vector<string> cells_;

        // number of lines read
        size_t line_ = 0;

        // optional columns
        bool has_observation_id = false;
        bool has_meta = false;

        // Get and parse the next cell of a CSV formatted stream.  Leading and
        // trailing whitespace is stripped unless within a quoted string.  The
        // maximum cell length is 1024 characters.  Returns cell data (without
        // quotes).
        string get_cell();

        // Get and parse the next line of a CSV formatted stream.  Blank lines
        // and lines starting with # are skipped.  Returns the number of lines
        // read and a vector of cells.
        std::vector<std::string> get_cells();
    };
}

#endif
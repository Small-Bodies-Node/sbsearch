#ifndef TABLE_H_
#define TABLE_H_

#define MAX_COLUMN_WIDTH 256

#include <algorithm>
#include <ostream>
#include <string>
#include <type_traits>
#include <vector>

#include "util.h"

using std::string;
using std::vector;

namespace sbsearch
{
    namespace table
    {
        class Table
        {
        public:
            typedef vector<string> Column;
            vector<Column> columns;
            const bool header;

            // Constructors.
            // Set h to true to print the headers.
            Table(const bool h) : header(h){};
            Table() : Table(true){};

            // Format data as a vector of strings, with or without the header,
            // and append to the table.  format is a printf-style string. Throws
            // std::range_error if column lengths do not agree.
            template <typename T>
            void add_column(const string name, const string format, const vector<T> &data);

            // The length of the table, including header rows, or -1 if there
            // are no columns.
            int length() const;

            // Get row i of the table.
            const string row(const size_t i) const;

            // Get all rows of the table.
            const vector<string> rows() const;

            friend std::ostream &operator<<(std::ostream &os, const Table &table);

        private:
            // Format a cell value as a string.
            template <typename T>
            const string format_cell(const string format, const T &value);
            const string format_cell(const string format, const bool &value);
            const string format_cell(const string format, const string &value);
        };

    }
}

#endif // TABLE_H_
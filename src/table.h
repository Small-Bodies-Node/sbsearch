#ifndef TABLE_H_
#define TABLE_H_

#define MAX_COLUMN_WIDTH 256

#include <algorithm>
#include <iostream>
#include <optional>
#include <ostream>
#include <string>
#include <type_traits>
#include <vector>

using std::optional;
using std::string;
using std::string_view;
using std::vector;

namespace sbsearch
{
    namespace table
    {
        class Column
        {
        public:
            // n is the column name, f is a printf-style string for the vector
            // d.  %s may be used to display bool as "true" or "false".
            template <typename T>
            Column(const string_view name,
                   const string_view format,
                   const vector<T> &data = {}) : format_(format)
            {
                cells.emplace_back(name);
                cells.emplace_back(""); // header divider
                for (auto const &datum : data)
                    cells.emplace_back(format_cell(datum));

                update_width();
            };

            // length of the column, including header
            string::size_type length() const { return cells.size(); };

            // get the ith cell; 0 and 1 are the header cells
            string_view operator[](const int i) const;

            // append a value to the vector
            template <typename T>
            void append(const T &value);

        private:
            string format_;
            vector<string> cells;

            string underline;

            unsigned short width = 0;

            // Format a cell value as a string.
            template <typename T>
            const string format_cell(const T &value);

            // Update the column width.
            void update_width();
        };

        class Table
        {
        public:
            vector<Column> columns;
            const bool header;

            // Constructors.
            // Set h to true to print the headers.
            Table(const bool h) : header(h) {};
            Table() : header(true) {};

            // Add a column.  Throws std::range_error if column lengths do not
            // agree.
            void add(Column column);

            // The length of the table, including header rows, or -1 if there
            // are no columns.
            int length() const;

            // Get row i of the table, 0 and 1 are always the header rows.
            string row(const size_t i) const;

            // Get all rows of the table, including header rows if header is true.
            vector<string> rows() const;

            friend std::ostream &operator<<(std::ostream &os, const Table &table);
        };

        /////////////////////////////////////////////////////////////////////////////////
        // implementation
        template <typename T>
        void Column::append(const T &value)
        {
            cells.push_back(format_cell(value));
            update_width();
        }

        template <typename T>
        const string Column::format_cell(const T &value)
        {
            // optional types without a value return "null"
            if constexpr ((std::is_same_v<T, optional<int>> == true) ||
                          (std::is_same_v<T, optional<int64_t>> == true) ||
                          (std::is_same_v<T, optional<double>> == true) ||
                          (std::is_same_v<T, optional<string>> == true) ||
                          (std::is_same_v<T, optional<bool>> == true))
                return value.has_value() ? format_cell(value.value()) : "null";

            char cell[MAX_COLUMN_WIDTH];

            if constexpr (std::is_same_v<T, bool> == true)
            {
                if (format_[format_.length() - 1] == 's')
                    return format_cell(value ? "true" : "false");

                sprintf(cell, format_.c_str(), value);
                return string(cell);
            }
            else if constexpr (std::is_same_v<T, string> == true)
            {
                std::cerr << value << std::endl;
                return format_cell(value.c_str());
            }

            sprintf(cell, format_.c_str(), value);
            return string(cell);
        }
    }
}

#endif // TABLE_H_
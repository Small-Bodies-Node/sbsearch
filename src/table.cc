#include <algorithm>
#include <ostream>
#include <string>
#include <type_traits>
#include <vector>

#include "table.h"
#include "util.h"

using std::string;
using std::vector;

namespace sbsearch
{
    namespace table
    {
        template <typename T>
        void Table::add_column(const string name, const string format, const vector<T> &data)
        {
            if ((length() >= 0) & ((data.size() + 2 * header) != length()))
                throw std::range_error("Refusing to create a table with inconsistent column lengths.");

            Column column;
            column.reserve(data.size() + 2);

            // first pass, get the string representation for every cell
            if (header)
            {
                column.push_back(name);
                column.push_back(string(name.size(), '-'));
            }

            for (auto value : data)
                column.push_back(format_cell(format, value));

            // no data and no header? then we're done
            if ((data.size() == 0) & !header)
                return;

            // second pass, get the maximum width of the column and resize
            // every cell to match
            unsigned short w = std::max_element(column.begin(),
                                                column.end(),
                                                [](const string &a, const string &b)
                                                { return a.size() < b.size(); })
                                   ->size();

            char fixed_width_format[10];
            sprintf(fixed_width_format, "%%%ds", w);
            for (size_t i = 0; i < column.size(); i++)
            {
                // special case for the header underline
                if ((i == 1) & header)
                {
                    column[i].resize(w, '-');
                    continue;
                }

                column[i] = format_cell(fixed_width_format, column[i]);
            }

            columns.push_back(column);
        }

        template void Table::add_column(const string, const string, const vector<bool> &);
        template void Table::add_column(const string, const string, const vector<int> &);
        template void Table::add_column(const string, const string, const vector<int64> &);
        template void Table::add_column(const string, const string, const vector<double> &);
        template void Table::add_column(const string, const string, const vector<string> &);

        int Table::length() const
        {
            if (columns.size() == 0)
                return -1;
            else
                return columns[0].size();
        }

        const string Table::row(const size_t i) const
        {
            vector<string> cells;
            cells.reserve(columns.size());

            for (const Column &col : columns)
                cells.push_back(col[i]);

            string r = join(cells, "  ");
            return r;
        }

        const vector<string> Table::rows() const
        {
            vector<string> r;

            // no data?  we're done
            if (columns.size() == 0)
                return r;

            r.reserve(length());
            for (int i = 0; i < length(); i++)
                r.push_back(row(i));
            return r;
        }

        std::ostream &operator<<(std::ostream &os, const Table &table)
        {
            for (const string &row : table.rows())
                os << row << "\n";
            return os;
        }

        template <typename T>
        const string Table::format_cell(const string format, const T &value)
        {
            char cell[MAX_COLUMN_WIDTH];
            sprintf(cell, format.c_str(), value);
            return string(cell);
        }

        template const string Table::format_cell(const string, const int &);
        template const string Table::format_cell(const string, const int64 &);
        template const string Table::format_cell(const string, const double &);

        const string Table::format_cell(const string format, const bool &value)
        {
            if (format[format.length() - 1] == 's')
                return format_cell(format, value ? "true" : "false");
            else
            {
                char cell[MAX_COLUMN_WIDTH];
                sprintf(cell, format.c_str(), value);
                return string(cell);
            }
        }

        const string Table::format_cell(const string format, const string &value)
        {
            return format_cell(format, value.c_str());
        }
    }
}

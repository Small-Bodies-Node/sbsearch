#include <algorithm>
#include <optional>
#include <ostream>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

#include "table.h"
#include "util/string.h"

using std::optional;
using std::string;
using std::string_view;
using std::vector;

namespace sbsearch
{
    namespace table
    {
        string_view Column::operator[](const int i) const
        {
            return cells.at(i);
        }

        void Column::update_width()
        {
            // get the maximum width of the column and resize
            // every cell to match
            unsigned short w = std::max_element(cells.begin(),
                                                cells.end(),
                                                [](string_view a, string_view b)
                                                { return a.size() < b.size(); })
                                   ->size();

            // if the width has not changed, then we're done!
            if (w == width)
                return;

            width = w;

            std::stringstream sstr;
            sstr << '%' << w << 's';
            string fixed_width_format(sstr.str());

            for (size_t i = 0; i < cells.size(); i++)
            {
                // special case for the header underline
                if (i == 1)
                {
                    cells[i].resize(w, '-');
                    continue;
                }

                char cell[MAX_COLUMN_WIDTH];
                sprintf(cell, fixed_width_format.c_str(), cells[i].c_str());
                cells[i] = {cell};
            }
        }

        void Table::add(Column column)
        {
            if ((columns.size() > 0) && (column.length() != columns[0].length()))
                throw std::range_error("Refusing to create a table with inconsistent column lengths.");
            columns.push_back(column);
        }

        int Table::length() const
        {
            if (columns.size() == 0)
                return -1;
            else
                return columns[0].length();
        }

        string Table::row(const size_t i) const
        {
            assert(i < length());

            auto column = columns.begin();
            std::stringstream s("");
            s << (*column)[i];
            std::advance(column, 1);

            while (column < columns.end())
            {
                s << "  " << (*column)[i];
                std::advance(column, 1);
            }

            return s.str();
        }

        vector<string> Table::rows() const
        {
            // no data?  we're done
            if (columns.size() == 0)
                return {};

            vector<string> rows;
            for (int i = 2 * !header; i < length(); i++)
                rows.push_back(row(i));

            return rows;
        }

        std::ostream &operator<<(std::ostream &os, const Table &table)
        {
            for (string_view row : table.rows())
                os << row << "\n";
            return os;
        }

    }
}

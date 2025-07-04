#include "config.h"

#include <optional>
#include <sstream>
#include <string>
#include <vector>
#include <gtest/gtest.h>

#include "table.h"

using std::optional;
using std::string;
using std::vector;

namespace sbsearch::table::testing
{
    TEST(TableTests, TestEmptyTable)
    {
        Table table;

        std::stringstream stream;
        stream << table;

        EXPECT_EQ(stream.str(), "");
    }

    TEST(TableTests, TestTableWithHeader)
    {
        vector<int> vi{0, 1, 2};
        vector<double> vd{0.0, 1.0, 2.000002};
        vector<string> vs{"a", "b", "csdf"};
        vector<optional<int64_t>> vo{1, std::nullopt, 3};

        Table table;
        table.add(Column("int", "%02d", vi));
        table.add(Column("double", "%lf", vd));
        table.add(Column("double.3", "%10.3f", vd));
        table.add(Column("string", "%s", vs));
        table.add(Column("string8", "%8s", vs));
        table.add(Column("optional", "%ld", vo));

        std::stringstream stream;
        stream << table;

        EXPECT_EQ(stream.str(),
                  "int    double    double.3  string   string8  optional\n"
                  "---  --------  ----------  ------  --------  --------\n"
                  " 00  0.000000       0.000       a         a         1\n"
                  " 01  1.000000       1.000       b         b      null\n"
                  " 02  2.000002       2.000    csdf      csdf         3\n");
    }

    TEST(TableTests, TestTableWithoutHeader)
    {
        vector<int> vi{0, 1, 2};
        vector<double> vd{0.0, 1.0, 2.000002};
        vector<string> vs{"a", "b", "csdf"};

        Table table(false);
        table.add(Column("int", "%02d", vi));
        table.add(Column("double", "%lf", vd));
        table.add(Column("double.3", "%10.3f", vd));
        table.add(Column("string", "%s", vs));
        table.add(Column("string8", "%8s", vs));

        std::stringstream stream;
        stream << table;

        EXPECT_EQ(stream.str(),
                  " 00  0.000000       0.000       a         a\n"
                  " 01  1.000000       1.000       b         b\n"
                  " 02  2.000002       2.000    csdf      csdf\n");
    }

    TEST(TableTests, TestTableWithoutData)
    {
        Table table;
        table.add(Column("int", "%02d", vector<int>()));
        table.add(Column("double", "%lf", vector<double>()));
        table.add(Column("double.3", "%10.3f", vector<double>()));
        table.add(Column("string", "%s", vector<string>()));
        table.add(Column("string8", "%8s", vector<string>()));

        std::stringstream stream;
        stream << table;

        EXPECT_EQ(stream.str(),
                  "int  double  double.3  string  string8\n"
                  "---  ------  --------  ------  -------\n");
    }

    TEST(TableTests, TestTableWithInconsistentColumnLengths)
    {
        vector<int> vi{0, 1, 2};
        vector<double> vd{0.0, 1.0};

        Table table;
        table.add(Column("int", "%02d", vi));
        EXPECT_THROW(table.add(Column("double", "%lf", vd)), std::range_error);
    }
}
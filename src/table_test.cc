#include "config.h"

#include <sstream>
#include <string>
#include <vector>
#include <gtest/gtest.h>

#include "table.h"

using std::string;
using std::vector;

namespace sbsearch
{
    namespace table
    {
        namespace testing
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

                Table table;
                table.add_column("int", "%02d", vi);
                table.add_column("double", "%lf", vd);
                table.add_column("double.3", "%10.3f", vd);
                table.add_column("string", "%s", vs);
                table.add_column("string8", "%8s", vs);

                std::stringstream stream;
                stream << table;

                EXPECT_EQ(stream.str(),
                          "int    double    double.3  string   string8\n"
                          "---  --------  ----------  ------  --------\n"
                          " 00  0.000000       0.000       a         a\n"
                          " 01  1.000000       1.000       b         b\n"
                          " 02  2.000002       2.000    csdf      csdf\n");
            }

            TEST(TableTests, TestTableWithoutHeader)
            {
                vector<int> vi{0, 1, 2};
                vector<double> vd{0.0, 1.0, 2.000002};
                vector<string> vs{"a", "b", "csdf"};

                Table table(false);
                table.add_column("int", "%02d", vi);
                table.add_column("double", "%lf", vd);
                table.add_column("double.3", "%10.3f", vd);
                table.add_column("string", "%s", vs);
                table.add_column("string8", "%8s", vs);

                std::stringstream stream;
                stream << table;

                EXPECT_EQ(stream.str(),
                          "00  0.000000       0.000     a         a\n"
                          "01  1.000000       1.000     b         b\n"
                          "02  2.000002       2.000  csdf      csdf\n");
            }

            TEST(TableTests, TestTableWithoutData)
            {
                Table table;
                table.add_column("int", "%02d", vector<int>());
                table.add_column("double", "%lf", vector<double>());
                table.add_column("double.3", "%10.3f", vector<double>());
                table.add_column("string", "%s", vector<string>());
                table.add_column("string8", "%8s", vector<string>());

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
                table.add_column("int", "%02d", vi);
                EXPECT_THROW(table.add_column("double", "%lf", vd), std::range_error);
            }
        }
    }
}
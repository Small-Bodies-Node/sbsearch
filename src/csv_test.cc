#include "config.h"

#include <sstream>
#include <string>
#include <gtest/gtest.h>

#include "csv.h"

using std::string;

namespace sbsearch::testing
{
    TEST(CsvTests, Read)
    {
        std::stringstream stream(R"(
# CSV file with test data
source,observatory,product_id,mjd_start,mjd_stop,fov
"Big Survey Project","I41","unique product ID1",60000.00,60000.01,"0:0, 1:0, 1:1, 0:1"
"Big Survey Project","I41","unique product ID2",60000.01,60000.02,"0:0, 1:0, 1:1, 0:1"
"Big Survey Project","I41","unique product ID3",60000.02,60000.03,"0:0, 1:0, 1:1, 0:1"
"Big Survey Project","I41","unique product ID4",60000.03,60000.04,"0:0, 1:0, 1:1, 0:1"
)");
        CsvStream csv(stream);
        Observation obs;

        csv >> obs;
        EXPECT_EQ(obs.source(), "Big Survey Project");
        EXPECT_EQ(obs.observatory(), "I41");
        EXPECT_EQ(obs.product_id(), "unique product ID1");
        EXPECT_EQ(obs.mjd_start(), 60000.00);
        EXPECT_EQ(obs.mjd_stop(), 60000.01);
        EXPECT_EQ(obs.fov(), "0:0, 1:0, 1:1, 0:1");

        csv >> obs;
        EXPECT_EQ(obs.source(), "Big Survey Project");
        EXPECT_EQ(obs.observatory(), "I41");
        EXPECT_EQ(obs.product_id(), "unique product ID2");
        EXPECT_EQ(obs.mjd_start(), 60000.01);
        EXPECT_EQ(obs.mjd_stop(), 60000.02);
        EXPECT_EQ(obs.fov(), "0:0, 1:0, 1:1, 0:1");

        csv >> obs;
        EXPECT_EQ(obs.source(), "Big Survey Project");
        EXPECT_EQ(obs.observatory(), "I41");
        EXPECT_EQ(obs.product_id(), "unique product ID3");
        EXPECT_EQ(obs.mjd_start(), 60000.02);
        EXPECT_EQ(obs.mjd_stop(), 60000.03);
        EXPECT_EQ(obs.fov(), "0:0, 1:0, 1:1, 0:1");

        csv >> obs;
        EXPECT_EQ(obs.source(), "Big Survey Project");
        EXPECT_EQ(obs.observatory(), "I41");
        EXPECT_EQ(obs.product_id(), "unique product ID4");
        EXPECT_EQ(obs.mjd_start(), 60000.03);
        EXPECT_EQ(obs.mjd_stop(), 60000.04);
        EXPECT_EQ(obs.fov(), "0:0, 1:0, 1:1, 0:1");

        csv >> obs;
        EXPECT_EQ(obs.source(), "");
        EXPECT_EQ(obs.observatory(), "");
        EXPECT_EQ(obs.product_id(), "");
        EXPECT_EQ(obs.mjd_start(), 0);
        EXPECT_EQ(obs.mjd_stop(), 0);
        EXPECT_EQ(obs.fov(), "");

        EXPECT_EQ(csv.good(), false);
        EXPECT_EQ(csv.fail(), true);
    }

    TEST(CsvTests, OptionalObservationID)
    {
        std::stringstream stream(R"(
# CSV file with test data and optional observation_id
source,observatory,product_id,mjd_start,mjd_stop,fov,observation_id
"Big Survey Project","I41","unique product ID1",60000.00,60000.01,"0:0, 1:0, 1:1, 0:1",1
"Big Survey Project","I41","unique product ID2",60000.01,60000.02,"0:0, 1:0, 1:1, 0:1",-1
)");
        CsvStream csv(stream);
        Observation obs;

        csv >> obs;
        EXPECT_EQ(obs.source(), "Big Survey Project");
        EXPECT_EQ(obs.observatory(), "I41");
        EXPECT_EQ(obs.product_id(), "unique product ID1");
        EXPECT_EQ(obs.mjd_start(), 60000.00);
        EXPECT_EQ(obs.mjd_stop(), 60000.01);
        EXPECT_EQ(obs.fov(), "0:0, 1:0, 1:1, 0:1");
        EXPECT_EQ(obs.observation_id(), 1);

        csv >> obs;
        EXPECT_EQ(obs.source(), "Big Survey Project");
        EXPECT_EQ(obs.observatory(), "I41");
        EXPECT_EQ(obs.product_id(), "unique product ID2");
        EXPECT_EQ(obs.mjd_start(), 60000.01);
        EXPECT_EQ(obs.mjd_stop(), 60000.02);
        EXPECT_EQ(obs.fov(), "0:0, 1:0, 1:1, 0:1");
        EXPECT_FALSE(obs.observation_id());
    }

    TEST(CsvTests, OptionalMeta)
    {
        std::stringstream stream(R"(
# CSV file with test data and optional observation_id
source,observatory,product_id,mjd_start,mjd_stop,fov,meta
"Big Survey Project","I41","unique product ID1",60000.00,60000.01,"0:0, 1:0, 1:1, 0:1","{\"maglim\": 25.1}"
"Big Survey Project","I41","unique product ID2",60000.01,60000.02,"0:0, 1:0, 1:1, 0:1",""
"Big Survey Project","I41","unique product ID3",60000.02,60000.03,"0:0, 1:0, 1:1, 0:1",,
)");
        CsvStream csv(stream);
        Observation obs;

        csv >> obs;
        EXPECT_EQ(obs.source(), "Big Survey Project");
        EXPECT_EQ(obs.observatory(), "I41");
        EXPECT_EQ(obs.product_id(), "unique product ID1");
        EXPECT_EQ(obs.mjd_start(), 60000.00);
        EXPECT_EQ(obs.mjd_stop(), 60000.01);
        EXPECT_EQ(obs.fov(), "0:0, 1:0, 1:1, 0:1");
        EXPECT_EQ(obs.meta(), R"({"maglim": 25.1})");

        csv >> obs;
        EXPECT_EQ(obs.source(), "Big Survey Project");
        EXPECT_EQ(obs.observatory(), "I41");
        EXPECT_EQ(obs.product_id(), "unique product ID2");
        EXPECT_EQ(obs.mjd_start(), 60000.01);
        EXPECT_EQ(obs.mjd_stop(), 60000.02);
        EXPECT_EQ(obs.fov(), "0:0, 1:0, 1:1, 0:1");
        EXPECT_FALSE(obs.meta());

        // trailing comma seems to work
        csv >> obs;
        EXPECT_EQ(obs.source(), "Big Survey Project");
        EXPECT_EQ(obs.observatory(), "I41");
        EXPECT_EQ(obs.product_id(), "unique product ID3");
        EXPECT_EQ(obs.mjd_start(), 60000.02);
        EXPECT_EQ(obs.mjd_stop(), 60000.03);
        EXPECT_EQ(obs.fov(), "0:0, 1:0, 1:1, 0:1");
        EXPECT_FALSE(obs.meta());
    }

    TEST(CsvTests, CellParsing)
    {
        std::stringstream stream(R"(
# Test cell parsing.  These lines should all be successful.
col1,col2,col3,col4
 a,bc , def ,ghijk
"asdf", "asdf","\"asdf\"" ," asdf "
a,b,"c
c",d
a,,c,d
)");

        CsvStream csv(stream);
        vector<string> cells;

        csv >> cells;
        EXPECT_EQ(cells, vector<string>({"a", "bc", "def", "ghijk"}));

        csv >> cells;
        EXPECT_EQ(cells, vector<string>({"asdf", "asdf", "\"asdf\"", " asdf "}));

        csv >> cells;
        EXPECT_EQ(cells, vector<string>({"a", "b", "c\nc", "d"}));

        csv >> cells;
        EXPECT_EQ(cells, vector<string>({"a", "", "c", "d"}));

        // no trailing newline after unquoted cell
        stream = std::stringstream("col1,col2,col3\naa,bb,cc");
        csv = CsvStream(stream);
        csv >> cells;
        EXPECT_EQ(cells, vector<string>({"aa", "bb", "cc"}));

        // no trailing newline after quoted cell
        stream = std::stringstream("col1,col2,col3\naa,bb,\"cc\"");
        csv = CsvStream(stream);
        csv >> cells;
        EXPECT_EQ(cells, vector<string>({"aa", "bb", "cc"}));

        // missing close quote at end of stream
        stream = std::stringstream("col1,col2,col3,col4\na,b,c,\"");
        csv = CsvStream(stream);
        csv >> cells;
        EXPECT_EQ(cells, vector<string>({"a", "b", "c", ""}));

        // too many characters in an unquoted cell
        stream = std::stringstream();
        stream << "A,B\na,";
        for (int i = 0; i < 1025; i++)
            stream << "b";
        csv = CsvStream(stream);
        EXPECT_THROW(csv >> cells, std::runtime_error);

        // too many characters in a quoted cell
        stream = std::stringstream();
        stream << "A,B\na,\"";
        for (int i = 0; i < 1025; i++)
            stream << "b";
        stream << '"';
        csv = CsvStream(stream);
        EXPECT_THROW(csv >> cells, std::runtime_error);
    }
}
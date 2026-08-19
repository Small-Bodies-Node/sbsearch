#include <cmath>
#include <sstream>
#include <string>
#include <exception>
#include <boost/json.hpp>
#include <s2/s2latlng.h>
#include <s2/s2latlng_rect.h>
#include <s2/s2metrics.h>
#include <s2/s2point.h>
#include <s2/s2polygon.h>
#include <s2/s2region_term_indexer.h>
#include <gtest/gtest.h>

#include "config.h"
#include "json.h"
#include "observation.h"
#include "polygons.h"
#include "sbsearch_testing.h"
#include "vertices.h"

namespace json = boost::json;

using sbsearch::Observation;
using std::string;
using std::vector;

namespace sbsearch::testing
{
    TEST(ObservationTests, ObservationInit)
    {
        Observation obs("test source", "G37", "product", 0, 0.1, "0:0, 0:1, 1:1");
        EXPECT_TRUE(obs.is_valid());
        EXPECT_FALSE(obs.observation_id());
        EXPECT_EQ(obs.mjd_start(), 0);
        EXPECT_EQ(obs.mjd_stop(), 0.1);

        obs = Observation("test source", "G37", "product", 1, 1.1, "0:0, 0:1, 1:1", {}, 1);
        EXPECT_EQ(obs.observation_id(), 1);

        // bad FOV: only two vertices
        EXPECT_THROW(Observation("test source", "G37", "product", 0, 0.1, "0:0, 0:1"), std::runtime_error);
        // bad FOV: not parsable as coordinates
        EXPECT_THROW(Observation("test source", "G37", "product", 0, 0.1, "asdf"), std::runtime_error);
        // stop is before start
        EXPECT_THROW(Observation("test source", "G37", "product", 0.1, 0, "0:0, 0:1, 1:1"), std::runtime_error);
    }

    TEST(ObservationTests, ObservationProperties)
    {
        EXPECT_TRUE(Observation() == Observation());

        Observation a("test source 1", "G37", "a", 1, 1.1, "0:0, 0:1, 1:1", {"asdf", "fdsa"}, 1);
        EXPECT_EQ(a.mjd_start(), 1);
        EXPECT_EQ(a.mjd_stop(), 1.1);
        EXPECT_EQ(a.mjd_mid(), 1.05);
        EXPECT_NEAR(a.exposure(), 8640, 1e-6);
        EXPECT_EQ(a.fov(), "0:0, 0:1, 1:1");
        EXPECT_EQ(a.terms(), vector<string>({"asdf", "fdsa"}));

        Observation b("test source 2", "G37", "b", 2, 2.1, "2:0, 2:1, 3:1", {"jkl;", ";lkj"}, 1);
        EXPECT_FALSE(a == b);
        EXPECT_TRUE(a != b);

        a.source("test source 2");
        EXPECT_FALSE(a == b);

        a.product_id("b");
        EXPECT_FALSE(a == b);

        a.mjd_start(2);
        EXPECT_FALSE(a == b);

        a.mjd_stop(2.1);
        EXPECT_FALSE(a == b);

        a.fov("2:0, 2:1, 3:1");
        EXPECT_TRUE(a == b); // terms are not compared
        EXPECT_FALSE(a != b);

        a.observation_id(2);
        EXPECT_FALSE(a == b);

        Observation c("test source", "G37", "c", 1, 1.1, "0:0, 0:1, 1:1");
        EXPECT_FALSE(c.observation_id());
        c.observation_id(2);
        EXPECT_EQ(c.observation_id(), 2);

        // update terms from a string
        a.terms("a b c");
        EXPECT_EQ(a.terms(), vector<string>({"a", "b", "c"}));
    }

    TEST(ObservationTests, ObservationTerms)
    {
        Observation obs("test source", "G37", "obs", 0, 0.1, "0:0, 0:1, 1:1", {"asdf", "fsda"});
        EXPECT_EQ(obs.terms(), vector<string>({"asdf", "fsda"}));
    }

    TEST(ObservationTests, ObservationToStream)
    {
        Observation obs("test source 2", "G37", "b", 2, 2.1, "2:0, 2:1, 3:1", {"jkl;", ";lkj"}, 1);

        std::stringstream stream;
        stream << obs;
        EXPECT_EQ(stream.str(), "1  \"test source 2\"  \"G37\"  \"b\"  2.000000  2.100000  8640");

        obs.format.show_fov = true;
        stream.str("");
        stream << obs;
        EXPECT_EQ(stream.str(), "1  \"test source 2\"  \"G37\"  \"b\"  2.000000  2.100000  8640  \"2:0, 2:1, 3:1\"");

        Observations observations({obs, obs});
        observations.format.show_fov = true;
        observations.format.date = Date::Format::Calendar;
        stream.str("");
        stream << observations;
        EXPECT_EQ(stream.str(),
                  "observation_id         source  product_id  observatory            mjd_start             mjd_stop  exposure            fov\n"
                  "--------------  -------------  ----------  -----------  -------------------  -------------------  --------  -------------\n"
                  "             1  test source 2           b          G37  1858-11-19 00:00:00  1858-11-19 02:24:00  8640.000  2:0, 2:1, 3:1\n"
                  "             1  test source 2           b          G37  1858-11-19 00:00:00  1858-11-19 02:24:00  8640.000  2:0, 2:1, 3:1\n");

        Observations no_observations;
        stream.str("");
        stream << no_observations;
        EXPECT_EQ(stream.str(),
                  "observation_id  source  product_id  observatory  mjd_start  mjd_stop  exposure\n"
                  "--------------  ------  ----------  -----------  ---------  --------  --------\n");

        no_observations.format.show_fov = true;
        stream.str("");
        stream << no_observations;
        EXPECT_EQ(stream.str(),
                  "observation_id  source  product_id  observatory  mjd_start  mjd_stop  exposure  fov\n"
                  "--------------  ------  ----------  -----------  ---------  --------  --------  ---\n");
    }

    TEST(ObservationTests, ObservationIsSameFov)
    {
        Observation a("test source", "G37", "a", 0, 0.1, "0:0, 0:1, 1:1");
        Observation b("test source", "G37", "b", 1, 1.1, "0:0, 0:1, 1:1");
        Observation c("test source", "G37", "b", 2, 2.1, "0.0:0.0, 0.0:1.0, 1.0:1.0");
        Observation d("test source", "G37", "b", 3, 3.1, "0:0, 0:1, 1:1.00001");
        EXPECT_TRUE(a != b);
        EXPECT_TRUE(a != c);
        EXPECT_TRUE(a != d);
        EXPECT_TRUE(a.is_same_fov(b));
        EXPECT_TRUE(a.is_same_fov(c));
        EXPECT_FALSE(a.is_same_fov(d));
    }

    TEST(ObservationTests, ObservationEquality)
    {
        Observation obs("test source", "G37", "a", 0, 0.1, "0:0, 0:1, 1:1", {}, 1);

        // same
        Observation other("test source", "G37", "a", 0, 0.1, "0:0, 0:1, 1:1", {}, 1);
        EXPECT_TRUE(obs == other);

        // different source
        other = Observation("another test source", "G37", "a", 0, 0.1, "0:0, 0:1, 1:1", {}, 1);
        EXPECT_FALSE(obs == other);

        // different product_id
        other = Observation("test source", "G37", "b", 0, 0.1, "0:0, 0:1, 1:1", {}, 1);
        EXPECT_FALSE(obs == other);

        // different mjd_start
        other = Observation("test source", "G37", "a", 0.05, 0.1, "0:0, 0:1, 1:1", {}, 1);
        EXPECT_FALSE(obs == other);

        // different mjd_stop
        other = Observation("test source", "G37", "a", 0, 0.15, "0:0, 0:1, 1:1", {}, 1);
        EXPECT_FALSE(obs == other);

        // different observation_id
        other = Observation("test source", "G37", "a", 0, 0.1, "0:0, 0:1, 1:1", {}, 2);
        EXPECT_FALSE(obs == other);

        // different terms
        other = Observation("test source", "G37", "a", 0, 0.1, "0:0, 0:1, 1:1", {"asdf"}, 1);
        EXPECT_TRUE(obs == other);

        // different FOV
        other = Observation("test source", "G37", "a", 0, 0.1, "0.05:0, 0:1, 1:1", {"asdf"}, 1);
        EXPECT_FALSE(obs == other);
    }

    TEST(ObservationTests, ObservationObservationId)
    {
        Observation obs("test source", "G37", "product", 0, 0.1, "0:0, 0:1, 1:1");
        obs.observation_id(1);
        EXPECT_EQ(obs.observation_id(), 1);
    }

    TEST(ObservationTests, ObservationAsPolygon)
    {
        Observation obs("test source", "G37", "product", 0, 1, "-1:-2,2:-2,2:2,-1:2");
        S2Polygon polygon;
        make_polygon(obs, polygon);
        S2Polygon expected;
        make_polygon(make_vertices("-1:-2,2:-2,2:2,-1:2"), expected);
        EXPECT_TRUE(polygon.Equals(expected));
    }

    TEST(ObservationTests, ObservationAsJSON)
    {
        Observation obs("test source", "G37", "product", 0, 1, "-1:-2,2:-2,2:2,-1:2");
        obs.format.show_fov = true;
        json::object obj(json::value_from(obs).as_object());
        EXPECT_EQ(obj["source"], "test source");
        EXPECT_EQ(obj["observatory"], "G37");
        EXPECT_EQ(obj["product_id"], "product");
        EXPECT_EQ(obj["mjd_start"], 0.);
        EXPECT_EQ(obj["mjd_stop"], 1.);
        EXPECT_EQ(obj["fov"], "-1:-2,2:-2,2:2,-1:2");
    }

    TEST(ObservationsTests, ObservationsGetters)
    {
        Observations observations({
            {"test source", "G37", "a", 0, 0.1, "0:0, 0:1, 1:1", {"1", "2"}, 10, "asdf", "none", 0},
            {"test source", "G37", "b", 1, 1.1, "0:0, 0:1, 1:1", {"2", "3"}, 11, "1a", "seeing 2", 1},
            {"test source", "G37", "c", 2, 2.2, "0.0:0.0, 0.0:1.0, 1.0:1.0", {"3", "4"}, 12, "2a", "bad", 2},
            {"another source", "F51", "d", 3, 3.1, "0:0, 0:1, 1:1.00001", {}, std::nullopt, std::nullopt, std::nullopt, 3},
        });

        EXPECT_EQ(observations.source(), vector<string>({"test source",
                                                         "test source",
                                                         "test source",
                                                         "another source"}));
        EXPECT_EQ(observations.observatory(), vector<string>({"G37", "G37", "G37", "F51"}));
        EXPECT_EQ(observations.product_id(), vector<string>({"a", "b", "c", "d"}));
        EXPECT_EQ(observations.observation_id(), vector<optional<int64_t>>({10, 11, 12, std::nullopt}));
        EXPECT_EQ(observations.mjd_start(), vector<double>({0, 1, 2, 3}));
        EXPECT_EQ(observations.mjd_mid(), vector<double>({0.05, 1.05, 2.1, 3.05}));
        EXPECT_EQ(observations.mjd_stop(), vector<double>({0.1, 1.1, 2.2, 3.1}));
        EXPECT_EQ(observations.exposure(), vector<double>({0.1 * 86400,
                                                           (1.1 - 1.0) * 86400,
                                                           (2.2 - 2.0) * 86400,
                                                           (3.1 - 3.0) * 86400}));
        EXPECT_EQ(observations.fov(), vector<string>({"0:0, 0:1, 1:1",
                                                      "0:0, 0:1, 1:1",
                                                      "0.0:0.0, 0.0:1.0, 1.0:1.0",
                                                      "0:0, 0:1, 1:1.00001"}));
        EXPECT_EQ(observations.center(), vector<optional<string>>({"asdf",
                                                                   "1a",
                                                                   "2a",
                                                                   std::nullopt}));
        EXPECT_EQ(observations.terms(), vector<Observation::Terms>({{"1", "2"},
                                                                    {"2", "3"},
                                                                    {"3", "4"},
                                                                    {}}));
        EXPECT_EQ(observations.meta(), vector<optional<string>>({"none",
                                                                 "seeing 2",
                                                                 "bad",
                                                                 std::nullopt}));
    }
}
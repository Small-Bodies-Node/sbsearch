#include "config.h"

#include <chrono>
#include <thread>
#include <string>
#include <boost/filesystem.hpp>
#include <gtest/gtest.h>

#include "date.h"
#include "ephemeris.h"
#include "files.h"
#include "horizons.h"
#include "moving_target.h"

using std::string;

namespace sbsearch
{
    namespace testing
    {
        TEST(HorizonsTests, TestFormatCommand)
        {
            string command;

            // Jupiter barycenter
            EXPECT_EQ(horizons::format_command("599", 0, false), "599");

            // Comet and ISO designations: expect NOFRAG and CAP
            EXPECT_EQ(horizons::format_command("2P"),
                      "DES=2P;NOFRAG;CAP;");

            EXPECT_EQ(horizons::format_command("2P", 60390),
                      "DES=2P;NOFRAG;CAP<2460390;");

            EXPECT_EQ(horizons::format_command("1I", 55390),
                      "DES=1I;NOFRAG;CAP<2455390;");

            EXPECT_EQ(horizons::format_command("3D", 55390),
                      "DES=3D;NOFRAG;CAP<2455390;");

            EXPECT_EQ(horizons::format_command("P/2001 YX127", 59990),
                      "DES=P/2001 YX127;NOFRAG;CAP<2459990;");

            EXPECT_EQ(horizons::format_command("C/1995 O1", 59990),
                      "DES=C/1995 O1;NOFRAG;CAP<2459990;");

            // asteroids
            EXPECT_EQ(horizons::format_command("AP"), "AP;"); // not a comet like 1P, etc.

            EXPECT_EQ(horizons::format_command("24"), "24;");

            EXPECT_EQ(horizons::format_command("europa"), "europa;");

            EXPECT_EQ(horizons::format_command("1999 JU3"), "DES=1999 JU3;");
        }

        TEST(HorizonsTests, TestFormatQuery)
        {
            string command = horizons::format_command("2P", Date("2024-01-01").mjd());
            EXPECT_EQ(
                horizons::format_query(command,
                                       "I41",
                                       Date("2024-01-01"),
                                       Date("2024-02-01"),
                                       "1d"),
                R"(
!$$SOF
MAKE_EPHEM=YES
COMMAND='DES=2P;NOFRAG;CAP<2460310;'
EPHEM_TYPE=OBSERVER
CENTER='I41'
START_TIME='2024-01-01'
STOP_TIME='2024-02-01'
STEP_SIZE='1d'
QUANTITIES='1,9,19,20,23,24,27,37,41'
REF_SYSTEM='ICRF'
CAL_FORMAT='JD'
CAL_TYPE='M'
TIME_DIGITS='MINUTES'
ANG_FORMAT='DEG'
APPARENT='AIRLESS'
RANGE_UNITS='AU'
SUPPRESS_RANGE_RATE='NO'
SKIP_DAYLT='NO'
SOLAR_ELONG='0,180'
EXTRA_PREC='YES'
R_T_S_ONLY='NO'
CSV_FORMAT='YES'
OBJ_DATA='YES'
)");
        }
    }

    TEST(HorizonsTests, TestQueryAndParse)
    {
        /////////////////////////////////////////
        // Test Ceres's approx position, no cache
        MovingTarget target("1");
        string command = horizons::format_command(target.designation());
        string parameters = horizons::format_query(command,
                                                   "I41",
                                                   Date("2005-07-01"),
                                                   Date("2005-07-02"),
                                                   "1d");
        fs::path fn = generate_cache_file_name(parameters);

        // clear previously cached data
        if (fs::exists(fn))
            fs::remove(fn);
        EXPECT_FALSE(fs::exists(fn));

        string table = horizons::query(parameters, false);
        Ephemeris eph(target, horizons::parse(table));

        // within a degree is fine for this test; reference values are from the Minor Planet Center
        EXPECT_EQ(eph.num_vertices(), 1);
        EXPECT_NEAR(eph.data(0).ra, 220.65, 1);
        EXPECT_NEAR(eph.data(0).dec, -10.52, 1);

        // verify the data was not cached
        EXPECT_FALSE(fs::exists(fn));

        ///////////////////////////////////
        // Run an Encke query with caching
        target = MovingTarget("2P");
        command = horizons::format_command(target.designation());
        parameters = horizons::format_query(command,
                                            "I41",
                                            Date("2005-07-01"),
                                            Date("2005-07-02"),
                                            "1d");

        // clear previously cached data
        fn = generate_cache_file_name(parameters);
        if (fs::exists(fn))
            fs::remove(fn);
        EXPECT_FALSE(fs::exists(fn));

        // query horizons and expect a new cache file
        table = horizons::query(parameters, true);
        eph = Ephemeris(target, horizons::parse(table));
        EXPECT_TRUE(fs::exists(fn));

        // the cached data has a timestamp, sleep for a bit so that a new query
        // would get a new timestamp
        std::this_thread::sleep_for(std::chrono::milliseconds(1500));

        // now, retrieve the cached data
        string cached_table = horizons::query(parameters, true);
        EXPECT_EQ(table, cached_table);

        // now, retrieve a fresh ephemeris
        string new_table = horizons::query(parameters, false);
        EXPECT_NE(table, new_table);
    }

    TEST(HorizonsTests, TestParse)
    {
        // missing $$EOE
        EXPECT_THROW(horizons::parse("$$SOE\n"), std::runtime_error);

        // missing $$SOE
        EXPECT_THROW(horizons::parse("$$EOE\n"), std::runtime_error);

        // a query for major body jupiter is ambiguous
        string parameters = horizons::format_query("jupiter",
                                                   "I41",
                                                   Date("2005-07-01"),
                                                   Date("2005-07-02"),
                                                   "1d");
        fs::path fn = generate_cache_file_name(parameters);

        // clear previously cached data
        if (fs::exists(fn))
            fs::remove(fn);
        EXPECT_FALSE(fs::exists(fn));

        string table = horizons::query(parameters, false);
        EXPECT_THROW(horizons::parse(table), std::runtime_error);
    }
}
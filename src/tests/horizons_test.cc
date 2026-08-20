#include "config.h"

#include <chrono>
#include <thread>
#include <string>
#include <boost/filesystem.hpp>
#include <gtest/gtest.h>

#include "date.h"
#include "files.h"
#include "horizons.h"
#include "moving_target.h"
#include "ephemeris/ephemeris.h"

using sbsearch::ephemeris::Ephemeris;
using std::string;

namespace sbsearch::testing
{
    TEST(HorizonsTests, FormatCommandAndParameters)
    {
        // Jupiter barycenter
        MovingTarget target("599", false);
        EXPECT_EQ(Horizons::format_command(target), "COMMAND='599'");

        // Comet and ISO designations: expect NOFRAG and CAP
        target = MovingTarget("2P");
        EXPECT_EQ(Horizons::format_command(target),
                  "COMMAND='DES=2P;NOFRAG;CAP;'");

        EXPECT_EQ(Horizons::format_command(target, 60390),
                  "COMMAND='DES=2P;NOFRAG;CAP<2460390;'");

        target = MovingTarget("1I");
        EXPECT_EQ(Horizons::format_command(target, 61000),
                  "COMMAND='DES=1I;NOFRAG;CAP<2461000;'");

        target = MovingTarget("3D");
        EXPECT_EQ(Horizons::format_command(target, 55390),
                  "COMMAND='DES=3D;NOFRAG;CAP<2455390;'");

        target = MovingTarget("P/2001 YX127");
        EXPECT_EQ(Horizons::format_command(target, 59990),
                  "COMMAND='DES=P/2001 YX127;NOFRAG;CAP<2459990;'");

        target = MovingTarget("C/1995 O1");
        EXPECT_EQ(Horizons::format_command(target, 59990),
                  "COMMAND='DES=C/1995 O1;NOFRAG;CAP<2459990;'");

        // asteroids
        target = MovingTarget("AP");
        EXPECT_EQ(Horizons::format_command(target), "COMMAND='AP;'"); // not a comet like 1P, etc.

        target = MovingTarget("24");
        EXPECT_EQ(Horizons::format_command(target), "COMMAND='24;'");

        target = MovingTarget("europa");
        EXPECT_EQ(Horizons::format_command(target), "COMMAND='europa;'");

        target = MovingTarget("1999 JU3");
        EXPECT_EQ(Horizons::format_command(target), "COMMAND='DES=1999 JU3;'");

        target = MovingTarget("2003 QD112");
        EXPECT_EQ(Horizons::format_command(target), "COMMAND='DES=2003 QD112;'");

        // orbits
        target.orbit(OrbitalElements{
            2450200.5l,            // epoch
            0.8241907231263196l,   // eccentricity
            0.532013766859137l,    // perihelion distance
            2450077.480966184235l, // perihelion time
            89.14262290335057l,    // longitude of the ascending node
            326.0591239257098l,    // argument of perihelion
            4.247821264821585l,    // inclination
            std::nullopt,          // mean anomaly
            std::nullopt,          // semi-major axis
            std::nullopt,          // mean motion
        });
        EXPECT_EQ(Horizons::format_command(target), "COMMAND=';'\nOBJECT='2003 QD112'\nEPOCH=2450200.5\nECLIP=J2000\nEC=0.8241907231263196\nOM=89.14262290335057\nW=326.0591239257098\nIN=4.247821264821585\nQR=0.532013766859137\nTP=2450077.480966184235");

        target.orbit(OrbitalElements{
            2450200.5l,          // epoch
            0.8241907231263196l, // eccentricity
            std::nullopt,        // perihelion distance
            std::nullopt,        // perihelion time
            89.14262290335057l,  // longitude of the ascending node
            326.0591239257098l,  // argument of perihelion
            4.247821264821585l,  // inclination
            1l,                  // mean anomaly
            2l,                  // semi-major axis
            std::nullopt,        // mean motion
        });
        EXPECT_EQ(Horizons::format_command(target), "COMMAND=';'\nOBJECT='2003 QD112'\nEPOCH=2450200.5\nECLIP=J2000\nEC=0.8241907231263196\nOM=89.14262290335057\nW=326.0591239257098\nIN=4.247821264821585\nMA=1\nA=2");

        target.orbit(OrbitalElements{
            2450200.5l,          // epoch
            0.8241907231263196l, // eccentricity
            std::nullopt,        // perihelion distance
            std::nullopt,        // perihelion time
            89.14262290335057l,  // longitude of the ascending node
            326.0591239257098l,  // argument of perihelion
            4.247821264821585l,  // inclination
            1l,                  // mean anomaly
            std::nullopt,        // semi-major axis
            0.012l,              // mean motion
        });
        EXPECT_EQ(Horizons::format_command(target), "COMMAND=';'\nOBJECT='2003 QD112'\nEPOCH=2450200.5\nECLIP=J2000\nEC=0.8241907231263196\nOM=89.14262290335057\nW=326.0591239257098\nIN=4.247821264821585\nMA=1\nN=0.012");
    }

    TEST(HorizonsTests, FormatQuery)
    {
        MovingTarget target("2P");
        string command = Horizons::format_command(target, Date("2024-01-01").mjd());
        EXPECT_EQ(
            Horizons::format_query(command,
                                   "I41",
                                   Date("2024-01-01"),
                                   Date("2024-02-01"),
                                   "1d"),
            R"(!$$SOF
MAKE_EPHEM=YES
COMMAND='DES=2P;NOFRAG;CAP<2460310;'
EPHEM_TYPE=OBSERVER
CENTER='I41'
START_TIME='2024-01-01 00:00:00'
STOP_TIME='2024-02-01 00:00:00'
STEP_SIZE='1d'
QUANTITIES='1,9,19,20,23,24,27,37,41,47'
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

    TEST(HorizonsTests, QueryAndParseRemote)
    {
        /////////////////////////////////////////
        // Test Ceres's approx position, no cache
        MovingTarget target("1");
        string center = "I41";
        Date start_date = Date("2005-07-01");
        Date stop_date = Date("2005-07-02");
        string time_step = "1d";

        Horizons horizons(target, center, start_date, stop_date, time_step, false);
        fs::path fn = generate_cache_file_name(horizons.parameters());

        // clear previously cached data
        if (fs::exists(fn))
            fs::remove(fn);
        EXPECT_FALSE(fs::exists(fn));

        // get the ephemeris data
        Ephemeris eph(target, horizons.get_ephemeris_data());

        // within a degree is fine for this test; reference values are from the Minor Planet Center
        EXPECT_EQ(eph.num_vertices(), 2);
        EXPECT_NEAR(eph.data(0).ra, 220.65, 1);
        EXPECT_NEAR(eph.data(0).dec, -10.52, 1);

        // verify the data was not cached
        EXPECT_FALSE(fs::exists(fn));

        ///////////////////////////////////
        // Run an Encke query with caching
        horizons.target(MovingTarget("2P"));
        horizons.cache(true);

        // clear previously cached data
        fn = generate_cache_file_name(horizons.parameters());
        if (fs::exists(fn))
            fs::remove(fn);
        EXPECT_FALSE(fs::exists(fn));

        // query horizons and expect a new cache file
        eph = Ephemeris(target, horizons.get_ephemeris_data());
        EXPECT_TRUE(fs::exists(fn));
        EXPECT_EQ(eph.num_vertices(), 2);

        // the cached data has a timestamp, sleep for a bit so that a new query
        // would get a new timestamp
        std::this_thread::sleep_for(std::chrono::milliseconds(1500));

        // now, retrieve the cached data
        string_view table = horizons.table();
        horizons.get_ephemeris_data();
        EXPECT_EQ(table, horizons.table());

        // now, retrieve a fresh ephemeris
        horizons.cache(false);
        horizons.get_ephemeris_data();
        string_view new_table = horizons.table();
        EXPECT_NE(table, new_table);

        // query with an orbit
        target = MovingTarget("Example");
        target.orbit(OrbitalElements{
            2450200.5l,            // epoch
            0.8241907231263196l,   // eccentricity
            0.532013766859137l,    // perihelion distance
            2450077.480966184235l, // perihelion time
            89.14262290335057l,    // longitude of the ascending node
            326.0591239257098l,    // argument of perihelion
            4.247821264821585l,    // inclination
            std::nullopt,          // mean anomaly
            std::nullopt,          // semi-major axis
            std::nullopt,          // mean motion
        });
        horizons.target(target);
        eph = Ephemeris(target, horizons.get_ephemeris_data());
        EXPECT_EQ(eph.data(1).ra, 250.568202312);
        EXPECT_EQ(eph.data(1).dec, -21.245005886);
    }

    TEST(HorizonsTests, ParseRemote)
    {
        // missing $$EOE
        EXPECT_THROW(Horizons::parse("$$SOE\n"), std::runtime_error);

        // missing $$SOE
        EXPECT_THROW(Horizons::parse("$$EOE\n"), std::runtime_error);

        // a query for major body jupiter is ambiguous
        string parameters = Horizons::format_query("COMMAND='jupiter'",
                                                   "I41",
                                                   Date("2005-07-01"),
                                                   Date("2005-07-02"),
                                                   "1d");
        EXPECT_THROW(Horizons::query(parameters, false), std::runtime_error);

        // A failed query should not be cached
        fs::path fn = generate_cache_file_name(parameters);
        EXPECT_FALSE(fs::exists(fn));
    }
}
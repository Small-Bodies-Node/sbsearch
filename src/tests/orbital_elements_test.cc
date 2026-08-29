#include <optional>
#include <string>
#include <boost/json.hpp>
#include <gtest/gtest.h>

#include "config.h"
#include "exceptions.h"
#include "json.h"
#include "orbital_elements.h"

using namespace sbsearch;
using std::optional;
using std::string;

namespace sbsearch::testing
{
    TEST(OrbitalElementsTests, IsComplete)
    {
        OrbitalElements a{
            2461200.5l,
            0.082l,
            std::nullopt,
            std::nullopt,
            128.075l,
            293.846l,
            2.401l,
            322.008l,
            3.108l,
            std::nullopt,
        };
        EXPECT_TRUE(a.is_complete());

        OrbitalElements b(a);
        b.ma = std::nullopt;
        EXPECT_FALSE(b.is_complete());

        b = a;
        b.a = std::nullopt;
        EXPECT_FALSE(b.is_complete());

        b.n = 0.180l;
        EXPECT_TRUE(b.is_complete());

        b.n = std::nullopt;
        b.ma = std::nullopt;
        b.Tp = 2461411.721l;
        EXPECT_FALSE(b.is_complete());

        b.qr = 2.853l;
        EXPECT_TRUE(b.is_complete());
    }

    TEST(OrbitalElementsTests, Comparison)
    {
        const OrbitalElements a{
            2461200.5l,
            0.082l,
            std::nullopt,
            std::nullopt,
            128.075l,
            293.846l,
            2.401l,
            322.008l,
            3.108l,
            std::nullopt,
        };

        using Member = long double OrbitalElements::*;
        Member required_parameters[] = {
            &OrbitalElements::epoch,
            &OrbitalElements::ec,
            &OrbitalElements::om,
            &OrbitalElements::w,
            &OrbitalElements::in,
        };

        for (Member parameter : required_parameters)
        {
            OrbitalElements b(a);
            EXPECT_EQ(a, b);
            b.*parameter = 5.;
            EXPECT_NE(a, b);
        }

        using OptionalMember = optional<long double> OrbitalElements::*;
        OptionalMember optional_parameters[] = {
            &OrbitalElements::qr,
            &OrbitalElements::ma,
            &OrbitalElements::a,
            &OrbitalElements::n,
        };

        for (OptionalMember parameter : optional_parameters)
        {
            OrbitalElements b(a);
            EXPECT_EQ(a, b);
            b.*parameter = 5.;
            EXPECT_NE(a, b);
        }
    }

    TEST(OrbitalElementsTests, Ostream)
    {
        std::stringstream stream;
        OrbitalElements orbit{
            2461200.5l,
            0.082l,
            2.853l,
            2461411.721l,
            128.075l,
            293.846l,
            2.401l,
            322.008l,
            3.108l,
            0.1234567890123456l,
        };
        stream << orbit;

        EXPECT_EQ(stream.str(),
                  "EPOCH=2461200.5\n"
                  "EC=0.082 QR=2.853 TP=2461411.721\n"
                  "OM=128.075 W=293.846 IN=2.401\n"
                  "MA=322.008 A=3.108 N=0.1234567890123456");
    }

    TEST(OrbitalElementsTests, Istream)
    {
        std::istringstream stream(R"(
# this is a comment, the next line is blank

actually, this is also ignored (no equal symbol)

# all caps
EPOCH=2461200.5

# lower case
ec=0.082

# commented parameter is ignored
# ec=1

# mixed case, and spaces around equals are OK
Qr = 2.853
Tp = 2461411.721

# spaces before parameter are OK
om = 128.075
 w = 293.846
in = 2.401

# lower-case null on optional parameters are OK
ma = null

# repeated parameters overwrite previous values
ma = 3.22008e+2

# multiple parameters separated by space are OK
a = 3.108 n = 1.80e-1
)");

        OrbitalElements a;
        stream >> a;
        OrbitalElements expected{
            2461200.5l,
            0.082l,
            2.853l,
            2461411.721l,
            128.075l,
            293.846l,
            2.401l,
            322.008l,
            3.108l,
            0.180l,
        };
        EXPECT_EQ(a, expected);

        // test aliases, new lines are optional
        stream = std::istringstream(R"(
epoch = 2461200.5
    e = 0.082      q = 2.853   Tp = 2461411.721
 node = 128.075 peri = 293.846  i = 2.401
    m = 322.008    a = 3.108    n = 0.180
)");

        a = OrbitalElements();
        EXPECT_NE(a, expected);
        stream >> a;
        EXPECT_EQ(a, expected);

        // test optional values and integers
        stream = std::istringstream(R"(
epoch = 2461200
   ec = 0.082
   om = 128.075
    w = 293.846
   in = 2.401
    m = null
)");
        a = OrbitalElements();
        stream >> a;
        EXPECT_EQ(a.epoch, 2461200l);
        EXPECT_FALSE(a.qr.has_value());
        EXPECT_FALSE(a.Tp.has_value());
        EXPECT_FALSE(a.ma.has_value());
        EXPECT_FALSE(a.a.has_value());
        EXPECT_FALSE(a.n.has_value());

        // error
        stream = std::istringstream(R"(
unknown = 1.0
)");
        EXPECT_THROW(stream >> a, OrbitError);
    }

    TEST(OrbitalElementsTests, JSON)
    {
        namespace json = boost::json;
        json::value jv = json::parse(R"(
{
  "epoch": "2461200.5",
  "ec": "0.082",
  "qr": "2.853",
  "Tp": "2461411.721",
  "om": "128.075",
  "w": "293.846",
  "in": "2.401",
  "ma": "322.008",
  "a": "3.108",
  "n": "0.180"
}
)");
        OrbitalElements orbit = OrbitalElements::from_json(jv.as_object());
        OrbitalElements expected{
            2461200.5l,
            0.082l,
            2.853l,
            2461411.721l,
            128.075l,
            293.846l,
            2.401l,
            322.008l,
            3.108l,
            0.180l,
        };
        EXPECT_EQ(orbit, expected);

        // repeat with alternate keys, name is ignored
        jv = json::parse(R"(
{
  "name": "6993",
  "epoch": "2461200.5",
  "e": "0.082",
  "q": "2.853",
  "Tp": "2461411.721",
  "node": "128.075",
  "peri": "293.846",
  "i": "2.401",
  "m": "322.008",
  "a": "3.108",
  "n": "0.180"
}
)");
        orbit = OrbitalElements::from_json(jv.as_object());
        EXPECT_EQ(orbit, expected);

        // repeat without nulls
        jv = json::parse(R"(
{
  "epoch": "2461200.5",
  "qr": null,
  "Tp": null,
  "e": "0.082",
  "node": "128.075",
  "peri": "293.846",
  "i": "2.401",
  "M": null
}
)");
        orbit = OrbitalElements::from_json(jv.as_object());
        expected = OrbitalElements{2461200.5l,
                                   0.082l,
                                   std::nullopt,
                                   std::nullopt,
                                   128.075l,
                                   293.846l,
                                   2.401l,
                                   std::nullopt,
                                   std::nullopt,
                                   std::nullopt};
        EXPECT_EQ(orbit, expected);

        // again, but without the optional keys
        jv = json::parse(R"(
{
  "epoch": "2461200.5",
  "e": "0.082",
  "node": "128.075",
  "peri": "293.846",
  "i": "2.401"
}
)");
        orbit = OrbitalElements::from_json(jv.as_object());
        expected = OrbitalElements{2461200.5l,
                                   0.082l,
                                   std::nullopt,
                                   std::nullopt,
                                   128.075l,
                                   293.846l,
                                   2.401l,
                                   std::nullopt,
                                   std::nullopt,
                                   std::nullopt};
        EXPECT_EQ(orbit, expected);
    }
}

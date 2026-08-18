#include <cmath>
#include <gtest/gtest.h>

#include "constants.h"
#include "util/spherical.h"

#define Deg(a) S1Angle::Degrees(a)

using std::ceil;
using std::floor;
using std::string;
using std::vector;

namespace sbsearch::util
{
    TEST(UtilSphericalTests, PositionAngle)
    {
        /*
        # verify with astropy
        >>> a = SkyCoord('0d 0d')
        >>> c = SkyCoord('0.01d 0.01d')
        >>> a.position_angle(c).deg
        44.999999563667686
        */
        const S2Point a = S2LatLng::FromDegrees(0, 0).ToPoint();
        const S2Point b = S2LatLng::FromDegrees(0, 0.01).ToPoint();
        const S2Point c = S2LatLng::FromDegrees(0.01, 0.01).ToPoint();
        const S2Point d = S2LatLng::FromDegrees(0.01, 0).ToPoint();
        const S2Point e = S2LatLng::FromDegrees(0, -0.01).ToPoint();
        const S2Point f = S2LatLng::FromDegrees(-0.01, 0).ToPoint();
        EXPECT_EQ(position_angle(a, b), 90 * DEG);
        EXPECT_EQ(position_angle(a, c), 44.999999563667686 * DEG);
        EXPECT_EQ(position_angle(a, d), 0);
        EXPECT_EQ(position_angle(a, e), -90 * DEG);
        EXPECT_EQ(position_angle(a, f), 180 * DEG);
    }

    // TEST(UtilSphericalTests, InnerAngle)
    // {
    //     EXPECT_NEAR(inner_angle(Deg(30), Deg(0)).radians(), Deg(30).radians(), 1e-8);
    //     EXPECT_NEAR(inner_angle(Deg(0), Deg(90)).radians(), Deg(90).radians(), 1e-8);
    //     EXPECT_NEAR(inner_angle(Deg(0), Deg(180)).radians(), Deg(180).radians(), 1e-8);
    //     EXPECT_NEAR(inner_angle(Deg(0), Deg(190)).radians(), Deg(170).radians(), 1e-8);
    //     EXPECT_NEAR(inner_angle(Deg(0), Deg(-270)).radians(), Deg(90).radians(), 1e-8);
    // }

    TEST(UtilSphericalTests, OffsetBy)
    {
        /*
        # verify with astropy
        >>> SkyCoord(0, 0, unit='rad').directional_offset_by(30 * u.deg, 2 * u.deg)
        <SkyCoord (ICRS): (ra, dec) in deg
            (1.00030471, 1.73196284)>
        >>> SkyCoord(0, 0, unit='rad').directional_offset_by(30 * u.deg, 2 * u.deg).ra.value
        1.0003047102322884
        >>> SkyCoord(0, 0, unit='rad').directional_offset_by(30 * u.deg, 2 * u.deg).dec.value
        1.7319628412802042
        */
        const S2LatLng coords = S2LatLng::FromDegrees(0, 0);

        const S2LatLng a = offset_by(coords, S1Angle::Degrees(0), S1Angle::Degrees(10));
        EXPECT_NEAR(a.lat().degrees(), 10, 1e-8);
        EXPECT_NEAR(a.lng().degrees(), 0, 1e-8);

        const S2LatLng b = offset_by(coords, S1Angle::Degrees(90), S1Angle::Degrees(10));
        EXPECT_NEAR(b.lat().degrees(), 0, 1e-8);
        EXPECT_NEAR(b.lng().degrees(), 10, 1e-8);

        const S2LatLng c = offset_by(coords, S1Angle::Degrees(180), S1Angle::Degrees(10));
        EXPECT_NEAR(c.lat().degrees(), -10, 1e-8);
        EXPECT_NEAR(c.lng().degrees(), 0, 1e-8);

        const S2LatLng d = offset_by(coords, S1Angle::Degrees(270), S1Angle::Degrees(10));
        EXPECT_NEAR(d.lat().degrees(), 0, 1e-8);
        EXPECT_NEAR(d.lng().degrees(), -10, 1e-8);

        const S2LatLng e = offset_by(coords, S1Angle::Degrees(-90), S1Angle::Degrees(10));
        EXPECT_NEAR(e.lat().degrees(), 0, 1e-8);
        EXPECT_NEAR(e.lng().degrees(), -10, 1e-8);

        const S2LatLng f = offset_by(coords, S1Angle::Degrees(30), S1Angle::Degrees(2));
        EXPECT_NEAR(f.lat().degrees(), 1.7319628412802042, 1e-8);
        EXPECT_NEAR(f.lng().degrees(), 1.0003047102322884, 1e-8);
    }

    TEST(UtilSphericalTests, Ellipse)
    {
        vector<S2LatLng> e = ellipse(4, S2LatLng::FromRadians(0, 0), 0.1, 0.05, 0);
        EXPECT_NEAR(e[0].lat().radians(), 0.1, 1e-8);
        EXPECT_NEAR(e[0].lng().radians(), 0, 1e-8);
        EXPECT_NEAR(e[1].lat().radians(), 0, 1e-8);
        EXPECT_NEAR(e[1].lng().radians(), 0.05, 1e-8);
        EXPECT_NEAR(e[2].lat().radians(), -0.1, 1e-8);
        EXPECT_NEAR(e[2].lng().radians(), 0, 1e-8);
        EXPECT_NEAR(e[3].lat().radians(), 0, 1e-8);
        EXPECT_NEAR(e[3].lng().radians(), -0.05, 1e-8);

        e = ellipse(4, S2LatLng::FromRadians(0, 0), 0.1, 0.05, PI_2);
        EXPECT_NEAR(e[0].lat().radians(), 0, 1e-8);
        EXPECT_NEAR(e[0].lng().radians(), 0.1, 1e-8);
        EXPECT_NEAR(e[1].lat().radians(), -0.05, 1e-8);
        EXPECT_NEAR(e[1].lng().radians(), 0, 1e-8);
        EXPECT_NEAR(e[2].lat().radians(), 0, 1e-8);
        EXPECT_NEAR(e[2].lng().radians(), -0.1, 1e-8);
        EXPECT_NEAR(e[3].lat().radians(), 0.05, 1e-8);
        EXPECT_NEAR(e[3].lng().radians(), 0, 1e-8);
    }
}
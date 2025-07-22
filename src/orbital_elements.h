#ifndef SBS_ORBITAL_ELEMENTS_H_
#define SBS_ORBITAL_ELEMENTS_H_

#include <optional>
#include <string>

using std::optional;
using std::string;

namespace sbsearch
{
    // Tied to Horizons format for heliocentric ecliptic elements.
    //
    //     EPOCH ..  Julian Day Number of osculating elements (TDB timescale, JDTDB)
    //     EC .....  Eccentricity
    //     QR .....  Perihelion distance (au)
    //     TP .....  Perihelion Julian Day Number
    //     OM .....  Longitude of ascending node (DEGREES) wrt IAU76/80 ecliptic
    //     W ......  Argument of perihelion (DEGREES) with respect to IAU76/80 ecliptic
    //     IN .....  Inclination (DEGREES) with respect to IAU76/80 ecliptic
    //     MA .....  Mean anomaly (DEGREES)
    //     A ......  Semi-major axis (au)
    //     N ......  Mean motion (DEG/DAY)
    //
    // Specify epoch, ec, om, w, in and one of {TP, QR}, {MA, A} or {MA, N}.
    //
    // Must be in J2000 reference frame: IAU76/80 J2000 obliquity of 84381.448
    // arcsec relative to the ICRF.
    struct OrbitalElements
    {
        long double epoch;        // Julian day, TDB
        long double ec;           // dimensionless
        optional<long double> qr; // au
        optional<long double> Tp; // Julian day
        long double om;           // degrees
        long double w;            // degrees
        long double in;           // degrees
        optional<long double> ma; // degrees
        optional<long double> a;  // au
        optional<long double> n;  // degrees/day
    };

    // Test for equality between two OrbitalElements objects.
    bool operator==(const OrbitalElements &a, const OrbitalElements &b);

    // Read orbital elements from a stream.  Each line is a single key-value
    // pair, separated by "=".  Key names are the same as Horizons format, case
    // insensitive.
    std::istream &operator>>(std::istream &is, OrbitalElements &orbit);
}

#endif

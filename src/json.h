#ifndef SBS_TO_JSON_H_
#define SBS_TO_JSON_H_

#include <initializer_list>
#include <string_view>
#include <boost/json.hpp>

#include "config.h"
#include "found.h"
#include "moving_target.h"
#include "observation.h"
#include "orbital_elements.h"
#include "query_info.h"
#include "ephemeris/ephemeris.h"

using sbsearch::ephemeris::Ephemeris;
using std::string_view;

namespace boost::json
{
    /** Observation to JSON object.
     *
     * Observation.format options are respected.
     *
     * {
     *   "source": string,
     *   "observatory": string,
     *   "product_id": string,
     *   "observation_id": number or null,
     *   "mjd_start": number,
     *   "mjd_stop": number,
     *   "fov": string,
     *   "center": string or null,
     *   "terms": array of strings,
     *   "meta": string or null,
     *   "mjd_added": number,
     * }
     *
     */
    void tag_invoke(const value_from_tag &, value &jv, const sbsearch::Observation &obs);

    /** Observations to JSON array.
     *
     * Observations.format options are respected.
     */
    void tag_invoke(const value_from_tag &, value &jv, const sbsearch::Observations &observations);

    /** OrbitalElements to JSON object.
     *
     * All values are converted to strings to preserve long double precision.
     *
     * {
     *   "description": string,
     *   "epoch": string,
     *   "ec": string,
     *   "qr": string or null,
     *   "Tp": string or null,
     *   "om": string,
     *   "w": string,
     *   "in": string,
     *   "ma": string or null,
     *   "a": string or null,
     *   "n": string or null,
     * }
     */
    void tag_invoke(const value_from_tag &, value &jv, const sbsearch::OrbitalElements &orbit);

    /** MovingTarget to JSON object.
     *
     * {
     *   "designation": string,
     *   "alternate_names": array of strings,
     *   "small_body": boolean,
     *   "orbit": OrbitalElements or null,
     *   "moving_target_id": number or null,
     * }
     */
    void tag_invoke(const value_from_tag &, value &jv, const sbsearch::MovingTarget &target);

    /** Ephemeris::Datum to JSON object.
     *
     * {
     *   "mjd": number or null,
     *   "tmtp": number or null,
     *   "ra": number or null,
     *   "dec": number or null,
     *   "mu": number or null,
     *   "mu_theta": number or null,
     *   "unc_a": number or null,
     *   "unc_b": number or null,
     *   "unc_theta": number or null,
     *   "rh": number or null,
     *   "delta": number or null,
     *   "phase": number or null,
     *   "selong": number or null,
     *   "true_anomaly": number or null,
     *   "sangle": number or null,
     *   "vangle": number or null,
     *   "vmag": number or null,
     * }
     */
    void tag_invoke(const value_from_tag &, value &jv, const Ephemeris::Datum &datum);

    /** Ephemeris::Data to JSON array. */
    void tag_invoke(const value_from_tag &, value &jv, const Ephemeris::Data &data);

    /** Ephemeris to JSON object.
     *
     * {
     *   "target": MovingTarget,
     *   "data": Ephemeris::Data,
     * }
     *
     */
    void tag_invoke(const value_from_tag &, value &jv, const Ephemeris &eph);

    /** Found to JSON object.
     *
     * All Observation, Ephemeris::target(), and Ephemeris::Data keys are used.
     */
    void tag_invoke(const value_from_tag &, value &jv, const sbsearch::Found &found);

    /** Founds to JSON array. */
    void tag_invoke(const value_from_tag &, value &jv, const sbsearch::Founds &founds);

    /** QueryInfo::Polygon to JSON array. */
    void tag_invoke(const value_from_tag &, value &jv, const sbsearch::QueryInfo::Polygon &polygon);
}

namespace sbsearch::json
{
    // Get a string value from obj[key] and convert it to a long double.
    //
    // Optionally try other keys if the first is not present.
    //
    // Throws sbsearch::KeyError if the keys are not present.  ValueError if a
    // key is present but not a string.
    long double get_string_as_long_double(boost::json::object &obj,
                                          std::initializer_list<string_view> keys);

    // Get a string or null value from obj[key] and convert it to an
    // optional<long double>.
    //
    // Optionally try other keys if the first is not present.
    //
    // Throws sbsearch::KeyError if the keys are not present.
    optional<long double> get_string_as_optional_long_double(boost::json::object &obj,
                                                             std::initializer_list<string_view> keys);
}

#endif
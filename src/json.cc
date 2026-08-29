#include <initializer_list>
#include <string>
#include <string_view>
#include <boost/json.hpp>

#include "config.h"
#include "exceptions.h"
#include "found.h"
#include "moving_target.h"
#include "observation.h"
#include "orbital_elements.h"
#include "query_info.h"
#include "ephemeris/ephemeris.h"
#include "ephemeris/interpolate.h"

using namespace sbsearch;
using sbsearch::ephemeris::Ephemeris;
using std::string;
using std::string_view;

namespace boost::json
{
    void tag_invoke(const value_from_tag &, value &jv, const Observation &obs)
    {
        auto format_date = [&obs](double mjd)
        { return obs.format.date == Date::Format::MJD
                     ? value_from(mjd)
                     : value_from(Date(mjd).iso()); };

        jv = {
            {"source", obs.source()},
            {"observatory", obs.observatory()},
            {"product_id", obs.product_id()},
            {"observation_id", value_from(obs.observation_id())},
            {"date_start", format_date(obs.mjd_start())},
            {"date_stop", format_date(obs.mjd_stop())},
        };

        if (obs.format.show_fov)
            jv.as_object().emplace("fov", obs.fov());

        jv.as_object().emplace("center", value_from(obs.center()));
        jv.as_object().emplace("terms", value_from(obs.terms()));

        if (obs.format.show_meta)
            jv.as_object().emplace("meta", value_from(obs.meta()));

        jv.as_object().emplace("date_added", format_date(obs.mjd_added()));
    }

    void tag_invoke(const value_from_tag &, value &jv, const Observations &observations)
    {
        jv.emplace_array();
        for (sbsearch::Observation obs : observations)
        {
            obs.format = observations.format;
            jv.as_array().emplace_back(value_from(obs));
        }
    }

    void tag_invoke(const value_from_tag &, value &jv, const OrbitalElements &orbit)
    {
        auto to_optional_string_value = [](optional<long double> x)
        {
            return value_from(x.has_value()
                                  ? std::make_optional(std::to_string(*x))
                                  : std::nullopt);
        };

        jv = {
            {"description", "Heliocentric ecliptic elements in J2000 reference frame, IAU76/80 obliquity.  Angles in degrees, distances in au, rates as per day, dates in TDB."},
            {"epoch", std::to_string(orbit.epoch)},
            {"ec", std::to_string(orbit.ec)},
            {"qr", to_optional_string_value(orbit.qr)},
            {"Tp", to_optional_string_value(orbit.Tp)},
            {"om", std::to_string(orbit.om)},
            {"w", std::to_string(orbit.w)},
            {"in", std::to_string(orbit.in)},
            {"ma", to_optional_string_value(orbit.ma)},
            {"a", to_optional_string_value(orbit.a)},
            {"n", to_optional_string_value(orbit.n)},
        };
    }

    void tag_invoke(const value_from_tag &, value &jv, const MovingTarget &target)
    {
        array alternates;
        for (auto const &name : target.alternate_names())
            alternates.emplace_back(name);

        jv = {
            {"designation", target.designation()},
            {"alternate_names", std::move(alternates)},
            {"small_body", target.small_body()},
            {"orbit", value_from(target.orbit())},
            {"moving_target_id", value_from(target.moving_target_id())},
        };
    }

    void tag_invoke(const value_from_tag &, value &jv, const Ephemeris::Datum &datum)
    {
        jv = {
            {"date", datum.mjd},
            {"tmtp", value_from(datum.tmtp)},
            {"ra", datum.ra},
            {"dec", datum.dec},
            {"mu", datum.mu},
            {"mu_theta", datum.mu_theta},
            {"unc_a", value_from(datum.unc_a)},
            {"unc_b", value_from(datum.unc_b)},
            {"unc_theta", value_from(datum.unc_theta)},
            {"rh", datum.rh},
            {"delta", datum.delta},
            {"phase", datum.phase},
            {"selong", value_from(datum.selong)},
            {"true_anomaly", value_from(datum.true_anomaly)},
            {"sangle", value_from(datum.sangle)},
            {"vangle", value_from(datum.vangle)},
            {"vmag", value_from(datum.vmag)},
        };
    }

    void tag_invoke(const value_from_tag &, value &jv, const Ephemeris::Data &data)
    {
        jv.emplace_array();
        for (auto const &datum : data)
            jv.as_array().emplace_back(value_from(datum));
    }

    void tag_invoke(const value_from_tag &, value &jv, const Ephemeris &eph)
    {
        jv.emplace_object();
        jv.as_object().emplace("target", value_from(eph.target()));
        jv.as_object().emplace("data", array{});

        for (auto const &datum : eph.data())
        {
            object obj = value_from(datum).as_object();
            jv.at("data").as_array().emplace_back(object{});
            for (auto const &[key, value] : obj)
            {
                if ((eph.format.date == Date::Format::Calendar) && (key == "date"))
                    jv.at("data").as_array().back().as_object().emplace("date", Date::MJDToCalendar(value.as_double()));
                else
                    jv.at("data").as_array().back().as_object().emplace(key, value);
            }
        }
    };

    void tag_invoke(const value_from_tag &, value &jv, const sbsearch::Found &found)
    {
        jv.emplace_object();

        json::object obj = value_from(found.observation).as_object();
        for (auto const &[key, value] : obj)
            jv.as_object().emplace(key, value);

        obj = value_from(found.ephemeris.target()).as_object();
        for (auto const &[key, value] : obj)
            jv.as_object().emplace(key, value);

        // if found.ephemeris is a segment, interpolate it to observation mid-time.
        Ephemeris eph;
        if (found.ephemeris.num_vertices() > 1)
            eph = ephemeris::interpolate(found.ephemeris, found.observation.mjd_mid());
        else
            eph = found.ephemeris;

        obj = value_from(eph).at("data").at(0).as_object();
        for (auto const &[key, value] : obj)
            jv.as_object().emplace(key, value);
    }

    void tag_invoke(const value_from_tag &, value &jv, const sbsearch::Founds &founds)
    {
        jv.emplace_array();
        for (sbsearch::Found found : founds.data)
        {
            found.observation.format = founds.observation_format;
            found.ephemeris.format = founds.ephemeris_format;
            jv.as_array().emplace_back(value_from(found));
        }
    }

    void tag_invoke(const value_from_tag &, value &jv, const QueryInfo::Polygon &polygon)
    {
        jv.emplace_array();
        for (int i = 0; i < 4; i++)
            jv.as_array().emplace_back(array{polygon[i][0], polygon[i][1]});
    }
}

namespace sbsearch::json
{
    long double get_string_as_long_double(boost::json::object &obj,
                                          std::initializer_list<string_view> keys)
    {
        for (string_view key : keys)
            if (obj.contains(key))
            {
                if (auto p = obj[key].if_string())
                {
                    char *end{};
                    return std::strtold(p->data(), &end);
                }
                throw ValueError(string(key) + " is not a string");
            }

        // no keys were present
        throw KeyError(util::join(vector<string_view>(keys), ", "));
    };

    optional<long double> get_string_as_optional_long_double(boost::json::object &obj,
                                                             std::initializer_list<string_view> keys)
    {
        for (string_view key : keys)
            if (obj.contains(key))
                if (obj[key].is_null())
                    return std::nullopt;
                else
                    return get_string_as_long_double(obj, {key});

        // no more keys to try
        return std::nullopt;
    };
}

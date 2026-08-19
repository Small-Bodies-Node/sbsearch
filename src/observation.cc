#include "config.h"

#include <algorithm>
#include <cinttypes>
#include <optional>
#include <ostream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <s2/s2error.h>
#include <s2/s2polygon.h>
#include <s2/s2latlng.h>

#include "date.h"
#include "observation.h"
#include "exceptions.h"
#include "polygons.h"
#include "table.h"
#include "vertices.h"
#include "util/optional.h"

using namespace sbsearch::table;
using std::string;
using std::vector;

namespace sbsearch
{
    Observation::Observation(string_view source,
                             string_view observatory,
                             string_view product_id,
                             const double mjd_start,
                             const double mjd_stop,
                             string_view fov,
                             const Terms &terms,
                             const optional<int64_t> &observation_id,
                             const optional<string> &center,
                             const optional<string> &meta,
                             const double mjd_added)
        : source_(source),
          observatory_(observatory),
          product_id_(product_id),
          observation_id_(observation_id),
          mjd_start_(mjd_start),
          mjd_stop_(mjd_stop),
          fov_(fov),
          terms_(terms),
          center_(center),
          meta_(meta),
          mjd_added_(mjd_added)
    {
        is_valid();
    }

    bool Observation::is_valid() const
    {
        // checks `fov`: must be parsable into at least 3 vertices
        // ensures that stop >= start
        vector<S2Point> vertices = make_vertices(fov_);
        if (vertices.size() < 3)
            throw std::runtime_error("FOV must be parsable into at least three vertices.");

        if (mjd_stop_ < mjd_start_)
            throw std::runtime_error("Observation stops before it starts.");

        return true;
    }

    std::ostream &operator<<(std::ostream &os, const Observation &observation)
    {
        auto format_date = [&observation](double mjd)
        { return observation.format.date == Date::Format::MJD
                     ? std::to_string(mjd)
                     : Date(mjd).iso(); };

        os
            << observation.observation_id().value_or(-1) << "  "
            << '"' << observation.source() << '"'
            << "  "
            << '"' << observation.observatory() << '"'
            << "  "
            << '"' << observation.product_id() << '"'
            << "  "
            << format_date(observation.mjd_start()) << "  "
            << format_date(observation.mjd_stop()) << "  "
            << (observation.mjd_stop() - observation.mjd_start()) * 86400;

        if (observation.format.show_fov)
            os << "  "
               << '"' << observation.fov() << '"';

        if (observation.format.show_meta)
            os << "  "
               << '"' << observation.meta().value_or("") << '"';

        return os;
    }

    bool Observation::is_same_fov(const Observation &other) const
    {
        // if the strings are the same, we're done
        if (fov_ == other.fov_)
            return true;

        S2Polygon this_polygon, other_polygon;
        make_polygon(*this, this_polygon);
        make_polygon(other, other_polygon);
        return this_polygon.BoundaryEquals(other_polygon);
    }

    bool Observation::operator==(const Observation &other) const
    {
        return (
            (source_ == other.source_) &
            (observatory_ == other.observatory_) &
            (product_id_ == other.product_id_) &
            (observation_id_ == other.observation_id_) &
            (mjd_start_ == other.mjd_start_) &
            (mjd_stop_ == other.mjd_stop_) &
            is_same_fov(other) &
            (meta_ == other.meta_));
    }

    vector<string> Observations::source() const
    {
        return get_data_vector<vector<string>>(&Observation::source);
    }

    vector<string> Observations::observatory() const
    {
        return get_data_vector<vector<string>>(&Observation::observatory);
    }

    vector<string> Observations::product_id() const
    {
        return get_data_vector<vector<string>>(&Observation::product_id);
    }

    vector<optional<int64_t>> Observations::observation_id() const
    {
        return get_data_vector<vector<optional<int64_t>>>(&Observation::observation_id);
    }

    vector<double> Observations::mjd_start() const
    {
        return get_data_vector<vector<double>>(&Observation::mjd_start);
    }

    vector<double> Observations::mjd_stop() const
    {
        return get_data_vector<vector<double>>(&Observation::mjd_stop);
    }

    vector<string> Observations::fov() const
    {
        return get_data_vector<vector<string>>(&Observation::fov);
    }

    vector<optional<string>> Observations::center() const
    {
        return get_data_vector<vector<optional<string>>>(&Observation::center);
    }

    vector<Observation::Terms> Observations::terms() const
    {
        return get_data_vector<vector<Observation::Terms>>(&Observation::terms);
    }

    vector<optional<string>> Observations::meta() const
    {
        return get_data_vector<vector<optional<string>>>(&Observation::meta);
    }

    vector<double> Observations::mjd_added() const
    {
        return get_data_vector<vector<double>>(&Observation::mjd_added);
    }

    vector<double> Observations::mjd_mid() const
    {
        return get_data_vector<vector<double>>(&Observation::mjd_mid);
    }

    vector<double> Observations::exposure() const
    {
        return get_data_vector<vector<double>>(&Observation::exposure);
    }

    std::ostream &operator<<(std::ostream &os, const Observations &observations)
    {
        auto format_dates = [](vector<double> mjds)
        {
            vector<string> dates(mjds.size());
            std::transform(mjds.begin(), mjds.end(), dates.begin(),
                           [](const double mjd)
                           { return Date(mjd).iso(); });
            return dates;
        };

        int n = observations.size();

        Table table;
        table.add(Column("observation_id", "%" PRId64,
                         util::optionals_to_values(observations.observation_id(),
                                                   (int64_t)-1)));
        table.add(Column("source", "%s", observations.source()));
        table.add(Column("product_id", "%s", observations.product_id()));
        table.add(Column("observatory", "%s", observations.observatory()));

        if (observations.format.date == Date::Format::MJD)
        {
            table.add(Column("mjd_start", "%.6lf", observations.mjd_start()));
            table.add(Column("mjd_stop", "%.6lf", observations.mjd_stop()));
        }
        else
        {
            table.add(Column("mjd_start", "%s", format_dates(observations.mjd_start())));
            table.add(Column("mjd_stop", "%s", format_dates(observations.mjd_stop())));
        }

        table.add(Column("exposure", "%.3lf", observations.exposure()));
        if (observations.format.show_fov)
            table.add(Column("fov", "%s", observations.fov()));
        if (observations.format.show_meta)
            table.add(Column("meta", "%s",
                             util::optionals_to_values(observations.meta(), string(""))));

        os << table;
        return os;
    }

    void Observations::remove_duplicate_observation_ids()
    {
        std::unordered_set<int64_t> unique_ids;
        unique_ids.reserve(data.size());

        auto is_duplicate = [&](auto const &observation)
        { return !unique_ids.insert(observation.observation_id().value_or(-1)).second; };

        data.erase(std::remove_if(data.begin(), data.end(), is_duplicate), data.end());
    }
}

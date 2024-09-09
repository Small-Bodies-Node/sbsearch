
#include "config.h"

#include <iterator>
#include <string>
#include <vector>
#include <gtest/gtest.h>

#include "ephemeris.h"
#include "found.h"
#include "observation.h"

using sbsearch::Ephemeris;
using sbsearch::Observation;
using std::string;
using std::vector;

namespace sbsearch
{
    namespace testing
    {
        TEST(FoundTests, StreamInsertOperator)
        {
            MovingTarget encke{"2P"};

            Observations observations = {
                Observation("test source", "I41", "a", 59252.01, 59252.019, "1:3, 2:3, 2:4, 1:4", "", 1),
                Observation("test source", "I41", "b", 59252.02, 59252.029, "2:3, 3:3, 3:4, 2:4", "", 2)};

            Ephemeris eph(encke, {{59252.01, 10.01, 0, 3.5, 0, 0, 0, 1, 1, 0},
                                  {59252.02, 10.02, 1.5, 3.5, 0, 0, 0, 1, 1, 0},
                                  {59252.03, 10.03, 2.5, 3.5, 0, 0, 0, 1, 1, 0},
                                  {59252.04, 10.04, 3.5, 3.5, 0, 0, 0, 1, 1, 0}});
            Founds founds;

            founds.append(Found(observations[0], eph.segment(0)));
            founds.append(Found(observations[1], eph.interpolate(observations[1].mjd_mid())));

            // Should be two found observations
            EXPECT_EQ(founds.size(), 2);

            std::stringstream stream;
            stream.str("");
            stream << founds;
            EXPECT_EQ(
                stream.str(),
                "observation_id       source  observatory     mjd_start      mjd_stop  exposure  moving_target_id  designation  small_body           mjd       tmtp        ra       dec      rh   delta  phase    selong  true_anomaly    sangle    vangle  unc_a  unc_b  unc_th    vmag\n"
                "--------------  -----------  -----------  ------------  ------------  --------  ----------------  -----------  ----------  ------------  ---------  --------  --------  ------  ------  -----  --------  ------------  --------  --------  -----  -----  ------  ------\n"
                "             1  test source          I41  59252.010000  59252.019000   777.600                -1           2P        true  59252.014500  10.014500  0.675000  3.500296  1.0000  1.0000  0.000  -999.000      -999.000  -999.000  -999.000  0.000  0.000   0.000  99.000\n"
                "             2  test source          I41  59252.020000  59252.029000   777.600                -1           2P        true  59252.024500  10.024500  1.950000  3.500132  1.0000  1.0000  0.000  -999.000      -999.000  -999.000  -999.000  0.000  0.000   0.000  99.000\n");

            stream.str("");
            founds.data[0].observation.format.show_fov = true;
            stream << founds;
            EXPECT_EQ(
                stream.str(),
                "observation_id       source  observatory     mjd_start      mjd_stop  exposure                 fov  moving_target_id  designation  small_body           mjd       tmtp        ra       dec      rh   delta  phase    selong  true_anomaly    sangle    vangle  unc_a  unc_b  unc_th    vmag\n"
                "--------------  -----------  -----------  ------------  ------------  --------  ------------------  ----------------  -----------  ----------  ------------  ---------  --------  --------  ------  ------  -----  --------  ------------  --------  --------  -----  -----  ------  ------\n"
                "             1  test source          I41  59252.010000  59252.019000   777.600  1:3, 2:3, 2:4, 1:4                -1           2P        true  59252.014500  10.014500  0.675000  3.500296  1.0000  1.0000  0.000  -999.000      -999.000  -999.000  -999.000  0.000  0.000   0.000  99.000\n"
                "             2  test source          I41  59252.020000  59252.029000   777.600  2:3, 3:3, 3:4, 2:4                -1           2P        true  59252.024500  10.024500  1.950000  3.500132  1.0000  1.0000  0.000  -999.000      -999.000  -999.000  -999.000  0.000  0.000   0.000  99.000\n");
        }

    }
}
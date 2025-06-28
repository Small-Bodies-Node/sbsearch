#include <gtest/gtest.h>

#include "./sbsdb_test.h"
#include "exceptions.h"
#include "observatory.h"
#include "sbsdb.h"

using namespace sbsearch;
using namespace sbsearch::sbsdb;
using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::vector;

namespace sbsearch::sbsdb::testing
{
    // SBSearchDatabaseTest fixture defined in sbsdb_test.cc

    TEST_F(SBSearchDatabaseTest, ObservatoryIO)
    {
        const Observatory ztf{243.14022, 0.836325, +0.546877, "I41"};
        const Observatory ldt{248.57749, 0.822887, 0.566916, "G37"};
        const Observatory maunakea{204.5278, 0.94171, +0.33725, "568"};
        const Observatory paranal{289.59569, 0.909943, -0.414336, "309"};

        add::observatory(&db, ztf);
        add::observatory(&db, ldt);
        add::observatory(&db, maunakea);
        add::observatory(&db, paranal);

        Observatory obs = get::observatory(&db, "I41");
        EXPECT_EQ(obs, ztf);

        obs = get::observatory(&db, "G37");
        EXPECT_EQ(obs, ldt);

        obs = get::observatory(&db, "568");
        EXPECT_EQ(obs, maunakea);

        obs = get::observatory(&db, "309");
        EXPECT_EQ(obs, paranal);

        Observatories observatories = get::all_observatories(&db);
        EXPECT_EQ(observatories["I41"], ztf);
        EXPECT_EQ(observatories["G37"], ldt);
        EXPECT_EQ(observatories["568"], maunakea);
        EXPECT_EQ(observatories["309"], paranal);

        remove::observatory(&db, "G37");
        EXPECT_THROW(get::observatory(&db, "G37"), ObservatoryError);

        EXPECT_THROW(add::observatory(&db, ztf), ObservatoryError);
    }
}
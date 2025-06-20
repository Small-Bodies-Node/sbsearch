#include "config.h"

#include <sstream>
#include <string>
#include <gtest/gtest.h>

#include "date.h"

using std::string;

namespace sbsearch
{
    namespace testing
    {
        TEST(DateTests, TestInitialization)
        {
            Date date("2024-03-21");
            EXPECT_EQ(date.mjd(), 60390.0);

            date = Date("2024-03-22 04");
            EXPECT_EQ(date.iso(), "2024-03-22 04:00:00");
            EXPECT_EQ(date.mjd(), 60391.0 + 4.0 / 24);

            date = Date("2024-03-22 ");
            EXPECT_EQ(date.iso(), "2024-03-22 00:00:00");
            EXPECT_EQ(date.mjd(), 60391.0);

            date = Date("2024-03-22 04:48");
            EXPECT_EQ(date.iso(), "2024-03-22 04:48:00");
            EXPECT_EQ(date.mjd(), 60391.2);

            date = Date("2024-03-22 17:16:48.0");
            EXPECT_EQ(date.iso(), "2024-03-22 17:16:48");
            EXPECT_EQ(date.mjd(), 60391.72);

            date = Date(60391.0);
            EXPECT_EQ(date.iso(), "2024-03-22 00:00:00");
            EXPECT_EQ(date.mjd(), 60391.0);

            date = Date("60391.2");
            EXPECT_EQ(date.iso(), "2024-03-22 04:48:00");
            EXPECT_EQ(date.mjd(), 60391.2);

            date = Date(60391.72);
            EXPECT_EQ(date.iso(), "2024-03-22 17:16:48");
            EXPECT_EQ(date.mjd(), 60391.72);

            date = Date(60391.999999);
            EXPECT_EQ(date.iso(), "2024-03-23 00:00:00");
            EXPECT_EQ(date.mjd(), 60391.999999);

            EXPECT_THROW(Date(2e9), std::range_error);
        }

        TEST(DateTests, TestStreamIO)
        {
            std::stringstream stream;
            Date date("2024-03-21");
            stream << date;
            EXPECT_EQ(stream.str(), "2024-03-21 00:00:00");

            stream.str("1999-12-31");
            stream >> date;
            EXPECT_EQ(date.iso(), "1999-12-31 00:00:00");
        }
    }
}
#include <regex>
#include <sstream>
#include <gtest/gtest.h>

#include "config.h"
#include "progress_widgets.h"

namespace sbsearch::testing
{
    TEST(ProgressWidgetsTest, ProgressPercent)
    {
        std::stringstream stream;
        ProgressPercent progress(5, stream);
        for (int i = 0; i < 5; i++)
        {
            progress.update();
            progress.status();
        }
        EXPECT_EQ(progress.count(), 5);

        progress.reset();
        EXPECT_EQ(progress.count(), 0);
        progress.update(3);
        EXPECT_EQ(progress.count(), 3);
        progress.status();
        EXPECT_EQ(stream.str(), "\r     20%\n\r     40%\n\r     60%\n\r     80%\n\r    100%\n\r     60%\n");
        stream.str("");

        progress.done();
        std::regex re("^[0-9]+e-[0-9]+ seconds elapsed.\n$");
        EXPECT_TRUE(std::regex_match(stream.str(), re));
    }

    TEST(ProgressWidgetsTest, ProgressTriangle)
    {
        std::stringstream stream;
        ProgressTriangle progress(stream);

        for (int i = 0; i < 15; i++)
            progress.update();
        progress.status();

        EXPECT_EQ(progress.count(), 15);
        EXPECT_EQ(stream.str(), ".\n..\n...\n15\n");

        stream.str("");
        progress.done();
        std::regex re("^[0-9]+e-[0-9]+ seconds elapsed.\n$");
        EXPECT_TRUE(std::regex_match(stream.str(), re));
    }
}
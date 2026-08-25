#include <regex>
#include <sstream>
#include <gtest/gtest.h>

#include "config.h"
#include "progress_widgets.h"

namespace sbsearch::testing
{
    TEST(ProgressWidgetsTest, ProgressFraction)
    {
        std::stringstream stream;
        ProgressFraction progress(5, stream);
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
        progress.total_count(50);
        progress.status();
        progress.prefix("asdf ");
        progress.suffix(" jkl;");
        progress.status();
        EXPECT_EQ(stream.str(), "1/5\n2/5\n3/5\n4/5\n5/5\n3/5\n3/50\nasdf 3/50 jkl;\n");
        stream.str("");

        progress.done();
        std::regex re("^[0-9].[0-9]+e-[0-9]+ seconds elapsed.\n$");
        EXPECT_TRUE(std::regex_match(stream.str(), re));
    }

    TEST(ProgressWidgetsTest, ProgressPercent)
    {
        std::stringstream stream;
        ProgressPercent progress(5, stream);
        progress.prefix(" ");
        progress.suffix("...");

        for (int i = 0; i < 5; i++)
        {
            progress.update();
            progress.status(false);
        }
        EXPECT_EQ(progress.count(), 5);

        progress.reset();
        EXPECT_EQ(progress.count(), 0);
        progress.update(3);
        EXPECT_EQ(progress.count(), 3);
        progress.status();
        EXPECT_EQ(stream.str(), " 20%... 40%... 60%... 80%... 100%... 60%...\n");
        stream.str("");
    }

    TEST(ProgressWidgetsTest, ProgressTriangle)
    {
        std::stringstream stream;
        ProgressTriangle progress(stream);

        for (int i = 0; i < 15; i++)
            progress.update();
        progress.status();

        EXPECT_EQ(progress.count(), 15);
        EXPECT_EQ(stream.str(), "    0 .\n    0 ..\n    0 ...\n15\n");
    }
}
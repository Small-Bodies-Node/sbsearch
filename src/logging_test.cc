#include "config.h"

#include <iostream>
#include <regex>
#include <sstream>
#include <unistd.h>
#include <gtest/gtest.h>

#include "logging.h"

#define DATE_PATTERN "20[0-9][0-9]-[01][0-9]-[0-3][0-9] [01][0-9]:[0-5][0-9]:[0-5][0-9]"

namespace sbsearch
{
    namespace testing
    {
        TEST(LoggingTest, LoggerBaseInit)
        {
            std::string messages;
            std::stringstream stream;

            char *filename = strdup("/tmp/tmpfileXXXXXX");
            mkstemp(filename);
            std::fstream file(filename, std::ios_base::in | std::ios_base::out);

            LoggerBase logger;
            logger.attach(&stream);
            logger.attach(&file);

            logger.info() << 1 << std::endl;
            std::regex re("^" DATE_PATTERN "::INFO::1\n$");
            messages = stream.str();
            EXPECT_TRUE(std::regex_match(messages, re));

            file.seekp(0);
            std::getline(file, messages);
            re = "^" DATE_PATTERN "::INFO::1$"; // no trailing new line expected
            EXPECT_TRUE(std::regex_match(messages, re));

            file.close();
        }

        TEST(LoggingTest, LoggerBaseLevels)
        {
            std::stringstream stream;
            LoggerBase logger;
            logger.attach(&stream);

            logger.log_level(DEBUG);
            logger.debug() << 1 << std::endl;
            logger.info() << 2 << std::endl;
            logger.warning() << 3 << std::endl;
            logger.error() << 4 << std::endl;
            std::string messages = stream.str();

            // Log format is
            // YYYY-MM-DD HH:mm:ss::LOGLEVEL::message
            std::regex re("^" DATE_PATTERN "::DEBUG::1\n" DATE_PATTERN "::INFO::2\n" DATE_PATTERN "::WARNING::3\n" DATE_PATTERN "::ERROR::4\n$");
            EXPECT_TRUE(std::regex_match(messages, re));

            stream.str("");
            logger.log_level(INFO);
            logger.debug() << 1 << std::endl;
            logger.info() << 2 << std::endl;
            logger.warning() << 3 << std::endl;
            logger.error() << 4 << std::endl;
            messages = stream.str();
            re = "^" DATE_PATTERN "::INFO::2\n" DATE_PATTERN "::WARNING::3\n" DATE_PATTERN "::ERROR::4\n$";
            EXPECT_TRUE(std::regex_match(messages, re));

            stream.str("");
            logger.log_level(WARNING);
            logger.debug() << 1 << std::endl;
            logger.info() << 2 << std::endl;
            logger.warning() << 3 << std::endl;
            logger.error() << 4 << std::endl;
            messages = stream.str();
            re = "^" DATE_PATTERN "::WARNING::3\n" DATE_PATTERN "::ERROR::4\n$";
            EXPECT_TRUE(std::regex_match(messages, re));

            stream.str("");
            logger.log_level(ERROR);
            logger.debug() << 1 << std::endl;
            logger.info() << 2 << std::endl;
            logger.warning() << 3 << std::endl;
            logger.error() << 4 << std::endl;
            messages = stream.str();
            re = "^" DATE_PATTERN "::ERROR::4\n$";
            EXPECT_TRUE(std::regex_match(messages, re));

            stream.str("");
            logger.log_level(ERROR + 1);
            logger.debug() << 1 << std::endl;
            logger.info() << 2 << std::endl;
            logger.warning() << 3 << std::endl;
            logger.error() << 4 << std::endl;
            messages = stream.str();
            re = "^$";
            EXPECT_TRUE(std::regex_match(messages, re));
        }

        TEST(LoggingTest, LoggerInit)
        {
            // The logger was likely initialized already, but just in case:
            Logger &logger = Logger::get_logger("/tmp/sbsearch_test_loggingtest.log");

            // best I think we can do is just exercise the code
            int level = logger.log_level();
            logger.log_level(INFO);
            Logger::debug() << 1 << std::endl;
            Logger::info() << 2 << std::endl;
            Logger::warning() << 3 << std::endl;
            Logger::error() << 4 << std::endl;
            logger.log_level(level);
        }

        TEST(LoggingTest, ProgressPercent)
        {
            std::stringstream stream;
            ProgressPercent progress(5, stream);
            for (int i = 0; i < 5; i++)
            {
                progress.update();
                progress.status();
            }

            progress.reset();
            progress.update(3);
            progress.status();
            EXPECT_EQ(stream.str(), "     20%\n     40%\n     60%\n     80%\n    100%\n     60%\n");
            stream.str("");

            progress.done();
            std::regex re("^[0-9]+e-[0-9]+ seconds elapsed.\n$");
            EXPECT_TRUE(std::regex_match(stream.str(), re));
        }
    }
}

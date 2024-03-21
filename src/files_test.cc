#include "config.h"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <regex>
#include <stdlib.h>
#include <string>
#include <boost/filesystem.hpp>
#include <gtest/gtest.h>

#include "files.h"
#include "logging.h"

#define DATE_PATTERN "20[0-9][0-9]-[01][0-9]-[0-3][0-9] [012][0-9]:[0-5][0-9]:[0-5][0-9]"

using std::string;
namespace fs = boost::filesystem;

namespace sbsearch
{
    namespace testing
    {
        TEST(FilesTests, TestReadFile)
        {
            std::ofstream out("sbsearch-test-file");
            out << "asdf";
            out.close();
            const string contents = read_file("sbsearch-test-file");
            EXPECT_EQ(contents, "asdf");
            fs::remove("sbsearch-test-file");
        }

        TEST(FilesTests, TestReadFileError)
        {
            EXPECT_THROW(read_file("sbsearch-test-file"), std::runtime_error);
        }

        TEST(FilesTests, TestCachedFileName)
        {
            const string home = std::getenv("HOME");
            setenv("HOME", "./", 1);
            EXPECT_EQ(generate_cache_file_name("test string").string(),
                      "./.cache/sbsearch/6f8db599de986fab7a21625b7916589c");

            unsetenv("HOME");
            EXPECT_EQ(generate_cache_file_name("test string").string(),
                      "/tmp/sbsearch/6f8db599de986fab7a21625b7916589c");

            setenv("HOME", home.c_str(), 1);
        }

        TEST(FilesTests, TestWriteToCache)
        {
            if (fs::exists("/tmp/sbsearch-testing"))
                fs::remove_all("/tmp/sbsearch-testing");

            struct stat sb;
            fs::create_directory("/tmp/sbsearch-testing");
            int error = stat("/tmp/sbsearch-testing", &sb);
            EXPECT_EQ(error, 0);
            chmod("/tmp/sbsearch-testing", S_IRUSR | S_IWUSR | S_IXUSR);

            write_to_cache(fs::path("/tmp/sbsearch-testing/cache-test"), "asdf");
            EXPECT_EQ(read_file("/tmp/sbsearch-testing/cache-test"), "asdf");
            fs::remove("/tmp/sbsearch-testing/cache-test");

            // make the directory unwritable and test for errors
            std::string messages;
            std::stringstream stream;

            Logger &logger = Logger::get_logger();
            logger.attach(&stream);

            // cannot create directory
            chmod("/tmp/sbsearch-testing", S_IRUSR | S_IXUSR);
            write_to_cache(fs::path("/tmp/sbsearch-testing/unwritable/cache-test"), "asdf");
            messages = stream.str();
            EXPECT_TRUE(
                std::regex_search(
                    messages,
                    std::regex(DATE_PATTERN "::ERROR::Could not write cache file /tmp/sbsearch-testing/unwritable/cache-test: boost::filesystem::create_directories: Permission denied.*")));
            stream.str("");

            // cannot create file
            chmod("/tmp/sbsearch-testing", S_IRUSR | S_IWUSR | S_IXUSR);
            fs::create_directory("/tmp/sbsearch-testing/unwritable");
            chmod("/tmp/sbsearch-testing/unwritable", S_IRUSR | S_IXUSR);
            write_to_cache(fs::path("/tmp/sbsearch-testing/unwritable/cache-test"), "asdf");
            messages = stream.str();
            EXPECT_TRUE(
                std::regex_search(
                    messages,
                    std::regex(DATE_PATTERN "::ERROR::Could not write cache file /tmp/sbsearch-testing/unwritable/cache-test: Unable to open file for writing.*")));
            stream.str("");

            fs::remove_all("/tmp/sbsearch-testing");

            logger.reset_buffer(); // i.e., detach string stream
        }
    }
}
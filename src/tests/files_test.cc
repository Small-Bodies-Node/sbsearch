#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <regex>
#include <string>
#include <gtest/gtest.h>

#include "config.h"
#include "files.h"
#include "logging.h"

#define DATE_PATTERN "20[0-9][0-9]-[01][0-9]-[0-3][0-9] [012][0-9]:[0-5][0-9]:[0-5][0-9]"

using std::string;
using namespace sbsearch;
namespace fs = std::filesystem;

namespace sbsearch::testing
{
    class FilesTests : public ::testing::Test
    {
    protected:
        static void SetUpTestSuite()
        {
            Logger::get_logger().attach(&capture);
        }

        static void TearDownTestSuite()
        {
            Logger::get_logger().reset_buffer();
        }

        void SetUp() override
        {
            capture.str("");
        }

        static std::stringstream capture;
    };

    std::stringstream FilesTests::capture;

    TEST_F(FilesTests, TestReadFile)
    {
        std::ofstream out("sbsearch-test-file");
        out << "asdf";
        out.close();
        const string contents = read_file("sbsearch-test-file");
        EXPECT_EQ(contents, "asdf");
        fs::remove("sbsearch-test-file");
    }

    TEST_F(FilesTests, TestReadFileError)
    {
        EXPECT_THROW(read_file("sbsearch-test-file"), std::runtime_error);
    }

    TEST_F(FilesTests, TestCachedFileName)
    {
        const string home = std::getenv("HOME");
        setenv("HOME", "./", 1);
        EXPECT_EQ(generate_cache_file_name("test string").string(),
                  "./.cache/sbsearch/6f8db599de986fab7a21625b7916589c");

        unsetenv("HOME");
        EXPECT_EQ(generate_cache_file_name("test string").string(),
                  "/tmp/sbsearch/6f8db599de986fab7a21625b7916589c");

        // make the directory unwritable and test for errors
        std::string messages;

        setenv("HOME", "/tmp/sbsearch-testing", 1);
        fs::create_directory("/tmp/sbsearch-testing");
        chmod("/tmp/sbsearch-testing", S_IRUSR | S_IXUSR);

        fs::path fn = generate_cache_file_name("test string");

        messages = capture.str();
        EXPECT_TRUE(
            std::regex_search(
                messages,
                std::regex(DATE_PATTERN
                           "::ERROR::"
                           "Could not create cache directory: "
                           "filesystem error: "
                           "cannot create directories: "
                           "Permission denied \\[/tmp/sbsearch-testing/.cache/sbsearch\\]")));

        chmod("/tmp/sbsearch-testing", S_IRUSR | S_IWUSR | S_IXUSR);
        fs::remove_all("/tmp/sbsearch-testing");
        setenv("HOME", home.c_str(), 1);
    }

    TEST_F(FilesTests, TestWriteToCache)
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

        // cannot create file
        chmod("/tmp/sbsearch-testing", S_IRUSR | S_IWUSR | S_IXUSR);
        fs::create_directory("/tmp/sbsearch-testing/unwritable");
        chmod("/tmp/sbsearch-testing/unwritable", S_IRUSR | S_IXUSR);
        write_to_cache(fs::path("/tmp/sbsearch-testing/unwritable/cache-test"), "asdf");

        string messages = capture.str();
        std::cerr << messages << "\n";
        EXPECT_TRUE(
            std::regex_search(
                messages,
                std::regex(DATE_PATTERN
                           "::ERROR::"
                           "Could not write cache file /tmp/sbsearch-testing/unwritable/cache-test: "
                           "Unable to open file for writing\\.")));
        capture.str("");

        fs::remove_all("/tmp/sbsearch-testing");
    }
}
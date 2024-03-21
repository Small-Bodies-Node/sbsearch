#include <fstream>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <streambuf>
#include <stdio.h>
#include <string>
#include <boost/filesystem.hpp>
#include <openssl/evp.h>

#include "files.h"
#include "logging.h"

using std::string;
namespace fs = boost::filesystem;

namespace sbsearch
{
    const string read_file(const string &file)
    {
        std::ifstream inf;
        try
        {
            inf.open(file);
        }
        catch (std::exception &e)
        {
            Logger::error() << e.what() << std::endl;
        }

        if (!inf.is_open())
            throw std::runtime_error("failed to open " + file);

        string contents;

        inf.seekg(0, std::ios::end);
        contents.reserve(inf.tellg());
        inf.seekg(0, std::ios::beg);

        contents.assign((std::istreambuf_iterator<char>(inf)),
                        std::istreambuf_iterator<char>());
        return contents;
    }

    size_t write_http_string_data(void *buffer, size_t size, size_t nmemb, void *data)
    {
        string *contents = (string *)data;
        size_t realsize = size * nmemb;

        try
        {
            if (contents->capacity() < (contents->size() + realsize))
                contents->reserve(contents->capacity() + realsize);
        }
        catch (std::exception &e)
        {
            cerr << "Cannot reserve string capacity " << contents->size() + realsize << ", capacity is already " << contents->capacity() << "\n";
            throw;
        }
        contents->append((char *)buffer, realsize);
        return realsize;
    }

    fs::path generate_cache_file_name(const string s)
    {
        fs::path path;
        const char *home = std::getenv("HOME");
        if (home == NULL)
            path = fs::path("/tmp/sbsearch");
        else
        {
            path = fs::path(home) / ".cache" / "sbsearch";
        }
        // hash the string to generate the file name
        EVP_MD_CTX *mdctx;
        const EVP_MD *md;

        md = EVP_get_digestbyname("MD5");
        mdctx = EVP_MD_CTX_new();
        unsigned char md_value[EVP_MAX_MD_SIZE];
        unsigned int md_len, i;
        EVP_DigestInit_ex(mdctx, md, NULL);
        EVP_DigestUpdate(mdctx, s.c_str(), s.size());
        EVP_DigestFinal_ex(mdctx, md_value, &md_len);
        EVP_MD_CTX_free(mdctx);

        char hash[2 * md_len + 1];
        for (int i = 0; i < md_len; i++)
            sprintf(hash + i * 2, "%02x", md_value[i]);
        hash[2 * md_len] = '\0';

        // append it to the path and return
        path /= string(hash);
        return path;
    }

    void write_to_cache(const fs::path filename, const string contents)
    {
        try
        {
            if (!fs::exists(filename.parent_path()))
            {
                Logger::info() << "Creating cache directories: " << filename.parent_path().string() << std::endl;
                fs::create_directories(filename.parent_path());
            }

            std::ofstream out(filename.string());
            if (!out)
                throw std::runtime_error("Unable to open file for writing.");
            out << contents;
            out.close();

            Logger::debug() << "Wrote " << contents.size() << " bytes to cache file: " << filename.string() << std::endl;
        }
        catch (std::exception &e)
        {
            Logger::error() << "Could not write cache file " << filename.string()
                            << ": " << e.what() << std::endl;
            return;
        }
    }
}
#ifndef FILES_H_
#define FILES_H_

#include "config.h"

#include <string>
#include <boost/filesystem.hpp>

using std::string;
namespace fs = boost::filesystem;

namespace sbsearch
{
    // Read a file and return the contents as a string
    const std::string read_file(const std::string &file);

    // Write a string data from CURL to a buffer.
    //
    // For example:
    //     curl_easly_setopt(handle, CURLOPT_WRITEFUNCTION, write_http_string_data)
    size_t write_http_string_data(void *buffer, size_t size, size_t nmemb, void *data);

    // Get a cached file name based on a hash of the given string.
    //
    // If the environment variable HOME exists, then the cache directory is
    // ${HOME}/.cache/sbearch, otherwise it will be /tmp/sbsearch
    fs::path generate_cache_file_name(const string s);

    // Write string data to the cache.
    void write_to_cache(const fs::path filename, const string contents);
}

#endif // FILES_H_
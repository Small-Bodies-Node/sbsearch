#ifndef FILES_H_
#define FILES_H_

#include "config.h"

#include <string>
#include <string_view>
#include <boost/filesystem.hpp>

using std::string;
using std::string_view;

namespace fs = boost::filesystem;

namespace sbsearch
{
    // Get the sbsearch ephemeris cache location.  The directories will be
    // created if they do not exist.
    fs::path get_cache_location();

    // Read a file and return the contents as a string
    string read_file(string_view file);

    // Write a string data from CURL to a buffer.
    //
    // For example:
    //     curl_easly_setopt(handle, CURLOPT_WRITEFUNCTION, write_http_string_data)
    size_t write_http_string_data(void *buffer, size_t size, size_t nmemb, void *data);

    // Get a cached file name based on a hash of the given string.
    //
    // If the environment variable HOME exists, then the cache directory is
    // ${HOME}/.cache/sbearch, otherwise it will be /tmp/sbsearch
    fs::path generate_cache_file_name(string_view s);

    // Write string data to the cache.
    void write_to_cache(const fs::path filename, string_view contents);
}

#endif // FILES_H_
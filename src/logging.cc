#include <chrono>
#include <iostream>
#include <string>

#include "config.h"
#include "logging.h"

using std::cout;

namespace sbsearch
{
    int LoggingBuffer::sync()
    {
        std::string s = str(); // get string in buffer
        if (s.length() == 0)
            return 0;

        str(""); // clear buffer

        std::time_t now = std::time(nullptr);
        for (std::ostream *os : streams)
        {
            (*os) << std::put_time(std::localtime(&now), "%F %T")
                  << "::" << s;
            os->flush(); // keep log files up to date
        }
        return 0;
    }

    std::ostream &operator<<(std::ostream &os, const LogLevel &level)
    {
        std::string s;
        switch (level)
        {
        case DEBUG:
            s = "DEBUG";
            break;
        case INFO:
            s = "INFO";
            break;
        case WARNING:
            s = "WARNING";
            break;
        case ERROR:
            s = "ERROR";
            break;
        }
        return os << s;
    }

    std::ostream &LoggerBase::log(LogLevel level, std::string label)
    {
        if (log_level() <= level)
            return (*this) << label << "::";
        else
            return NullStream::get();
    }

    Logger &Logger::get_logger(const std::string &filename)
    {
        static Logger logger(filename);
        return logger;
    }

    std::ostream &Logger::log(LogLevel level)
    {
        Logger &logger = Logger::get_logger();
        if (logger.log_level() <= level)
        {
            logger << level << "::";
            return logger;
        }
        else
            return NullStream::get();
    }
}

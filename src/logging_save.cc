#include "config.h"
#include "logging.h"

#include <chrono>
#include <iomanip>
#include <iostream>
#include <s2/base/integral_types.h>

using std::cout;

namespace sbsearch
{
    std::ostream &LoggerBase::log(LogLevel level, std::string label)
    {
        if (log_level() <= level)
        {
            (*this) << label << "::";
            return (*this);
        }
        else
            return NullStream::get();
    }

    Logger &Logger::get_logger()
    {
        static Logger logger;
        return logger;
    }

    /*
    int LoggerStringBuf::sync()
    {
        std::string s = str();
        if (s.length() == 0)
            return 0;

        str("");

        std::time_t now = std::time(nullptr);
        os << std::put_time(std::localtime(&now), "%F %T")
           << "::" << s;
        return 0;
    }

    std::ostream &LoggerBase::log(LogLevel level, std::string label)
    {
        if (log_level_ <= level)
        {
            os << label << "::";
            return os;
        }
        else
            return NullStream::get();
    }

    Logger &Logger::get_logger()
    {
        static Logger logger("sbsearch.log");
        return logger;
    }

    std::ostream &Logger::log(LogLevel level, std::string label)
    {
        Logger &logger = Logger::get_logger();
        if (logger.log_level() <= level)
        {
            logger << label << "::";
            return logger;
        }
        else
            return NullStream::get();
    }
  */
    void ProgressWidget::reset()
    {
        count = 0;
        t0 = std::chrono::steady_clock::now();
    }

    void ProgressWidget::elapsed()
    {
        std::chrono::duration<double> diff = std::chrono::steady_clock::now() - t0;
        log << std::setprecision(0) << " seconds elapsed.\n";
    }

    void ProgressWidget::done()
    {
        elapsed();
    }

    void ProgressPercent::status()
    {
        log << "\r" << std::setprecision(3) << std::setw(7) << float(count) / total_count * 100 << "%" << std::flush;
    }

    void ProgressPercent::update(int64 increment)
    {
        count += increment;
        status();
    }
}
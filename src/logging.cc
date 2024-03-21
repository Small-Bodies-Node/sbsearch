#include "config.h"
#include "logging.h"

#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>

using std::cout;

namespace sbsearch
{
    // A string buffer that prepends log info with text and writes to multiple streams
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

    int64 ProgressWidget::count()
    {
        return count_;
    }

    void ProgressWidget::reset()
    {
        count_ = 0;
        t0 = std::chrono::steady_clock::now();
    }

    ProgressWidget &ProgressWidget::operator+=(const int64 increment)
    {
        update(increment);
        return *this;
    }

    ProgressWidget &ProgressWidget::operator++()
    {
        update(1);
        return *this;
    }

    double ProgressWidget::elapsed()
    {
        std::chrono::duration<double> diff = std::chrono::steady_clock::now() - t0;
        return diff.count();
    }

    void ProgressWidget::done()
    {
        log << std::setprecision(1) << elapsed() << " seconds elapsed." << std::endl;
    }

    void ProgressPercent::status()
    {
        log << std::setprecision(3) << std::setw(7) << float(count_) / total_count * 100 << "%" << std::endl;
    }

    void ProgressPercent::update(int64 increment)
    {
        count_ += increment;
    }

    void ProgressTriangle::status()
    {
        log << count_ << std::endl;
    }

    void ProgressTriangle::update(int64 increment)
    {
        count_ += increment;
        float x = std::log10(count_);
        while (x >= next_update)
        {
            log << std::string(next_update, '.') << " " << std::pow(10, next_update) << std::endl;
            next_update += 1;
        }
    }
}
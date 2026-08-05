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

    void ProgressPercent::status(const bool end_line)
    {
        log << "\r" << std::setprecision(3) << std::setw(7)
            << float(count_) / total_count * 100 << "%"
            << std::flush;
        if (end_line)
            log << std::endl;
    }

    void ProgressPercent::update(int64 increment)
    {
        count_ += increment;
    }

    void ProgressTriangle::status(const bool end_line)
    {
        log << count_;
        if (end_line)
            log << std::endl;
    }

    void ProgressTriangle::update(int64 increment)
    {
        count_ += increment;
        float x = (base_ == 2) ? std::log2(count_) : std::log10(count_);
        while (x >= next_update)
        {
            log << std::string(next_update, '.') << std::endl;
            next_update += 1;
        }
    }

    void ProgressTriangle::base(int b)
    {
        switch (b)
        {
        case 2:
        case 10:
            base_ = b;
	    break;
	default:
            throw std::invalid_argument("Base must be 2 or 10.");
	}
    }
}

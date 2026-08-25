
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <string_view>

#include "config.h"
#include "progress_widgets.h"

namespace sbsearch
{
    int64 ProgressWidget::total_count()
    {
        return total_count_;
    };

    void ProgressWidget::total_count(int64 n)
    {
        total_count_ = n;
    };

    int64 ProgressWidget::count()
    {
        return count_;
    }

    void ProgressWidget::reset()
    {
        count_ = 0;
        t0 = std::chrono::steady_clock::now();
    }

    void ProgressWidget::prefix(std::string_view s)
    {
        prefix_ = s;
    }

    void ProgressWidget::suffix(std::string_view s)
    {
        suffix_ = s;
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

    double ProgressWidget::rate()
    {
        return count_ / elapsed();
    }

    void ProgressWidget::done()
    {
        log << std::setprecision(4) << elapsed() << " seconds elapsed." << std::endl;
    }

    void ProgressFraction::status(const bool end_line)
    {
        log << prefix_
            << count_ << "/" << total_count_
            << suffix_;
        if (end_line)
            log << "\n";
        log << std::flush;
    }

    void ProgressFraction::update(int64 increment)
    {
        count_ += increment;
    }

    void ProgressPercent::status(const bool end_line)
    {
        log << prefix_
            << std::setprecision(4)
            << float(count_) / total_count_ * 100 << "%"
            << suffix_;
        if (end_line)
            log << "\n";
        log << std::flush;
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
            log << std::setw(5) << static_cast<int>(elapsed()) << " "
                << std::string(next_update, '.')
                << std::endl;
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
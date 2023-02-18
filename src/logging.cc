#include "config.h"
#include "logging.h"

#include <chrono>
#include <iomanip>
#include <iostream>
#include <s2/base/integral_types.h>

using std::cout;

namespace sbsearch
{
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
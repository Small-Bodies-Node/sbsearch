#ifndef SBS_PROGRESS_WIDGETS_H_
#define SBS_PROGRESS_WIDGETS_H_

#include <chrono>
#include <iostream>
#include <string>
#include <string_view>
#include <vector>
#include <s2/base/integral_types.h>

using std::cerr;

namespace sbsearch
{
    // Abstract base class to visualize task progress
    class ProgressWidget
    {
    public:
        ProgressWidget(int64 n, std::ostream &stream = std::cout)
            : total_count_(n), t0(std::chrono::steady_clock::now()), log(stream) {};

        // total count
        int64 total_count();

        // set total_count
        void total_count(int64 n);

        // counter count
        int64 count();

        // reset the counter
        void reset();

        // Set a prefix for status().
        void prefix(std::string_view s);

        // Set a suffix for status().
        void suffix(std::string_view s);

        // log current status
        virtual void status(const bool end_line = true) = 0;

        // update counter
        virtual void update(const int64 increment) = 0;

        ProgressWidget &operator+=(const int64 increment);
        ProgressWidget &operator++();

        // elapsed time, seconds
        double elapsed();

        // number of updates per second
        double rate();

        // log elapsed time
        void done();

    protected:
        int64 total_count_;
        int64 count_ = 0;
        std::chrono::time_point<std::chrono::steady_clock> t0;
        std::ostream &log;
        std::string prefix_, suffix_;
    };

    // Show progress as a fraction.
    class ProgressFraction : public ProgressWidget
    {
    public:
        ProgressFraction(int64 n, std::ostream &stream = std::cerr) : ProgressWidget(n, stream) {};
        void status(const bool end_line = true) override;
        void update(int64 increment = 1) override;
    };

    // Show progress with a percentage value.
    class ProgressPercent : public ProgressWidget
    {
    public:
        ProgressPercent(int64 n, std::ostream &stream = std::cerr) : ProgressWidget(n, stream) {};
        void status(const bool end_line = true) override;
        void update(int64 increment = 1) override;
    };

    // Show progress with a growing set of dots.
    class ProgressTriangle : public ProgressWidget
    {
    public:
        ProgressTriangle(std::ostream &stream = std::cerr) : ProgressWidget(0, stream) {};
        void status(const bool end_line = true) override;
        void update(int64 increment = 1) override;

        // Set the logarithm base for output steps: 2 (default) or 10
        void base(int b);

    private:
        int next_update = 1;
        int base_ = 2;
    };
}

#endif
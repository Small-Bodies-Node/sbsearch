#ifndef SBS_LOGGING_H_
#define SBS_LOGGING_H_

#include <chrono>
#include <iostream>
#include <s2/base/integral_types.h>

namespace sbsearch
{
    class ProgressWidget
    {
    public:
        ProgressWidget(int64 n, std::ostream &stream = std::cout)
            : total_count(n), t0(std::chrono::steady_clock::now()), log(stream){};

        // reset the counter
        void reset();

        // log current status
        virtual void status() = 0;

        // update counter
        virtual void update(int64 increment) = 0;

        // log elapsed time
        void elapsed();

        // log elapsed time
        void done();

    protected:
        int64 total_count;
        int64 count = 0;
        std::chrono::time_point<std::chrono::steady_clock> t0;
        std::ostream &log;
    };

    class ProgressPercent : public ProgressWidget
    {
    public:
        ProgressPercent(int64 n, std::ostream &stream = std::cout) : ProgressWidget(n, stream){};
        void status() override;
        void update(int64 increment = 1) override;
    };
}

#endif // SBS_LOGGING_H_

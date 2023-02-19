#ifndef SBS_LOGGING_H_
#define SBS_LOGGING_H_

#include <chrono>
#include <iostream>
#include <sstream>
#include <s2/base/integral_types.h>

using std::cerr;

// following singleton design pattern example:
// https://stackoverflow.com/questions/1008019/c-singleton-design-pattern

// and the prefix example at:
// https://stackoverflow.com/questions/37490881/overloading-operator-in-c-with-a-prefix

namespace sbsearch
{
    struct LoggerStringBuf : public std::stringbuf
    {
    public:
        LoggerStringBuf(std::ostream &stream) : std::stringbuf(std::ios_base::out), os(stream) {}
        ~LoggerStringBuf() { sync(); }

    protected:
        virtual int sync()
        {
            std::string s = str();
            str("");
            os << "[timestamp] " << s << std::endl;
            return 0;
        };

    private:
        std::ostream &os;
    };

    class Logger : public std::ostream
    {
    public:
        Logger(Logger const &) = delete;
        void operator=(Logger const &) = delete;

        static Logger &get_logger()
        {
            static Logger logger;
            return logger;
        }

        enum LogLevel
        {
            debug = 10,
            info = 20,
            error = 30,
            warning = 40
        };

        // get/set log level
        int log_level() { return log_level_; };
        void log_level(int level) { log_level_ = level; };

    private:
        // Constructor
        Logger() : std::ostream(0), sbuf(cerr)
        {
            init(&sbuf);
        }

        int log_level_ = 20;
        LoggerStringBuf sbuf;
    };

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

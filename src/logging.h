#ifndef SBS_LOGGING_H_
#define SBS_LOGGING_H_

#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <streambuf>
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
        virtual int sync();

    private:
        std::ostream &os;
    };

    // a no-op ostream singleton
    // get a copy with NullStream::get()
    class NullStream : public std::ostream
    {
    public:
        NullStream(NullStream const &) = delete;
        void operator=(NullStream const &) = delete;

        static NullStream &get()
        {
            static NullStream null;
            return null;
        }

    private:
        NullStream() : std::ostream(&m_sb) {}

        class NullBuffer : public std::streambuf
        {
        public:
            int overflow(int c) { return c; }
        } m_sb;
    };

    // An application-wide logger.
    //
    // Always logs to "sbsearch.log".
    //
    // Use the log-level interfaces:
    //   Logger::debug() << message << std::endl;
    //   Logger::info() << message << std::endl;
    //   Logger::warning() << message << std::endl;
    //   Logger::error() << message << std::endl;
    class Logger : public std::ostream
    {
    public:
        Logger(Logger const &) = delete;
        void operator=(Logger const &) = delete;

        static Logger &get_logger();

        template <typename T>
        std::ostream &operator<<(const T &data)
        {
            return static_cast<std::ostream &>(*this) << data;
        }

        enum LogLevel
        {
            DEBUG = 10,
            INFO = 20,
            WARNING = 30,
            ERROR = 40
        };

        // get/set log level
        int log_level() { return log_level_; };
        void log_level(int level) { log_level_ = level; };

        // log to specific log levels
        static std::ostream &log(LogLevel level, std::string label);
        static std::ostream &debug() { return Logger::log(DEBUG, "DEBUG"); }
        static std::ostream &info() { return Logger::log(INFO, "INFO"); }
        static std::ostream &warning() { return Logger::log(WARNING, "WARNING"); }
        static std::ostream &error() { return Logger::log(ERROR, "ERROR"); }

    private:
        // Constructor
        Logger(const std::string &filename) : std::ostream(0), fstream(filename, std::ios_base::app), sbuf(fstream)
        {
            init(&sbuf);
        }

        int log_level_ = DEBUG;
        std::ofstream fstream;
        LoggerStringBuf sbuf;
    };

    // Abstract base class to visualize task progress
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

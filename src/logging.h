#ifndef SBS_LOGGING_H_
#define SBS_LOGGING_H_

#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <streambuf>
#include <string>
#include <vector>
#include <s2/base/integral_types.h>

using std::cerr;

// following singleton design pattern example:
// https://stackoverflow.com/questions/1008019/c-singleton-design-pattern

// and the prefix example at:
// https://stackoverflow.com/questions/37490881/overloading-operator-in-c-with-a-prefix

namespace sbsearch
{
    struct LoggingBuffer : public std::stringbuf
    {
    public:
        LoggingBuffer() : std::stringbuf(std::ios_base::out) {}
        ~LoggingBuffer() { sync(); }

        void attach(std::ostream *stream) { streams.push_back(stream); }

    protected:
        virtual int sync();

    private:
        std::vector<std::ostream *> streams;
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

    enum LogLevel
    {
        DEBUG = 10,
        INFO = 20,
        WARNING = 30,
        ERROR = 40
    };

    // A logging class that writes to multiple streams.
    // LoggerBase is separated out from Logger to facilitate testing
    class LoggerBase : public std::ostream
    {
    public:
        // Constructor
        LoggerBase() : std::ostream(0)
        {
            // buffer.attach(&fstream);
            // buffer.attach(&std::clog);
            init(&buffer);
        }

        LoggerBase(LoggerBase const &) = delete;
        void operator=(LoggerBase const &) = delete;

        void attach(std::ostream *stream) { buffer.attach(stream); }

        template <typename T>
        std::ostream &operator<<(const T &data)
        {
            return static_cast<std::ostream &>(*this) << data;
        }

        // get/set log level
        int log_level() { return log_level_; };
        void log_level(int level) { log_level_ = level; };

        std::ostream &log(LogLevel level, std::string label);
        std::ostream &debug() { return log(DEBUG, "DEBUG"); }
        std::ostream &info() { return log(INFO, "INFO"); }
        std::ostream &warning() { return log(WARNING, "WARNING"); }
        std::ostream &error() { return log(ERROR, "ERROR"); }

    private:
        int log_level_ = DEBUG;
        std::ofstream fstream;
        LoggingBuffer buffer;
    };

    // An application-wide logger.
    //
    // Use the log-level interfaces:
    //   Logger::debug() << message << std::endl;
    //   Logger::info() << message << std::endl;
    //   Logger::warning() << message << std::endl;
    //   Logger::error() << message << std::endl;
    //
    // Use std::endl to indicate the end of the message, and syncs the buffer with the device (console or file).
    class Logger : public std::ostream
    {
    public:
        Logger(Logger const &) = delete;
        void operator=(Logger const &) = delete;

        static Logger &get_logger(const std::string &filename = "/dev/null");

        template <typename T>
        std::ostream &operator<<(const T &data)
        {
            return static_cast<std::ostream &>(*this) << data;
        }

        // get/set log level
        int log_level() { return log_level_; };
        void log_level(int level) { log_level_ = level; };

        static std::ostream &log(LogLevel level, std::string label);
        static std::ostream &debug() { return Logger::log(DEBUG, "DEBUG"); }
        static std::ostream &info() { return Logger::log(INFO, "INFO"); }
        static std::ostream &warning() { return Logger::log(WARNING, "WARNING"); }
        static std::ostream &error() { return Logger::log(ERROR, "ERROR"); }

    private:
        // Constructor
        Logger(const std::string &filename) : std::ostream(0), fstream(filename, std::ios_base::app)
        {
            buffer.attach(&fstream);
            buffer.attach(&std::cerr);
            init(&buffer);
        }

        int log_level_ = DEBUG;
        std::ofstream fstream;
        LoggingBuffer buffer;
    };

    // Abstract base class to visualize task progress
    class ProgressWidget
    {
    public:
        ProgressWidget(int64 n, std::ostream &stream = std::cout)
            : total_count(n), t0(std::chrono::steady_clock::now()), log(stream){};

        // counter count
        int64 count();

        // reset the counter
        void reset();

        // log current status
        virtual void status() = 0;

        // update counter
        virtual void update(const int64 increment) = 0;

        ProgressWidget &operator+=(const int64 increment);
        ProgressWidget &operator++();

        // elapsed time, seconds
        double elapsed();

        // log elapsed time
        void done();

    protected:
        int64 total_count;
        int64 count_ = 0;
        std::chrono::time_point<std::chrono::steady_clock> t0;
        std::ostream &log;
    };

    class ProgressPercent : public ProgressWidget
    {
    public:
        ProgressPercent(int64 n, std::ostream &stream = std::cerr) : ProgressWidget(n, stream){};
        void status() override;
        void update(int64 increment = 1) override;
    };

    class ProgressTriangle : public ProgressWidget
    {
    public:
        ProgressTriangle(std::ostream &stream = std::cerr) : ProgressWidget(0, stream){};
        void status() override;
        void update(int64 increment = 1) override;

    private:
        int next_update = 1;
    };
}

#endif // SBS_LOGGING_H_

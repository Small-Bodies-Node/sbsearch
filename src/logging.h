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
        virtual int sync()
        {
            std::string s = str();
            if (s.length() == 0)
                return 0;

            str("");

            std::time_t now = std::time(nullptr);
            os << std::put_time(std::localtime(&now), "%F %T")
               << "::" << s;
            return 0;
        };

    private:
        std::ostream &os;
    };

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

    class Logger : public std::ostream
    {
    public:
        Logger(Logger const &) = delete;
        void operator=(Logger const &) = delete;

        static Logger &get_logger()
        {
            static Logger logger("sbsearch.log");
            return logger;
        }

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
        static std::ostream &debug()
        {
            Logger &logger = Logger::get_logger();
            if (logger.log_level() <= DEBUG)
            {
                logger << "DEBUG::";
                return logger;
            }
            else
                return NullStream::get();
        }

        static std::ostream &info()
        {
            Logger &logger = Logger::get_logger();
            if (logger.log_level() <= INFO)
            {
                logger << "INFO::";
                return logger;
            }
            else
                return NullStream::get();
        }

        static std::ostream &warning()
        {
            Logger &logger = Logger::get_logger();
            if (logger.log_level() <= WARNING)
            {
                logger << "WARNING::";
                return logger;
            }
            else
                return NullStream::get();
        }

        static std::ostream &error()
        {
            Logger &logger = Logger::get_logger();
            if (logger.log_level() <= ERROR)
            {
                logger << "ERROR::";
                return logger;
            }
            else
                return NullStream::get();
        }

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

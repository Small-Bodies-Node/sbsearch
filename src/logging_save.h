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

// and the multiple streams example at:
// https://stackoverflow.com/questions/16135370/class-to-output-to-several-ostreams-file-and-console

namespace sbsearch
{
    // A stringbuf that will enable us to prepend log info with text, writes to multiple streams
    struct LoggingBuffer : public std::stringbuf
    {
    public:
        ~LoggingBuffer() { sync(); }

        void attach(std::ostream *stream)
        {
            ostreams.push_back(stream);
        }

    protected:
        // When syncing, prepend lines with time of message
        virtual int sync()
        {
            std::string s = str(); // get string in buffer
            if (s.length() == 0)
                return 0;
            str(""); // clear buffer

            std::time_t now = std::time(nullptr);
            for (std::ostream *os : ostreams)
                (*os) << std::put_time(std::localtime(&now), "%F %T")
                      << "::" << s;
            return 0;
        }

    private:
        std::vector<std::ostream *> ostreams;
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

    class LoggerBase : public std::ostream
    {
    public:
        // LoggerBase(const std::string &filename) : file(filename, std::ios_base::app), console(std::clog.rdbuf())
        // {
        //     buf.attach(&console);
        //     buf.attach(&file);
        //     init(&buf);
        // }

        void attach(std::ostream *stream) { buf.attach(stream); };

        // get/set log level
        int log_level() { return log_level_; };
        void log_level(int level) { log_level_ = level; };

        // log to specific log levels
        std::ostream &log(LogLevel level, std::string label);
        std::ostream &debug() { return log(DEBUG, "DEBUG"); }
        std::ostream &info() { return log(INFO, "INFO"); }
        std::ostream &warning() { return log(WARNING, "WARNING"); }
        std::ostream &error() { return log(ERROR, "ERROR"); }

    private:
        int log_level_ = DEBUG;
        std::ofstream file;
        std::ostream console;
        LoggingBuffer buf;
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
    class Logger : public LoggerBase
    {
    public:
        Logger(Logger const &) = delete;
        void operator=(Logger const &) = delete;

        static Logger &get_logger();
        static void attach(std::ostream *stream)
        {
            Logger &logger = Logger::get_logger();
            logger.attach(stream);
        };

        // log to specific log levels
        static std::ostream &log(LogLevel level, std::string label);
        static std::ostream &debug() { return log(DEBUG, "DEBUG"); }
        static std::ostream &info() { return log(INFO, "INFO"); }
        static std::ostream &warning() { return log(WARNING, "WARNING"); }
        static std::ostream &error() { return log(ERROR, "ERROR"); }

    private:
        // Constructor
        Logger();
        LoggingBuffer buf;
    };

    /*
    class LoggerBase : public std::ostream
    {
    public:
        LoggerBase(std::ostream &stream, const std::string &filename) : os(stream.rdbuf()), fstream(filename, std::ios_base::app),
        {

            streams.push_back(&stream);
            if (filename != "")
            {
                fstream = std::ofstream(filename, std::ios_base::app);
            }
            streams.push_back(&fstream);
        }

    private:
        int log_level_ = DEBUG;
        std::vector<std::ostream *> streams;
        std::ofstream fstream;
        std::ostream os;
        LoggerStringBuf fbuf, cbuf; // file and console buffers
    }

        class LoggerBase : public std::ostream
        {
        public:
            LoggerBase(std::ostream &stream, const std::string &filename = "")
            {
                std::vector<std::ostream *> streams;
                streams.push_back(&stream);
                if (filename != "")
                {
                    fstream = std::ofstream(filename, std::ios_base::app);
                }
                streams.push_back(&fstream);
            }

            template <typename T>
            std::ostream &operator<<(const T &data)
            {
                return (std::ostream &)(*this) << data;
            }

            // get/set log level
            int log_level() { return log_level_; };
            void log_level(int level) { log_level_ = level; };

            // log to specific log levels
            std::ostream &log(LogLevel level, std::string label);
            std::ostream &debug() { return log(DEBUG, "DEBUG"); }
            std::ostream &info() { return log(INFO, "INFO"); }
            std::ostream &warning() { return log(WARNING, "WARNING"); }
            std::ostream &error() { return log(ERROR, "ERROR"); }

        private:
            int log_level_ = INFO;
            LoggerStringBuf sbuf;
            std::ofstream fstream;

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
    */
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

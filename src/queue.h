#ifndef SBS_QUEUE_H_
#define SBS_QUEUE_H_

#include <condition_variable>
#include <mutex>
#include <optional>
#include <queue>

#include "exceptions.h"

using std::optional;

namespace sbsearch
{
    // Thread-safe queue.
    //
    // Add items with put().  Get the next item (FIFO) with next().  Indicate
    // that no more items will be added with finish().  If the queue is ever
    // empty, std::nullopt will be returned.
    template <class T>
    class Queue
    {
    public:
        // Append a new item to the queue.
        void put(T const &item);

        // Remove and return the next item in the queue.  If there isn't another
        // item, next() will wait for one.  In some circumstances, e.g., the
        // queue is empty and no more will be added (finished), a std::nullopt
        // will be returned.
        optional<T> next();

        // True if the task queue is empty.
        bool empty() { return items.empty(); };

        // Indicate to the queue that no more tasks will be added.
        void finish();

        // Returns true if the queue is empty and no more items will be added.
        bool finished() { return empty() && finish_; };

    private:
        std::mutex access;
        std::queue<T> items;
        std::condition_variable condition;
        bool finish_ = false;
    };

    ////////////////////////////////////////////////////////////////////////////
    // implementation

    template <typename T>
    void Queue<T>::put(T const &item)
    {
        if (finish_)
            throw SBSException("Cannot add items to a finished queue.");

        std::lock_guard lock{access};
        items.push(std::move(item));
        condition.notify_one();
    }

    template <typename T>
    optional<T> Queue<T>::next()
    {
        std::unique_lock lock{access};

        // wait until the queue is not empty or is empty and finished
        condition.wait(lock, [this]()
                       { return !empty() || finished(); });

        if (finished())
            return std::nullopt;

        T item = items.front();
        items.pop();

        condition.notify_one();

        return item;
    }

    template <typename T>
    void Queue<T>::finish()
    {
        finish_ = true;
        condition.notify_all();
    }
}

#endif
#ifndef SBS_QUEUE_H_
#define SBS_QUEUE_H_

#include <condition_variable>
#include <mutex>
#include <optional>
#include <queue>
#include <set>

#include "exceptions.h"

using std::optional;

namespace sbsearch
{
    // Thread-safe queue.
    //
    // Add items with put().  Get the next item (FIFO) with next().  Indicate
    // that no more items will be added with finish().  If the queue is ever
    // empty, std::nullopt will be returned.
    template <typename T>
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

        // Number of items currently in the queue.
        size_t size() { return items.size(); };

        // Number of items placed in the queue in total.
        size_t enqueued() { return enqueued_; };

        // True if the task queue is empty.
        bool empty() { return items.empty(); };

        // Indicate to the queue that no more tasks will be added.
        void finish();

        // Returns true if the queue is empty and no more items will be added.
        bool finished() { return empty() && finish_; };

    protected:
        std::mutex access;
        std::queue<T> items;
        std::condition_variable condition;
        size_t enqueued_ = 0;
        bool finish_ = false;
    };

    // Thread-safe queue that only keeps unique items based on a key.
    template <typename K, typename T>
    class UniqueQueue : public Queue<T>
    {
    public:
        // Append a new item to the queue, but only if `key` wasn't already
        // added before.  Returns true if the item was added.
        bool put(K const &key, T const &item);

        // Number of items that were attempted to be put in the queue.
        size_t total_puts() { return puts_; };

        using Queue<T>::empty;
        using Queue<T>::enqueued;
        using Queue<T>::finish;
        using Queue<T>::finished;
        using Queue<T>::next;
        using Queue<T>::size;

    private:
        std::set<K> keys;
        size_t puts_ = 0;

        using Queue<T>::access;
        using Queue<T>::items;
        using Queue<T>::enqueued_;
        using Queue<T>::condition;
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
        enqueued_++;
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

    template <typename K, typename T>
    bool UniqueQueue<K, T>::put(K const &key, T const &item)
    {
        if (finished())
            throw SBSException("Cannot add items to a finished queue.");

        std::lock_guard lock{access};

        puts_++;
        if (keys.find(key) != keys.end())
            return false; // already added

        keys.insert(key);
        items.push(std::move(item));
        enqueued_++;
        condition.notify_one();

        return true;
    }
}

#endif
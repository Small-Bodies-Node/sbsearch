#include <chrono>
#include <future>
#include <thread>
#include <gtest/gtest.h>
#include "queue.h"

namespace sbsearch::testing
{
    void producer(Queue<int> &queue)
    {
        int i = 5;
        while (i > 0)
        {
            int sleep = i * 10;
            queue.put(i--);
            std::this_thread::sleep_for(std::chrono::milliseconds(sleep));
        }
    }

    int consumer(Queue<int> &queue, int scale)
    {
        int consumed = 0;
        bool running = true;
        while (running)
        {
            if (queue.finished())
            {
                running = false;
                continue;
            }
            auto item = queue.next();
            if (!item)
                continue;

            int sleep = item.value() * scale;
            std::this_thread::sleep_for(std::chrono::milliseconds(sleep));
            consumed += item.value();
        }
        return consumed;
    }

    TEST(QueueTests, FastProduction)
    {
        Queue<int> queue;
        std::thread production_thread(producer, std::ref(queue));

        std::packaged_task consumer_task(consumer);
        std::future<int> consumed = consumer_task.get_future();
        std::thread consumption_thread(std::move(consumer_task), std::ref(queue), 3);

        production_thread.join();
        queue.finish();

        consumption_thread.join();
        EXPECT_EQ(consumed.get(), 15);
    }

    TEST(QueueTests, FastConsumption)
    {
        Queue<int> queue;
        std::thread production_thread(producer, std::ref(queue));

        std::packaged_task consumer_task(consumer);
        std::future<int> consumed = consumer_task.get_future();
        std::thread consumption_thread(std::move(consumer_task), std::ref(queue), 3);

        production_thread.join();
        queue.finish();

        consumption_thread.join();
        EXPECT_EQ(consumed.get(), 15);
    }

    TEST(QueueTests, ConsumeAfterFinish)
    {
        Queue<int> queue;
        std::thread production_thread(producer, std::ref(queue));

        production_thread.join();
        queue.finish();

        std::packaged_task consumer_task(consumer);
        std::future<int> consumed = consumer_task.get_future();
        std::thread consumption_thread(std::move(consumer_task), std::ref(queue), 3);
        consumption_thread.join();
        EXPECT_EQ(consumed.get(), 15);
    }
}

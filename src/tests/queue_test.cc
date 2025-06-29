#include <chrono>
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
            // std::cerr << "producing " << sleep << " ms..." << std::endl;
            queue.put(i--);
            std::this_thread::sleep_for(std::chrono::milliseconds(sleep));
        }
    }

    void consumer(Queue<int> &queue, int scale)
    {
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

            // std::cerr << "consuming " << sleep << " ms..." << std::endl;
            std::this_thread::sleep_for(std::chrono::milliseconds(sleep));
        }
    }

    TEST(QueueTests, FastProduction)
    {
        Queue<int> queue;
        std::thread query_thread(producer, std::ref(queue));
        std::thread find_matches_thread(consumer, std::ref(queue), 30);

        query_thread.join();
        queue.finish();

        EXPECT_NO_THROW(find_matches_thread.join());
    }

    TEST(QueueTests, FastConsumption)
    {
        Queue<int> queue;
        std::thread query_thread(producer, std::ref(queue));
        std::thread find_matches_thread(consumer, std::ref(queue), 3);

        query_thread.join();
        queue.finish();

        EXPECT_NO_THROW(find_matches_thread.join());
    }
}

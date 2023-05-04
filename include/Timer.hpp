#pragma once

#include <chrono>
#include <iostream>
class Timer
{
    public:
    Timer()
    {
        start = std::chrono::high_resolution_clock::now();
    }

    ~Timer()
    {
        end = std::chrono::high_resolution_clock::now();
        duration = end - start;

        double microseconds = duration.count() * 1e6;
        std::cout << "\tElapsed time : " << microseconds << " us.\n";
    }

    private:

    std::chrono::time_point<std::chrono::steady_clock> start, end;
    std::chrono::duration<double> duration;
};
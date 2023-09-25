#include <iostream>
#include <array>
#include <fmt/ranges.h>
#include "dist.h"
int main()
{
    Distribution::TruncatedNormal dist(0.5, 0.1, 0.375, 0.625);
    std::array<float, 5> arr;
    for(auto& elem : arr)
    {
        elem = dist();
    }
    std::cout<<fmt::format("Samples: {}", arr)<<std::endl;
}

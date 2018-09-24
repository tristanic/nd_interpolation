/**
 * @Author: Tristan Croll <tic20>
 * @Date:   18-Apr-2018
 * @Email:  tic20@cam.ac.uk
 * @Last modified by:   tic20
 * @Last modified time: 26-Apr-2018
 * @License: Free for non-commercial use (see license.pdf)
 * @Copyright: 2017-2018 Tristan Croll
 */

#include <vector>
#include <iostream>
#include <random>
#include <chrono>

#include "nd_interp_mgr.h"


#ifdef _WIN32
# define EXPORT __declspec(dllexport)
#else
# define EXPORT __attribute__((__visibility__("default")))
#endif

int main()
{
    //3x3x3 array with value 1 at [1,1,1], 2 at [2,2,2], 0 everywhere else
    std::vector<float> data {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,2.0};

    Interpolation_Mgr<float> mgr;
    size_t dim=3;
    size_t axis_lengths[3] = {3,3,3};
    float minv[3] = {0,0,0};
    float maxv[3] = {2,2,2};
    mgr.add_interpolator("test1", dim, axis_lengths, minv, maxv, data.data());
    mgr.add_interpolator("test2", dim, axis_lengths, minv, maxv, data.data());

    std::random_device rd{};
    std::mt19937 gen{rd()};

    std::normal_distribution<float> d(1, 0.5);
    float single_test_vals[3];
    for (size_t i=0; i<10; ++i)
    {
        for (size_t j=0; j<3; ++j)
            single_test_vals[j] = d(gen);
        try {
            auto result = mgr.interpolate("test1", single_test_vals);
            std::cout << "Interpolated value for {" <<single_test_vals[0]<<","<<single_test_vals[1]<<","<<single_test_vals[2]
                << "} is " << result << std::endl << std::flush;
        } catch(std::range_error) {
            std::cout << "{"<<single_test_vals[0]<<","<<single_test_vals[1]<<","<<single_test_vals[2]<<"} is outside of the interpolation range!" << std::endl << std::flush;
        }
    }

    // Create a large array of random values, all guaranteed to be in the
    // range of the data
    std::vector<size_t> array_sizes {100, 500, 1000, 5000, 10000, 50000, 100000, 500000, 1000000};
    std::vector<float> speedups;
    for (const auto& large_array_size: array_sizes)
    {
        std::uniform_real_distribution<float> d_u(0.1, 1.9);
        std::vector<float> many_test_vals;
        for (size_t i=0; i<large_array_size; ++i)
            for (size_t j=0; j<3; ++j)
                many_test_vals.push_back(d_u(gen));
        std::vector<float> return_vals(large_array_size);

        std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
        mgr.interpolate("test1", many_test_vals.data(), large_array_size, return_vals.data());
        std::chrono::steady_clock::time_point t2 = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();
        std::cout << "Interpolating "<<large_array_size<<" 3D points took " << duration << " microseconds" << std::endl;

        const char* names[2] = {"test1", "test2"};

        auto many_test_vals_2 = many_test_vals;
        auto return_vals_2 = return_vals;
        float* all_vals[2] = {many_test_vals.data(), many_test_vals_2.data()};
        size_t counts[2] = {large_array_size, large_array_size};
        float* returns[2] = {return_vals.data(), return_vals_2.data()};
    
        t1 = std::chrono::steady_clock::now();
        mgr.interpolate(2, names, counts, all_vals, returns, false);
        t2 = std::chrono::steady_clock::now();
        auto duration_single = std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();
        std::cout << "Interpolating 2 sets of "<<large_array_size<<" 3D points in serial took " << duration_single << " microseconds" << std::endl;

        t1 = std::chrono::steady_clock::now();
        mgr.interpolate(2, names, counts, all_vals, returns, true);
        t2 = std::chrono::steady_clock::now();
        auto duration_thread = std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count();
        std::cout << "Interpolating 2 sets of "<<large_array_size<<" 3D points in parallel took " << duration_thread << " microseconds" << std::endl;
        std::cout << "First set: first value: " << returns[0][0] << " last value: " << returns[0][large_array_size-1] << std::endl;
        std::cout << "Second set: first value: " << returns[1][0] << " last value: " << returns[1][large_array_size-1] << std::endl;
        speedups.push_back((float)duration_single/duration_thread);
    }
    for (size_t i=0; i<array_sizes.size(); ++i)
        std::cout << "Array size: " << array_sizes[i] <<"\t Thread speedup: " << speedups[i] << std::endl;
    return 0;
}

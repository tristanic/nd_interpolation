#include <unordered_map>
#include <string>
#include <future>

#include "nd_interp.h"

template<typename T>
class Interpolation_Mgr_Impl
{
public:
    void add_interpolator(const char* name, size_t dim, size_t* n,
        T* minv, T* maxv, T* data);
    T interpolate(const char* name, T* axis_vals) const;
    void interpolate(const char* name, T* axis_vals, size_t n, T* return_values) const;
    void interpolate(size_t n_cases, const char** names, size_t* case_counts,
        T** axis_vals, T** return_values, bool threaded=true) const;
    const char** names();
    size_t num_maps() const;
    size_t dim(const char* name) const;


private:
    std::unordered_map<std::string, RegularGridInterpolator<T>> interpolators_;
    std::vector<const char*> names_;
}; //Interpolation_Mgr_Impl

// IMPLEMENTATIONS

template<typename T> void
Interpolation_Mgr_Impl<T>::add_interpolator(const char* name, size_t dim, size_t* n,
    T* minv, T* maxv, T* data)
{
    interpolators_[std::string(name)] = RegularGridInterpolator<T>(dim, n, minv, maxv, data);
}

template<typename T> T
Interpolation_Mgr_Impl<T>::interpolate(const char* name, T* axis_vals) const
{
    return interpolators_.at(std::string(name)).interpolate(axis_vals);
}

template<typename T> void
Interpolation_Mgr_Impl<T>::interpolate(const char* name, T* axis_vals, size_t n, T* return_values) const
{
    interpolators_.at(std::string(name)).interpolate(axis_vals, n, return_values);
}

template<typename T> void
Interpolation_Mgr_Impl<T>::interpolate(size_t n_cases, const char** names, size_t* case_counts,
    T** axis_vals, T** return_values, bool threaded) const
{
    if (!threaded)
    {
        for (size_t i=0; i<n_cases; ++i)
            interpolate(names[i], axis_vals[i], case_counts[i], return_values[i]);
        return;
    } else {
        std::vector<std::future<void>> threads;
        for (size_t i=0; i<n_cases; ++i)
        {
            threads.push_back(std::async(std::launch::async,
                static_cast<void (Interpolation_Mgr_Impl<T>::*)(const char*, T*, size_t, T*) const>(&Interpolation_Mgr_Impl<T>::interpolate),
                this, names[i], axis_vals[i], case_counts[i], return_values[i]
            ));
        }
        // Make sure all threads have completed before returning
        for (auto& t: threads)
            t.wait();
    }
}

template<typename T> const char**
Interpolation_Mgr_Impl<T>::names()
{
    names_.clear();
    for (auto kv : interpolators_)
        names_.push_back(kv.first.c_str());
    return names_.data();
}

template<typename T> size_t
Interpolation_Mgr_Impl<T>::num_maps() const
{
    return interpolators_.size();
}

template<typename T> size_t
Interpolation_Mgr_Impl<T>::dim(const char* name) const
{
    return interpolators_.at(std::string(name)).dim();
}

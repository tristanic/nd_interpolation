/**
 * @Author: Tristan Croll <tic20>
 * @Date:   18-Apr-2018
 * @Email:  tic20@cam.ac.uk
 * @Last modified by:   tic20
 * @Last modified time: 26-Apr-2018
 * @License: Free for non-commercial use (see license.pdf)
 * @Copyright: 2017-2018 Tristan Croll
 */



//! N-dimensional regular grid interpolator

#pragma once

#include <stdint.h>
#include <vector>
#include <math.h>
#include <iostream>


template <typename T>
class RegularGridInterpolator
{

public:
    RegularGridInterpolator() {} // null constructor
    ~RegularGridInterpolator() {} // destructor
    //! Construct a RegularGridInterpolator object for the given data
    /*!
     * This implementation requires one data value for every grid point
     * dim:  the number of dimensions
     * n:    the number of points in each dimension
     * minv: the minimum axis value for each dimension
     * maxv: the maximum axis value for each dimension
     * data: the actual data to be interpolated, as a single array in C-style
     *       row order (must match the dimensions defined by the previous
     *       arguments)
     */
    RegularGridInterpolator(const size_t &dim, size_t* n, T* minv, T* maxv, T* data);

    //! Interpolate a single point
    T interpolate(T *axis_vals) const;
    //! Interpolate a single point
    T interpolate(std::vector<T> axis_vals) const;
    //! Interpolate n points
    void interpolate(T* axis_vals, const size_t &n, T* values) const;
    const std::vector<T> &min() const {return min_;}
    const std::vector<T> &max() const {return max_;}
    const size_t &dim() const {return dim_;}
    const std::vector<T> &data() const {return data_;}
    const std::vector<size_t> &length() const {return n_;}

private:
    void corner_values(const size_t &lb_indices, std::vector<T> &corners) const;
    void lb_index_and_offsets(T *axis_vals, size_t &lb_index,
        std::vector<std::pair<T, T> > &offsets) const;
    void _interpolate(const size_t &dim, std::vector<T> &corners, size_t size,
    const std::vector<std::pair<T, T> > &offsets, T *value) const;
    void _interpolate1d(const std::pair<T, T> &offset, const T& lower, const T& upper, T *val) const;
    void corner_offsets();

    size_t dim_;
    size_t n_corners_;
    std::vector<size_t> n_;
    std::vector<T> min_;
    std::vector<T> max_;
    std::vector<T> step_;
    std::vector<std::vector<T> > axes_;

    //TODO: Replace data_ with a std::unordered_map sparse array implementation
    //      to minimise memory use for higher dimensions
    std::vector<T> data_;
    std::vector<size_t> corner_offsets_;
    std::vector<size_t> jump_;

}; //RegularGridInterpolator

// IMPLEMENTATIONS

template <typename T>
RegularGridInterpolator<T>::RegularGridInterpolator(const size_t& dim,
        size_t* n, T* minv, T* maxv, T* data)
{
    dim_ = dim;
    size_t this_n, d_count=1;
    T this_min, this_max;
    T step, dval;
    for (size_t i=0; i<dim; ++i) {
        this_n = n[i];
        this_min = minv[i];
        this_max = maxv[i];
        n_.push_back(this_n);
        min_.push_back(this_min);
        max_.push_back(this_max);
        step = (this_max-this_min)/(T)(this_n-1);
        step_.push_back(step);
        jump_.push_back((size_t)pow(2,i));
        dval = this_min;
        std::vector<T> axis;

        for (;dval<this_max+0.5*step;) {
            axis.push_back(dval);
            dval+=step;
        }
        axes_.push_back(axis);
        d_count *= this_n;
    }

    for (size_t i=0; i<d_count; ++i) {
        data_.push_back(data[i]);
    }
    n_corners_ = (size_t)pow(2.0, (T)dim);
    corner_offsets();
} //RegularGridInterpolator

template<typename T>
void
RegularGridInterpolator<T>::lb_index_and_offsets(T *axis_vals, size_t &lb_index,
    std::vector<std::pair<T, T> > &offsets) const
{
    size_t axis_prod = 1;
    for (int axis=dim_-1; axis>=0; --axis) {
        const T &max = max_[axis];
        const T &min = min_[axis];
        const T &value = axis_vals[axis];
        if (value <= min || value >= max) {
            std::cerr << "Value " << value << " is outside of the range " << min << ".." << max << std::endl;
            throw std::range_error("Value outside of interpolation range!");
        }
        size_t li = (size_t)floor((value-min)/step_[axis]);
        lb_index +=axis_prod*li;
        axis_prod*=n_[axis];
        auto &this_axis = axes_[axis];
        const T &low = this_axis[li++];
        const T &high = this_axis[li];
        T offset = (value-low)/(high-low);
        offsets[axis]=(std::pair<T, T> (offset, 1-offset));
    }
}

/*
 * This comes out looking a little like black magic, so requires a bit
 * of explanation. We want to get the values at all the corners
 * in a well-defined order. Using the 3D case as an example, if our
 * lower bound is (0,0,0), we want the corners in the order:
 * ((0,0,0),(0,0,1),(0,1,0),(0,1,1),(1,0,0),(1,0,1),(1,1,0),(1,1,1))
 * ... which is 0 to 7 in binary. The logic below simply extends this
 * to n dimensions.
 */
template<typename T>
void
RegularGridInterpolator<T>::corner_offsets()
{
    for (size_t i=0; i < n_corners_; ++i) {
        size_t corner = 0;
        size_t dim_prod = 1;
        for (size_t j=0; j<dim_; ++j) {
            corner += dim_prod * ((i & (1<<j))>>j);
            dim_prod *= n_[dim_-j-1];
        }
        corner_offsets_.push_back(corner);
    }
}

template<typename T>
void
RegularGridInterpolator<T>::corner_values(const size_t &lb_index, std::vector<T> &corners) const
{
    for (size_t i=0; i<corner_offsets_.size(); i++) {
        corners[i]=(data_[lb_index + corner_offsets_[i]]);
    }
}

// Reduces the vector of corners in-place for efficiency
template<typename T>
void
RegularGridInterpolator<T>::_interpolate(const size_t &dim, std::vector<T> &corners, size_t size,
    const std::vector<std::pair<T, T> > &offsets, T* value) const
{
    for (size_t i=0; i<dim; ++i) {
        const std::pair<T, T> &this_offset=offsets[dim-i-1];
        for (size_t ind=0, j=0; j<size; ind++, j+=2) {
            _interpolate1d(this_offset, corners[j], corners[j+1], &corners[ind]);
        }
        size/=2;
    }
    *value=corners[0];
}

template<typename T>
void
RegularGridInterpolator<T>::_interpolate1d(const std::pair<T, T> &offset, const T &lower, const T &upper, T *val) const
{
    *val= offset.first*upper + offset.second*lower;
}

template<typename T>
T
RegularGridInterpolator<T>::interpolate(T *axis_vals) const
{
    T value[1];
    interpolate(axis_vals, 1, value);
    return value[0];
}

template<typename T>
T
RegularGridInterpolator<T>::interpolate(std::vector<T> axis_vals) const
{
    return interpolate(axis_vals.data());
}


template<typename T>
void
RegularGridInterpolator<T>::interpolate (T* axis_vals, const size_t &n, T* values) const
{
    std::vector<std::pair<T, T>> offsets(dim_);
    std::vector<T> corners(n_corners_);
    for (size_t i=0; i<n; ++i) {
        size_t lb_index = 0;
        // find the minimum corner of the hypercube, and the offsets
        // along each axis
        lb_index_and_offsets(axis_vals+i*dim_, lb_index, offsets);
        // ... and get values at all the corners surrounding the target
        // position.
        corner_values(lb_index, corners);
        _interpolate(dim_, corners, corners.size(), offsets, values++);
    }
}

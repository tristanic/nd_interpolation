#pragma once

#include <memory>
#include "nd_interp_mgr_impl.h"

//! Central manager for interpolation of values on N-dimensional grids
/*! Maintains a <name: interpolator> mapping (e.g. for residue to Ramachandran
 *  contour data) and methods to obtain interpolated values for either a single
 *  lookup or an array of values.
 */
template <typename T>
class Interpolation_Mgr
{
public:
    Interpolation_Mgr()
    {
        impl_ = std::unique_ptr<Interpolation_Mgr_Impl<T>> (new Interpolation_Mgr_Impl<T>());
    }
    //! Create a new interpolator for a given N-dimensional dataset
    /*! Points are expected to be evenly spaced along each axis, but each axis
     *  is allowed to have its own grid spacing. Values must be provided for all
     *  grid points (including zeros).
     *  Arguments:
     *    name: a unique name for this interpolator. Any existing interpolator
     *          with the same name will be replaced.
     *    dim:  the number of dimensions in the grid
     *    n:    the number of points in each dimension.
     *    minv:  the minimum axis value for each dimension
     *    maxv:  the maximum axis value for each dimension
     *    data: the gridded data, as a single array in C-style row order (must
     *          match the dimensions defined by the previous arguments)
    */
    void add_interpolator(const char* name, size_t dim, size_t* n,
        T* minv, T* maxv, T* data)
    {
        impl_->add_interpolator(name, dim, n, minv, maxv, data);
    }
    //! Return the interpolated value for a single instance.
    T interpolate(const char* name, T* axis_vals) const
    {
        return impl_->interpolate(name, axis_vals);
    }
    //! Interpolate an arbitrary number of instances of the same type
    /*! Arguments:
            name:       The name of the interpolator to use. Will throw an error
                        if no such interpolator exists.
            axis_vals:  The values for each look-up, as a one-dimensional
                        c-style array of length (n * {interpolator dimension}).
            n:          The number of interpolations to perform
            return_values: An array of length n, to be filled with the
                           interpolated values.

    */
    void interpolate(const char* name, T* axis_vals, size_t n, T* return_values) const
    {
        impl_->interpolate(name, axis_vals, n, return_values);
    }
    //! Interpolate an arbitrary number of instances for multiple types
    /*! Arguments:
            n_cases:        The number of different cases to be interpolated
            names:          An array of strings of length n_cases, defining the
                            interpolators to be used
            case_counts:    An array of length n_cases, defining the size of
                            the interpolation task for each case
            axis_vals:      The arrays of axis values for each case
            return_values:  The arrays to be filled with interpolated values
            threaded:       If true, each interpolation case will be performed
                            in its own thread.
    */
    void interpolate(size_t n_cases, const char** names, size_t* case_counts,
        T** axis_vals, T** return_values, bool threaded=false) const
    {
        impl_->interpolate(n_cases, names, case_counts, axis_vals, return_values, threaded);
    }

    //! Number of maps managed by this interpolator
    size_t num_maps() const { return impl_->num_maps(); }
    //! Maximum size in bytes for the name of a map.
    size_t max_name_length() const { return impl_->max_name_length(); }
    //! Fills buf with the name of the map at position i.
    void map_name(size_t i, char* buf ) const { impl_->name(i, buf); }
    //! Dimensionality of the interpolator matching the given name
    size_t dim(const char* name) const { return impl_->dim(name); }
private:
    std::unique_ptr<Interpolation_Mgr_Impl<T>> impl_;
}; // Interpolation_Mgr

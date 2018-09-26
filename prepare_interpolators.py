'''
Loads data arrays from standard MolProbity text files, and pre-compiles them
into a C++ object with fast interpolation functions accessible from Python.
'''

import os
from math import degrees, radians, pi


_rama_files = {
    'CISPRO': 'rama8000-cispro.data',
    'TRANSPRO': 'rama8000-transpro.data',
    'GLY':'rama8000-gly-sym.data',
    'ILEVAL': 'rama8000-ileval-nopreP.data',
    'PREPRO': 'rama8000-prepro-noGP.data',
    'GENERAL': 'rama8000-general-noGPIVpreP.data',
}

_rota_files = {
    'ARG': 'rota8000-arg.data',
    'ASN': 'rota8000-asn.data',
    'ASP': 'rota8000-asp.data',
    'CYS': 'rota8000-cys.data',
    'GLN': 'rota8000-gln.data',
    'GLU': 'rota8000-glu.data',
    'HIS': 'rota8000-his.data',
    'ILE': 'rota8000-ile.data',
    'LEU': 'rota8000-leu.data',
    'LYS': 'rota8000-lys.data',
    'MET': 'rota8000-met.data',
    'PHE': 'rota8000-phetyr.data',
    'PRO': 'rota8000-pro.data',
    'SER': 'rota8000-ser.data',
    'THR': 'rota8000-thr.data',
    'TRP': 'rota8000-trp.data',
    'TYR': 'rota8000-phetyr.data',
    'VAL': 'rota8000-val.data',
}

_packages = {
    'ramachandran': ['rama.cpp',],
    'rotamer': ['rota.cpp',],
}

def prepare_constructor(file_dict):
    constructor_text = ''
    for name, filename in file_dict.items():
        ndim, axis_lengths, min_vals, max_vals, grid_data = generate_interpolator_data(
            filename, wrap_axes=True
        )
        constructor_text += '''
    size_t {}_ndim = {};
    size_t {}_axis_lengths[{}] = {{{}}};
    fp_type {}_min_vals[{}] = {{{}}};
    fp_type {}_max_vals[{}] = {{{}}};
    std::vector<fp_type> {}_data {{{}}};
    add_interpolator("{}", {}_ndim, {}_axis_lengths, {}_min_vals, {}_max_vals, {}_data.data());

        '''.format(
            name, ndim,
            name, ndim, ','.join([str(al) for al in axis_lengths]),
            name, ndim, ','.join([str(mv) for mv in min_vals]),
            name, ndim, ','.join([str(mv) for mv in max_vals]),
            name, ','.join((str(d) for d in grid_data.ravel())),
            name, name, name, name, name, name
            )
    return constructor_text


def prepare_ramachandran(rama_files):
    '''
    Prepares a ready-to-compile C++ file containing the data and interpolation
    methods for Ramachandran contours.
    '''

    constructor_text = prepare_constructor(rama_files)

    module_docstring = '''Fast C++ interpolator for Ramachandran validation.'''

    interpolate_single_docstring = ('Interpolate a single residue\'s phi and psi'
        ' angles on the Ramachandran plot.\\\n'
        ' Arguments: \\\n'
        '   class: the class of Ramachandran contour for this residue (e.g. \'GENERAL\')\\\n'
        '   angles: an iterable (list, 1D NumPy array, etc.) of two angles in radians'
        )

    interpolate_multiple_docstring = ('Find the Ramachandran scores for a set of '
        'residues of the same class in a single call. Returns a NumPy array of '
        'scores. \\\n'
        'Arguments: \\\n'
        '   class: the class of Ramachandran contour for these residues(e.g. \'GENERAL\')\\\n'
        '   angles: a (n x 2) NumPy array containing (phi, psi) for each residue in radians.')

    extra_defs = '''
    .def("interpolate_single",
        [](const $C_CLASS_NAME& self, const char* name,
            const fp_type& phi, const fp_type& psi)
        {
            fp_type data[2] {phi, psi};
            return self.interpolate(name, data);
        })
    '''

    outfile = open('rama.cpp', 'wt')
    with open('mgr_base.cpp.in', 'rt') as infile:
        for line in infile:
            line = line.replace('$EXTRA_DEFS', extra_defs)
            line = line.replace('$C_CLASS_NAME', 'Ramachandran_Mgr')
            line = line.replace('$PY_MODULE_NAME', 'ramachandran')
            line = line.replace('$PY_CLASS_NAME', 'Rama_Mgr')
            line = line.replace('$PY_DOCSTRING', module_docstring)
            line = line.replace('$INTERPOLATE_SINGLE_DOCSTRING', interpolate_single_docstring)
            line = line.replace('$INTERPOLATE_MULTIPLE_DOCSTRING', interpolate_multiple_docstring)
            line = line.replace('$CONSTRUCTOR', constructor_text)
            outfile.write(line)

    outfile.close()

def prepare_rotamers(rota_files):
    '''
    Prepares a ready-to-compile C++ file containing the data and interpolation
    methods for amino acid rotamers.
    '''
    constructor_text = prepare_constructor(rota_files)

    module_docstring = '''Fast C++ interpolator for rotamer validation.'''

    interpolate_single_docstring = ('Calculate the rotamer probability for a single residue\\\n'
        ' Arguments: \\\n'
        '   resname: the 3-letter name of this residue (e.g. \'LEU\')\\\n'
        '   chi_angles: an iterable (list, 1D NumPy array, etc.) of the residue\'s chi angles in radians'
        )

    interpolate_multiple_docstring = ('Find the rotamer scores for a set of '
        'residues of the same type in a single call. Returns a NumPy array of '
        'scores. \\\n'
        'Arguments: \\\n'
        '   residue_name: the 3-letter residue code for these residues(e.g. \'LEU\')\\\n'
        '   angles: a (n x num_chi) NumPy array containing chi angles for each residue in radians.')

    extra_defs = ''

    outfile = open('rota.cpp', 'wt')
    with open('mgr_base.cpp.in', 'rt') as infile:
        for line in infile:
            line = line.replace('$EXTRA_DEFS', extra_defs)
            line = line.replace('$C_CLASS_NAME', 'Rotamer_Mgr')
            line = line.replace('$PY_MODULE_NAME', 'rotamer')
            line = line.replace('$PY_CLASS_NAME', 'Rota_Mgr')
            line = line.replace('$PY_DOCSTRING', module_docstring)
            line = line.replace('$INTERPOLATE_SINGLE_DOCSTRING', interpolate_single_docstring)
            line = line.replace('$INTERPOLATE_MULTIPLE_DOCSTRING', interpolate_multiple_docstring)
            line = line.replace('$CONSTRUCTOR', constructor_text)
            outfile.write(line)

    outfile.close()


def generate_interpolator_data(file_name, wrap_axes = True):
    '''
    Load a MolProbity data set and format it into a form ready for generation of
    a RegularGridInterpolator object for later fast interpolation of values.
    '''
    import numpy, pickle
    infile = open(file_name, 'r')
    # Throw away the first line - we don't need it
    infile.readline()
    # Get number of dimensions
    ndim = int(infile.readline().split()[-1])
    # Throw away the next line - it's just headers
    infile.readline()
    lower_bounds = []
    upper_bounds= []
    axis_lengths = []

    step_sizes = []
    min_vals = []
    max_vals = []
    #~ axes = []

    # Read in the header to get the dimensions and step size for each
    # axis, and initialise the axis arrays
    for i in range(ndim):
        line = infile.readline().split()
        lb = float(line[2])
        lower_bounds.append(lb)
        ub = float(line[3])
        upper_bounds.append(ub)
        nb = int(line[4])
        axis_lengths.append(nb)

        ss = (ub - lb)/nb
        step_sizes.append(ss)
        # Values are at the midpoint of each bin
        fs = lb + ss/2
        min_vals.append(fs)
        ls = ub - ss/2
        max_vals.append(ls)
        #~ axis.append(numpy.linspace(fs,ls,nb))

    infile.close()

    grid_data = numpy.zeros(axis_lengths)

    # Slurp in the actual numerical data as a numpy array
    data = numpy.loadtxt(file_name)

    # Convert each coordinate to an integral number of steps along each
    # axis
    axes = []
    for i in range(ndim):
        axes.append([])

    for i in range(ndim):
        ss = step_sizes[i]
        fs = min_vals[i]
        lb = lower_bounds[i]
        axis_vals = data[:,i]
        axes[i]=(((axis_vals - ss/2 - lb) / ss).astype(int))

    grid_data[axes] = data[:,ndim]

    '''
    At this point we should have the full n-dimensional matrix, with
    all values not present in the text file present as zeros.
    Now we have to consider periodicity. Since we're just going to
    be doing linear interpretation, the easiest approach is to simply
    pad the array on all sides with the values from the opposite extreme
    of the relevant matrix. This can be handily done with numpy.pad
    '''
    if wrap_axes:
        grid_data = numpy.pad(grid_data, 1, 'wrap')

        # ... and we need to extend each of the axes by one step to match
        for i, a in enumerate(axes):
            min_v = min_vals[i]
            max_v = max_vals[i]
            ss = step_sizes[i]
            min_vals[i] = min_v - ss
            max_vals[i] = max_v + ss
            axis_lengths[i] += 2
    # Replace all zero or negative values with the minimum positive non-zero
    # value, so that we can use logs
    grid_data[grid_data<=0] = numpy.min(grid_data[grid_data > 0])
    # Finally, convert the axes to radians
    min_vals = numpy.radians(numpy.array(min_vals, numpy.double))
    max_vals = numpy.radians(numpy.array(max_vals, numpy.double))

    # Pickle a tuple containing the axes and grid for fast loading in
    # future runs

    axis_lengths = numpy.array(axis_lengths, numpy.int)
    return (ndim, axis_lengths, min_vals, max_vals, grid_data)

def _extensions():
    import sys
    from setuptools import Extension
    extra_compile_args=['-std=c++11',]
    libraries = []
    inc_dirs=[]
    lib_dirs=[]
    extra_link_args=[]
    if sys.platform in ('darwin', 'linux'):
        extra_compile_args.append('-fvisibility=hidden')
        extra_link_args.append('-s')
    from numpy.distutils.misc_util import get_numpy_include_dirs
    inc_dirs.extend(get_numpy_include_dirs())

    extensions = []
    for package_name, package_files in _packages.items():
        extensions.append(Extension(package_name,
            define_macros=None,
            extra_compile_args=extra_compile_args,
            include_dirs=inc_dirs,
            library_dirs=lib_dirs,
            libraries = libraries,
            extra_link_args=extra_link_args,
            sources=package_files)
            )
    return extensions





if __name__ == '__main__':
    prepare_ramachandran(_rama_files)
    prepare_rotamers(_rota_files)
    extensions = _extensions()
    from distutils.core import setup
    setup(name="MolProbity_Interpolators",
        version="1.0",
        description="Ramachandran, rotamer and CABLAM interpolator classes",
        ext_modules = extensions)

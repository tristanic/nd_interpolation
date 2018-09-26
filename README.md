# nd_interpolation
## Fast C++ N-dimensional interpolation on pre-defined data

Originally designed for fast validation of molecular models against Ramachandran, rotamer and CABLAM (and future other...?) contours from MolProbity. The core classes are defined as pure header-only templates and should be readily applicable to a range of tasks with minor modification. 

### Dependencies

- a Python installation with (as a minimum) setuptools, wheel and numpy modules
- PyBind11
- The MolProbity [Ramachandran](https://github.com/rlabduke/reference_data/tree/master/Top8000/Top8000_ramachandran_pct_contour_grids), [rotamer](https://github.com/rlabduke/reference_data/tree/master/Top8000/Top8000_rotamer_pct_contour_grids) and [CaBLAM](https://github.com/rlabduke/reference_data/tree/master/Top8000/Top8000_cablam_pct_contour_grids) contours from the Richardson lab at Duke University.
  

### Getting started

- `test_nd_interp.cpp` contains a brief example of usage at the C++ API level. Compile it to an executable with:
```
gcc test_nd_interp.cpp -o test_nd_interp -pthread -std=c++11 -lstdc++ -lm -g -O2
```
- `python prepare_interpolators.py bdist_wheel` will parse the rama*.data and rota*.data files into C++ source and compile them into Python-compatible libraries. The resulting library files in build/lib... can be copied to the destination of your choice. Put
them in the same directory as `test_interpolators.py` and run it to test.


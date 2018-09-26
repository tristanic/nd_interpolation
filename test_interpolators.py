import sys, os, numpy
from math import pi
from time import time

import timeit

sys.path.insert(0, os.curdir)

print ('Importing Ramachandran interpolator')
from ramachandran import Rama_Mgr

t = timeit.timeit(
    "rmgr = Rama_Mgr()",
    setup="from ramachandran import Rama_Mgr",
    number=100000
)

print ("Importing Rama_Mgr 100k times took {} seconds".format(t))

rmgr = Rama_Mgr()
num_interpolations = 10000
angles = ((numpy.random.rand(num_interpolations,2)-1)*pi)

print ('Ramachandran case names: {}'.format(rmgr.keys()))

t = timeit.timeit(
    "rmgr.interpolate_single('GENERAL', 1., 1.)",
    setup='''
from ramachandran import Rama_Mgr
rmgr = Rama_Mgr()
    ''',
    number=100000
)
print("Performing 100,000 single Rama evaluations took {} seconds".format(t))

# Check against some known values. Tolerances are relatively loose since
# known values were calculated in double precision but current implementation is
# float
abs_tol=1e-4
from math import isclose
known_rama_scores = (
    [['GENERAL', -1.30604503,  2.68515261], 0.38533466754824536],
    [['TRANSPRO', -1.02151863, -0.69621326], 0.7914812407411016],
    [['CISPRO', -1.20270843,  2.79918665], 0.8549866164172167],
    [['ILEVAL', -2.44296282, -0.35592498], 0.0036246624729488735],
    [['GLY', 1.1895601 , 0.61036634], 0.823064915185082]
)

for r in known_rama_scores:
    try:
        assert isclose(rmgr.interpolate_single(*r[0]), r[1], abs_tol=abs_tol)
    except AssertionError:
        calc_score = rmgr.interpolate_single(*r[0])
        print('WARNING: Calculated score of {} for {} case phi={}, psi={} does not match stored value of {}.'
            .format(calc_score, *r[0], r[1]))



t = timeit.timeit(
    "rmgr.interpolate('GENERAL', angles)",
    setup='''
from ramachandran import Rama_Mgr
rmgr = Rama_Mgr()
from __main__ import angles
    ''',
    number=10
)
print('Performing 10x{} concerted Rama evaluations on random (phi, psi) double-precision pairs took {} seconds'.format(num_interpolations, t))

angles_float = angles.astype(numpy.float32)

t = timeit.timeit(
    "rmgr.interpolate('GENERAL', angles)",
    setup='''
from ramachandran import Rama_Mgr
rmgr = Rama_Mgr()
from __main__ import angles
    ''',
    number=10
)
print('Performing 10x{} concerted Rama evaluations on random (phi, psi) single-precision pairs took {} seconds'.format(num_interpolations, t))

dp_results = rmgr.interpolate('GENERAL', angles)
sp_results = rmgr.interpolate('GENERAL', angles_float)

assert numpy.allclose(dp_results, sp_results)

print('Importing rotamer interpolator')
from rotamer import Rota_Mgr

rota_mgr = Rota_Mgr()

print ('Rotamer names: {}'.format(rota_mgr.keys()))

known_rota_scores = (
    [['THR', [-1.18680692]], 0.27169472688690244],
    [['LYS', [-0.97795798, -3.01192304,  2.85388646,  1.4183923]], 0.2325102334891604],
    [['LEU', [-2.70087757, -2.41489383]], 8.858258998e-06],
    [['ILE', [-0.93842671, -1.05547921]], 0.42725676541200613],
    [['TYR', [2.94456112, 1.56040726]], 0.3246952257152924],
    [['PHE', [2.94456112, 1.56040726]], 0.3246952257152924],
    [['SER', [1.14593017]], 0.9879297437362637],
    [['PRO', [-0.49725932]], 0.9593997327648075],
    [['ASN', [1.32791294, -0.29426418]], 0.23133084762617592],
    [['ASP', [-1.33717312,  2.77402037]], 0.7552585586278451],
    [['VAL', [3.02299915]], 0.8141136782245025],
    [['HIS', [-0.93634878,  1.79943171]], 0.33809833158450403],
    [['ARG', [-2.50029818,  0.86957468, -2.81427486,  2.00635224]], 0.003460325535481471],
    [['GLN', [-1.31204367,  1.53333072, -0.27079344]], 0.2286698477020941],
    [['GLU', [-1.30131147,  1.43252533,  3.11070798]], 0.5837979452557274],
    [['CYS', [3.12039768]], 0.4830497003515302],
    [['MET', [-2.51420927,  2.69937058,  2.71132163]], 0.006540101745710986],
    [['TRP', [-1.14757094,  1.75275024]], 0.9677868018396398],
)

for r in known_rota_scores:
    try:
        assert isclose(rota_mgr.interpolate_single(*r[0]), r[1], abs_tol=abs_tol)
    except AssertionError:
        calc_score = rota_mgr.interpolate_single(*r[0])
        print('WARNING: Calculated rotamer score of {} for {} with chi angles {} does not match stored value of {}.'
            .format(calc_score, r[0][0], ','.join((str(c) for c in r[0][1])), r[1]))



lys_angles = ((numpy.random.rand(num_interpolations,4)-1)*pi)

t = timeit.timeit(
    "rota_mgr.interpolate('LYS', lys_angles)",
    setup='''
from rotamer import Rota_Mgr
rota_mgr = Rota_Mgr()
from __main__ import lys_angles
    ''',
    number=10
)
print('Performing 10x{} concerted lysine rotamer evaluations on random chi angles took {} seconds'.format(num_interpolations, t))

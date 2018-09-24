import sys, os, numpy
from math import pi

sys.path.insert(0, os.curdir)

print ('Importing Ramachandran interpolator')
from ramachandran import Rama_Mgr

rmgr = Rama_Mgr()
num_interpolations = 10000
angles = ((numpy.random.rand(num_interpolations,2)-1)*pi)

print ('Ramachandran case names: {}'.format(rmgr.keys()))


print('Testing Ramachandran with {} double-precision phi/psi pairs: '.format(num_interpolations))
print(rmgr.interpolate('GENERAL', angles))

angles = angles.astype(numpy.float32)
print('Testing Ramachandran with {} single-precision phi/psi pairs: '.format(num_interpolations))
print(rmgr.interpolate('GENERAL', angles))

print('Testing single Ramachandran interpolation')
print(rmgr.interpolate_single('GENERAL', [0.1, 0.1]))


print('Importing rotamer interpolator')
from rotamer import Rota_Mgr

rota_mgr = Rota_Mgr()

print ('Rotamer names: {}'.format(rota_mgr.keys()))

lys_angles = ((numpy.random.rand(num_interpolations,4)-1)*pi)

print ('Testing rotamer validation for {} random double precision LYS rotamers: '.format(num_interpolations))
print(rota_mgr.interpolate('LYS', lys_angles))

lys_angles = lys_angles.astype(numpy.float32)
print ('Testing rotamer validation for {} random single precision LYS rotamers: '.format(num_interpolations))
print(rota_mgr.interpolate('LYS', lys_angles))

print('Testing single lysine rotamer interpolation')
print(rota_mgr.interpolate_single('LYS', [0.1,0.1,0.1,0.1]))

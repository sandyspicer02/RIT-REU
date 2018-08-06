
from __future__ import division, print_function

import illustris_api as ia

import matplotlib.pyplot as pyplot

# Illustris API Configuration
#ia.set_apikey('')
ia.set_apikey('4e2c93455bbc8e991a91a61ef5e11daa')

h = 0.704 # Apparently

z = 0.
lower = 100.#40.
upper = 120.#50.

def test_plot(fname):
    f = h5py.File(fname,'r')
    x = f['PartType0']['Coordinates'][:,0]
    z = f['PartType0']['Coordinates'][:,1]
    dens = numpy.log10(f['PartType0']['Masses'][:])
    import matplotlib.pyplot as pyplot
    pyplot.figure()
    pyplot.hist2d(z,x,weights=dens,bins=[150,100])
    pyplot.show()


files = None
#files = ['cutout_41094.hdf5.orig'] # uncomment and set this to suppress to downloading

if files is None: 

    halo_data_prelim = ia.obtain_halos_by_mass(lower, upper, z)

    c = halo_data_prelim['count']
    if ia.DEBUG: print(c)

    halo_ids, halo_data = ia.obtain_subhalos_data(halo_data_prelim, z)

    halo_ids, halo_data = ia.filter_haloids(halo_ids, halo_data, filters=[ia.filter_non_zero_mass, ia.filter_star_mass, ia.filter_fluid_mass], args=[None, {'lower': 1., 'upper': 100}, {'lower': 1., 'upper': 100.}])

    if ia.DEBUG: print(halo_ids[:])

    files = ia.obtain_cutouts(halo_ids, 0., gascols=['Coordinates', 'Density', 'GFM_Metallicity', 'Masses', 'NeutralHydrogenAbundance','SmoothingLength', 'Velocities'], starcols=['Coordinates', 'Masses', 'Velocities','GFM_Metallicity'], count=3)

    ia.detilt_files(files, halo_data)

ia.process_cutouts_for_splash(files)

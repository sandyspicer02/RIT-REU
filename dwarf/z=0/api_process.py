
from __future__ import division, print_function

import illustris_api as ia

import matplotlib.pyplot as pyplot

# Illustris API Configuration
ia.set_apikey('')
#ia.set_apikey('a013b499584f72d63fe4c1aee4960c57')

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

    ia.store_metadata(halo_data)
    
    # This is an example of how to unpickle, not a real bit of code :-)
    # You probably don't want to uncomment this!S
    #new_halos = ia.load_metadata(['data_9.pickle', 'data_66085.pickle', 'data_183685.pickle'])
    #print(len(new_halos))

    files = ia.obtain_cutouts(halo_data, gascols=['Coordinates', 'Density', 'GFM_Metallicity', 'Masses', 'SmoothingLength', 'Velocities'], starcols=['Coordinates', 'Masses', 'Velocities','GFM_Metallicity'], count=1)

    Ldata = []
    ia.detilt_files(files, halo_data, L_outputs=Ldata)
    if ia.DEBUG: print(Ldata)

ia.process_cutouts_for_splash(files)

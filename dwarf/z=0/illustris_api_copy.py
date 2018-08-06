"""
This is the main library of API and related functions for illustris.

Before using this you will need to obtain an API key, and then after doing
    import illustris_api as ia
set it using
    ia.set_apikey(YOURKEY)
    
Debugging is enabled by calling ia.enable_debugging(). 

The recommended namespace (because Ben is a lazy typist) is ia. *Not* il since that is used for illustris_python which we would like to support too.

Requirements
------------
requests, h5py, and numpy (all obtainable in PIP).
"""
from __future__ import division, print_function

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning) # Silence that annoying warning about h5py being antiquated.

import requests, os, h5py, sys, numpy, math, shutil

__author__ = "Benjamin Lewis, Sukanya Chakrabarti"
__license__ = "GPLv2 (without an 'or later versions' clause)"
__maintainer__ = "Benjamin Lewis" # For his sins.
__email__ = "ben.lewis@benl.co.uk"

apikey = None
baseUrl = 'http://www.illustris-project.org/api/'

DEBUG = False # Set this to True to get lots of terminal spam

def config(debug=None, key=None):
    """
    Little helper function to set globals properly. Use the functions in 'See Also' instead.

    Parameters
    ----------
    debug : bool
    key : str
        Your API key

    Returns
    -------
    None
    
    See Also
    --------
    enable_debug()
    disable_debug()
    set_apikey()
    """
    global DEBUG
    global apikey
    if debug is not None: DEBUG = debug
    if key is not None: apikey = key
    
def enable_debug():
    """
    Enable debugging. Produces copious output.

    Parameters
    ----------
    None

    Returns
    -------
    None
    
    See Also
    --------
    enable_debug()
    disable_debug()
    set_apikey()
    config()
    """
    config(debug=True)

def disable_debug():
    """
    Disable debugging. This is the default state.

    Parameters
    ----------
    None

    Returns
    -------
    None
    
    See Also
    --------
    enable_debug()
    disable_debug()
    set_apikey()
    config()
    """
    config(debug=False)

def set_apikey(key):
    """
    Use this to set the api key.

    Parameters
    ----------
    key : str
        Your API key

    Returns
    -------
    None
    
    See Also
    --------
    enable_debug()
    disable_debug()
    set_apikey()
    config()
    """
    if key is None or key == '': return
    print('Setting API Key to %s' % key)
    config(key=key)

def get(path, params=None):
    """
    Function to get things from the api.

    A deliberatly unhandled exception will be raised if something goes wrong.

    Parameters
    ----------
    path : str
        Url to get
    params : list of str
        Not presently used for anything. Do not use without consulting blewis.

    Returns
    -------
    JSON data
        JSON data returned by the API.
                  OR
    str
        filename of downloaded data.
    """
    if DEBUG: print('apikey:', apikey)
    if apikey is None or apikey == '':
        print(' ERROR: no api key has been set. You need to call set_apikey().')
        sys.exit(666)
    headers = {"api-key":apikey}
    # make HTTP GET request to path
    r = requests.get(path, params=params, headers=headers)
    # raise exception if response code is not HTTP SUCCESS (200)
    r.raise_for_status()
    if r.headers['content-type'] == 'application/json':
        return r.json() # parse json responses automatically


    if 'content-disposition' in r.headers:
        filename = r.headers['content-disposition'].split("filename=")[1]
        with open(filename, 'wb') as f:
            f.write(r.content)
        return filename # return the filename string

    return r

def obtain_halos_by_mass(lower, upper, z, limit=10000000):
    """
    Get a list of halo ids within a certain range of virial masses.

    Parameters
    ----------
    lower : float
        Lower bound on M200 (the virial mass)
    upper : float
        Upper bound on M200 (the virial mass)
    z : float
        The redshift
    limit : int
        (optional, default=1000) Limit API call to this many results. Use for debugging only.

    Returns
    -------
    JSON data
        JSON data returned by the API, formatted as a Python list of lists.
        

    See Also
    --------
    obtain_subhalos_data : Same idea, but for the full meta data.
    """
    lower = lower
    upper = upper
    
    search_query = "?limit=" + str(limit) + "&mass__gt=" + str(lower) + "&mass__lt=" + str(upper)
    url = "http://www.illustris-project.org/api/Illustris-1/snapshots/z=%i/subhalos/%s" % (z, search_query)
    resp = get(url)
    
    return resp

def obtain_subhalo_data(idn, z):
    """
    Get metadata for a particular subhalo

    (generally called in a loop by obtain_subhalos_data)

    Parameters
    ----------
    idn : int
        Subhalo id number
    z : float
        The redshift

    Returns
    -------
    JSON data
        JSON data returned by the API, formatted as a Python list.
        

    See Also
    --------
    obtain_subhalos_data : Obtain many subhalos worth of metadata.
    """
    url = "http://www.illustris-project.org/api/Illustris-1/snapshots/z=%i/subhalos/%s" % (z, idn)
    json = get(url)
    return json

def obtain_subhalos_data(resp, z, limit=None):
    """
    Get metadata for a particular subhalo

    (generally called in a loop by obtain_subhalos_data)

    Parameters
    ----------
    resp : JSON data
        JSON data in the format returned by obtain_halos_by_mass
    z : float
        The redshift
    limit : int
        (optional, default=None) Limit number of API queries. Useful for debugging only.

    Returns
    -------
    JSON data
        JSON data returned by the API, formatted as a Python list of lists.

    See Also
    --------
    obtain_subhalo_data : Obtain data for one subhalo. This function calls this one in a loop.
    """
    #print(resp['count'])
    if limit is None: limit = int(resp['count'])
    if DEBUG: print('resp - count %i' % int(resp['count']))
    ids = extract_haloids(resp, limit=limit)
    dat = []
    print(' Obtaining metadata for %i subhalos' % limit)
    thresh = int(resp['count']) // 10.
    i = 0
    for idn in ids:
        dat.append(obtain_subhalo_data(idn, z))
        if i % thresh == 0: print('     %i %% done' % (i//thresh * 10))
        i = i + 1
    if DEBUG: print(len(dat))
    return ids, dat

def extract_haloids(resp, limit=None):
    """
    Extract the halo id numbers from the JSON data produced by obtain_halos_by_mass

    Parameters
    ----------
    resp : JSON data
        JSON data in the format returned by obtain_halos_by_mass
    limit : int
        (optional, default=None) Limit number of API queries. Useful for debugging only.

    Returns
    -------
    list of ints
        List of integers, each corresponding to a halo id number.

    See Also
    --------
    obtain_subhalos_data : Obtain many subhalos worth of metadata.
    """
    #print(int(resp['count']))
    if limit is None: limit = int(resp['count'])
    ids = [ resp['results'][i]['id'] for i in range(limit) ]
    return ids

def filter_eliminate_halos(halos_data, ids, to_keep):
    """
    Given an array of metadata and a list of ids, eliminate all unwanted halos.

    Parameters
    ----------
    halos_data : list of JSON data
        Python list containing the JSON data for each subhalo
    ids : list of ints
        Python list of integers corresponding to each halo id
    to_keep : list of ints
        Python list of indexes into halos_data and ids of halos *to keep*.

    Returns
    -------
    list of JSON data
        Python list containing the JSON data for each subhalo
    list of ints
        Python list of integers corresponding to each halo id
    """
    halos_data2 = []
    ids2 = []
    for i in to_keep:
        halos_data2.append(halos_data[i])
        ids2.append(ids[i])
    return ids2, halos_data2

def filter_non_zero_mass(halos_data, ids, arg=None):
    """
    Filter function to remove defective halos.
    
    This filter should *always* be used since halos with empty arrays will crash the code later.

    Parameters
    ----------
    halos_data : list of JSON data
        Python list containing the JSON data for each subhalo
    ids : list of ints
        Python list of integers corresponding to each halo id
    arg : dict
        Dictionary of arguments: only 'epsilon' is accepted, defaults to 1e-5.

    Returns
    -------
    list of JSON data
        Python list containing the JSON data for each subhalo
    list of ints
        Python list of integers corresponding to each halo id
    """
    if arg == None:
        epsilon = 1e-5
    else:
        epsilon = arg['epsilon']
    
    to_keep = []
    
    for i in range(len(ids)):
        halo = halos_data[i]
        if halo['mass_gas'] < epsilon or halo['mass_stars'] < epsilon:
            print(' %i eliminated due to less than epsilon mass' % i)
        else:
            to_keep.append(i)
            
    ids2, halos_data2 = filter_eliminate_halos(halos_data, ids, to_keep)
    return ids2, halos_data2

def filter_mass_range(halos_data, ids, lower, upper, phase):
    """
    Filter function to remove halos with masses out of range.
    
    NOTE: this is _not_ to be put into the filter stack. It is an internal function used by other filters only.

    Parameters
    ----------
    halos_data : list of JSON data
        Python list containing the JSON data for each subhalo
    ids : list of ints
        Python list of integers corresponding to each halo id
    lower : float
    upper : float
        Lower and upper mass limits
    phase : str
        Which particle phase, e.g. 'gas' or 'stars'

    Returns
    -------
    list of JSON data
        Python list containing the JSON data for each subhalo
    list of ints
        Python list of integers corresponding to each halo id
    """
    to_keep = []

    if DEBUG: print('len halos - ', len(halos_data))

    for i in range(len(ids)):
        halo = halos_data[i]
        mstar = halo['mass_%s' % phase]
        #print('len halos 2- ', len(halos_data))
        if mstar < upper and mstar > lower:
            to_keep.append(i)
        else:
            if DEBUG: print('limit exceeded for %i with m%s %f' % (i, phase, mstar))

    ids2, halos_data2 = filter_eliminate_halos(halos_data, ids, to_keep)
    if DEBUG: print('len halos 2- ', len(halos_data2))
    return ids2, halos_data2

def filter_star_mass(halos_data, ids, arg=None):
    """
    Filter function to remove halos with stellar masses out of range.

    Parameters
    ----------
    halos_data : list of JSON data
        Python list containing the JSON data for each subhalo
    ids : list of ints
        Python list of integers corresponding to each halo id
    arg : dict
        Dictionary of arguments: only 'lower' and 'upper' are accepted, corresponding to the lower and upper mass limits
        If set to None the filter becomes a no-op.

    Returns
    -------
    list of JSON data
        Python list containing the JSON data for each subhalo
    list of ints
        Python list of integers corresponding to each halo id
    """
    if arg is None: return ids, halos_data

    lower = arg['lower']
    upper = arg['upper']

    ids2, halos_data2 = filter_mass_range(halos_data, ids, lower, upper, 'stars')
    
    return ids2, halos_data2

def filter_fluid_mass(halos_data, ids, arg=None):
    """
    Filter function to remove halos with fluid masses out of range.

    Parameters
    ----------
    halos_data : list of JSON data
        Python list containing the JSON data for each subhalo
    ids : list of ints
        Python list of integers corresponding to each halo id
    arg : dict
        Dictionary of arguments: only 'lower' and 'upper' are accepted, corresponding to the lower and upper mass limits
        If set to None the filter becomes a no-op.

    Returns
    -------
    list of JSON data
        Python list containing the JSON data for each subhalo
    list of ints
        Python list of integers corresponding to each halo id
    """
    if arg is None: return ids, halos_data

    lower = arg['lower']
    upper = arg['upper']

    ids2, halos_data2 = filter_mass_range(halos_data, ids, lower, upper, 'gas')
    
    return ids2, halos_data2

def filter_haloids(ids, halos_data, filters=None, args=None):
    """
    Main filter routine which runs the full filtering stack.

    Parameters
    ----------
    ids : list of ints
        Python list of integers corresponding to each halo id
    halos_data : list of JSON data
        Python list containing the JSON data for each subhalo
    filters : list of function pointers
        Pointers to each filter function
    arg : list of dicts
        List of argument dictionaries. You *must* use None when a filter needs no args.

    Returns
    -------
    list of ints
        Python list of integers corresponding to each halo id
    list of JSON data
        Python list containing the JSON data for each subhalo
    """
    initial_idn = ids
    initial_data = halos_data
    print('Initial idn set = %i' % len(initial_idn))
    if filters == None: return initial_idn, initial_data #Rejug this to eliminate the copy.

    partial_i = initial_idn
    print(len(partial_i))
    for i in range(len(filters)):
        ftr = filters[i]
        arg = args[i]
        print('Applying filter %s' % str(ftr.__name__)) # I love function pointers.
        partial_i, halos_data = ftr(halos_data, partial_i, arg)
        print('idn set now %i ' %len(partial_i))
        
    return partial_i, halos_data

def calc_com(name, ptypes=(0,)):
    """
    Calculate the centre of mass for all particles in a file.

    Parameters
    ----------
    name : str
        Filename of the hdf5 file in question
    pytypes : tuple of ints
        (optional, default=(0,)) Calculate the CoM over these particle phases only.

    Returns
    -------
    numpy array
        vector representing the location of the CoM
    """
    f = h5py.File(name,'r')
    comtot = numpy.zeros(3)
    mtot = 0.
    for p in ptypes:
        dataset = 'PartType%i' % p
        dat = f[dataset]
        r = dat['Coordinates'][:,:]
        m = dat['Masses'][:]
        rm = r.copy() # numpy being a pain
        rm[:,0] = m[:]*r[:,0]
        rm[:,1] = m[:]*r[:,1]
        rm[:,2] = m[:]*r[:,2]
        comtot[0] += rm[:,0].sum()
        comtot[1] += rm[:,1].sum()
        comtot[2] += rm[:,2].sum()
        mtot += m.sum()
    com = comtot/mtot
    if DEBUG: print('CoM: ', com)
    f.close()
    return com

def calc_com_vel(name, ptypes=(0,)):
    """
    Calculate the centre of mass velocity for all particles in a file.

    Parameters
    ----------
    name : str
        Filename of the hdf5 file in question
    pytypes : tuple of ints
        (optional, default=(0,)) Calculate the CoM velocity over these particle phases only.

    Returns
    -------
    numpy array
        vector representing the CoM velocity
    """
    f = h5py.File(name,'r')
    comtot = numpy.zeros(3)
    mtot = 0.
    for p in ptypes:
        dataset = 'PartType%i' % p
        dat = f[dataset]
        v = dat['Velocities'][:,:]
        m = dat['Masses'][:]
        vm = v.copy() # numpy being a pain
        vm[:,0] = m[:]*v[:,0]
        vm[:,1] = m[:]*v[:,1]
        vm[:,2] = m[:]*v[:,2]
        comtot[0] += vm[:,0].sum()
        comtot[1] += vm[:,1].sum()
        comtot[2] += vm[:,2].sum()
        mtot += m.sum()
    com = comtot/mtot
    if DEBUG: print('CoM vel: ', com)
    f.close()
    return com

def transpose_cutout(outname, idn, ptypes=(0,), comptypes=None):
    """
    Transpose a cutout into the zero momentum CoM frame.
    
    NOTE: function edits HDF5 files in place.

    Parameters
    ----------
    outname : str
        Filename of the hdf5 file in question
    idn : int
        No-op. Retained for backwards compatibility. Ignored.
    ptypes : tuple of ints
        (optional, default=(0,)) Transpose these phases. You should set this to include all phases you download otherwise the results will be weird.
    comptypes : tuple of ints
        (optional, default=None, i.e. use the defaults in the CoM functions) Calculate the CoM velocity over these particle phases only.

    Returns
    -------
    str
        Filename, identical to outname, for backwards compatibility.
    """
    
    if comptypes == None: comptypes=ptypes
    
    com = calc_com(outname, ptypes=ptypes)
    comvel = calc_com_vel(outname, ptypes=ptypes)
    
    f = h5py.File(outname,'r+')
    for p in ptypes:
        dataset = 'PartType%i' % p
        dat = f[dataset]
        r = dat['Coordinates'][:,:]
        v = dat['Velocities'][:,:]
        if DEBUG: print('mean r = ',r[:,0].mean())
        if DEBUG: print('mean x = ',dat['Coordinates'][:,0].mean())
        if DEBUG: print('mean vx = ',dat['Velocities'][:,0].mean())
        #print(r[:,0].mean())
        r[:,0] = r[:,0] - com[0]#posx
        r[:,1] = r[:,1] - com[1]#posy
        r[:,2] = r[:,2] - com[2]#posz
        v[:,0] = v[:,0] - comvel[0]#velx
        v[:,1] = v[:,1] - comvel[1]#vely
        v[:,2] = v[:,2] - comvel[2]#velz
        dat['Coordinates'][:,:] = r[:,:] # poor life choices early in h5py's development require this.
        dat['Velocities'][:,:] = v[:,:]
        if DEBUG: print('mean r = ',r[:,0].mean())
        if DEBUG: print('mean x = ',dat['Coordinates'][:,0].mean())
        if DEBUG: print('mean vx = ',dat['Velocities'][:,0].mean())
    
    #outs = outname.split('.')
    #outfile = '%s_transpose.%s' %(outs[0], outs[1])
    #create_new_file(outfile, f)
    f.close()
    outfile = outname
    return outfile

def obtain_cutout(idn, z, dmcols=None, starcols=None, gascols=None, bhcols=None):
    """
    Download cutout files. Normally called in a loop by obtain_cutouts()

    Parameters
    ----------
    idn : int
        ID number of the halo
    z : float
        Redshift
    dmcols, starcols, gascols, bhcols : lists of strs
        Column names to download (see illustris documentation)

    Returns
    -------
    str
        Filename HDF5 file has been saved to.
    """
    print(' Obtaining cutout of halo %i at redshift %f' % (idn, z))
    print(' Columns requested:', end='') # In the ideal world we would use the * operator here but None is not iterable
    cutout_query = ''
    if dmcols is not None:
        cutout_query += 'dm=%s&' % ','.join(dmcols)
        print(*dmcols, end=' ')
    if starcols is not None:
        cutout_query += 'stars=%s&' % ','.join(starcols)
        print(*starcols, end=' ')
    if gascols is not None:
        cutout_query += 'gas=%s&' % ','.join(gascols)
        print(*gascols, end=' ')
    if bhcols is not None:
        cutout_query += 'bhs=%s&' % ','.join(bhcols)
        print(*bhcols, end=' ')
    print('')
    
    url = "http://www.illustris-project.org/api/Illustris-1/snapshots/z=%i/subhalos" % z
    full_url = '%s/%i/cutout.hdf5?%s' % (url, idn, cutout_query)
    
    outname = 'cutout_%i.hdf5' % idn
    ret = os.system('wget "%s" --header="api-key: %s" -O %s -o wget.log' % (full_url, apikey, outname))
    if ret != 0:
        print(' Wget returned a non-zero value. Check the wget.log file!')
        sys.exit(ret)
    copyname = 'cutout_%i.hdf5.orig' % idn
    shutil.copy2(outname, copyname)
    
    tname = transpose_cutout(outname, idn, ptypes=(0,4), comptypes=(4,))
    
    return tname
    
    if DEBUG: print(full_url)
    #resp = get(full_url)
    
    #return resp

def obtain_cutouts(ids, z, dmcols=None, starcols=None, gascols=None, bhcols=None, count=None):
    """
    Download cutout files.

    Parameters
    ----------
    ids : list of ints
        Halo ID numbers
    z : float
        Redshift
    dmcols, starcols, gascols, bhcols : lists of strs
        Column names to download (see illustris documentation)

    Returns
    -------
    list of strs
        List of filenames the HDF5 files have been saved to.
    """
    r = []
    if count is None:
        count = len(ids)
        print('Cutouts will be obtained for %i halos' % count)
    for i in range(count):
        r.append(obtain_cutout(ids[i], z, dmcols=dmcols, starcols=starcols, gascols=gascols, bhcols=bhcols))
    return r

def process_cutout_for_splash(cutout, outfile, ptypes=(0,4)):
    """
    Process a cutout file so that it can be plotted in SPLASH. I.e. add a header to it.
    
    Normally called in a loop by process_cutouts_for_splash()

    Parameters
    ----------
    cutout : str
        Input filename
    outfile : str
        Output filename (conventionally foo_splash.hdf5)
    ptypes : tuple of ints
        (optional, default=(0,4)) Which phases to store in the output file. Probably should be the same as in cutout but included for corner cases.

    Returns
    -------
    None.
    """
    print(' "splashifying" %s, output will be saved to %s' % (cutout, outfile))
    f = h5py.File(cutout, 'r')
    np = [0, 0, 0, 0, 0, 0]
    for p in ptypes:
        x = f['PartType%i'%p]['Coordinates'][:,0]
        np[p] = x.shape[0]
    
    head = f['Header'].attrs
    ###head['NumFiles']
    # recyclable headers
    Massarr = head['MassTable']
    Time = head['Time']
    Redshift = head['Redshift']
    Flag_Sfr = head['Flag_Sfr']
    Flag_Feedback = head['Flag_Feedback']
    FlagCooling = head['Flag_Cooling'] # NOT A TYPO
    BoxSize = head['BoxSize']
    Omega0 = head['Omega0']
    OmegaLambda = head['OmegaLambda']
    HubbleParam = head['HubbleParam']
    FlagMetals = head['Flag_Metals'] # NOT A TYPO
    
    # Not recyclable
    Npart = np
    Nall = Npart
    NallHW = [0, 0, 0, 0, 0, 0]
    NumFiles = 1
    
    fnew = h5py.File(outfile, 'w')
    
    g = fnew.create_group('Header')
    ga = g.attrs
    ga['MassTable'] = Massarr
    ga['Time'] = Time
    ga['Redshift'] = Redshift
    ga['Flag_Sfr'] = Flag_Sfr
    ga['Flag_Feedback'] = Flag_Feedback
    ga['FlagCooling'] = FlagCooling
    ga['BoxSize'] = BoxSize
    ga['Omega0'] = Omega0
    ga['OmegaLambda'] = OmegaLambda
    ga['HubbleParam'] = HubbleParam
    ga['FlagMetals'] = FlagMetals
    ga['Npart'] = Npart
    ga['NumPart_Total'] = Nall
    ga['NumPart_Total_HighWord'] = NallHW
    ga['NumFiles'] = NumFiles
    ga['NumFilesPerSnapshot'] = NumFiles
    ga['NumPart_ThisFile'] = Npart
    
    #fnew.create_group('PartType0')
    #for d in f['PartType0'].keys():
    #    fnew['PartType0'][d] = f['PartType0'][d].copy()
    for p in ptypes:
        f.copy('PartType%i'%p, fnew)
    
    if DEBUG: print(Npart)
    
    f.close()
    fnew.close()

def process_cutouts_for_splash(cutouts):
    """
    Process cutouts file so that they can be plotted in SPLASH.
    
    You should call this after all destructive modifications to the cutout have been done.

    Parameters
    ----------
    cutouts : list of strs
        Input filenames

    Returns
    -------
    None.
    """
    for c in cutouts:
        outs = c.split('.')
        outfile = '%s_splash.%s' %(outs[0], outs[1])
        process_cutout_for_splash(c, outfile)


def cross_product(a, b):
    """
    Compute the vector product of a and b in situations where the numpy built-in is inappropriate.

    Parameters
    ----------
    a : numpy array
    b : numpy array
        Two length 3 arrays

    Returns
    -------
    numpy array
        Vector product of a and b
    """
    i = a[1]*b[2] - a[2]*b[1]
    j = a[0]*b[2] - a[2]*b[0]
    k = a[0]*b[1] - a[1]*b[0]
    out = numpy.zeros(3)
    out[0] = i
    out[1] = j
    out[2] = k
    return out

def calc_angular_momentum(f, ptypes=(0,), rlim=None):
    """
    Compute the angular momentum of a cutout.

    Parameters
    ----------
    f : h5py.File object
        cutout file handle
    ptypes : tuple of ints
        (optional, default=(0,)) Which particle phases to compute L over.
    rlim : float
        (optional, default=30.) Limit the angular momentum calculation to particles within this radius. Beware the default -- it is probably not what you want. Set to None to include everything.

    Returns
    -------
    numpy array
        Angular momentum
    float
        Modulus of the angular momentum
    """
    if rlim == None: rlim = 1000.
    out = []
    for p in ptypes:
        dataset = 'PartType%i' % p
        dat = f[dataset]
        r = dat['Coordinates'][:,:]
        m = dat['Masses'][:].reshape(-1,1) # Felix told me to add this
        v = dat['Velocities'][:,:]
        linmomenta = m*v
        if DEBUG: print(linmomenta)
        L = numpy.copy(linmomenta)
        L[:,:] = 0.
        r2 = r[:,0]**2 + r[:,1]**2 + r[:,2]**2
        radius = numpy.sqrt(r2)
        for i in range(linmomenta.shape[0]): # FIXME: redo this using numpy slices
            if radius[i] < rlim:
                L[i,:] = cross_product(r[i,:], linmomenta[i,:])
        totL = L.sum(axis=0)
        if DEBUG: print(totL.shape, totL)
        modL = math.sqrt(totL[0]*totL[0] +  totL[1]*totL[1] +  totL[2]*totL[2])
        if DEBUG: print(modL)
        out.append([totL, modL])
    if len(out) == 1: # Handle the simple case 'simply'
        return out[0][0], out[0][1]
    else:
        return out
    
def create_new_file(outname, f, ptypes=(0, 4)):
    fnew = h5py.File(outname, 'w')
    f.copy('Header', fnew)
    for p in ptypes:
        dataset = 'PartType%i' % p
        f.copy(dataset, fnew)
    fnew.close()
    
def skew_sym_vecp(v):
    """
    Try not to think too hard about this.

    Parameters
    ----------
    v : numpy array

    Returns
    -------
    numpy array
    """
    U = numpy.zeros((3,3))
    U[0,1] = -v[2]
    U[0,2] = v[1]
    U[1,0] = v[2]
    U[1,2] = v[0]
    U[2,0] = v[1]
    U[2,1] = v[0]
    return U


def detilt_orig(fname, ptypes=(0, 4), angmomptypes=(4,)):
    """
    The original now deprecated detilting algorithm. DO NOT USE.

    Parameters
    ----------
    fname : str
        filename
    ptypes : tuple of ints
        (optional, default=(0,4)) Which particle phases to detilt -- probably should be all phases in the file to avoid wierdness
    angmomptypes : tuple of ints
        (optional, default=(4,)) Which phases to compute the angular momentum over. Need not match ptypes.

    Returns
    -------
    None
    """
    print(' Detilting %s' % fname)
    f = h5py.File(fname,'r+')
    totL, modL = calc_angular_momentum(f, ptypes=angmomptypes)
    Lhat = totL/modL # this is the direction the total angular momentum vector points
    if DEBUG: print(Lhat)
    # Change the coordinate basis so that Lhat is the new x,y,z axis
    zhat = numpy.array([0,0,1])
    v = cross_product(zhat, Lhat)
    if DEBUG: print(v)
    s = math.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])
    c = numpy.dot(zhat, Lhat)
    vstar = skew_sym_vecp(v)
    
    R = numpy.identity(3) + vstar + numpy.matmul(vstar,vstar) * (1 - c)/(s*s)
    if DEBUG: print(R)
    
    for p in ptypes:
        # I havent yet divined how to not do this as a loop
        dataset = 'PartType%i' % p
        dat = f[dataset]
        r = dat['Coordinates'][:,:]
        v = dat['Velocities'][:,:]
        # NOTE: no other vector quantites (are there any?) are rotated.
        n = r.shape[0]
        for i in range(n):
            r1 = r[i,:]
            v1 = v[i,:]
            rnew = numpy.matmul(R, r1)
            vnew = numpy.matmul(R, v1)
            r[i,:] = rnew
            v[i,:] = vnew
        dat['Coordinates'][:,:] = r
        dat['Velocities'][:,:] = v
        
    totL, modL = calc_angular_momentum(f, ptypes=angmomptypes)
    if DEBUG: print(totL, modL)
        
    
    #outs = fname.split('.')
    #outfile = '%s_rotated.%s' %(outs[0], outs[1])
    #create_new_file(outfile, f)
    f.close()
    
    

def dot_product(a, b):
    return  a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

def get_stellar_halfmass(halo_data):
    # This lever should be used if you want to fudge things
    ratio = 1.
    return ratio * float(halo_data['halfmassrad_stars'])

# def get_stellar_halfmass2(fname,halo_data,angmomptypes=(4,)):
#     # This lever should be used if you want to fudge things
#     ratio = 1.
#     f = h5py.File(fname,'r+')
#     rlim = get_stellar_halfmass2(fname,halo_data)
#     totL, modL = calc_angular_momentum(f, ptypes=angmomptypes, rlim=rlim)
#     Lhat = totL/modL
#     return ratio * float(halo_data['halfmassrad_stars']),Lhat
    

def detilt(fname, halo_data, ptypes=(0, 4), angmomptypes=(4,), rlim_func=get_stellar_halfmass):
    """
    The new, rotation matrix based, detilt algorithm. USE THIS ONE.

    Parameters
    ----------angmomptypes=(4,)
    fname : str
        filename of halo.
    halo_data : JSON data
        JSON data corresponding to this halo
    ptypes : tuple of ints
        (optional, default=(0,4)) Which particle phases to detilt -- probably should be all phases in the file to avoid wierdness
    angmomptypes : tuple of ints
        (optional, default=(4,)) Which phases to compute the angular momentum over. Need not match ptypes.
    rlim_func : function pointer
        (optional, default=get_stellar_halfmass) Set this to compute change the radius cut in the angular momentum calculation. The default is probably adequate always. Set to None to include everything.

    Returns
    -------
    None
    """
    if DEBUG: print(halo_data['id'])
    
    if rlim_func is not None: 
        rlim = rlim_func(halo_data)
        print('Limiting angular momentum calculation to r < %g' % rlim)

    else:
        rlim = None
        print(' NO rlim applied.')
    
    f = h5py.File(fname,'r+')
    
    totL, modL = calc_angular_momentum(f, ptypes=angmomptypes, rlim=rlim)
    Lhat = totL/modL # this is the direction the total angular momentum vector points
    
    print(Lhat)
   
    
    
    xhat = [1., 0., 0.]
    yhat = [0., 1., 0.]
    zhat = [0., 0., 1.]
    cosa = dot_product(Lhat, xhat)
    cosb = dot_product(Lhat, yhat)
    cosc = dot_product(Lhat, zhat)
    
    alpha = math.acos(cosc)
    beta = math.acos(cosa)
    gamma = math.acos(cosb)
    
    if DEBUG: print('Angles: ', math.degrees(math.acos(cosa)), math.degrees(math.acos(cosb)), math.degrees(math.acos(cosc)))
    xmat = [[1., 0., 0.], [0., math.cos(alpha), -math.sin(alpha)], [0., math.sin(alpha), math.cos(alpha)]]
    ymat = [[math.cos(beta), 0., math.sin(0)], [0., 1., 0.], [-math.sin(beta), 0., math.cos(beta)]]
    zmat = [[math.cos(gamma), -math.sin(gamma), 0.], [math.sin(gamma), math.cos(gamma), 0.], [0., 0., 1.]]
    
    Rp = numpy.matmul(xmat, ymat)
    R = numpy.matmul(Rp, zmat)
    #R = numpy.matmul(zmat, ymat)
    #R = Rp
    # I havent yet divined how to not do this as a loop
    for p in ptypes:
        dataset = 'PartType%i' % p
        dat = f[dataset]
        r = dat['Coordinates'][:,:]
        v = dat['Velocities'][:,:]
        # NOTE: no other vector quantites (are there any?) are rotated.
        n = r.shape[0]
        for i in range(n):
            r1 = r[i,:]
            v1 = v[i,:]
            rnew = numpy.matmul(R, r1)
            vnew = numpy.matmul(R, v1)
            r[i,:] = rnew
            v[i,:] = vnew
        dat['Coordinates'][:,:] = r
        dat['Velocities'][:,:] = v
        
    totL, modL = calc_angular_momentum(f, ptypes=angmomptypes, rlim=rlim)
    if DEBUG: print(totL, modL)
        
    
    #outs = fname.split('.')
    #outfile = '%s_rotated.%s' %(outs[0], outs[1])
    #create_new_file(outfile, f)
    f.close()
    
    
    return Lhat



def detilt_files(resps, halos_data, ptypes=(0, 4), detilt_func=detilt):
    """
    The new, rotation matrix based, detilt algorithm. USE THIS ONE.

    Parameters
    ----------
    resps : list of strs
        List of filenames of halos.
    halos_data : JSON data
        List of JSON data corresponding to these halos
    ptypes : tuple of ints
        (optional, default=(0,4)) Which particle phases to detilt -- probably should be all phases in the file to avoid wierdness
    detilt_func : function pointer
        (optional, default=detilt) Which detilting algorithm to use. *Use the default*

    Returns
    -------
    None
    """

    LHAT = []
    for i, f in enumerate(resps):
        Lhat = detilt_func(f, halos_data[i], ptypes=ptypes)
        LHAT.append(Lhat)
        
    return LHAT
        

    
# def return_Lhat(fname,halo_data,angmomptypes=(4,),rlim_func2=get_stellar_halfmass2):
#     f = h5py.File(fname,'r+')
    
#     if rlim_func2 is not None: 
#         rlim = rlim_func2(halo_data)

#     else:
#         rlim = None
#         print(' NO rlim applied.')
        
#     totL, modL = calc_angular_momentum(f, ptypes=angmomptypes, rlim=rlim)
#     Lhat = totL/modL 
#     return Lhat
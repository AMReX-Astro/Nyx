import numpy as np
import sys
import h5py

class AMReXParticleHeader(object):
    '''

    This class is designed to parse and store the information 
    contained in an AMReX particle header file. 

    Usage:

        header = AMReXParticleHeader("plt00000/particle0/Header")
        print(header.num_particles)
        print(header.version_string)

    etc...

    '''

    def __init__(self, header_filename):

        self.real_component_names = []
        self.int_component_names = []
        with open(header_filename, "r") as f:
            self.version_string = f.readline().strip()

            particle_real_type = self.version_string.split('_')[-1]
            particle_real_type = self.version_string.split('_')[-1]
            if particle_real_type == 'double':
                self.real_type = np.float64
            elif particle_real_type == 'single':
                self.real_type = np.float32
            else:
                raise RuntimeError("Did not recognize particle real type.")
            self.int_type = np.int32

            self.dim = int(f.readline().strip())
            self.num_int_base = 2
            self.num_real_base = self.dim
            self.num_real_extra = int(f.readline().strip())
            for i in range(self.num_real_extra):
                self.real_component_names.append(f.readline().strip())
            self.num_int_extra = int(f.readline().strip())
            for i in range(self.num_int_extra):
                self.int_component_names.append(f.readline().strip())
            self.num_int = self.num_int_base + self.num_int_extra
            self.num_real = self.num_real_base + self.num_real_extra
            self.is_checkpoint = bool(int(f.readline().strip()))
            self.num_particles = int(f.readline().strip())
            self.max_next_id = int(f.readline().strip())
            self.finest_level = int(f.readline().strip())
            self.num_levels = self.finest_level + 1

            if not self.is_checkpoint:
                self.num_int_base = 0
                self.num_int_extra = 0
                self.num_int = 0

            self.grids_per_level = np.zeros(self.num_levels, dtype='int64')
            for level_num in range(self.num_levels):
                self.grids_per_level[level_num] = int(f.readline().strip())

            self.grids = [[] for _ in range(self.num_levels)]
            for level_num in range(self.num_levels):
                for grid_num in range(self.grids_per_level[level_num]):
                    entry = [int(val) for val in f.readline().strip().split()]
                    self.grids[level_num].append(tuple(entry))

                    
def read_amrex_binary_particle_file(fn, ptype="particle0"):
    '''

    This function returns the particle data stored in a particular 
    plot file and particle type. It returns two numpy arrays, the
    first containing the particle integer data, and the second the
    particle real data. For example, if a dataset has 3000 particles,
    which have two integer and five real components, this function will
    return two numpy arrays, one with the shape (3000, 2) and the other
    with the shape (3000, 5).

    Usage:
    
        idata, rdata = read_particle_data("plt00000", "particle0")

    '''
    base_fn = fn + "/" + ptype
    header = AMReXParticleHeader(base_fn + "/Header")
    
    idtype = "(%d,)i4" % header.num_int    
    if header.real_type == np.float64:
        fdtype = "(%d,)f8" % header.num_real
    elif header.real_type == np.float32:
        fdtype = "(%d,)f4" % header.num_real
    
    ip = 0
    for lvl, level_grids in enumerate(header.grids):
        for (which, count, where) in level_grids:
            print(which,count,where)

            if count == 0: continue
            fn = base_fn + "/Level_%d/DATA_%05d" % (lvl, which)

            with open(fn, 'rb') as f:
                f.seek(where)
                if header.is_checkpoint:
                    ints   = np.fromfile(f, dtype = idtype, count=count)
                    #idata[ip:ip+count] = ints

                floats = np.fromfile(f, dtype = fdtype, count=count)
                yield ints,floats,count


if __name__ == "__main__":
    import glob

    if len(sys.argv) < 3:
        print("Usage: python plotfile2HDF5_particles_amr.py <plotfile> <outfile>")

    fn = sys.argv[1]
    ptype = 'DM' #particle type

    base_fn = fn + "/" + ptype
    header = AMReXParticleHeader(base_fn + "/Header")
    nopart = header.num_particles
    print("num_int:", header.num_int)
    print("num_real:", header.num_real)
    print("num_particles:", nopart)
    print("grids:", sum(len(g) for g in header.grids))
        
    # Prepare output
    outfn = sys.argv[2]
    f     = h5py.File(outfn, 'w')
    variables = ['x','y','z','m','vx','vy','vz']
    dset = {}
    for v in variables:
        dset[v] = f.create_dataset(v, (nopart,), dtype='f4')

    ip = 0
    for i, (ints, floats, count) in enumerate(read_amrex_binary_particle_file(fn, ptype)):
        print("Grid %d, filling %d:%d" % (i,ip,ip+count))
        for it,v in enumerate(variables):
           dset[v][ip:ip+count] = floats[:,it]
        ip += count

        # Report memory usage
        import os
        import psutil
        process = psutil.Process(os.getpid())
        print("Memory usage: %d MB" % (process.memory_info().rss / 1024**2))


    print ('x:  ', np.mean(dset['x']),np.amin(dset['x']),np.amax(dset['x']))
    print ('y:  ', np.mean(dset['y']),np.amin(dset['y']),np.amax(dset['y']))
    print ('z:  ', np.mean(dset['z']),np.amin(dset['z']),np.amax(dset['z']))
    print ('m:  ', np.mean(dset['m']))
    print ('vx: ', np.mean(dset['vx']))
    print ('vy: ', np.mean(dset['vy']))
    print ('vz: ', np.mean(dset['vz']))

    print ('Problem particles: ', (np.where(dset['m']==0)))
    f.close()

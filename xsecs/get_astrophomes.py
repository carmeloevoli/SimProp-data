import sys
sys.path.append('/Users/carmeloevoli/Work/libs/AstroPhoMes/')
from config import *
sys.path.append(global_path)
from photomeson_lib.photomeson_models import *

def make_np_table():
    spm = SingleParticleModel()

    # neutron
    cs = spm.cs_nonel(100)
    E, n = cs[0], cs[1]

    # proton
    cs = spm.cs_nonel(101)
    assert(np.array_equal(E, cs[0]))
    p = cs[1]

    # print
    f = open('tables/xsecs_photopion_inelastic_astrophomes.txt', 'w')
    f.write("#epsilon_r [GeV] - (n,g) - (p,g) [mb]\n")
    for i,j,k in zip(E, n, p):
        f.write('{:.5e} {:.5e} {:.5e}\n'.format(i,j,k))
    f.close()
    
if __name__== "__main__":
    print ("remind: workon phomes")
    make_np_table()



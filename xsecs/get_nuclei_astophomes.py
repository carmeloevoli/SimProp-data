import sys
sys.path.append('/Users/carmeloevoli/Work/libs/AstroPhoMes/')
from config import *
sys.path.append(global_path)
from photomeson_lib.photomeson_models import *

#    dict = {
#        "neutron" : 100,
#        "proton" : 101,
#        "He4" : 402,
#        "C12" : 1206,
#        "Ca40": 4020,
#        "Fe56": 5626,
#    }

def make_nucleus_table(id, filename):
    em = EmpiricalModel()
    cs_em = em.cs_nonel(id)

    spm = SingleParticleModel()
    cs_spm = spm.cs_nonel(id)

    rdm = ResidualDecayModel()
    cs_rdm = rdm.cs_nonel(id)

    filename = 'tables/' + filename
    # print
    f = open(filename, 'w')
    f.write("#epsilon_r [GeV] - SPM - Empirical Model - Residual Decay Model [mb]\n")
    for i,j,k,z in zip(cs_em[0], cs_spm[1], cs_em[1], cs_rdm[1]):
        f.write('{:.5e} {:.5e} {:.5e} {:.5e}\n'.format(i,j,k,z))
    f.close()
    
def make_pion_table(filename):
    em = EmpiricalModel()
    cs_piplus = em.cs_incl(101, 2)
    cs_piminus = em.cs_incl(101, 3)
    cs_pi0 = em.cs_incl(101, 4)
    cs_inel = em.cs_nonel(101)
    
    filename = 'tables/' + filename
    f = open(filename, 'w')
    f.write("#epsilon_r [GeV] - pi0 - pi+ - pi- [mb]\n")
    for i,j,k,z,t in zip(cs_pi0[0], cs_pi0[1], cs_piplus[1], cs_piminus[1], cs_inel[1]):
        f.write('{:.5e} {:.5e} {:.5e} {:.5e} {:.5e}\n'.format(i,j,k,z,t))
    f.close()

if __name__== "__main__":
    print ("remind: workon phomes")
#    make_nucleus_table(1206, 'xsecs_photopion_C12_astrophomes.txt')
#    make_nucleus_table(4020, 'xsecs_photopion_Ca40_astrophomes.txt')
#    make_nucleus_table(5626, 'xsecs_photopion_Fe56_astrophomes.txt')
    make_pion_table('xsecs_photopion_pions_astrophomes.txt')

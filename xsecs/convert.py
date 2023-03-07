import numpy as np
from scipy import integrate

protonmass = 0.93827208816 # GeV
pionmass = 0.1349768 # GeV

def eps2s(eps):
    return protonmass * protonmass + 2. * protonmass * eps

def load_photopion_sophia(filename):
    eps, sigma = np.loadtxt(filename, usecols=(0,1), unpack=True, skiprows=2)
    sigma /= 1e3 # microbarn -> mbarn
    s = eps2s(eps)
    return s, sigma
    
def load_photopion_astrophomes(filename, doProton):
    eps, neutron, proton = np.loadtxt(filename, usecols=(0,1,2), unpack=True, skiprows=1)
    s = eps2s(eps)
    sigma = proton if doProton else neutron
    sigma /= 1e3 # microbarn -> mbarn
    return s, sigma

def load_photopion_v2r4():
    filename = 'tables/xsecs_photopion_v2r4_lows.txt'
    s_low, sigma_low = np.loadtxt(filename, usecols=(0,1), unpack=True, skiprows=1)

    filename = 'tables/xsecs_photopion_v2r4_highs.txt'
    logsqrts, logsigma = np.loadtxt(filename, usecols=(0,1), unpack=True, skiprows=1)
    s_high = (np.power(10., logsqrts))**2.0
    sigma_high = np.power(10., logsigma)
    
    return np.concatenate((s_low, s_high)), np.concatenate((sigma_low, sigma_high))

def compute_sigma(E, t, h1, x1, w1, c):
    size = len(E)
    epsilon1 = 30. # MeV
    sigma = np.zeros(size)
    for i in range(size):
        if E[i] > epsilon1:
            sigma[i] = c
        elif E[i] > t:
            sigma[i] = h1 * np.exp(-(E[i] - x1)**2.0 / w1)
    return sigma

def dump_v2r4(filename, doAlpha):
    url = 'tables/xsect_Gauss2_TALYS-restored.txt'
    
    if doAlpha:
        A, Z, t, h1, x1, w1, c = np.loadtxt(url, skiprows=1, unpack=True, usecols=(0,1,7,8,9,10,11))
    else:
        A, Z, t, h1, x1, w1, c = np.loadtxt(url, skiprows=1, unpack=True, usecols=(0,1,2,3,4,5,6))
        
    print("#t range : ", min(t), max(t))
    
    E = np.logspace(0, 2, 100) # MeV
    
    size = len(A)
    
    f = open(filename, 'w')
    f.write("#A - Z - sigma [mb]\n")
    
    for i in range(size):
        sigma = compute_sigma(E, t[i], h1[i], x1[i], w1[i], c[i])
        f.write('{},{},'.format(int(A[i]),int(Z[i])))
        for s in sigma[:-1]:
            f.write('{:.5e},'.format(s))
        f.write('{:.5e}'.format(sigma[-1]))
        f.write('\n')
    f.close()
    print ('dumped xsecs table on ', filename)

def dump_xsecs_file(s, sigma, filename):
    assert(s.size == sigma.size)
    file = open(filename, "w")
    file.write(f'# size in energy {s.size}\n')
    file.write(f'# s [GeV^2] - sigma [mbarn] - phi(s) [GeV^4 mbarn]\n')
    for i in range(s.size):
        if i > 0:
            phi = integrate.simpson(s[0:i] * (s[0:i] - protonmass**2.) * sigma[0:i], np.log(s[0:i]), even='first')
        else:
            phi = 0.
        file.write(f'{s[i]:10.5e}, {sigma[i]:10.5e}, {phi:10.5e}\n')
    file.close()
    print ('dumped xsecs table on ', filename)
    
if __name__== "__main__":
    s, sigma = load_photopion_sophia('tables/xs_proton.txt')
    dump_xsecs_file(s, sigma, 'xsecs_photopion_proton_sophia.txt')

    s, sigma = load_photopion_sophia('tables/xs_neutron.txt')
    dump_xsecs_file(s, sigma, 'xsecs_photopion_neutron_sophia.txt')

    dump_v2r4('xsecs_photodisintegration_v2r4_alpha.txt', True)
    dump_v2r4('xsecs_photodisintegration_v2r4_singlenucleon.txt', False)

#    s, sigma = load_photopion_astrophomes('tables/xsecs_photopion_inelastic_astrophomes.txt', doProton=True)
#    dump_xsecs_file(s, sigma, 'xsecs_photopion_proton_astrophomes.txt')
#
#    s, sigma = load_photopion_astrophomes('tables/xsecs_photopion_inelastic_astrophomes.txt', doProton=False)
#    dump_xsecs_file(s, sigma, 'xsecs_photopion_neutron_astrophomes.txt')
#
#    s, sigma = load_photopion_v2r4()
#    dump_xsecs_file(s, sigma, 'xsecs_photopion_proton_v2r4.txt')

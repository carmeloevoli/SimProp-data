import numpy as np
import math
from scipy import integrate

def wavelength2energy(l_mu):
    hc = 1.23984198 # [eV mu]
    E = hc / l_mu # [eV]
    return E

def brightness2density(b, l):
    E = wavelength2energy(l)
    cLight = 299792458.0 # [m / s]
    n = 1e10 / 1.60218 * 4. * math.pi / cLight * b / E / E # [eV^-1 m^-3]
    return n
        
def load_Dominguez2011_file(z, original):
    filename = 'tables/' + original
    w, b = np.loadtxt(filename, usecols=(0,1), unpack=True, skiprows=10)
    eps = wavelength2energy(w)
    eps = np.array(eps[::-1])
    density = np.zeros((z.size, eps.size))
    for i in range(z.size):
        w, b = np.loadtxt(filename, usecols=(0, i+1), unpack=True, skiprows=10)
        n = brightness2density(b, w)
        density[i] = np.power(1. + z[i], 3.) * n[::-1]
    return eps, density
    
def load_Saldana2021_error_file(z, doLower):
    fiducial = 'tables/ebl_saldana2021_fiducial.txt'
    error = 'tables/ebl_saldana2021_uncertainties.txt'
    w, b = np.loadtxt(fiducial, usecols=(0,1), unpack=True, skiprows=10)
    eps = wavelength2energy(w)
    eps = np.array(eps[::-1])
    density = np.zeros((z.size, eps.size))
    for i in range(z.size):
        w, b_fiducial = np.loadtxt(fiducial, usecols=(0, i+1), unpack=True, skiprows=10)
        w, b_error = np.loadtxt(error, usecols=(0, i+1), unpack=True, skiprows=10)
        n_fiducial = brightness2density(b_fiducial, w)
        n_error = brightness2density(b_error, w)
        if doLower:
            density[i] = np.power(1. + z[i], 3.) * (n_fiducial[::-1] - n_error[::-1])
        else:
            density[i] = np.power(1. + z[i], 3.) * (n_fiducial[::-1] + n_error[::-1])
    return eps, density

def load_Gilmore2012_file(z, original):
    def wavelength2energy(w):
        hc = 1.23984198e4 # eV A
        return hc / w # eV
    def brightness2density(b, w):
        cLight = 29970254700. # cm s^-1
        b *= w # erg s^-1 cm^-2 sr^-1
        b *= 4. * math.pi / cLight # erg cm^-3
        E = wavelength2energy(w)
        b /= E * E # eV^-2 erg cm^-3
        b *= 1e6 * 6.242e11 # eV^-1 m^-3
        return b

    filename = 'tables/' + original
    w, b = np.loadtxt(filename, usecols=(0,1), unpack=True, skiprows=2)
    eps = wavelength2energy(w)
    eps = np.array(eps[::-1])
    density = np.zeros((z.size, eps.size))
    for i in range(z.size):
        w, b = np.loadtxt(filename, usecols=(0, i+1), unpack=True, skiprows=2)
        n = brightness2density(b, w)
        density[i] = n[::-1] # / (1. + z[i])**3.
    return eps, density
    
def compute_Igamma_integral(eps, n):
    size = len(n)
    assert(size == len(n))
    I = []
    for i in range(size):
        I.append(integrate.simpson(np.array(n[i:] / eps[i:]), np.log(eps[i:]), even='first'))
    return np.array(I)
    
def dump_ebl_file(eps, z, density, filename):
    assert(density.size == eps.size * z.size)
    file = open(filename, "w")
    file.write(f'# size in redshift {z.size} - size in energy {eps.size}\n')
    file.write(f'# z - eps [eV] - n [1/eV/m3] - I_gamma [1/eV2/m3]\n')
    for i in range(z.size):
        I_gamma = compute_Igamma_integral(eps, density[i])
        for j in range(eps.size):
            file.write("%6.2f, %10.6e, %10.6e, %10.6e\n" % (z[i], eps[j], density[i][j], I_gamma[j]))
    file.close()
    print ('dumped ebl table on ', filename)

if __name__== "__main__":
    # Dominguez2011 fiducial
    z = np.array([0, 0.01, 0.03, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0,
                  1.25, 1.5, 2.0, 2.5, 3.0, 3.9])
    eps, density = load_Dominguez2011_file(z, 'ebl_dominguez2011_fiducial.txt')
    dump_ebl_file(eps, z, density, 'ebl_Dominguez2011_fiducial.txt')

    # Dominguez2011 lower
    eps, density = load_Dominguez2011_file(z, 'ebl_dominguez2011_lower.txt')
    dump_ebl_file(eps, z, density, 'ebl_Dominguez2011_lower.txt')

    # Dominguez2011 upper
    eps, density = load_Dominguez2011_file(z, 'ebl_dominguez2011_upper.txt')
    dump_ebl_file(eps, z, density, 'ebl_Dominguez2011_upper.txt')

    # Saldana2021 fiducial
    z = np.array([0., 0.01, 0.03, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0,
                  1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4,
                  3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.])
    eps, density = load_Dominguez2011_file(z, 'ebl_saldana2021_fiducial.txt')
    dump_ebl_file(eps, z, density, 'ebl_Saldana2021_fiducial.txt')
    
    # Saldana2021 lower
    eps, density = load_Saldana2021_error_file(z, True)
    dump_ebl_file(eps, z, density, 'ebl_Saldana2021_lower.txt')

    # Saldana2021 upper
    eps, density = load_Saldana2021_error_file(z, False)
    dump_ebl_file(eps, z, density, 'ebl_Saldana2021_upper.txt')

    # Gilmore2012 fiducial
    z = np.array([0., 0.015, 0.025, 0.044, 0.05, 0.2, 0.4, 0.5, 0.6, 0.8, 1.0,
                  1.25, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0])
    eps, density = load_Gilmore2012_file(z, 'ebl_gilmore2012_fiducial.txt')
    dump_ebl_file(eps, z, density, 'ebl_Gilmore2012_fiducial.txt')

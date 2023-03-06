import numpy as np

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
        
    E = np.logspace(-1, 3, 100) # MeV
    
    size = len(A)
    
    f = open('tables/' + filename, 'w')
    f.write("#A - Z - sigma [mb]\n")
    
    for i in range(size):
        sigma = compute_sigma(E, t[i], h1[i], x1[i], w1[i], c[i])
        f.write('{},{},'.format(int(A[i]),int(Z[i])))
        for s in sigma:
            f.write('{:.5e},'.format(s))
        f.write('\n')
    
    f.close()
    
if __name__== "__main__":
    dump_v2r4('xsecs_photodisintegration_v2r4_alpha.txt', True)
    dump_v2r4('xsecs_photodisintegration_v2r4_singlenucleon.txt', False)

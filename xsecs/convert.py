import numpy as np
from scipy import interpolate

#def load_losses(filename):
#    file = open(filename, 'r')
#    lines = file.readlines()
#
#    xstr = lines[0].split(',')
#    #ystr = list (map (lambda x: float(x),xstr))
#    x = [float(i) for i in xstr]
#
#    ystr = lines[1].split(',')
#    y = [float(i) for i in ystr]
#
#    return np.array(x), np.array(y)
#
#def convert_losses():
#    filename = 'original/pairlossesInSimprop.txt'
#    x, y = load_losses(filename)
#
#    expected_length = 501
#    file = open('losses_pair_BGG2002.txt', 'w')
#    file.write('# size in energy %d\n' % expected_length)
#
#    counter = 0
#    for i in range(0, len(x), 2):
#        counter += 1
#        file.write("%10.5e, %10.5e\n" % (x[i], y[i]))
#
#    assert(counter == expected_length)
#
#    file.close()
#
#def interpolate_ppp_sigma(s):
#    s_low = [1.160,1.165,1.195,1.250,1.305,1.375,1.455,1.470,1.485,1.500,1.625,1.680,1.750,1.795,1.820,1.875,1.950,1.980,2.000,2.050,2.100,2.150,2.200,2.250,2.300,2.350,2.400,2.450,2.500,2.625,2.750,2.875,3.000,3.250,3.500,4.000,4.500,5.000]
#    sigma_low = [.0000,.0003,.0502,.1279,.1952,.3173,.4970,.5229,.5414,.5317,.2799,.2233,.1851,.1719,.1676,.1811,.2014,.2115,.2189,.2253,.2343,.2491,.2690,.2874,.2915,.2752,.2490,.2257,.2112,.2051,.2166,.2109,.1873,.1620,.1564,.1493,.1395,.1359]
#
#    logsqrts_high = [.3495,.5044,1.,2.,3.,4.,4.548]
#    logsigma_high = [-.867,-.915,-.941,-.851,-.684,-.514,-.423]
#
#    if (s < s_low[0]):
#        return 0.
#    elif (s < s_low[-1]):
#        tck = interpolate.splrep(s_low, sigma_low, s=0)
#        return interpolate.splev(s, tck, der=0)
#    else:
#        s = min(s, 1e9)
#        tck = interpolate.splrep(logsqrts_high, logsigma_high, s=0)
#        sqrts = np.sqrt(s)
#        return np.power(10., interpolate.splev(np.log10(sqrts), tck, der=0))

#def compute_ppp_sigma():
#    expected_length = 9999
#    file = open('xsec_ppp.txt', 'w')
#    file.write('# s [GeV^2] - sigma [mbarn] - phi [GeV^4 mbarn]\n')
#
#    s = np.logspace(0, 10, expected_length)
#    phi_s = 0.
#    lnf = np.log(s[1] / s[0])
#    protonmass = 0.93827208816 # GeV
#    for s_i in s:
#        sigma_i = interpolate_ppp_sigma(s_i)
#        phi_s += lnf * s_i * (s_i - protonmass * protonmass) * sigma_i
#        file.write("%10.6e, %10.6e, %10.6e\n" % (s_i, sigma_i, phi_s))
#
#
#def interpolate_ppp_sigma_proton(s_sophia, sigma_sophia, s):
#    if (s < s_sophia[0]):
#        return 0.
#    if (s > s_sophia[-1]):
#        return sigma_sophia[-1]
#    tck = interpolate.splrep(s_sophia, sigma_sophia, s=0)
#    return interpolate.splev(s, tck, der=0)
#    file = open(outputfile, 'w')
#    file.write('# s [GeV^2] - sigma [mbarn] - phi [GeV^4 mbarn]\n')
#
#    file.write("#\n")
#    for i in range(s.size):
#        file.write("%10.6e, %10.6e, %10.6e\n" % (s[i], sigma[i], phi[i]))
#    file.close()


#    compute_ppp_sigma_sophia('xs_proton.txt', 'xsec_ppp_sophia_proton.txt')
#    compute_ppp_sigma_sophia('xs_neutron.txt', 'xsec_ppp_sophia_neutron.txt')



def eps2s(eps):
    protonmass = 0.93827208816 # GeV
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

def dump_xsecs_file(s, sigma, filename):
    assert(s.size == sigma.size)
    file = open(filename, "w")
    file.write(f'# size in energy {s.size}\n')
    file.write(f'# s [GeV^2] - sigma [mbarn] - phi(s) [GeV^4 mbarn]\n')
    for i in range(s.size):
        phi = 0.
        file.write(f'{s[i]:10.5e}, {sigma[i]:10.5e}, {phi:10.5e}\n')
    file.close()
    print ('dumped xsecs table on ', filename)
    
if __name__== "__main__":
    s, sigma = load_photopion_sophia('tables/xs_proton.txt')
    dump_xsecs_file(s, sigma, 'xsecs_photopion_proton_sophia.txt')

    s, sigma = load_photopion_sophia('tables/xs_neutron.txt')
    dump_xsecs_file(s, sigma, 'xsecs_photopion_neutron_sophia.txt')

    s, sigma = load_photopion_astrophomes('tables/xsecs_photopion_inelastic_astrophomes.txt', doProton=True)
    dump_xsecs_file(s, sigma, 'xsecs_photopion_proton_astrophomes.txt')

    s, sigma = load_photopion_astrophomes('tables/xsecs_photopion_inelastic_astrophomes.txt', doProton=False)
    dump_xsecs_file(s, sigma, 'xsecs_photopion_neutron_astrophomes.txt')

    s, sigma = load_photopion_v2r4()
    dump_xsecs_file(s, sigma, 'xsecs_photopion_proton_v2r4.txt')

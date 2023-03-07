import numpy as np

def load_losses(filename):
    file = open(filename, 'r')
    lines = file.readlines()
    
    xstr = lines[0].split(',')
    x = [float(i) for i in xstr]

    ystr = lines[1].split(',')
    y = [float(i) for i in ystr]

    return np.array(x), np.array(y)
    
def dump_losses(outname):
    filename = 'tables/pairlossesInSimprop.txt'
    x, y = load_losses(filename)
    
    expected_length = 1001
    file = open(outname, 'w')
    file.write('# size in energy %d\n' % expected_length)
    file.write('# log10(E/eV) - log10(b/year^-1) \n')
    
    counter = 0
    for i in range(0, len(x)):
        counter += 1
        file.write("%10.5e, %10.5e\n" % (x[i], y[i]))
    
    assert(counter == expected_length)
    file.close()
    print('dumped losses table on ', outname)

if __name__== "__main__":
    dump_losses('losses_pair_BGG2002.txt')

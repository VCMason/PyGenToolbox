def mean(lst):
    return sum(lst) / len(lst)


def plot_histogram(x, lowerlimit=1000):
    import matplotlib.pyplot as plt
    print('Mean of fragment size: %f' % (mean(x)))
    # x = [value1, value2, value3, ....]
    plt.hist(x, bins=100)
    plt.show()
    plt.close()

    x = [i for i in x if i >= lowerlimit]
    print('Mean of mean fragment size > %d: %f' % (lowerlimit, mean(x)))
    plt.hist(x, bins=100)
    plt.show()
    plt.close()

    x = [i for i in x if i >= lowerlimit*2]
    print('Mean of mean fragment size > %d: %f' % (lowerlimit*2, mean(x)))
    plt.hist(x, bins=100)
    plt.show()
    plt.close()


def calculate_distribution_of_re_cut_fragments(pattern, seqs, efficiency=0.9):
    distances = []
    gate = 0
    import re
    import random
    p = re.compile(pattern)
    print('Cuttting efficiency = %f' % efficiency)
    for k, v in seqs.items():
        gate = 0
        for m in p.finditer(v):
            # print(m.start(), m.group()) # the starting position of each match, and then the string (.group()) matched
            if gate == 0:
                current_pos = m.start()
                distances.append(current_pos - 0)
                previous_pos = m.start()
                gate = 1
            elif (gate == 1) and (random.randint(1, 10) <= 10 * efficiency):
                current_pos = m.start()
                distances.append(current_pos - previous_pos)
                previous_pos = m.start()
    print(distances[:100])
    return distances


def CalcExpectedNumberOfRESites(l, p):
    k = len(p) # length of kmer (pattern)
    NumExp = l / 4 ** k  # 4 because four possible nucleotides
    print('Expected number of RE recognition sites for sequence of length %d and pattern of length %d = %d' % (l, k, NumExp))


def CountPatternMatchesInSeqs(seqs, pattern):
    ''' seqs is a dictionary, key is line with carrot '>', value is sequence concatenated to one line (no CRLF or LF) '''
    t = 0 # total number of RE recognition sites
    l = 0  # total length of all fasta sequences summed together
    lpat = len(pattern)
    for k, v in seqs.items():
        t += v.count(pattern)
        l += len(v)
    return(t, l)


def ParseFasta(f):
    print(f)
    seqs = {}
    FILE = open(f, 'r')
    l = FILE.readline()
    while l:
        if l[0] == '>':
            n = l[1:].strip()
            seqs[n] = ''
        else:
            seqs[n] += l.strip().upper()
        l = FILE.readline()
    FILE.close()
    return(seqs)


def main(filepath='', repattern=''):
    #import pkg_resources
    #DATA_PATH = pkg_resources.resource_filename('pygentoolbox', 'data/seqs/ptetraurelia_mac_51.fa')

    if filepath == '':
        filepath = input('Enter full path to file: ')
    if repattern == '':
        repattern = input('Enter restriction enzyme pattern, ex. GATC: ')
    #f = 'D:\\LinuxShare\\Ciliates\\Genomes\\Seqs\\ptetraurelia_mac_51.fa'
    #f = 'D:\\LinuxShare\\Ciliates\\Genomes\\Seqs\\ptetraurelia_mac_51_with_ies.fa'
    #f = 'D:\\LinuxShare\\Ciliates\\Genomes\\Seqs\\ptetraurelia_mic2.fa'
    #f = 'D:\\LinuxShare\\Ciliates\\Genomes\\Seqs\\contigs_ABK_COSP_best_k51_no_scaf.fa'
    
    seqs = ParseFasta(filepath)
    print('Number of sequences: %d' % (len(seqs.keys())))
    total, LenFastaSeqs = CountPatternMatchesInSeqs(seqs, repattern)
    print('Number of RE sites matching %s: %d' % (repattern, total))
    CalcExpectedNumberOfRESites(LenFastaSeqs, repattern)
    distances = calculate_distribution_of_re_cut_fragments(repattern, seqs, efficiency=0.9)
    plot_histogram(distances)

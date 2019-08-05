
def CalcExpectedNumberOfRESites(l, p):
    k = len(p) # length of kmer (pattern)
    NumExp = l / 4 ** k
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
    
def main_CountRESites():
    #import pkg_resources
    #DATA_PATH = pkg_resources.resource_filename('pygentoolbox', 'data/seqs/ptetraurelia_mac_51.fa')

    f = input('Enter full path to file: ')
    pattern = input('Enter restriction enzyme pattern, ex. GATC: ')
    #f = 'D:\\LinuxShare\\Ciliates\\Genomes\\Seqs\\ptetraurelia_mac_51.fa'
    #f = 'D:\\LinuxShare\\Ciliates\\Genomes\\Seqs\\ptetraurelia_mac_51_with_ies.fa'
    #f = 'D:\\LinuxShare\\Ciliates\\Genomes\\Seqs\\ptetraurelia_mic2.fa'
    #f = 'D:\\LinuxShare\\Ciliates\\Genomes\\Seqs\\contigs_ABK_COSP_best_k51_no_scaf.fa'
    
    seqs = ParseFasta(f)
    print('Number of sequences: %d' % (len(seqs.keys())))
    total, LenFastaSeqs = CountPatternMatchesInSeqs(seqs, pattern)
    print('Number of RE sites matching %s: %d' % (pattern, total))
    CalcExpectedNumberOfRESites(LenFastaSeqs, pattern)
    
main_CountRESites()
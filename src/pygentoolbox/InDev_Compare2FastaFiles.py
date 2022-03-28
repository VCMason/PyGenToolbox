
def CompareSeqLength(s1, s2):


def ParseFasta(f):
    print (f)
    seqs = {}
    FILE = open(f, 'r')
    l = FILE.readline()
    while l:
        if l[0] == '>':
            n = l[1:].strip()
            seqs[n] = ''
        else:
            seqs[n] += l.strip()
        l = FILE.readline()
    FILE.close()
    return(seqs)

def Compare2FastaFiles():
    f1 = input('Enter first fasta file: ')
    f2 = input('Enter second fasta file: ')

    s1 = ParseFasta(f1)
    s2 = ParseFasta(f2)

    CompareSeqLength()

Compare2FastaFiles()
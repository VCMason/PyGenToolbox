
def query_len(cigar_string):
    """
    Given a CIGAR string, return the number of bases consumed from the
    query sequence.
    """
    from itertools import groupby

    read_consuming_ops = ("M", "I", "S", "=", "X")
    seqlength = 0
    cig_iter = groupby(cigar_string, lambda chr: chr.isdigit())
    for _, length_digits in cig_iter:
        length = int(''.join(length_digits))
        op = next(next(cig_iter)[1])
        if op in read_consuming_ops:
            seqlength += length

    return seqlength


#####  Calculations  #####

def sturges_rule(n):
    import math
    # n is number of observations
    numbins = round(1 + 3.322 * math.log(n, 10))

    print('Number if bins claculated by Sturge\'s Rule: %d' % binnumber)

    return numbins


def pearson_correlation(x, y):
    from scipy.stats.stats import pearsonr
    print('Pearson Correlation:')
    print(pearsonr(x, y))


def ttest_two_sided_independent(group1, group2):
    from scipy.stats import ttest_ind

    ttest_res = ttest_ind(group1, group2)  # groups are lists  # ttest_res is a tuple with two entries
    #ttest_res = (The calculated t-statistic, and The two-tailed p-value)

    return(ttest_res)

def tpm(df):
    """ calculate TPM """
    """ 
    #1. Divide the read counts by the length of each gene in kilobases. This gives you reads per kilobase (RPK).
    #2. Count up all the RPK values in a sample and divide this number by 1,000,000. This is your “per million” scaling factor.
    #3. Divide the RPK values by the “per million” scaling factor. This gives you TPM.
    """

    df['RPK'] = df.iloc[:, 9] / (df.iloc[:, 4] - df.iloc[:, 3] + 1)
    # print(df['RPK'].dtype)
    # print(df['RPK'].size)
    # print(df['RPK'].sum())
    # print('If the above line says \'inf\' there is a problem')
    # for val in df['RPK'].iteritems():
    #     if val[1] > 1:
    #         print('BIG')
    #         print(val)
    #     elif val[1] < 0.0:
    #         print('small')
    #         print(val)
    df['ScaleFactor'] = df['RPK'].sum() / 1000000.0
    df['TPM'] = df['RPK'] / df['ScaleFactor']

    print('Sum of all TPM values = %f' % df['TPM'].sum())

    return df


def median_of_ratios(df, verbose):
    """ Normalization method that imitates DESeq2 normalization. Give pandas DF of raw counts. """
    """ assumes unique index present to identify, gene/scaffold/position """
    """ following: https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html """
    import math
    import statistics

    # for each gene calculate the geometric mean sqrt(a*b*c*d...)
    # EXAMPLE OF APPLY: product = df.apply(lambda x: x**2) how to square all values
    # multiply all values in each row by one another, take sqrt of each, save Geometric mean as new column
    df['GeoMean'] = df.prod(axis=1).apply(math.sqrt)

    # divide each sample's gene count by the geometric mean for each gene (step 2)
    dfratios = df.iloc[:, :-1].div(df.GeoMean, axis=0)

    # Normalization factor: calculate the median of all ratios (from step 2). One median value per sample.
    NFs = dfratios.apply(statistics.median, axis=0)  # axis 0 is default.. # one median value per column (sample)

    # For each sample, divide the all raw gene counts by the one normalization factor for that sample
    dfNorm = df.iloc[:, :-1].div(NFs, axis=1)
    # dfNorm.insert(loc=0, column='GeneNames', value=df.iloc[:, 1])  # add names (first column) to normalized output df

    if verbose == True:
        print('Geometric Means:')
        print(df.GeoMean.head())
        print('Ratios:')
        print(dfratios.head())
        print('Normalization Factors:')
        print(NFs)

    return dfNorm


def l2fc(controls, cases):
    import math
    res = math.log2(mean(cases) / mean(controls))
    return(res)


def log2_fold_change(df, a, b, verbose=False):
    """ Accepts pandas df, a is index of numerator, b is denominator """
    """ Returns pandas series of log2(a/b) fold changes """
    import math

    #divide columns
    div = df.iloc[:, a] / df.iloc[:, b]
    log2FCseries = div.apply(math.log2)

    if verbose == 'True':
        print(len(div))
        print(div.head())
        print(len(div > 0))

    return(log2FCseries)


def mean(lst):
    return sum(lst) / len(lst)


def all_combinations_of_size(l, n):
    """l is list, n is number. 2 gives all combinations of size 2 of values from list"""
    import itertools

    return(list(itertools.combinations(l, n)))


#####  OPERATIONS  #####


def run_gzip_compress(fullpath):
    import subprocess

    # pack files (using full path)
    print('Start packing files:\n%s' % fullpath)
    cmd = "gzip -f -k %s" % fullpath
    print(cmd)
    cmdlist = cmd.split()
    subprocess.call(cmdlist)
    filegz = fullpath + '.gz'
    print('Finished packing:\n%s\n' % filegz)

    return filegz


def run_gzip_decompress(fullpath):
    import subprocess
    # unpack fastp trimmed file (using full path)
    print('Start unpacking trimmed.gz file:\n%s' % fullpath)
    cmd = "gzip -f -d -k %s" % fullpath
    print(cmd)
    cmdlist = cmd.split()
    subprocess.call(cmdlist)
    filetrim = fullpath[:-len('.gz')]
    print('Finished unpacking:\n%s\n' % filetrim)

    return filetrim


def make_directory(dirName):
    # dirName is directory in cwd or full path to directory
    import os

    if not os.path.exists(dirName):
        os.mkdir(dirName)
        print("Directory ", dirName, " Created ")
    else:
        print("Directory ", dirName, " already exists")


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
    k = len(p)  # length of kmer (pattern)
    NumExp = l / 4 ** k  # 4 because four possible nucleotides
    print('Expected number of RE recognition sites for sequence of length %d and pattern of length %d = %d' % (
    l, k, NumExp))


def CountPatternMatchesInSeqs(seqs, pattern):
    ''' seqs is a dictionary, key is line with carrot '>', value is sequence concatenated to one line (no CRLF or LF) '''
    t = 0  # total number of RE recognition sites
    l = 0  # total length of all fasta sequences summed together
    lpat = len(pattern)
    for k, v in seqs.items():
        t += v.count(pattern)
        l += len(v)
    return (t, l)


def subsample_dictionary(d):
    import random
    subsampled_dictionary = {k: v for (k, v) in random.sample(d.items(), k=numreads)}
    return subsampled_dictionary


def translate_AA_from_dna_human_terminate_at_stop(seq):
    #allpossiblecodons = {'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A', 'AAC': 'B', 'AAT': 'B', 'GAC': 'B', 'GAT': 'B',
    #           'TGC': 'C', 'TGT': 'C', 'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E', 'TTC': 'F', 'TTT': 'F',
    #           'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G', 'CAC': 'H', 'CAT': 'H', 'ATA': 'I', 'ATC': 'I',
    #           'ATT': 'I', 'AAA': 'K', 'AAG': 'K', 'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L', 'TTA': 'L',
    #           'TTG': 'L', 'ATG': 'M', 'AAC': 'N', 'AAT': 'N', 'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    #           'CAA': 'Q', 'CAG': 'Q', 'AGA': 'R', 'AGG': 'R', 'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    #           'AGC': 'S', 'AGT': 'S', 'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S', 'ACA': 'T', 'ACC': 'T',
    #           'ACG': 'T', 'ACT': 'T', 'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V', 'TGG': 'W', 'NNN': 'X',
    #           'TAC': 'Y', 'TAT': 'Y', 'CAA': 'Z', 'CAG': 'Z', 'GAA': 'Z', 'GAG': 'Z', 'TAA': '*', 'TAG': '*',
    #           'TGA': '*'}
    noambigcodons = {'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
               'TGC': 'C', 'TGT': 'C', 'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E', 'TTC': 'F', 'TTT': 'F',
               'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G', 'CAC': 'H', 'CAT': 'H', 'ATA': 'I', 'ATC': 'I',
               'ATT': 'I', 'AAA': 'K', 'AAG': 'K', 'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L', 'TTA': 'L',
               'TTG': 'L', 'ATG': 'M', 'AAC': 'N', 'AAT': 'N', 'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
               'CAA': 'Q', 'CAG': 'Q', 'AGA': 'R', 'AGG': 'R', 'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
               'AGC': 'S', 'AGT': 'S', 'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S', 'ACA': 'T', 'ACC': 'T',
               'ACG': 'T', 'ACT': 'T', 'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V', 'TGG': 'W', 'NNN': 'X',
               'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
               'TGA': '*'}

    AAseq = ''
    count = 0
    while count < len(seq):
        sub = seq[count:count + 3]
        if noambigcodons[sub] == '*':
            break
        else:
            AAseq += noambigcodons[sub]
        count += 3

    return AAseq


def translate_AA_from_dna_human(seq):
    allpossiblecodons = {'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A', 'AAC': 'B', 'AAT': 'B', 'GAC': 'B', 'GAT': 'B',
               'TGC': 'C', 'TGT': 'C', 'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E', 'TTC': 'F', 'TTT': 'F',
               'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G', 'CAC': 'H', 'CAT': 'H', 'ATA': 'I', 'ATC': 'I',
               'ATT': 'I', 'AAA': 'K', 'AAG': 'K', 'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L', 'TTA': 'L',
               'TTG': 'L', 'ATG': 'M', 'AAC': 'N', 'AAT': 'N', 'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
               'CAA': 'Q', 'CAG': 'Q', 'AGA': 'R', 'AGG': 'R', 'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
               'AGC': 'S', 'AGT': 'S', 'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S', 'ACA': 'T', 'ACC': 'T',
               'ACG': 'T', 'ACT': 'T', 'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V', 'TGG': 'W', 'NNN': 'X',
               'TAC': 'Y', 'TAT': 'Y', 'CAA': 'Z', 'CAG': 'Z', 'GAA': 'Z', 'GAG': 'Z', 'TAA': '*', 'TAG': '*',
               'TGA': '*'}
    noambigcodons = {'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
               'TGC': 'C', 'TGT': 'C', 'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E', 'TTC': 'F', 'TTT': 'F',
               'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G', 'CAC': 'H', 'CAT': 'H', 'ATA': 'I', 'ATC': 'I',
               'ATT': 'I', 'AAA': 'K', 'AAG': 'K', 'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L', 'TTA': 'L',
               'TTG': 'L', 'ATG': 'M', 'AAC': 'N', 'AAT': 'N', 'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
               'CAA': 'Q', 'CAG': 'Q', 'AGA': 'R', 'AGG': 'R', 'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
               'AGC': 'S', 'AGT': 'S', 'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S', 'ACA': 'T', 'ACC': 'T',
               'ACG': 'T', 'ACT': 'T', 'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V', 'TGG': 'W', 'NNN': 'X',
               'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
               'TGA': '*'}

    AAseq = ''
    count = 0
    while count < len(seq):
        sub = seq[count:count + 3]
        AAseq += noambigcodons[sub]
        count += 3

    return AAseq


def reverse_complement_to_dna(seq):
    revcomp = ''
    rev = seq[::-1]
    rev = rev.upper()
    for nucl in rev:
        if nucl == 'A':
            revcomp += 'T'
        #if nucl == 'A':
            #revcomp += 'U'
        elif nucl == 'T':
            revcomp += 'A'
        #elif nucl == 'U':
            #revcomp += 'A'
        elif nucl == 'G':
            revcomp += 'C'
        elif nucl == 'C':
            revcomp += 'G'
        elif nucl == 'N':
            revcomp += 'N'
    return revcomp, rev
    # revcomp, revline = ReverseComplement(line)


def reverse_complement_to_rna(seq):
    revcomp = ''
    rev = seq[::-1]
    rev = rev.upper()
    for nucl in rev:
        #if nucl == 'A':
            #revcomp += 'T'
        if nucl == 'A':
            revcomp += 'U'
        #elif nucl == 'T':
            #revcomp += 'A'
        elif nucl == 'U':
            revcomp += 'A'
        elif nucl == 'G':
            revcomp += 'C'
        elif nucl == 'C':
            revcomp += 'G'
        elif nucl == 'N':
            revcomp += 'N'
    return revcomp, rev
    # revcomp, revline = ReverseComplement(line)


def Merge(dict1, dict2):
    ''' Python code to merge dict using update() method '''
    ''' update dict1 with values/keys of dict2 '''
    #return dict2.update(dict1)
    return {**dict1, **dict2}


def window(seq, n=2):
    """ Returns a sliding window (of width n) over data from the iterable """
    """   s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...                   """
    from itertools import islice

    it = iter(seq)
    result = tuple(islice(it, n))
    if len(result) == n:
        yield result
    for elem in it:
        result = result[1:] + (elem,)
        yield result


def window_2(iterable, size=2):
    i = iter(iterable)
    win = []
    for e in range(0, size):
        win.append(next(i))
    yield win
    for e in i:
        win = win[1:] + [e]
        yield win


def extract_seqs_atcg_above_length(seq):
    ''' collect_all_continuous_strings_of_ATCG '''
    import re
    matchObjs = re.finditer('[ATCG]+', seq, re.I) # record all stretches of ATCG regardless of case
    #for m in matchObjs:
    #    print(m.group())
    #    print(m.span()[0])
    #    print(m.span()[1])
    return matchObjs


def extract_spans_atcg(seq):
    ''' collect_all_continuous_strings_of_ATCG '''
    import re
    spans = {} # record spans of all indels from all taxa
    matchObjs = re.finditer('[ATCG]+', seq, re.I) # record all stretches of ATCG regardless of case
    for m in matchObjs:
        s = '%d_%d' % (m.span()[0], m.span()[1])
        spans[s] = spans.get(s, []) + [n]
    return spans


def get_spans_of_char(char, d):
    ''' collect_all_continuous_strings_of_char_from_dict_Span '''
    import re
    spans = {}
    for n in d.keys(): # record spans of all indels from all taxa
        matchObjs = [m for m in re.finditer(r'%s+' % (char), d[n])] # record all indel events
        for m in matchObjs:
            s = '%d_%d' % (m.span()[0], m.span()[1])
            spans[s] = spans.get(s, []) + [n]
    return spans


def clone_features_up_down():
    ''' Will give improper coordinates if cloned features extend before or after start of contig '''
    f = input('Full path to .gff3 file: ')
    # D:\LinuxShare\Ciliates\Genomes\Annotations\internal_eliminated_sequence_PGM_ParTIES.pt_51_with_ies.gff3
    e = input('Designate cloned feature base name (ex: IES_clone: ')

    with open(f, 'r') as FILE:
        header = []
        out = []
        countnegstart = 0
        line = FILE.readline()
        while line:
            if line[0] == '#':
                header.append(line.strip())
            else:  # clone each feature (one upstream, one downstream). Name cloned features with e + 'Up' or + 'Down'.
                x = line.strip().split('\t')
                start, end, featurelength = int(x[3]), int(x[4]), int(x[4]) - int(x[3])
                up, down = x[:], x[:]  # [:] clones the list so i can change the three independently
                up[2], up[3], up[4] = e + '_Up', str(start - featurelength - 1), str(start - 1)  # subtract one to not overlap feature
                down[2], down[3], down[4] = e + '_Down', str(end + 1), str(end + featurelength + 1)
                if int(up[3]) < 0:  # if cloned start coordinate is < 0, then set it to 1
                    up[3] = str(1)
                    out.append(up)
                    out.append(x)
                    out.append(down)
                    countnegstart += 1
                else:  # am not controlling for coordinates > contig length if int(down[4]) > len(contig)
                    out.append(up)
                    out.append(x)
                    out.append(down)

            line = FILE.readline()

    outfile = '.'.join(f.split('.')[:-1] + ['clone'] + [f.split('.')[-1]])
    with open(outfile, 'w') as OUT:
        OUT.write('\n'.join(header) + '\n')
    with open(outfile, 'a') as OUT:
        OUT.write('\n'.join(['\t'.join(line) for line in out]))
    print('Finished cloning features to output file: %s' % (outfile))
    print('Number of cloned features with negative start coordinates adjusted to 1 = %d' % (countnegstart))


#####  FORMATING  #####


def make_circos_karyotype_file(d, outpath, sn):
    ''' input dictionary (key = chromosome name, value = chrom sequence (non-interleaved)), and full output path'''
    # usually would run read_interleaved_fasta_as_noninterleaved(filename) before
    # d is dictionary of all "chromosomes" for karyotype file
    # sn is species name Genus species
    print('begin formatting for karyotype file')
    print('Species: %s' % sn)
    from natsort import natsorted, ns
    # orderedkeys = natsorted(d.keys(), alg=ns.IGNORECASE)  # using natural sort
    karylist = [' '.join(['chr', '-', sn[0].lower() + sn.split(' ')[1][0] + k.split('_')[1], k.split('_')[1], '0', str(len(d[k])), 'chr' + k.split('_')[1]]) for k in natsorted(d.keys(), alg=ns.IGNORECASE)]  # chr - pt_1 1 0 len(scaff) chr1
    with open(outpath, 'w') as OUT:
        OUT.write('\n'.join(karylist))
    print('output karyotype file to: %s' % outpath)
    return


def ParseGff3SkipComments(f):
    print(f)
    seqs = {}
    with open(f, 'r') as FILE:
        for l in FILE:
            if l[0] == '#':
                pass
            else:
                n = '_'.join([l.strip().split('\t')[0]] + l.strip().split('\t')[3:5])  # n == scaffold_start_end
                seqs[n] += l.strip().split('\t')
    return (seqs)


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


def parse_bed_to_dict(bedfile):
    ''' read bedfile in as dictionary. Not Memory efficient '''
    d = {}
    with open(bedfile) as BED:
        for line in BED:
            d[line.strip().split('\t')[0]] = d.get(line.strip().split('\t')[0], []) + [line.strip()]
    return d


def write_dict_to_fasta(names, d, outfile):
    output = []
    for n in names:
        output.append('%s\n%s' % (n, d[n]))
    with open(outfile, 'w') as OUT:
        OUT.write('\n'.join(output))

def read_interleaved_fasta_as_noninterleaved_human_proteins(filename):

    ''' filename (full path) '''
    count = 0
    print('Counting total lines')
    with open(filename, 'r') as GENOME:
        for line in GENOME:
            count += 1
    print('Number of lines in genome file: %d' % count)

    print('Reading in interleaved fasta as noninterleaved fasta dictionary')
    with open(filename, 'r') as GENOME:
        names = []  # record names to maintain order
        d = {}  # dictionary of sequences, key is >line.strip() and in names, value is noninterleaved sequence
        count = 0
        for line in GENOME:
            if '>' != line[0]:
                d[names[-1]] += [line.strip()]
            elif '>' == line[0]:
                if len(names) > 0:
                    d[names[-1]] = ''.join(d[names[-1]])  # join the list of all lines (optional)
                d[line.strip().split('|')[1]] = []
                names.append(line.strip().split('|')[1])
            else:
                print('Problem!!!')
            count += 1
            if count % 1000000 == 0:
                print('Current line: %d' % count)
    print('Finished reading file')
    return names, d


def read_interleaved_fasta_as_noninterleaved(filename):

    ''' filename (full path) '''
    count = 0
    print('Counting total lines')
    with open(filename, 'r') as GENOME:
        for line in GENOME:
            count += 1
    print('Number of lines in genome file: %d' % count)

    print('Reading in interleaved fasta as noninterleaved fasta dictionary')
    with open(filename, 'r') as GENOME:
        names = []  # record names to maintain order
        d = {}  # dictionary of sequences, key is >line.strip() and in names, value is noninterleaved sequence
        count = 0
        for line in GENOME:
            if '>' != line[0]:
                d[names[-1]] += [line.strip()]
            elif '>' == line[0]:
                if len(names) > 0:
                    d[names[-1]] = ''.join(d[names[-1]])  # join the list of all lines (optional)
                d[line.strip()[1:].split()[0]] = []
                names.append(line.strip()[1:].split()[0])
            else:
                print('Problem!!!')
            count += 1
            if count % 5000000 == 0:
                print('Current line: %d' % count)
    print('Finished reading in interleaved fasta as non-interleaved')

    return names, d


def interlevaed_fasta_to_noninterleaved_output(filename):
    ''' filename (full path) '''

    HANDLE = open(filename, 'r')
    with open(filename, 'w') as OUT:  # clear contents of file before appending
        OUT.write('')
    OUT = open(filename + '.reformat', 'a')

    gate = 'closed'
    line = HANDLE.readline()
    while line:
        if '>' not in line:
            OUT.write(line.strip())
        elif '>' in line and gate == 'closed':
            OUT.write(line.strip() + '\n')
            gate = 'open'
        elif '>' in line and gate == 'open':
            OUT.write('\n' + line.strip() + '\n')
        else:
            print('Problem!!!')
        line = HANDLE.readline()

    HANDLE.close()
    OUT.close()


def intersect_dfs_by_common_indices(df1, df2):
    """ provide pandas df """
    import pandas as pd
    inter = pd.Index.intersection(df1.index, df2.index)  # inter is an index object for rows indices common to both dfs

    # To find duplicated indices:
    # print(df1[df1.index.duplicated()])
    # print(df2[df2.index.duplicated()])
    # To remove rows with duplicated indices, use: df = df[~df.index.duplicated()]

    df1trim, df2trim = df1.filter(inter, axis=0), df2.filter(inter, axis=0)  # df.filter(axis=0) filters with index object on rows
    return df1trim, df2trim


def extract_columns_from_files(flist, delim, n):  # n is column index, flist is list of all files in order
    out = []
    for f in flist:
        col = []
        with open(f, 'r') as FILE:
            line = FILE.readline()
            col.append(line.strip().split(delim)[n])
        out.append(col)
    return(out)  # out is a 2d list of lists (inner lists each represent column (of index n) from each file in flist)


def pandas_output_to_matrix(p, f, df, delim):
    import os

    fullpath = os.path.join(p, f)
    print('Output File: %s' % (fullpath))
    df.to_csv(fullpath, sep=delim)


def fastq_to_fasta(fullpath):  # fullpath = string, full path to file
    """ Converts a .fastq file to fasta file """
    print('converting fastq to fasta')
    import os

    path, f = os.path.split(fullpath)
    fout = os.path.join(path, '.'.join(f.split('.')[:-1]) + '.fa')
    if os.path.exists(fout):  # if the output file exists delete it before appending
        os.remove(fout)

    FILE = open(fullpath, 'r')
    OUT = open(fout, 'a')

    line = FILE.readline()
    count = 0
    while line:
        if count % 4 == 0:  # if the entire line is only a +
            o = '>' + line[1::].strip() + '\n'
            OUT.write(o)
        elif count % 4 == 1:
            o = line.strip() + '\n'
            OUT.write(o)
        line = FILE.readline()
        count += 1

    FILE.close()
    OUT.close()

    print('Converted fastq to fasta:\n%s\n%s\n' % (fullpath, fout))

    return fout


def read_fastq(fullpath):  # fullpath = string, full path to file
    """ Converts a .fastq file to fasta file """
    import os

    path, f = os.path.split(fullpath)
    fout = os.path.join(path, '.'.join(f.split('.')[:-1]) + '.fa')
    if os.path.exists(fout):  # if the output file exists delete it before appending
        os.remove(fout)

    FILE = open(fullpath, 'r')
    OUT = open(fout, 'a')

    line = FILE.readline()
    count = 0
    d = {}
    while line:
        if count % 4 == 0:  # if the entire line is only a +

            d[line[1:].strip()] = ''
        elif count % 4 == 1:
            o = line.strip() + '\n'
            OUT.write(o)
        line = FILE.readline()
        count += 1

    FILE.close()
    OUT.close()

    print('Converted fastq to fasta:\n%s\n%s' % (fullpath, fout))


def example_pandas_df():
    import pandas as pd

    # Create a DataFrame
    d = {
        'Name': ['Alisa', 'Bobby', 'Cathrine', 'Madonna', 'Rocky', 'Sebastian', 'Jaqluine',
                 'Rahul', 'David', 'Andrew', 'Ajay', 'Teresa'],
        'Score1': [62, 47, 55, 74, 31, 77, 85, 63, 42, 32, 71, 57],
        'Score2': [89, 87, 67, 55, 47, 72, 76, 79, 44, 92, 99, 69],
        'Score3': [56, 86, 77, 45, 73, 62, 74, 89, 71, 67, 97, 68]}

    df = pd.DataFrame(d)
    return(df)

#########################
#####  Plotting     #####
#########################


def bar_plot_two_sided(labels, up, low, outpath, figwidth=10, figheight=5):  # pass in all lists
    import matplotlib.pyplot as plt
    import numpy as np

    up = np.array(up)
    low = np.array(low) * -1

    fig, ax = plt.subplots(figsize=(figwidth, figheight))
    ax.bar(labels, up)
    ax.bar(labels, low)

    # Formatting x labels
    plt.xticks(rotation=90)
    plt.tight_layout()
    # Use absolute value for y-ticks
    ticks = ax.get_yticks()
    ax.set_yticklabels([int(abs(tick)) for tick in ticks])
    ax.legend(['Sense', 'Antisense'], fontsize=14)

    plt.savefig(outpath)
    plt.show()
    plt.close()


#####  Simple Commands  #####


def longest_dict_value():
    longest = max(map(len, d.values()))  # length of longest dict value
    return longest


def keep_rows_df_greater_than(df, cutoff, n):  # dataframe, cutoff value, column index to filter
    dftrim = df[df.iloc[:, n] > cutoff]

    return dftrim
    # df[df['column name'].map(len) < 2]
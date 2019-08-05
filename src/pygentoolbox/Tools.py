
#####  Calculations  #####


def pearson_correlation(x, y):
    from scipy.stats.stats import pearsonr
    print('Pearson Correlation:')
    print(pearsonr(x, y))


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


def Merge(dict1, dict2):
    ''' Python code to merge dict using update() method '''
    ''' update dict1 with values/keys of dict2 '''
    #return dict2.update(dict1)
    return {**dict1, **dict2}


def window(seq, n=2):
    "Returns a sliding window (of width n) over data from the iterable"
    "   s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...                   "
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


def read_interlevaed_fasta_as_noninterleaved(filename):
    ''' filename (full path) '''

    HANDLE = open(filename, 'r')
    names = []  # record names to maintain order
    d = {}  # dictionary of sequences, key is >line.strip() and in names, value is noninterleaved sequence

    gate = 'closed'
    line = HANDLE.readline()
    while line:
        if '>' not in line:
            d[names[-1]] += line.strip()
        elif '>' in line and gate == 'closed':
            d[line.strip()] = ''
            names.append(line.strip())
        else:
            print('Problem!!!')
        line = HANDLE.readline()
    HANDLE.close()

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


#####  Simple Commands  #####


def longest_dict_value():
    longest = max(map(len, d.values()))  # length of longest dict value
    return longest


def keep_rows_df_greater_than(df, cutoff, n):  # dataframe, cutoff value, column index to filter
    dftrim = df[df.iloc[:, n] > cutoff]

    return dftrim
    # df[df['column name'].map(len) < 2]
######################################################################
##########   Assumes .sam files as input (w/ extension)     ##########
######################################################################


def reverse_complement(x, reverse=True, complement=True):
    if reverse is True:
        x = x[::-1]  # reverse string
    if complement is True:
        trans = str.maketrans('ATGC', 'TACG')
        x = x.translate(trans)  # complement string
    return x

def plot_scatter(x, y, outpath, title_analysis, exitplotcount, countplot, legenlocation='upper right', prefix='NaN', xlab='Distance From Read End', ylab='Depth', figlength=10, figheight=5):
    ''' x is list of x values, y is list of y values, path is full path to output file '''
    import matplotlib.pyplot as plt
    # import math
    # import os
    # font = {'family': 'normal',
    #         'weight': 'normal',
    #         'size': 20}
    if countplot == 1:
        plt.rcParams.update({'font.size': 16})
        plt.figure(figsize=(figlength, figheight))

    # y = [math.log(i, 2) for i in y]
    plt.plot(x, y, label=prefix)
    plt.grid(True)
    plt.title(title_analysis)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.legend(loc=legenlocation, fontsize=10)
    plt.savefig(outpath)
    if countplot == exitplotcount:  # 2 = len(infiles), infiles passed to phasing function
        plt.show()
        plt.close()


def trend_and_periodicity(df, analysis, season=23):
    ''' Calcutate trend, subtract trend from data, determine periodicity '''
    # import os
    # import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    # import seaborn as sns

    # plot trends
    rollmeans = [df.iloc[:, i].rolling(season).mean() for i in range(df.shape[1])]
    df_rm = pd.concat(rollmeans, axis=1)
    df_rm.plot(figsize=(10, 5), linewidth=2, fontsize=16)
    plt.grid(True)
    plt.title('Trend of %s' % analysis, fontsize=16)
    plt.ylabel('Depth', fontsize=16)
    plt.xlabel('5\' -> 5\' distance on same strand', fontsize=16)
    plt.show()
    plt.close()

    # plot "seasonality" of 5' 5' sRNA
    offset = 3
    asymptotes = [i + offset for i in range(df.shape[0]) if i % season == 0 ]
    print(asymptotes)
    df.diff().plot(figsize=(10, 5), linewidth=2, fontsize=16)
    for i in asymptotes:
        plt.axvline(x=i, linewidth=0.5)
    #plt.grid(True)
    plt.title('Periodicity of %s' % analysis, fontsize=16)
    plt.ylabel('Depth', fontsize=16)
    plt.xlabel('5\' -> 5\' distance on same strand', fontsize=16)
    plt.show()
    plt.close()
    print(df.diff().corr())


def make_logo(df):
    # do imports
    import matplotlib.pyplot as plt
    import logomaker as logomaker

    # load crp energy matrix
    # crp_df = -logomaker.get_example_matrix('crp_energy_matrix', print_description=False)

    # create Logo object
    fig_logo = logomaker.Logo(df, shade_below=.5, fade_below=.5)  # , font_name='Arial Rounded MT Bold'

    # style using Logo methods
    fig_logo.style_spines(visible=False)
    fig_logo.style_spines(spines=['left', 'bottom'], visible=True)
    fig_logo.style_xticks(rotation=90, fmt='%d', anchor=0)

    # style using Axes methods
    fig_logo.ax.set_ylabel("Base Frequency", labelpad=-1)
    fig_logo.ax.xaxis.set_ticks_position('none')
    fig_logo.ax.xaxis.set_tick_params(pad=-1)

    # style and show figure
    # fig_logo.fig.show()  # commented this out to remove error message in jupyter lab


def calculate_depth_of_coverage(sam, maxdepth=0, mindepth=0):
    ''' if depth > maxdepth then set depth to maxdepth '''
    ''' if depth < mindepth then set depth to 0 '''
    print('Calculating depth of coverage for: %s' % sam)
    d = {}
    names = []
    numberseqs = 0
    with open(sam, 'r') as FILE:
        for line in FILE:
            if line[0] != '@':
                for i in range(len(line.strip().split('\t')[9])):  # for each base in sequence
                    # subtract one from position because sam file left-most positions are 1-based
                    position = int(line.strip().split('\t')[3]) + i - 1  # calculate coordinate in scaffold/chromosome
                    key = '_'.join([line.strip().split('\t')[2], str(position)])
                    # for each position in each scaffold + 1 if base present
                    d[key] = d.setdefault(key, 0) + 1  # for adding use this structure
                    names.append(key)
                    numberseqs += 1
        for n in names:
            if (maxdepth != 0) and (d[n] > maxdepth):
                d[n] = maxdepth
            if d[n] < mindepth:
                d[n] = 0
    return d, names, numberseqs


def single_stranded_rna(samfilelist_f, samfilelist_r, maxdepth):
    ''' Compares depth of coverage in forward vs reverse files. If imbalanced, maybe source was single stranded RNA '''
    from natsort import natsorted, ns
    import os

    for forward, reverse in zip(samfilelist_f, samfilelist_r):
        d_for, names_for, nseqs_for = calculate_depth_of_coverage(forward, maxdepth)  # input is sam file # output does not include positions of 0 coverage
        d_rev, names_rev, nseqs_rev = calculate_depth_of_coverage(reverse, maxdepth)  # key is scafffold_position # value is depth
        # print('Total number of sequences in %s = %d' % forward, nseqs_for)
        # print('Total number of sequences in %s = %d' % reverse, nseqs_rev)
        # keys will not always be present between forward and reverse seqs. Make them match
        for name in names_for:
            try:
                d_rev[name]
            except:  # there is not key name in d_rev
                d_rev[name] = 0
            else:
                pass
        for name in names_rev:
            try:
                d_for[name]
            except:
                d_for[name] = 0
            else:
                pass

        allnames = list(set().union(names_for, names_rev))  # union two name lists
        # sort key values naturally # scaffolds may not be in same order as sam but coordinates within each scaffold should be ordered
        sortnames = natsorted(allnames, key=lambda y: y.lower())  # or # natsorted(x, alg=ns.IGNORECASE)  # or alg=ns.IC

        outfileprefix = 'SingleStrandRNA.FvsR.Depth'
        analysis = 'Single Stranded RNA Analysis'
        path, file = os.path.split(forward)

        with open(os.path.join(path, '.'.join(file.split('.')[:-1] + [outfileprefix, 'tsv'])), 'w') as OUT:
            OUT.write('\n'.join(['\t'.join([n, str(d_for[n]), str(d_rev[n])]) for n in sortnames]))

        plot_scatter(list(range(len(sortnames))), [d_for[n] for n in sortnames],
                         os.path.join(path, '.'.join(file.split('.')[:-1] + [outfileprefix, 'png'])), analysis,
                         2, 1, 'upper right', '_'.join(file.split('.')[:4]), 'Genomic Position', figlength=20,
                         figheight=10)

        path, file = os.path.split(reverse)
        plot_scatter(list(range(len(sortnames))), [-1*d_rev[n] for n in sortnames],
                         os.path.join(path, '.'.join(file.split('.')[:-1] + [outfileprefix, 'png'])), analysis,
                         2, 2, 'upper right', '_'.join(file.split('.')[:4]), 'Genomic Position', figlength=20,
                         figheight=10)  # plot reverse reads with negative value


def calculate_base_composition(bases, verbose=False, basekey=['A', 'T', 'G', 'C']):
    ''' calculate percentage of each capital letter nucleotide in list bases '''
    # import math
    bases = [b.upper() for b in bases]  # make sure all strings are uppercase
    percentbases = [bases.count(i) / len(bases) for i in basekey]
    totalbases = len(bases)
    if verbose == True:
        for base, percent in zip(basekey, percentbases):
            print('Nucleotide: %s, Percent: %.2f, total: %d' % (base, percent, totalbases))
    # bits = [frequency * math.log2(frequency) for frequency in percentbases ]
    return basekey, percentbases, totalbases  # , bits


def plus_one_u(samfiles, infiles, x=1, orientation='forward'):
    '''  assumes Phasing35 analysis gave positive result @ +1 '''
    ''' find % nucleotide composition @ every 5' base + x bases away from all 3' ends on same strand  '''
    import os
    import pandas as pd

    end1 = 3  # access number of 3' ends
    end2 = 2  # access number of 5' ends
    if orientation == 'reverse':
        x = x * -1

    allbases = []  # list of all + n bases away from 3' ends of reads
    for samfile, f in zip(samfiles, infiles):
        with open(samfile, 'r') as FILE:
            # record first (5') base in sequence from sam file. Reverse reads may need reverse complementation
            # skip header
            #dsam = {'_'.join(line.strip().split('\t')[2:4]) : line.strip().split('\t')[9][0] for line in FILE if line[0] != '@'}
            # keep entire sequence # this only records one sequence per position in scaffold.
            #dsam = {'_'.join(line.strip().split('\t')[2:4]) : line.strip().split('\t')[9] for line in FILE if line[0] != '@'}
            dsam = {}
            for line in FILE:
                if line[0] != '@':
                    #dict_x.setdefault(key, []).append(value)
                    if orientation == 'reverse':
                        # seqs in my sam files are oriented 3'<----|5', so to get 5' coordinate need to all len(seq) - 1 to 3' position coordinate
                        dsam.setdefault('_'.join([line.strip().split('\t')[2], str(int(line.strip().split('\t')[3]) + len(line.strip().split('\t')[9]) - 1)]), []).append(line.strip().split('\t')[9])
                    else:
                        dsam.setdefault('_'.join(line.strip().split('\t')[2:4]), []).append(line.strip().split('\t')[9])

        d = {}
        names = []
        with open(f, 'r') as FILE:
            for line in FILE:
                s = line.strip().split('\t')
                d['_'.join(s[0:2])] = s
                names.append('_'.join(s[0:2]))
        bases = []
        seqs = []
        for name in names:  # basically for each split line
            s = d[name]  # list of elements of line previously called s
            if int(s[end1]) != 0:  # for every position with a 3' end of read
                try:
                    d['_'.join([s[0], str(int(s[1]) + x)])]  # is the + x 5' position present on same strand?
                except:  # trying to pass key errors, because i don't have numbers for every position
                    pass  # skip of no 5' end at + x position on same strand
                else:
                    if int(d['_'.join([s[0], str(int(s[1]) + x)])][end2]) != 0:  # is + x 5' end != 0 from infile
                        try:
                            dsam['_'.join([s[0], str(int(s[1]) + x)])]  # if 5' position exists
                        except:
                            pass
                        else:
                            for seq in dsam['_'.join([s[0], str(int(s[1]) + x)])]:
                                if orientation == 'reverse':
                                    bases.append(reverse_complement(seq[-1], reverse=True, complement=True))  # append(base)
                                    seqs.append(reverse_complement(seq, reverse=True, complement=True))  # append(sequence) # seqs in my sam files are oriented 3'<----|5'
                                else:
                                    bases.append(seq[0])  # append(base)
                                    seqs.append(seq)  # append(sequence)
        allbases.append(bases)

        # output sequences
        outfileprefix = 'Plus%dU.35' % x
        path, file = os.path.split(f)
        outfile = os.path.join(path, '.'.join(file.split('.')[:-1] + [outfileprefix, 'out']))
        print('Writing sequences that have 5\' base + x from 3\' end on same strand to file: %s' % outfile)
        with open(outfile, 'w') as OUT:
            OUT.write('\n'.join(seqs))

        print('Nucleotide composition (frequency) of all 5\' nucleotides + %d bases away from all 3\' ends' % (x))
        print('In file: %s' % (f))
        # calculate base composition of + x position from 3' end
        basekey, basesfreq, totalbases = calculate_base_composition(bases, True)

        # lets make a sequence logo
        allbasefrequencies = []
        tseqs = list(map(list, zip(*seqs)))  # transpose list so each row represents one of 23 bases of pi RNA
        for nucleotides in tseqs:
            basekey, basesfreq, totalbases = calculate_base_composition(nucleotides)
            allbasefrequencies.append(basesfreq)
        dflogo = pd.DataFrame(allbasefrequencies, columns=basekey)

        print('Logo for sequences in file: %s' % (outfile))
        make_logo(dflogo)


def ping_pong(infiles_f, infiles_r, lower=0, upper=20, analysis='PingPong55'):
    ''' Assumes first forward and reverse files are from the same sam/bam file '''
    ''' Main difference from phasing analysis is comparing forward to reverse reads '''
    import os
    from natsort import natsorted, ns

    end1 = 2  # access 5' end information
    end2 = 2  # access 5' end information
    print('lower: %d, upper: %d, orientation not applicable' % (lower, upper))

    allsignals = []
    alldicer = []
    for forward, reverse in zip(infiles_f, infiles_r):
        print(forward)
        print(reverse)
        df, dr = {}, {}
        namesf, namesr = [], []
        with open(forward, 'r') as FILE:
            for line in FILE:
                s = line.strip().split('\t')
                df['_'.join(s[0:2])] = s
                namesf.append('_'.join(s[0:2]))
        with open(reverse, 'r') as FILE:
            for line in FILE:
                s = line.strip().split('\t')
                dr['_'.join(s[0:2])] = s
                namesr.append('_'.join(s[0:2]))

        dicer_signal_location = {}  # record locations and strength of dicer signal at position upper - 3 (should be - 2 if count was 1 based)
        # only need to compare forward to reverse, not reverse to forward. This would be redundant
        signals = [0] * abs(upper - lower)
        count = 0  # keep track of element in list signals that i should add too
        for x in range(lower, upper):  # range excludes the last value of 50 so use 51 # distance from 3' end
            minimums = []
            for name in namesf:  # basically for each name in forward reads
                s = df[name]  # list of elements of line previously called s
                if int(s[end1]) != 0:  # for every position with a 5' end in forward read orientation
                    # now try to compare to 5' ends in reverse file
                    Mi = int(s[end1])  # get 5' counts from current position for ping pong
                    try:
                        dr['_'.join([s[0], str(int(s[1]) + x)])]  # is same scaffold_position is present in reverse?
                    except:  # trying to pass key errors, because i don't have numbers for every position
                        Nix = 0  # pass 0 if key not present in reverse file
                    else:  # get 5' counts from i+x position
                        Nix = int(dr['_'.join([s[0], str(int(s[1]) + x)])][end2])  # get 5' counts from i+x position
                    minimums.append(min(Mi, Nix))  # min(Mi, Nix) # Mi * Nix # assume that reads must be paired 1/1 ratio
                    if (count == upper - 3) and (min(Mi, Nix) != 0):  # count is zero based so it should be -2 normally for dicer products
                        dicer_signal_location[name] = dicer_signal_location.setdefault(name, 0) + min(Mi, Nix)
            signals[count] = str(sum(minimums))  # sum of all minimum pair counts for distance x from position
            count += 1  # add one before going to next distance value
        allsignals.append(signals)
        alldicer.append(dicer_signal_location)

    outfilename = 'PingPong.FvsR.55'
    # add name of files to output
    output = []
    countplot = 1
    for fullpath, signal in zip(infiles_f, allsignals):
        path, file = os.path.split(fullpath)
        output.append('\t'.join([file] + signal))
        plot_scatter(list(range(lower, upper)), [int(sig) for sig in signal], os.path.join(path, '.'.join(file.split('.')[:-1] + [outfilename, 'pdf'])), analysis, len(allsignals), countplot, 'upper left', '_'.join(file.split('.')[:4]))
        countplot += 1

    # output data
    path, file = os.path.split(infiles_f[0])
    outfile = os.path.join(path, '.'.join(file.split('.')[:-1] + [outfilename, 'out']))
    with open(outfile, 'w') as OUT:
        OUT.write('\n'.join(output))

    for fullpath, dicer_dictionary in zip(infiles_f, alldicer):
        path, file = os.path.split(fullpath)
        # sort dictionary by value and output. output is position info as key and sorted by value # of 5-5
        # {k: v for k, v in sorted(x.items(), key=lambda item: item[1])}
        output = [f'{k}\t{v}' for k, v in sorted(dicer_dictionary.items(), key=lambda item: item[1], reverse=True)]
        outfile = os.path.join(path, '.'.join(file.split('.')[:-1] + [outfilename, 'DicerSignal.SortStrength.out']))
        print(f'Writing to output file: {outfile}')
        with open(outfile, 'w') as OUT:
            OUT.write('\n'.join(output))
        output = [f'{k}\t{v}' for k, v in natsorted(dicer_dictionary.items(), key=lambda item: item[0])]
        outfile = os.path.join(path, '.'.join(file.split('.')[:-1] + [outfilename, 'DicerSignal.SortPosition.out']))
        print(f'Writing to output file: {outfile}')
        with open(outfile, 'w') as OUT:
            OUT.write('\n'.join(output))


def phasing(infiles, lower=-10, upper=50, seqlength=23, analysis='Phasing35', orientation='forward'):
    ''' infiles, list of full paths to files output from count53 function '''
    ''' lower int(), upper int(), analysis == 'Phasing35' or 'PeriodicPeaks55' '''
    import os
    import pandas as pd

    if analysis == 'Phasing35':
        end1 = 3  # access 3' end
        end2 = 2  # access 5' end
        print('Performing 3\' to 5\' phasing anaylsis on same strand')
        print('lower = %d, upper = %d, %s' % (lower, upper, orientation))
    elif analysis == 'PeriodicPeaks55':
        end1 = 2  # access 5' end
        end2 = 2  # access 5' end
        print('Performing 5\' to 5\' periodic peaks anaylsis on same strand')
        print('lower = %d, upper = %d, %s' % (lower, upper, orientation))
        # lower, upper = 20, 200
    else:
        print('analysis specification incorrect')
        print('Please specify analysis = \'Phasing35\' or \'PeriodicPeaks55\' when calling function')

    allsignals = []
    for f in infiles:
        d = {}
        names = []
        with open(f, 'r') as FILE:
            for line in FILE:
                s = line.strip().split('\t')
                d['_'.join(s[0:2])] = s
                names.append('_'.join(s[0:2]))
            # for each distance from 3' end sum all minimum values of Mi and Ni+x (Mi = number of 3' bases at position i, Ni+x = number of 5' bases at position i+x

        signals = [0] * abs(upper - lower)
        count = 0  # keep track of element in list signals that i should add too
        for x in range(lower, upper):  # range excludes the last value of 50 so use 51 # distance from 3' end
            if orientation == 'reverse':
                x = x * -1  # for reverse we need +10 -> -50 not -10 -> 50
            minimums = []
            for name in names:  # basically for each split line
                s = d[name]  # list of elements of line previously called s
                if int(s[end1]) != 0:  # for every position with a 3' end of read # or 5' end for periodic peaks
                    Mi = int(s[end1])  # get 3' counts from current position # or 5' for periodic peaks
                    try:
                        d['_'.join([s[0], str(int(s[1]) + x)])]
                    except:  # trying to pass key errors, because i don't have numbers for every position
                        Nix = 0  # pass 0 if key not present in dict
                    else:  # get 5' counts from i+x position
                        Nix = int(d['_'.join([s[0], str(int(s[1]) + x)])][end2])  # get 5' counts from i+x position
                    minimums.append(min(Mi, Nix))
            signals[count] = str(sum(minimums))  # sum of all minimum pair counts (of 3' or 5') for distance x from position
            count += 1  # add one before going to next distance value
        allsignals.append(signals)

    if analysis == 'Phasing35':
        outfilename = 'Phasing.35'
    elif analysis == 'PeriodicPeaks55':
        outfilename = 'PeriodicPeaks.55'
    else:
        raise ValueError('Analysis specification incorrect when calling phasing function')
    # add name of files to output
    output = []
    countplot = 1
    for fullpath, signal in zip(infiles, allsignals):
        path, file = os.path.split(fullpath)
        output.append('\t'.join([file] + signal))
        plot_scatter(list(range(lower, upper)), [int(sig) for sig in signal],
                     os.path.join(path, '.'.join(file.split('.')[:-1] + [outfilename, 'pdf'])),
                     analysis, len(allsignals), countplot, 'upper right', '_'.join(file.split('.')[:4]))
        countplot += 1

    # output data
    path, file = os.path.split(infiles[0])
    outfile = os.path.join(path, '.'.join(file.split('.')[:-1] + [outfilename, 'out']))
    with open(outfile, 'w') as OUT:
        OUT.write('\n'.join(output))

    if analysis == 'PeriodicPeaks55':
        colnames = []
        for fullpath in infiles:
            path, file = os.path.split(fullpath)
            colnames.append('.'.join(file.split('.')[:4]))
        allsignals = [[int(signal) for signal in signals] for signals in allsignals]
        tallsignals = list(map(list, zip(*allsignals)))
        df = pd.DataFrame(tallsignals, columns=colnames)
        trend_and_periodicity(df, analysis, seqlength)


def get_rightmost_reference_based_alignment_coordinate(CIGAR, leftmost_coordinate):
    import re
    cigar = re.findall(r'\d+[A-Z]', CIGAR)
    if cigar == []:  # if there was a match # sometimes CIGAR string == * # this skips unmapped reads
        print(f'Provided CIGAR string: {CIGAR} does not match CIGAR pattern \\d+[A-Z]')
        rightmost_position = 0  # assumes unmapped read
    else:  # then read should be mapped
        rightmost_position = leftmost_coordinate - 1  # subtract 1 because leftmost base is 1-based
        for i in cigar:
            if i[-1] in ['M', 'N', 'D', 'X', '=']:
                rightmost_position += int(i[:-1])
            elif i[-1] in ['I', 'S', 'H', 'P']:
                pass
            else:
                pass

    return rightmost_position


def count53(filelist, maxdepth, orientation='forward'):
    from natsort import natsorted, ns

    outnames = []
    for f in filelist:
        print('Working on: %s' % (f))
        # filelist = ['D:\\LinuxShare\\Projects\\Theresa\\Hisat2\\FlagHA_Pt08\\Pt_51_MacAndIES\\52_Late.Uniq.23M.F.sort.sam', 'D:\\LinuxShare\\Projects\\Theresa\\Hisat2\\FlagHA_Pt08\\Pt_51_MacAndIES\\52_Late.Uniq.23M.R.sort.sam', 'D:\\LinuxShare\\Projects\\Theresa\\Hisat2\\FlagHA_Pt08\\Pt_51_MacAndIES\\53_Later.Uniq.23M.F.sort.sam', 'D:\\LinuxShare\\Projects\\Theresa\\Hisat2\\FlagHA_Pt08\\Pt_51_MacAndIES\\53_Later.Uniq.23M.R.sort.sam']
        if maxdepth > 0:  # if maxdepth specified chance output file name
            outfile = f[:-len('.sam')] + '.53.MaxDepth.tsv'
        else:
            outfile = f[:-len('.sam')] + '.53.tsv'
        outnames.append(outfile)

        d5 = {}
        d3 = {}
        with open(f, 'r') as FILE:
            for line in FILE:
                line = line.strip()
                if line[0] == '@':
                    pass
                else:
                    leftmost_coordinate = int(line.split('\t')[3])
                    CIGAR = line.strip().split('\t')[5]
                    if orientation == 'forward':
                        scaff_pos5 = '_'.join(line.split('\t')[2:4])  # scaffold name + 5' position
                        # subtract one from coordinate below because 3' end # reasoning: start=11, end = 15, len=5, but 11+5 = 16 not 15
                        # scaff_pos3 = '_'.join([line.split('\t')[2]] + [str(int(line.split('\t')[3]) + len(line.split('\t')[9]) - 1)])  # scaffold name + 3' position
                        scaff_pos3 = '_'.join([line.split('\t')[2]] + [str(get_rightmost_reference_based_alignment_coordinate(CIGAR, leftmost_coordinate))])  # scaffold name + 3' position
                    elif orientation == 'reverse':  # switch 5' and 3' values for reverse reads
                        scaff_pos3 = '_'.join(line.split('\t')[2:4])  # scaffold name + 3' position
                        # subtract one from coordinate below because 3' end # reasoning: start=11, end = 15, len=5, but 11+5 = 16 not 15
                        # scaff_pos5 = '_'.join([line.split('\t')[2]] + [str(int(line.split('\t')[3]) + len(line.split('\t')[9]) - 1)])  # scaffold name + 5' position
                        scaff_pos5 = '_'.join([line.split('\t')[2]] + [str(get_rightmost_reference_based_alignment_coordinate(CIGAR, leftmost_coordinate))])
                    # could have done below two lines, but doesn't keep order
                    # d5[scaff_pos5].get(scaff_pos5, 0) + 1  # kinda slow # counts number of 5' positions for all scaffolds
                    # d3[scaff_pos3].get(scaff_pos3, 0) + 1  # kinda slow # counts number of 3' positions for all scaffolds
                    # count 5' ends at position
                    try:
                        d5[scaff_pos5]
                    except:
                        d5[scaff_pos5] = 1  # initiate with 1 because we want to count the first
                    else:
                        d5[scaff_pos5] += 1
                    try:  # to maintain the same keys between d5 and d3
                        d5[scaff_pos3]
                    except:  # if there is no scaff_pos3 entry make it zero
                        d5[scaff_pos3] = 0

                    # count 3' ends at each position
                    try:
                        d3[scaff_pos3]
                    except:
                        d3[scaff_pos3] = 1  # initiate  with 1 because we want to count the first
                    else:
                        d3[scaff_pos3] += 1
                    try:  # to maintain the same keys between d5 and d3
                        d3[scaff_pos5]
                    except:  # if there is no scaff_pos5 entry in d3 make it zero
                        d3[scaff_pos5] = 0
        # sort key values naturally # scaffolds may not be in same order as sam but coordinates within each scaffold should be ordered
        names = list(d5.keys())  # assumes d5 and d3 have same key values
        sortnames = natsorted(names, key=lambda y: y.lower())  # or # natsorted(x, alg=ns.IGNORECASE)  # or alg=ns.IC
        # print('First sorted names:')
        # print(sortnames[:4])

        with open(outfile, 'w') as OUT:
            output = []
            for name in sortnames:
                scaffold, position = '_'.join(name.split('_')[:-1]), name.split('_')[-1]
                if maxdepth == 0:  # then do not control depth # proceed without changing depth values
                    output.append('\t'.join([scaffold, position, str(d5[name]), str(d3[name])]))
                else:  # enforce specified maxdepth value
                    if d5[name] > maxdepth:
                        d5[name] = maxdepth
                    if d3[name] > maxdepth:
                        d3[name] = maxdepth
                    output.append('\t'.join([scaffold, position, str(d5[name]), str(d3[name])]))
            OUT.write('\n'.join(output))
        print('Output file of 5\' and 3\' counts: %s\n' % (outfile))
    return(outnames)


def main(samfilelist_f, samfilelist_r, seqlength=23, maxdepth=0):  # 0 maxdepth will not limit depth > 0 will limit to max depth
    infiles_f = count53(samfilelist_f, maxdepth, orientation='forward')  # infiles = scaffold\tposition\t#of5'ends\t#of3'ends
    infiles_r = count53(samfilelist_r, maxdepth, orientation='reverse')
    print('##################################################')
    print('Start phasing analysis')
    phasing(infiles_f, -10, 50, seqlength, analysis='Phasing35', orientation='forward')
    phasing(infiles_r, -10, 50, seqlength, analysis='Phasing35', orientation='reverse')
    print('Done with phasing analysis')
    print('##################################################')
    print('Starting periodic peaks analysis')
    phasing(infiles_f, 20, 200, seqlength, analysis='PeriodicPeaks55', orientation='forward')
    phasing(infiles_r, 20, 200, seqlength, analysis='PeriodicPeaks55', orientation='reverse')
    print('Done with periodic peaks analysis')
    print('##################################################')
    print('Starting ping pong analysis')
    ping_pong(infiles_f, infiles_r, 0, seqlength)  # compare 5' forward to 5' reverse reads and vice versa
    print('Done with ping pong analysis')
    print('##################################################')
    print('Starting: Is my RNA single stranded? analysis')
    single_stranded_rna(samfilelist_f, samfilelist_r, maxdepth)
    print('Done with Is my RNA single stranded? analysis')
    print('##################################################')
    print('Starting +1 U analysis')
    plus_one_u(samfilelist_f, infiles_f, 1)
    plus_one_u(samfilelist_r, infiles_r, 1, orientation='reverse')
    print('Done with +1 U analysis')
    print('##################################################')

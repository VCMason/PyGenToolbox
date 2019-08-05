######################################################################
##########   Assumes .sam files as input (w/ extension)     ##########
######################################################################

from natsort import natsorted, ns


def plot_scatter(x, y, outpath, prefix='NaN'):
    ''' x is list of x values, y is list of y values, path is full path to output file '''
    import matplotlib.pyplot as plt
    # import os
    plt.plot(x, y, label=prefix)
    plt.grid(True)
    # plt.title(os.path.split(outpath)[1])
    plt.xlabel('All Scaffold Positions')
    plt.ylabel('Depth')
    plt.legend()
    plt.savefig(outpath)
    plt.close()


def periodic_peaks_forward(infiles, lower=15, upper=200):
    ''' infiles, list of full paths to files output from count53 function '''
    ''' lower int(), upper int() '''
    import os

    allsignals = []
    for f in infiles:
        d = {}
        names = []
        with open(f, 'r') as FILE:
            for line in FILE:
                s = line.strip().split('\t')
                d['_'.join(s[0:2])] = s
                names.append('_'.join(s[0:2]))
            # for each distance from 5' end sum all minimum values of Mi and Ni+x (Mi = number of 5' bases at position i, Ni+x = number of 5' bases at position i+x

        signals = [0] * abs(upper - lower)
        count = 0  # keep track of element in list signals that i should add too
        for x in range(lower, upper):  # range excludes the last value of 50, could use 51 # distance from 5' end
            minimums = []
            for name in names:  # basically for each split line
                s = d[name]  # list of elements of line previously called s
                if s[2] != 0:  # for every position with a 5' end of read
                    Mi = int(s[2])  # get 5' counts from current position
                    try:
                        d['_'.join([s[0], str(int(s[1]) + x)])]
                    except:  # trying to pass key errors, because i don't have numbers for every position
                        Nix = 0  # pass 0 if key not present in dict
                    else: # get 5' counts from i+x position
                        Nix = int(d['_'.join([s[0], str(int(s[1]) + x)])][2])  # get 5' counts from i+x position
                    minimums.append(min(Mi, Nix))
            signals[count] = str(sum(minimums))  # sum of all minimum pair counts (of 3' or 5') for distance x from position
            count += 1  # add one before going to next distance value
        allsignals.append(signals)

    # add name of files to output
    output = []
    for fullpath, signal in zip(infiles, allsignals):
        path, file = os.path.split(fullpath)
        output.append('\t'.join([file] + signal))
        plot_scatter(list(range(lower, upper)), [int(sig) for sig in signal], os.path.join(path, 'Phasing.35.pdf'), '_'.join(file.split('.')[:4]))

    # output data
    path, file = os.path.split(infiles[0])
    outfile = os.path.join(path, 'phasing.35.out')
    with open(outfile, 'w') as OUT:
        OUT.write('\n'.join(output))



def phasing_reverse(infiles, lower, upper):
    ''' infiles, list of full paths to files output from count53 function '''
    ''' lower int(), upper int() '''
    import os

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
            x = x*-1  # for reverse we need +10 -> -50 not -10 -> 50
            minimums = []
            for name in names:  # basically for each split line
                s = d[name]  # list of elements of line previously called s
                if s[3] != 0:  # for every position with a 3' end of read
                    Mi = int(s[3])  # get 3' counts from current position
                    try:
                        d['_'.join([s[0], str(int(s[1]) + x)])]
                    except:  # trying to pass key errors, because i don't have numbers for every position
                        Nix = 0  # pass 0 if key not present in dict
                    else: # get 5' counts from i+x position
                        Nix = int(d['_'.join([s[0], str(int(s[1]) + x)])][2])  # get 5' counts from i+x position
                    minimums.append(min(Mi, Nix))
            signals[count] = str(sum(minimums))  # sum of all minimum pair counts (of 3' or 5') for distance x from position
            count += 1  # add one before going to next distance value
        allsignals.append(signals)
    # add name of files to output
    output = []
    for fullpath, signal in zip(infiles, allsignals):
        path, file = os.path.split(fullpath)
        output.append('\t'.join([file] + signal))
        plot_scatter(list(range(lower, upper)), [int(sig) for sig in signal], os.path.join(path, 'Phasing.35.pdf'), '_'.join(file.split('.')[:4]))

    # output data
    path, file = os.path.split(infiles[0])
    outfile = os.path.join(path, 'phasing.35.out')
    with open(outfile, 'w') as OUT:
        OUT.write('\n'.join(output))

def phasing_forward_copy(infiles, lower, upper, analysis='Phasing35', orientation='forward'):
    ''' infiles, list of full paths to files output from count53 function '''
    ''' lower int(), upper int(), analysis == 'Phasing35' or 'PeriodicPeaks55' '''
    import os

    if analysis == 'Phasing35':
        end1 = 3  # access 3' end
        end2 = 2  # access 5' end
        print('Performing 3\' to 5\' phasing anaylsis on same strand')
    elif analysis == 'PeriodicPeaks55':
        end1 = 2  # access 5' end
        end2 = 2  # access 5' end
        print('Performing 5\' to 5\' periodic peaks anaylsis on same strand')
        print('setting lower = 20, and upper = 200')
        lower, upper = 20, 200
    else:
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
                x = x * -1
            minimums = []
            for name in names:  # basically for each split line
                s = d[name]  # list of elements of line previously called s
                if s[3] != 0:  # for every position with a 3' end of read
                    Mi = int(s[end1])  # get 3' counts from current position
                    try:
                        d['_'.join([s[0], str(int(s[1]) + x)])]
                    except:  # trying to pass key errors, because i don't have numbers for every position
                        Nix = 0  # pass 0 if key not present in dict
                    else: # get 5' counts from i+x position
                        Nix = int(d['_'.join([s[0], str(int(s[1]) + x)])][end2])  # get 5' counts from i+x position
                    minimums.append(min(Mi, Nix))
            signals[count] = str(sum(minimums))  # sum of all minimum pair counts (of 3' or 5') for distance x from position
            count += 1  # add one before going to next distance value
        allsignals.append(signals)

    # add name of files to output
    output = []
    for fullpath, signal in zip(infiles, allsignals):
        path, file = os.path.split(fullpath)
        output.append('\t'.join([file] + signal))
        plot_scatter(list(range(lower, upper)), [int(sig) for sig in signal], os.path.join(path, 'Phasing.35.pdf'), '_'.join(file.split('.')[:4]))

    # output data
    path, file = os.path.split(infiles[0])
    outfile = os.path.join(path, 'phasing.35.out')
    with open(outfile, 'w') as OUT:
        OUT.write('\n'.join(output))


def phasing_forward(infiles, lower, upper):
    ''' infiles, list of full paths to files output from count53 function '''
    ''' lower int(), upper int() '''
    import os

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
            minimums = []
            for name in names:  # basically for each split line
                s = d[name]  # list of elements of line previously called s
                if s[3] != 0:  # for every position with a 3' end of read
                    Mi = int(s[3])  # get 3' counts from current position
                    try:
                        d['_'.join([s[0], str(int(s[1]) + x)])]
                    except:  # trying to pass key errors, because i don't have numbers for every position
                        Nix = 0  # pass 0 if key not present in dict
                    else: # get 5' counts from i+x position
                        Nix = int(d['_'.join([s[0], str(int(s[1]) + x)])][2])  # get 5' counts from i+x position
                    minimums.append(min(Mi, Nix))
            signals[count] = str(sum(minimums))  # sum of all minimum pair counts (of 3' or 5') for distance x from position
            count += 1  # add one before going to next distance value
        allsignals.append(signals)

    # add name of files to output
    output = []
    for fullpath, signal in zip(infiles, allsignals):
        path, file = os.path.split(fullpath)
        output.append('\t'.join([file] + signal))
        plot_scatter(list(range(lower, upper)), [int(sig) for sig in signal], os.path.join(path, 'Phasing.35.pdf'), '_'.join(file.split('.')[:4]))

    # output data
    path, file = os.path.split(infiles[0])
    outfile = os.path.join(path, 'phasing.35.out')
    with open(outfile, 'w') as OUT:
        OUT.write('\n'.join(output))



def count53_reverse(filelist):
    outnames = []
    for f in filelist:
        print('Working on: %s' % (f))
        # filelist = ['D:\\LinuxShare\\Projects\\Theresa\\Hisat2\\FlagHA_Pt08\\Pt_51_MacAndIES\\52_Late.Uniq.23M.F.sort.sam', 'D:\\LinuxShare\\Projects\\Theresa\\Hisat2\\FlagHA_Pt08\\Pt_51_MacAndIES\\52_Late.Uniq.23M.R.sort.sam', 'D:\\LinuxShare\\Projects\\Theresa\\Hisat2\\FlagHA_Pt08\\Pt_51_MacAndIES\\53_Later.Uniq.23M.F.sort.sam', 'D:\\LinuxShare\\Projects\\Theresa\\Hisat2\\FlagHA_Pt08\\Pt_51_MacAndIES\\53_Later.Uniq.23M.R.sort.sam']
        outfile = f[:-len('.sam')] + '.53.tsv'
        outnames.append(outfile)

        d5 = {}
        d3 = {}
        with open(f, 'r') as FILE:
            for line in FILE:
                line = line.strip()
                if line[0] == '@':
                    pass
                else:  # switch 5' and 3' values for reverse reads
                    scaff_pos3 = '_'.join(line.split('\t')[2:4])  # scaffold name + 3' position
                    # subtract one from coordinate below because 3' end # reasoning: start=11, end = 15, len=5, but 11+5 = 16 not 15
                    scaff_pos5 = '_'.join([line.split('\t')[2]] + [str(int(line.split('\t')[3]) + len(line.split('\t')[9]) - 1)])  # scaffold name + 5' position
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
        print('First sorted names:')
        print(sortnames[:4])

        with open(outfile, 'w') as OUT:
            output = []
            for name in sortnames:
                scaffold, position = '_'.join(name.split('_')[:-1]), name.split('_')[-1]
                output.append('\t'.join([scaffold, position, str(d5[name]), str(d3[name])]))
            OUT.write('\n'.join(output))
        print('Output file of 5\' and 3\' counts: %s\n' % (outfile))
    return(outnames)



def count53_forward(filelist):
    outnames = []
    for f in filelist:
        print('Working on: %s' % (f))
        # filelist = ['D:\\LinuxShare\\Projects\\Theresa\\Hisat2\\FlagHA_Pt08\\Pt_51_MacAndIES\\52_Late.Uniq.23M.F.sort.sam', 'D:\\LinuxShare\\Projects\\Theresa\\Hisat2\\FlagHA_Pt08\\Pt_51_MacAndIES\\52_Late.Uniq.23M.R.sort.sam', 'D:\\LinuxShare\\Projects\\Theresa\\Hisat2\\FlagHA_Pt08\\Pt_51_MacAndIES\\53_Later.Uniq.23M.F.sort.sam', 'D:\\LinuxShare\\Projects\\Theresa\\Hisat2\\FlagHA_Pt08\\Pt_51_MacAndIES\\53_Later.Uniq.23M.R.sort.sam']
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
                    scaff_pos5 = '_'.join(line.split('\t')[2:4])  # scaffold name + 5' position
                    # subtract one from coordinate below because 3' end # reasoning: start=11, end = 15, len=5, but 11+5 = 16 not 15
                    scaff_pos3 = '_'.join([line.split('\t')[2]] + [str(int(line.split('\t')[3]) + len(line.split('\t')[9]) - 1)])  # scaffold name + 3' position
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
        print('First sorted names:')
        print(sortnames[:4])

        with open(outfile, 'w') as OUT:
            output = []
            for name in sortnames:
                scaffold, position = '_'.join(name.split('_')[:-1]), name.split('_')[-1]
                output.append('\t'.join([scaffold, position, str(d5[name]), str(d3[name])]))
            OUT.write('\n'.join(output))
        print('Output file of 5\' and 3\' counts: %s\n' % (outfile))
    return(outnames)


def main(filelist_f, filelist_r, lower=-10, upper=50):
    infiles_f = count53_forward(filelist_f)
    infiles_r = count53_reverse(filelist_r)
    print('Start phasing analysis')
    phasing_forward(infiles_f, lower, upper)
    phasing_reverse(infiles_r, lower, upper)
    periodic peaks_forward(infiles_f):
    # periodic_peaks_reverse(infiles_r)
    print('Done with phasing analysis')

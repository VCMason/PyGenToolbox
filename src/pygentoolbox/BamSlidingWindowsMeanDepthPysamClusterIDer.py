
def plot_per_base_strandedness(dclusters, clusternames, prefix, dbed, bednames, numplots=4):
    import matplotlib.pyplot as plt

    for i in range(0, len(clusternames), numplots):
        numsubplots = int(numplots / 2)
        fig, ax = plt.subplots(numsubplots, numsubplots, figsize=(10, 10))  # sharex=False
        a, b = 0, 0
        for j in range(i, i + numplots):
            if j == len(clusternames):
                break  # break to avoid out of range clusternames[j]
            if j == i:
                pass
            elif j == i + numplots:
                b += 1
            elif j % 2 == 1:
                b += 1
            elif j % 2 == 0:
                a += 1
                b -= 1
            name = clusternames[j]
            scaff = name.split('|')[0]
            start = name.split('|')[1]
            end = name.split('|')[2]
            # srcids is list of all src cluster ids that overlap window
            srcids = src_bed_overlap(dbed, bednames, scaff, int(start), int(end))
            ids = ','.join(srcids)
            if ids == '':
                ids = 'NewCluster'
            x = dclusters[name][0]  # genomic positions, x axis coordinates
            y = dclusters[name][1]  # number of forward reads per x axis coordinate
            z = dclusters[name][2]  # number of reverse reads per x axis coordinate
            ax[a, b].plot(x, y, color='blue', label='F')
            ax[a, b].plot(x, z, color='orange', label='R')
            ax[a, b].set_title(ids + ': ' + name[len('scaffold51_'):])
            ax[a, b].set_xlim([int(start), int(end)])
            ax[a, b].legend(loc='upper left')
            for k, label in enumerate(ax[a, b].get_xticklabels()):
                if k % 2 == 1:
                    label.set_visible(True)
                elif k % 2 == 0:
                    label.set_visible(False)
        outpath = prefix + 'perbase.strandedness.plots' + str(i) + '.pdf'
        print(outpath)
        plt.savefig(outpath)
        plt.show()
        plt.close()

    return


def per_base_strandedness(bamfile, windowfile):
    import pysam

    with open(windowfile) as FILE:
        windows = [line.strip().split('\t') for line in FILE]

    samfile = pysam.AlignmentFile(bamfile, "rb")

    dclusters = {}
    clusternames = []
    for window in windows:
        scaff = window[0]
        start = int(window[1])
        end = int(window[2])
        name = '|'.join([scaff, window[1], window[2]])
        clusternames.append(name)
        coordinates, f, r = [], [], []  # f = sum of all forward reads per base, r = sum of all reverse reads per base
        # positions = [i for i in range(start, end)]
        # for pos in range(start, end):
        for pileupcolumn in samfile.pileup(scaff, start, end):  # pileup gives depth values outside region
            # if pileupcolumn.pos == pos:  # it's a bitch but need to check that the position is correct
            if (pileupcolumn.pos <= end) and (pileupcolumn.pos >= start):  # it's a bitch but need to check that the position is correct
                forward, reverse = 0, 0
                for pileupread in pileupcolumn.pileups:
                    if not pileupread.is_del and not pileupread.is_refskip:
                        if pileupread.alignment.is_reverse == True:
                            reverse -= 1
                        elif pileupread.alignment.is_reverse == False:
                            forward += 1
                coordinates.append(pileupcolumn.pos)
                f.append(forward)
                r.append(reverse)
        dclusters[name] = [coordinates, f, r]

    return dclusters, clusternames


def gff3_gene_overlap(dgff3, gff3names, scaff, start, end):
    geneids = []
    orientations = []
    for k in gff3names:
        genescaff = dgff3[k][0]
        genestart = int(dgff3[k][3])
        geneend = int(dgff3[k][4])
        geneorientation = dgff3[k][6]
        id = dgff3[k][8].split(';')[0][3:]  # isolate gene id, from metadata column
        if (scaff == genescaff) and (start <= geneend) and (end >= genestart):
            geneids.append(id)
            orientations.append(geneorientation)

    return geneids, orientations


def src_bed_overlap(dbed, bednames, scaff, start, end):
    srcids = []
    for k in bednames:
        bedscaff = dbed[k][0]
        bedstart = int(dbed[k][1])
        bedend = int(dbed[k][2])
        id = dbed[k][3]
        if (scaff == bedscaff) and (start <= bedend) and (end >= bedstart):
            srcids.append(id)

    return srcids


def fasta_seq_gc_percent(d, scaff, start, end):
    GCcount = 0
    seq = d[scaff][start:end+1]
    for base in seq:
        if (base == 'G') or (base == 'C'):
            GCcount += 1
    refGC = round((GCcount / len(seq)) * 100, 2)

    return refGC, seq


def strandedness(bamfile, windowfile, dbed, bednames, dgff3, gff3names, dgenome):
    import pysam

    with open(windowfile) as FILE:
        windows = [line.strip().split('\t') for line in FILE]

    newwindows = []
    samfile = pysam.AlignmentFile(bamfile, "rb")  # pysam.AlignmentFile("ex1.bam", "rb")
    for window in windows:
        forward, reverse = 0, 0
        scaff = window[0]
        start = int(window[1])
        end = int(window[2])
        length = end - start + 1
        refseqGCpercent, seq = fasta_seq_gc_percent(dgenome, scaff, start, end)
        for read in samfile.fetch(scaff, start, end):
            if read.is_reverse == True:
                reverse += 1
            elif read.is_reverse == False:
                forward += 1

        totalreads = forward + reverse
        strandfor = round((forward / totalreads) * 100, 2)
        strandrev = round((reverse / totalreads) * 100, 2)
        srcids = src_bed_overlap(dbed, bednames, scaff, start, end)  # srcids is list of all src cluster ids that overlap window
        ids = ','.join(srcids)
        if ids == '':
            ids = 'NewCluster'

        geneids, geneorientations = gff3_gene_overlap(dgff3, gff3names, scaff, start, end)
        polarity = []
        gids = ','.join(geneids)
        orientations = ','.join(geneorientations)
        if gids == '':
            gids = 'NA'  # no gene overlap
            orientations = 'NA'
            polarity.append('X')  # polarity will be S, A, N, -, or X
        else:
            for gorientation in geneorientations:
                if strandfor >= 75.0:
                    if gorientation == '+':
                        polarity.append('S')  # S == sense
                    elif gorientation == '-':
                        polarity.append('A')  # A == antisense
                    elif gorientation == '.':
                        polarity.append('N')  # N == no gene orientation information
                elif strandrev >= 75.0:
                    if gorientation == '+':
                        polarity.append('A')  # S == sense
                    elif gorientation == '-':
                        polarity.append('S')  # A == antisense
                    elif gorientation == '.':
                        polarity.append('N')  # N == no gene orientation information
                else:  # if if strandfor and strandrev are both < 75.0 (no majority of reads pointing in same direction)
                    polarity.append('-')  # - == no strong strand bias for cluster (can't say sense or antisense)
        pols = ','.join(polarity)

        newwindows.append('\t'.join([ids] + window + [str(refseqGCpercent), str(length), str(forward), str(reverse), str(strandfor), str(strandrev), gids, orientations, pols, seq]))

    header = '\t'.join(['srcID', 'scaffold', 'start', 'end', 'average_depth', 'GCPercent_reads', 'GCPercent_reference', 'length', 'forward_reads', 'reverse_reads', 'forward_percent', 'reverse_percent', 'geneIDs', 'gene_orientations', 'clusterVSgene_polarity', 'sequence'])
    samfile.close()

    outwindows = [header] + newwindows

    return outwindows


def refine_window_coordinates(sortwindows, bamfile):
    import pysam
    import numpy as np

    refinedwindows = []
    positions = set()
    # d = {}  # i don't keep all read depth values
    samfile = pysam.AlignmentFile(bamfile, "rb")  # pysam.AlignmentFile("ex1.bam", "rb")
    for win in sortwindows:
        scaff = win.split('\t')[0]
        winstart = int(win.split('\t')[1])
        winend = int(win.split('\t')[2])
        print(f'Refining Window Coordinates.')
        # array.arrays A, C, G, T. bases, per coordinate
        coverage = samfile.count_coverage(scaff, winstart, winend, read_callback='nofilter')
        numbases = np.sum(coverage)
        depth = list(np.sum(coverage, axis=0))  # sum four arrays vertically with np.sum axis=0
        # numATs = np.sum(coverage[0] + coverage[3])  # number of AT bases in window
        numGCs = np.sum(coverage[1] + coverage[2])  # total number of GC bases for all aligned reads in window
        countstart, countend = 0, 0
        while depth[0] == 0:
            countstart += 1
            del depth[0]
        while depth[-1] == 0:
            countend += 1
            del depth[-1]
        refinedstart = winstart + countstart
        refinedend = winend - countend
        refinedwinlength = refinedend - refinedstart + 1
        average_depth = mean_list(depth, denominator=refinedwinlength)
        GCpercent = round((numGCs / numbases) * 100, 2)  # GC percentage of all aligned reads
        refinedwindows.append([scaff, str(refinedstart), str(refinedend), str(average_depth), str(GCpercent)])

    samfile.close()

    return refinedwindows, positions


def write_to_bed(windows, outpath):
    print('Writing windows to: %s' % outpath)
    with open(outpath, 'w') as OUT:
        OUT.write('\n'.join(windows))


def collapse_windows(sortwindows):
    # this function is janky
    # sortwindows is a naturally sorted list of strings
    # the strings inside the list are formatted as scaffold\tstartposition\tendposition
    collapsedsortwindows = []
    collapsewindepth = []
    collapseGC = []
    overlapgate = 1  # record start position of first window in overlap
    for count, win in enumerate(sortwindows):
        scaff = win.split('\t')[0]
        winstart = int(win.split('\t')[1])
        winend = int(win.split('\t')[2])
        windepth = float(win.split('\t')[3])
        GCpercent = float(win.split('\t')[4])
        if count == 0:
            pass
        elif (scaff == prevscaff) and (winstart <= prevwinend) and (prevwinstart <= winend):
            # some type of overlap # also equivalent to (i think) a[1] > b[0] and a[0] < b[1]
            # determine if there is overlap between previous window and current window, then combine them
            if overlapgate == 1:
                collapsewindow = '%s\t%s' % (scaff, prevwinstart)  # will add end coordinate when i find end of overlap
                overlapgate = 0
                collapsewindepth.append(prevwindepth)
                collapsewindepth.append(windepth)
                collapseGC.append(prevwinGC)
                collapseGC.append(GCpercent)
        elif overlapgate == 0:
            # if it makes it to this point, there is no overlap with previous window
            # but if overlapgate == 0 then there was overlap between at least two windows and we combine coordinates
            collapsedsortwindows.append('%s\t%s\t%.2f\t%.2f' % (collapsewindow, prevwinend, sum(collapsewindepth) / len(collapsewindepth), sum(collapseGC) / len(collapseGC)))
            collapsewindepth = []  # reset for next series of overlapping windows
            collapseGC = []
            overlapgate = 1  # reset overlapgate to allow for first start coordinates to be recorded
        elif overlapgate == 1:
            # if no overlap with previous, and overlapgate == 1 then no windows were overlapping
            # check if this window WILL overlap with the NEXT window and add it to collapsedsortwindows if NO
            if win != sortwindows[-1]:  # need to have this to avoid list indexing error
                if (scaff == sortwindows[count+1].split('\t')[0]) and (winstart <= int(sortwindows[count+1].split('\t')[2])) and (int(sortwindows[count+1].split('\t')[1]) <= winend):
                    pass
                else:
                    # no overlap with NEXT window
                    collapsedsortwindows.append(win)
                    collapsewindepth = []
                    collapseGC = []
            else:  # if it is the last window here output the window
                collapsedsortwindows.append(win)
                collapsewindepth = []
                collapseGC = []

        prevscaff = win.split('\t')[0]
        prevwinstart = int(win.split('\t')[1])
        prevwinend = int(win.split('\t')[2])
        prevwindepth = float(win.split('\t')[3])
        prevwinGC = float(win.split('\t')[4])

    return collapsedsortwindows


def mean_list(lst, denominator):
    mean = sum(lst) / denominator
    return mean


def mean_with_denominator(value, denominator):
    # mean = sum(lst) / denominator
    mean = value / denominator
    return mean


def calculate_depth_with_pysam_sliding_windows(bamfile, dscaff_max, step=25, winsize=100, cutoff=100):
    import pysam
    import numpy as np

    windows = []
    positions = set()
    # d = {}  # i don't keep all read depth values
    samfile = pysam.AlignmentFile(bamfile, "rb")  #pysam.AlignmentFile("ex1.bam", "rb")
    for scaff, maxcoord in dscaff_max.items():
        print(scaff, end=', ')
        for i in range(0, maxcoord, step):
            ## tempdepth = []  # temppos = []
            ## .pileup skips bases when depth = 0  # need to figure out mpileup -a option in pysam
            ## tempdepth = [pileupcolumn.nsegments for pileupcolumn in samfile.pileup(scaff, i, i + winsize)]
            #tempdepth = []
            #for pos in range(i, i + winsize):
            #    depth = 0
            #    for pileupcolumn in samfile.pileup(scaff, pos, pos+1):  # pileup gives depth values outside region
            #        if pileupcolumn.pos == pos:  # it's a bitch but need to check that the position is correct
            #            #for pileupread in pileupcolumn.pileups:
            #            #    if not pileupread.is_del and not pileupread.is_refskip:
            #            #        depth += 1
            #            #tempdepth.append(depth)
            #            ##tempdepth.append(pileupcolumn.nsegments)  # pileupcolumn.nsegments # n depricated
            coverage = samfile.count_coverage(scaff, i, i + winsize, read_callback='nofilter')  # array.arrays A, C, G, T. bases, per coordinate
            numbases = np.sum(coverage)
            # numATs = np.sum(coverage[0] + coverage[3])  # number of AT bases in window
            numGCs = np.sum(coverage[1] + coverage[2])  # number of GC bases in window

            if numbases >= 10:
                average_depth = mean_with_denominator(numbases, denominator=winsize)
                GCpercent = round((numGCs / numbases) * 100, 2)
            else:
                average_depth = 0
                GCpercent = 0
            if average_depth >= cutoff:
                # print(f'{average_depth}, {scaff}, {i+1}, {i+winsize}, {tempdepth}')
                # print(str(average_depth), end=', ')
                # scaff, start, end, averagedepth
                windows.append([scaff, str(i + 1), str(i + winsize), str(average_depth), str(GCpercent)])
                for j in range(i + 1, i + winsize):
                    if j <= maxcoord:
                        positions.add('_'.join([scaff, str(j)]))
    samfile.close()

    return windows, positions


def read_gff3(f):
    # assumes gff3 is trimmed to 'gene' lines only (no exon/CDS/whatever else), one line per gene
    dgff3 = {}
    gff3names = []
    with open(f, 'r') as FILE:
        for line in FILE:
            if line[0] != '#':
                s = line.strip().split('\t')
                dgff3['|'.join([s[0], s[3], s[4]])] = line.strip().split('\t')
                gff3names.append('|'.join([s[0], s[3], s[4]]))

    return dgff3, gff3names


def read_bed(f):
    dbed = {}
    bednames = []
    with open(f, 'r') as FILE:
        header = FILE.readline()
        for line in FILE:
            dbed['|'.join(line.strip().split('\t')[:3])] = line.strip().split('\t')
            bednames.append('|'.join(line.strip().split('\t')[:3]))

    return dbed, bednames


def read_fasta_as_dict(f):
    d = {}  # fasta names are key, values are sequences as string
    namelist = []
    with open(f, 'r') as FILE:
        for line in FILE:
            if line[0] == '>':
                if ' ' in line:
                    name = line.strip().split()[0][1:-len('_with_IES')]
                    namelist.append(name)
                    d[name] = []
                else:
                    name = line.strip()[1:]
                    namelist.append(name)
                    d[name] = []
            elif line.strip() != '':  # else: # trying to prevent some crap happening on the last line
                d[name].append(line.strip())
    for name in namelist:
        d[name] = ''.join(d[name])  # join list of partial sequences. Useful if interleaved fasta

    return d, namelist


def length_of_fasta_sequences(genomefile):
    import os

    print('Counting lengths of all scaffolds')
    path, f = os.path.split(genomefile)
    dgenome, names = read_fasta_as_dict(genomefile)
    d = {k: len(v) for k, v in dgenome.items()}

    return d, dgenome, names


def max_coords_per_scaffold(d):
    dmax_coords = {}
    dscaff_max = {}
    for n in list(d.keys()):
        # collect all coordinates in list, for each scaffold
        dmax_coords.setdefault('_'.join(n.split('_')[:-1]), []).append(int(n.split('_')[-1]))
    for scaff in list(dmax_coords.keys()):
        # get maximum coordinate per scaffold
        dscaff_max[scaff] = max(dmax_coords[scaff])
    return dscaff_max


def main(bamfiles, genomefile, srcbedfile, gff3file, step=25, winsize=100, cutoff=100, numplots=4):
    # input is sam file # output does not include positions of 0 coverage
    # key is scafffold_position # value is depth
    import os
    from natsort import natsorted, ns

    for bam in bamfiles:
        path, file = os.path.split(bam)

        prefix = os.path.join(path, '.'.join(file.split('.')[:-1]))
        print('step: %d, window size: %d, cutoff: %d' % (step, winsize, cutoff))

        dscafflengths, dgenome, names = length_of_fasta_sequences(genomefile)
        windows, positions = calculate_depth_with_pysam_sliding_windows(bam, dscafflengths, step, winsize, cutoff)

        sortwindows = natsorted(['\t'.join(w) for w in windows], key=lambda y: y.lower())
        write_to_bed(sortwindows, '%s.windows.s%d.w%d.d%d.pileup.bed' % (prefix, step, winsize, cutoff))
        collapsedsortwindows = collapse_windows(sortwindows)
        #collapsedwindowfilename = '%s.collapsedwindows.s%d.w%d.d%d.pileup.bed' % (prefix, step, winsize, cutoff)
        refinedsortwindows = refine_window_coordinates(collapsedsortwindows, bam)
        refinedwindowfilename = '%s.refinedwindows.s%d.w%d.d%d.pileup.bed' % (prefix, step, winsize, cutoff)
        write_to_bed(refinedsortwindows, refinedwindowfilename)

        dbed, bednames = read_bed(srcbedfile)
        dgff3, gff3names = read_gff3(gff3file)

        strandwindows = strandedness(bam, refinedwindowfilename, dbed, bednames, dgff3, gff3names, dgenome)
        strandedfilename = '%s.refinedwindows.stranded.s%d.w%d.d%d.pileup.tsv' % (prefix, step, winsize, cutoff)
        write_to_bed(strandwindows, strandedfilename)

        dclusters, clusternames = per_base_strandedness(bam, refinedwindowfilename)
        plot_per_base_strandedness(dclusters, clusternames, prefix, dbed, bednames, numplots)

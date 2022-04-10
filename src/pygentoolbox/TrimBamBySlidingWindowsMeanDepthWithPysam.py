

def write_out_sam(header, outsam, outpath):
    print('Writing trimmed sam to: %s' % outpath)
    with open(outpath, 'w') as OUT:
        OUT.write('\n'.join(header + outsam))


def trim_sam_by_windows(sam, positions):
    header = []
    output = []
    with open(sam, 'r') as FILE:
        for line in FILE:
            if line[0] == '@':
                header.append(line.strip())
            else:
                if '_'.join(line.strip().split('\t')[2:4]) in positions:
                    output.append(line.strip())
    return header, output


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
    overlapgate = 1  # record start position of first window in overlap
    for count, win in enumerate(sortwindows):
        scaff = win.split('\t')[0]
        winstart = int(win.split('\t')[1])
        winend = int(win.split('\t')[2])
        windepth = float(win.split('\t')[3])
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
        elif overlapgate == 0:
            # if it makes it to this point, there is no overlap with previous window
            # but if overlapgate == 0 then there was overlap between at least two windows and we combine coordinates
            collapsedsortwindows.append('%s\t%s\t%.2f' % (collapsewindow, prevwinend, sum(collapsewindepth) / len(collapsewindepth)))
            collapsewindepth = []  # reset for next series of overlapping windows
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
            else:  # if it is the last window here output the window
                collapsedsortwindows.append(win)
                collapsewindepth = []

        prevscaff = win.split('\t')[0]
        prevwinstart = int(win.split('\t')[1])
        prevwinend = int(win.split('\t')[2])
        prevwindepth = float(win.split('\t')[3])

    return collapsedsortwindows


def mean_with_denominator(lst, denominator):
    mean = sum(lst) / denominator
    return mean


def calculate_depth_with_pysam_sliding_windows(bamfile, dscaff_max, step=25, winsize=100, cutoff=100):
    import pysam

    windows = []
    positions = set()
    # d = {}  # i don't keep all read depth values
    samfile = pysam.AlignmentFile(bamfile, "rb")  #pysam.AlignmentFile("ex1.bam", "rb")
    for scaff, maxcoord in dscaff_max.items():
        print(scaff, end=', ')
        for i in range(0, maxcoord, step):
            # tempdepth = []  # temppos = []
            # .pileup skips bases when depth = 0  # need to figure out mpileup -a option in pysam
            tempdepth = [pileupcolumn.nsegments for pileupcolumn in samfile.pileup(scaff, i, i + winsize)]
            #for pileupcolumn in samfile.pileup(scaff, i, i + winsize):
            #    #print("\ncoverage at base %s = %s" % (pileupcolumn.pos, pileupcolumn.n))
            #    #temppos.append(pileupcolumn.reference_pos + 1)  # pos depricated
            #    tempdepth.append(pileupcolumn.nsegments)  # n depricated
            if len(tempdepth) != 0:
                average_depth = mean_with_denominator(tempdepth, denominator=winsize)
            else:
                average_depth = 0
            if average_depth >= cutoff:
                # print(str(average_depth), end=', ')
                # scaff, start, end, averagedepth
                windows.append([scaff, str(i + 1), str(i + winsize), str(average_depth)])
                for j in range(i + 1, i + winsize):
                    if j <= maxcoord:
                        positions.add('_'.join([scaff, str(j)]))
    samfile.close()

    return windows, positions


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

    return d, names


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


def main(bamfiles, genomefile, step=25, winsize=100, cutoff=100, maxdepth=0):
    # input is sam file # output does not include positions of 0 coverage
    # key is scafffold_position # value is depth
    import os
    from natsort import natsorted, ns

    for bam in bamfiles:
        path, file = os.path.split(bam)

        prefix = os.path.join(path, '.'.join(file.split('.')[:-1]))
        print('step: %d, window size: %d, cutoff: %d' % (step, winsize, cutoff))
        #d, names, nseqs = calculate_depth_of_coverage(sam, maxdepth)
        # dscaff_max = max_coords_per_scaffold(d)
        dscafflengths, names = length_of_fasta_sequences(genomefile)
        windows, positions = calculate_depth_with_pysam_sliding_windows(bam, dscafflengths, step, winsize, cutoff)
        #d = fill_in_zero_depth(d, dscaff_max)
        #windows, positions = sliding_windows(d, dscaff_max, interval, cutoff)

        sortwindows = natsorted(['\t'.join(w) for w in windows], key=lambda y: y.lower())
        write_to_bed(sortwindows, '%s.windows.s%d.w%d.d%d.pileup.bed' % (prefix, step, winsize, cutoff))
        collapsedsortwindows = collapse_windows(sortwindows)
        write_to_bed(collapsedsortwindows, '%s.collapsedwindows.s%d.w%d.d%d.pileup.bed' % (prefix, step, winsize, cutoff))
        #header, outsam = trim_sam_by_windows(sam, positions)
        #write_out_sam(header, outsam, '%s.windows.w%d.d%d.sam' % (prefix, interval, cutoff))

        # d_for, names_for, nseqs_for = calculate_depth_of_coverage(forward, maxdepth)
        # d_rev, names_rev, nseqs_rev = calculate_depth_of_coverage(reverse, maxdepth)






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


def mean(lst):
    return sum(lst) / len(lst)


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


def sliding_windows(d, dscaff_max, interval, cutoff):
    print('Starting sliding window analysis')
    windows = []
    positions = set()
    for scaff in list(dscaff_max.keys()):
        print(scaff)
        maxcoord = dscaff_max[scaff]
        # for each scaffold
        # iterate through each position (starting at 1) and find mean depth per window
        for i in range(1, maxcoord + 1):
            average_depth = mean([d['_'.join([scaff, str(j)])] for j in range(i, i + interval) if j <= maxcoord])
            if average_depth >= cutoff:
                windows.append([scaff, str(i), str(i + interval), str(average_depth)])
                for j in range(i, i + interval):
                    if j <= maxcoord:
                        positions.add('_'.join([scaff, str(j)]))
    print('Finished with sliding windows')
    return windows, positions


def fill_in_zero_depth(d, dscaff_max):
    for scaff in list(dscaff_max.keys()):
        maxcoord = dscaff_max[scaff]
        for i in range(1, maxcoord + 1):
            try:
                d['_'.join([scaff, str(i)])]
            except:
                # there is no position recorded on this scaffold with depth
                d['_'.join([scaff, str(i)])] = 0
            else:
                # there is a depth value for this position on this scaffold
                pass
    return d


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
                    position = int(line.strip().split('\t')[3]) + i  # calculate coordinate in scaffold/chromosome
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


def main(samfiles, interval=100, cutoff=1000, maxdepth=0):
    # input is sam file # output does not include positions of 0 coverage
    # key is scafffold_position # value is depth
    import os
    from natsort import natsorted, ns

    for sam in samfiles:
        path, file = os.path.split(sam)

        prefix = os.path.join(path, '.'.join(file.split('.')[:-1]))
        print('interval: %d, cutoff: %d' % (interval, cutoff))
        d, names, nseqs = calculate_depth_of_coverage(sam, maxdepth)
        dscaff_max = max_coords_per_scaffold(d)
        d = fill_in_zero_depth(d, dscaff_max)
        windows, positions = sliding_windows(d, dscaff_max, interval, cutoff)
        sortwindows = natsorted(['\t'.join(w) for w in windows], key=lambda y: y.lower())
        write_to_bed(sortwindows, '%s.windows.w%d.d%d.bed' % (prefix, interval, cutoff))
        collapsedsortwindows = collapse_windows(sortwindows)
        write_to_bed(collapsedsortwindows, '%s.collapsedwindows.w%d.d%d.bed' % (prefix, interval, cutoff))
        header, outsam = trim_sam_by_windows(sam, positions)
        write_out_sam(header, outsam, '%s.windows.w%d.d%d.sam' % (prefix, interval, cutoff))

        # d_for, names_for, nseqs_for = calculate_depth_of_coverage(forward, maxdepth)
        # d_rev, names_rev, nseqs_rev = calculate_depth_of_coverage(reverse, maxdepth)




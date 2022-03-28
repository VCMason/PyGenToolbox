def read_filter_contact_file_by_thickness(f, cutoff):
    # f is full path to circos contact file
    # cutoff is the lower limit of contact thickness allowed for it to be kept in file

    # example contact file lines, space delimited file
    # pt1 1591 1592 pt1 376344 376345 thickness=2.321928094887362p
    # pt1 1591 1592 pt511 2589 2590 thickness=1.0p
    # pt1 1705 1706 pt29 188715 188716 thickness=1.0p
    # pt1 2066 2067 pt96 177358 177359 thickness=2.0p
    import os

    outlist = []
    outarrowlist = []
    count = 0
    countout = 0
    with open(f, 'r') as FILE:
        for line in FILE:
            # line.strip().split()[-1][10:-1] == isolates 1.0 from: pt1 1591 1592 pt511 2589 2590 thickness=1.0p
            if float(line.strip().split()[-1][10:-1]) >= cutoff:
                outlist.append(line.strip())
                outarrowlist.append(' '.join(line.strip().split()[3:-1] + ['0']))
                countout += 1
            count += 1
    pathtodir, name = os.path.split(f)
    fout = '.'.join(['DNAArrows', 'scatter', 'Cutoff%.2f' % cutoff, 'txt'])
    fout = os.path.join(pathtodir, fout)
    with open(fout, 'w') as OUT:
        OUT.write('\n'.join(outarrowlist))

    fout = '.'.join(f.split('.')[:-1] + ['Cutoff%.2f' % cutoff, 'txt'])
    with open(fout, 'w') as OUT:
        OUT.write('\n'.join(outlist))

    print('Number of lines in input file: %s = %d' % (f, count))
    print('Number of lines in output file: %s = %d' % (fout, countout))


def main(f, cutoff):
    # f is full path to circos contact file
    # cutoff is the lower limit of contact thickness allowed for it to be kept in file

    # example contact file lines, space delimited file
    # pt1 1591 1592 pt1 376344 376345 thickness=2.321928094887362p
    # pt1 1591 1592 pt511 2589 2590 thickness=1.0p
    # pt1 1705 1706 pt29 188715 188716 thickness=1.0p
    # pt1 2066 2067 pt96 177358 177359 thickness=2.0p
    read_filter_contact_file_by_thickness(f, cutoff)
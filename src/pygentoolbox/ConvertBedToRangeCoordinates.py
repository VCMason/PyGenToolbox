def main(bedfile):
    with open(bedfile, 'r') as FILE:
        output = []
        for line in FILE:
            scaffold = line.strip().split('\t')[0]
            start = int(line.strip().split('\t')[1])
            end = int(line.strip().split('\t')[2])
            for i in range(start, end+1):  # need to add 1 to end to get final position recorded
                # each line of output will be: scaffold51_8_14740	0	0
                output.append('\t'.join(['%s_%d' % (scaffold, i), '0', '0']))

    outfile = bedfile + '.range'
    with open(outfile, 'w') as OUT:
        OUT.write('\n'.join(output))
    print('Fisnished converting Bed file to range of coordinates with 0 values for counts of forward and reverse reads')
    print('output file: %s' % outfile)


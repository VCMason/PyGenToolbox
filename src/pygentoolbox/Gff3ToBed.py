
def main(gff):
    bed = '.'.join(gff.split('.')[:-1] + ['bed'])
    print('Reading input: %s' % gff)
    outlist = []
    with open(gff, 'r') as FILE:
        for line in FILE:
            if line[0] == '#':
                pass
            else:
                outlist.append('\t'.join([line.strip().split('\t')[0], line.strip().split('\t')[3], line.strip().split('\t')[4], line.strip().split('\t')[8].split(';')[0][3:]]))

    with open(bed, 'w') as OUT:
        OUT.write('\n'.join(outlist))

    print('output file: %s' % bed)
    print('finished converting gff3 to bed file')
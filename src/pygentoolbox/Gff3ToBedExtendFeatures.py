
def main(gff, window):
    extension = int(window/2)
    bed = '.'.join(gff.split('.')[:-1] + ['%dwin' % window, 'bed'])
    print('Reading input: %s' % gff)
    outlist = []
    with open(gff, 'r') as FILE:
        for line in FILE:
            if line[0] == '#':
                pass
            else:
                extendstart = int(line.strip().split('\t')[3]) - extension
                if extendstart < 0:
                    extendstart = 0
                outlist.append('\t'.join([line.strip().split('\t')[0], str(extendstart), str(int(line.strip().split('\t')[4]) + extension), line.strip().split('\t')[8].split(';')[0][3:]]))

    with open(bed, 'w') as OUT:
        OUT.write('\n'.join(outlist))

    print('output file: %s' % bed)
    print('finished converting gff3 to bed file')
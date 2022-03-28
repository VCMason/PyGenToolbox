def filter_variants(variants):
    # variants is a list of lists
    keepvars = []
    for line in variants:
        # line is list
        if line[-3].split(';')[0] == 'INDEL':
            # process indel variants
            if (line[-1].split(':')[0] == '1/1') and (int(line[-3].split(';')[1][4:]) >= 5) and (int(line[-3].split(';')[-1][3:]) >= 30):
                keepvars.append(line)
        elif (line[-1].split(':')[0] == '1/1') and (int(line[-3].split(';')[0][3:]) >= 5) and (int(line[-3].split(';')[-1][3:]) >= 30):
            keepvars.append(line)

    return(keepvars)


def parse_vcf(f):
    # assumes file compressed as .gz
    # assumes all lines are variants (no reference matching genotypes)
    import gzip

    header = []
    with gzip.open(f, 'rt') as FILE:
        for line in FILE:
            if line[0] == '#':
                header.append(line.strip())
            else:
                break
    with gzip.open(f, 'rt') as FILE:
        variants = [line.strip().split('\t') for line in FILE if line[0] != '#']

    return header, variants


def main(vcffile):
    import gzip
    # vcffile is full path to vcf file. vcf file is gzipped .gz
    print('Working on file %s: ' % vcffile)
    print('This program requires variants to be: \n'
          'homozygous non-reference (1/1)\n'
          'overall depth of >= 5 (DP & IDV)\n'
          'MQ >= 30\n')
          #'DP4 non-reference ratio >= 0.90%. Ex: DP4=1,1,10,8)')

    # with gzip.open('file.gz', 'wt') as f:
    #    f.write('Hello world!')
    # with gzip.open('file.gz', 'wb') as f:
    #    f.write('Hello world!'.encode())
    header, variants = parse_vcf(vcffile)
    keepvars = filter_variants(variants)
    with gzip.open('.'.join(vcffile.split('.')[:-1] + ['filtered', 'gz']), 'wt') as OUT:
        OUT.write('\n'.join(header + ['\t'.join(l) for l in keepvars]))

    print('output file: %s' % '.'.join(vcffile.split('.')[:-1] + ['filtered', 'gz']))


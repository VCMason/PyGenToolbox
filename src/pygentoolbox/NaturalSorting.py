def main(f):
    import os
    from natsort import natsorted

    path, filename = os.path.split(f)
    outfile = os.path.join(path, filename + '.natsorted')

    with open(f, 'r') as FILE:
        alllines = [line.strip().split('\t') for line in FILE]

    # sort each column by natsorted(col1)
    output = natsorted(alllines)

    with open(outfile, 'w') as OUT:
        OUT.write('\n'.join(['\t'.join(l) for l in output]))

    print('Finished Sorting .tsv file by first column')
    print(f'output file: {outfile}')

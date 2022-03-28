def trim_gff3(f, list_of_scaffolds):
    import os
    outlist = []
    with open(f, 'r') as FILE:
        for line in FILE:
            if line.strip().split('\t')[0] in list_of_scaffolds:
                outlist.append(line.strip())
    path, file = os.path.split(f)
    outfile = os.path.join(path, f[:-len('.gff3')] + '.TrimScaffolds.gff3')
    with open(outfile, 'w') as OUT:
        OUT.write('\n'.join(outlist))
    print('Number of lines output = %d' % len(outlist))
    print('Output file: %s' % outfile)
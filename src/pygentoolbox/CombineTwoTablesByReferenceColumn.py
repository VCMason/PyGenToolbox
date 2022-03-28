# assumes files tab delimited with extension .tsv, all paths to file must be full paths
# f1 is the base file that will have f2 added to the end of each row with matching first column
# f1 names in reference column must be a subset of names in file f2
# f1 = 'D:\\LinuxShare\\Projects\\Theresa\\Limma\\EBayes\\mRNA_EVE12_vs_PtCAF1E12\\VoomEBayes.fit.DiffExpRes.TopTable.EVE12_vs_PtCAF1E12.tsv'
# f2 is the additional file that will be added to the end of rows of f1
# f2 = 'D:\\LinuxShare\\Ciliates\\Genomes\\Annotations\\Ptetraurelia51_GeneNames.GOTerms.AutogamyExpression.tsv'
# outfile is the output filename


def combine_two_tables_by_first_column_id(f1, f2, refcol1, refcol2):
    with open(f2, 'r') as F2:
        header2 = F2.readline().strip()
        print('Number of columns in file2: %d' % len(header2.split('\t')))
        d = {line.strip().split('\t')[refcol2]: line.rstrip().split('\t') for line in F2}

    with open(f1, 'r') as F1:
        outlist = []
        header1 = F1.readline().strip()
        for line in F1:
            try:
                d[line.strip().split('\t')[refcol1]]
            except:
                outlist.append('\t'.join(line.strip().split() + ['NA']*len(header2.split('\t'))))
            else:
                outlist.append('\t'.join(line.strip().split() + d[line.rstrip().split('\t')[refcol1]]))
        outlines = '\n'.join(outlist)
        output = header1 + '\t' + header2 + '\n' + outlines

    outfile = '.'.join(f1.split('.')[:-1] + ['Combined.tsv'])
    with open(outfile, 'w') as OUT:
        OUT.write(output)


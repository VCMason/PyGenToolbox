
def output_file(outfile, d, names):
    output = [d[n] for n in names]
    with open(outfile, 'w') as OUT:
        OUT.write('\n'.join(output))
    print(f'Output file: {outfile}')


def trim_dictionary(dgff3, dies, gff3ids, iesids, cutoff):
    # trims dgff3 by values in dies
    dtrimup = {}  # {k:v for k, v in d.items() if }
    dtrimdown = {}
    upnames = []
    downnames = []
    for gff3id in gff3ids:
        if gff3id in iesids:
            if dies[gff3id] >= cutoff:
                dtrimup[gff3id] = dgff3[gff3id]
                upnames.append(gff3id)
            elif dies[gff3id] < cutoff:
                dtrimdown[gff3id] = dgff3[gff3id]
                downnames.append(gff3id)
    return dtrimup, dtrimdown, upnames, downnames


def read_table(tablefile, n):
    # assumes ID column of tablefile is 0
    # n is column number
    # requires header row
    print('Reading table file: %s' % tablefile)
    d = {}
    names = []
    with open(tablefile, 'r') as FILE:
        for count, line in enumerate(FILE):
            if count == 0:
                head = line.strip().split('\t')[n]
            else:
                key = line.strip().split('\t')[0]
                value = line.strip().split('\t')[n]
                if value != 'NA':
                    d[key] = float(value)
                    names.append(key)
    print(f'Number of rows in table: {len(names)}')
    print('NA column values excluded')
    return d, names, head


def read_gff3(GFF3file, features=['all']):
    print('Reading Gff3 file: %s' % GFF3file)
    d = {}
    names = []
    with open(GFF3file, 'r') as FILE:
        for line in FILE:
            if (line[0] == '#') or (line.strip() == ''):
                pass
            elif features == ['all']:  # keep all lines
                # key = line.strip().split('\t')[0] + ':' +  line.strip().split('\t')[3] + '_' + line.strip().split('\t')[4]
                key = line.strip().split('\t')[8].split(';')[0][len('ID='):]
                d[key] = line.strip()  # .split('\t')
                names.append(key)
            elif line.strip().split('\t')[2] in features:  # keep lines only if in features
                # key = line.strip().split('\t')[0] + ':' + line.strip().split('\t')[3] + '_' + line.strip().split('\t')[4]
                key = line.strip().split('\t')[8].split(';')[0][len('ID='):]
                d[key] = line.strip()  # .split('\t')
                names.append(key)
    print('Number of features: %d' % len(list(d.keys())))
    # print(list(d.keys())[:10])
    return d, names


def main(gff3, table, column=2, cutoff=0.2):
    # gff3 is full path to annotation file
    # table is full path to tab delimited table file, columns are float values
    # column is column number (0-based), assumes ID ist first column i.e. 0

    import os

    dgff3, gff3ids = read_gff3(gff3, features=['all'])
    dies, iesids, head = read_table(table, column)
    dgff3trimup, dgff3trimdown, upnames, downnames = trim_dictionary(dgff3, dies, gff3ids, iesids, cutoff)

    path, f = os.path.split(gff3)

    out = '.'.join(f.split('.')[:-1] + ['trim', head, str(cutoff), 'up', 'gff3'])
    outfile = os.path.join(path, out)
    output_file(outfile, dgff3trimup, upnames)

    out = '.'.join(f.split('.')[:-1] + ['trim', head, str(cutoff), 'down', 'gff3'])
    outfile = os.path.join(path, out)
    output_file(outfile, dgff3trimdown, downnames)

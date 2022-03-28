
def write_out(outlines, outfile):
    with open(outfile, 'w') as OUT:
        OUT.write('\n'.join(outlines))
    print(f'Output file written to: {outfile}')


def create_intron_annotations(dexons, ids):
    # input is dictionary, key is transctipt ID
    # value is list of exon lines from gff3 file for mRNA with two or more exons
    # exon should be written as CDS
    introns = []
    exonsandintrons = []
    for id in ids:
        for count, exon in enumerate(dexons[id]):
            if count == 0:
                intronstart = int(exon[4]) + 1  # start of introns is one base passed the end of previous exon
                exonsandintrons.append('\t'.join(exon))
            elif count >= 1:
                intronend = int(exon[3]) - 1  # end of intron is one base before the start of next exon
                metadata = exon[8].split(';')
                #intronid = 'ID=' + '.'.join(id.split('.')[:-1] + ['I' + id.split('.')[-1][1:]])
                for data in metadata:
                    if data[:len('ID=')] == 'ID=':
                        intronid = '.'.join(data.split('.')[:3] + ['I' + data.split('.')[3][1:].split(':')[0] + f':{intronstart}..{intronend}'])  # data.split('.')[3][1:]] + data.split('.')[4:]
                    if data[:len('Name=')] == 'Name=':
                        intronname = '.'.join(data.split('.')[:3] + ['I' + data.split('.')[3][1:].split(':')[0] + f':{intronstart}..{intronend}'])  # data.split('.')[4:]
                    if data[:len('Parent=')] == 'Parent=':
                        parent = data
                intron = exon[:2] + ['intron', str(intronstart), str(intronend)] + exon[5:8] + [';'.join([intronid, intronname, parent])]
                introns.append('\t'.join(intron))
                exonsandintrons.append('\t'.join(intron))
                exonsandintrons.append('\t'.join(exon))
                intronstart = int(exon[4]) + 1  # start of introns is one base passed the end of previous exon

    return introns, exonsandintrons


def read_cds_gff3(gff3file):
    print('Reading Gff3 file: %s' % gff3file)
    d = {}
    mrnaids = []  # list of transcript IDs of mRNA features
    with open(gff3file, 'r') as FILE:
        for line in FILE:
            if (line[0] == '#') or (line.strip() == ''):
                pass
            elif line.strip().split('\t')[2] == 'mRNA':  # keep all lines
                mrnaids.append(line.strip().split('\t')[8].split(';')[0][len('ID='):])
            elif line.strip().split('\t')[2] == 'CDS':  # keep lines only if in features # can use exon
                for metadata in line.strip().split('\t')[8].split(';'):
                    if metadata[:len('Parent=')] == 'Parent=':
                        key = metadata[len('Parent='):]  # isolate parent transcript ID
                    else:
                        pass
                # key is parent transcript ID, value is list of lists where internal list is gff3line for exon
                d.setdefault(key, []).append(line.strip().split('\t'))
    dtrim = {id: d[id] for id in mrnaids}
    print('Number of mRNA features: %d' % len(list(dtrim.keys())))
    dtrim2 = {id: d[id] for id in mrnaids if len(d[id]) >= 2}
    mrnaidstrim = [id for id in mrnaids if id in dtrim2.keys()]
    print('Number of mRNA features with 2 or more exons: %d' % len(list(dtrim2.keys())))
    # returns exons of mRNA features that have two or more exons
    return dtrim2, mrnaidstrim


def main(gff3file):
    # dexons actually isolates CDS features...
    dexons, ids = read_cds_gff3(gff3file)
    introns, exonsandintrons = create_intron_annotations(dexons, ids)

    outfile = '.'.join(gff3file.split('.')[:-1] + ['introns', 'gff3'])
    write_out(introns, outfile)
    outfile = '.'.join(gff3file.split('.')[:-1] + ['introns', 'CDS', 'gff3'])
    write_out(exonsandintrons, outfile)
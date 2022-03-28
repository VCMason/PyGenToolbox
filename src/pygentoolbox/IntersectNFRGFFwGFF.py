
def write_out(outpath, output):
    print('Writing out file: %s' % outpath)
    with open(outpath, 'w') as OUT:
        OUT.write('\n'.join(output))


def sort_dictionary_by_value(d):
    ''' returns two sorted lists '''
    tkeys = list(d.keys())
    tvalues = [d[n] for n in tkeys]
    # *sorted would be increasing
    sortvalues, sortkeys = (list(t) for t in zip(*sorted(zip(tvalues, tkeys), reverse=True)))
    return sortvalues, sortkeys


def read_sam(sam, positions=[]):
    ''' key is scaffold, value is line.strip() '''
    print('Reading sam file: %s' % sam)
    d = {}
    dhead = {}
    countreads = 0
    with open(sam, 'r') as FILE:
        for line in FILE:
            if line[0] != '@':
                if int(line.strip().split('\t')[3]) not in positions:
                    key = line.strip().split('\t')[2]
                    d.setdefault(key, []).append(line.strip().split('\t'))
                    countreads += 1
            elif line[:3] == '@SQ':
                dhead[line.strip().split('\t')[1][3:]] = int(line.strip().split('\t')[2].split(':')[1])
    print('Number of scaffolds (in header): %d' % len(list(dhead.keys())))
    print('Total number of lines in sam file (excluding header): %d' % countreads)
    return d, dhead, countreads


def read_gff3(GFF3file, features=['all'], silent=False):
    if silent == False:
        print('Reading Gff3 file: %s' % GFF3file)
    d = {}
    dhead = {}
    with open(GFF3file, 'r') as FILE:
        for line in FILE:
            if line[:5] == '##gff':
                pass
            elif line[:17] == '##sequence-region':
                dhead[line.strip().split()[1]] = int(line.strip().split()[3])
            elif line[0] == '#':
                pass
            elif features == ['all']:  # keep all lines
                # key is scaffold name (zeroth column)
                key = line.strip().split('\t')[0]
                d.setdefault(key, []).append(line.strip().split('\t'))
            elif line.strip().split('\t')[2] in features:  # keep lines only if in features
                # key is scaffold name (zeroth column)
                key = line.strip().split('\t')[0]
                d.setdefault(key, []).append(line.strip().split('\t'))
    if silent == False:
        print('Number of scaffolds: %d' % len(list(d.keys())))  # dhead.keys()
    # print(list(d.values())[:2])
    return d, dhead


def mean(lst):
    return sum(lst) / len(lst)


def intersect(dgff3, dNFR, t1='gff', t2='gff'):

    if t1 == 'gff' and t2 == 'gff':
        dfeatures = {}  # number of reads that overlap each feature(key == ID) (could even be one base)
        names = []  # every ID= for each feature with reads overlapping
        dfeaturetypes = {}  # number of reads that overlap each feature type. gene, mRNA, exon, CDS
        intersectfeatures = []
        # overlap represent the percent of each sequence covered by each feature
        overlaptype = {}  # for each feature type, record percent of overlap of each aligned seq to each feature
        completeoverlap = []
        dcompleteoverlapcounts = {}
        dbinary = {}
        dfeaturesrightandleft = {}  # number of reads that overlap each feature(key == ID) (could even be one base)

        # allfeaturelengths = []  # not used
        for scaffold in list(dgff3.keys()):
            for gff3line in dgff3[scaffold]:
                start_gff = int(gff3line[3])
                end_gff = int(gff3line[4])
                try:
                    dNFR[scaffold]
                except:
                    pass  # nothing aligned to scaffold in dNFR
                else:
                    for feature in dNFR[scaffold]:
                        ncRNAtype = ''  # will be equal to tRNA, rRNA, or snoRNA if present in ncRNA line
                        start_NFR = int(feature[3])
                        end_NFR = int(feature[4])
                        # if we assume that  ranges are well-formed (so that x1 <= x2 and y1 <= y2) then
                        # x1 <= y2 && y1 <= x2
                        # meaning the aligned sequence overlaps in some way to the gff3 feature
                        if (start_gff <= end_NFR) and (start_NFR <= end_gff):
                            # gff3ID = gff3line[8].split(';')[0]
                            names.append(gff3line[8].split(';')[0][3:])  # [3:] would cut out 'ID=...'
                            dfeatures[gff3line[8].split(';')[0][3:]] = dfeatures.setdefault(gff3line[8].split(';')[0][3:], 0) + 1
                            dfeaturetypes[gff3line[2]] = dfeaturetypes.setdefault(gff3line[2], 0) + 1
                            intersectfeatures.append('\t'.join(gff3line + feature))
                            if (start_NFR <= start_gff) and (end_NFR >= end_gff):
                                # if read aligned is larger than feature aligned to...
                                #    although read completely covers feature, the feature does not completely cover read
                                completeoverlap.append('\t'.join(gff3line + feature))
                                dcompleteoverlapcounts[gff3line[8].split(';')[0][3:]] = dcompleteoverlapcounts.setdefault(gff3line[8].split(';')[0][3:], 0) + 1
                                dbinary[gff3line[8].split(';')[0][3:]] = dbinary.setdefault(gff3line[8].split(';')[0][3:], 0) + 1
                            else:
                                pass
                # dbinary[gff3line[8].split(';')[0][3:]] = dbinary.setdefault(gff3line[8].split(';')[0][3:], 0) + 0

    return dfeatures, dfeaturetypes, names, intersectfeatures, completeoverlap, dcompleteoverlapcounts, dbinary


def main(GFF3file, NFRgfflist):
    import os

    # dTE, dTEhead = read_sam(TEfile, positions=[0])  # sam file of aligned TEs
    # keys of dgff3 and dsRNA is scaffold name
    dgff3, dgff3head = read_gff3(GFF3file)  # , features=['mRNA']
    allnames = []
    dallfeatures = {}
    dallfeaturetypes = {}
    allintersectfeatures = []
    allcompleteoverlap = []
    dallcompleteoverlapcounts = {}
    # dall binary will only contain IESs with value 1 when NFR feature COMPLETELY overlaps Gff3 feature
    dallbinary = {}
    # dbinarybackground makes a dictionary of all IES ID's with value == 0
    # then update this value with
    dbinarybackground = {gff3line[8].split(';')[0][3:]: 0 for scaffold, gff3lines in dgff3.items() for gff3line in gff3lines}
    allbinarynames = [gff3line[8].split(';')[0][3:] for scaffold, gff3lines in dgff3.items() for gff3line in gff3lines]

    for count, NFRgff in enumerate(NFRgfflist, start=1):  # for NFR.gff in NFR.gff
        # dsRNA, dsRNAhead, countreads = read_sam(sam, positions=[0])  # dsRNA stands for dictionary small RNA
        dNFR, dNFRhead = read_gff3(NFRgff, ['all'], silent=True)

        dfeatures, dfeaturetypes, names, intersectfeatures, completeoverlap, dcompleteoverlapcounts, dbinary = intersect(dgff3, dNFR, 'gff', 'gff')

        allnames = allnames + names
        dallfeatures.update(dfeatures)
        dallfeaturetypes.update(dfeaturetypes)
        allintersectfeatures = allintersectfeatures + intersectfeatures
        allcompleteoverlap = allcompleteoverlap + completeoverlap
        dallcompleteoverlapcounts.update(dcompleteoverlapcounts)
        dallbinary.update(dbinary)

        # print(allintersectfeatures[-3:])
        print(count, end=', ')

    dbinarybackground.update(dallbinary)

    path, file = os.path.split(NFRgfflist[0])

    print('Number of NFR features overlapping Gff3 features: %d' % len(allintersectfeatures))
    output = allintersectfeatures
    outpath = os.path.join(path, '.'.join(file.split('.')[:-1] + ['Gff3', 'NFR', 'features', 'tsv']))
    write_out(outpath, output)

    print('Number of NFR features Completely overlapping Gff3 features: %d' % len(allintersectfeatures))
    output = allcompleteoverlap
    outpath = os.path.join(path, '.'.join(file.split('.')[:-1] + ['Gff3', 'NFR', 'features', 'CompleteOverlap', 'tsv']))
    write_out(outpath, output)

    # dbinary should only be 1 if a NFR feature completely overlaps IES
    # output = ['\t'.join(n.split('_') + [str(dallfeatures[n])]) for n in allnames]
    output = ['\t'.join([n, str(dbinarybackground[n])]) for n in allbinarynames]
    outpath = os.path.join(path, '.'.join(file.split('.')[:-1] + ['allfeatures', 'CompleteOverlap', 'binary', 'tsv']))
    write_out(outpath, output)

    completenames = [elem.split('\t')[8].split(';')[0][3:] for elem in allcompleteoverlap]
    # output = ['\t'.join(n.split('_') + [str(dallfeatures[n])]) for n in completenames]
    output = ['\t'.join([n, str(dallcompleteoverlapcounts[n])]) for n in completenames]
    outpath = os.path.join(path, '.'.join(file.split('.')[:-1] + ['features', 'CompleteOverlap', 'tsv']))
    write_out(outpath, output)

    # output = ['\t'.join(n.split('_') + [str(dallfeatures[n])]) for n in allnames]
    output = ['\t'.join([n, str(dallfeatures[n])]) for n in allnames]
    outpath = os.path.join(path, '.'.join(file.split('.')[:-1] + ['features', 'tsv']))
    write_out(outpath, output)



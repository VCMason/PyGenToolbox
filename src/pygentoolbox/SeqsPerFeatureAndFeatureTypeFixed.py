
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


def get_rightmost_reference_based_alignment_coordinate(CIGAR, leftmost_coordinate):
    import re
    cigar = re.findall(r'\d+[A-Z]', CIGAR)
    if cigar == []:  # if there was a match # sometimes CIGAR string == * # this skips unmapped reads
        print(f'Provided CIGAR string: {CIGAR} does not match CIGAR pattern \\d+[A-Z]')
        rightmost_position = 0  # assumes unmapped read
    else:  # then read should be mapped
        rightmost_position = leftmost_coordinate - 1  # subtract 1 because leftmost base is 1-based
        for i in cigar:
            if i[-1] in ['M', 'N', 'D', 'X', '=']:
                rightmost_position += int(i[:-1])
            elif i[-1] in ['I', 'S', 'H', 'P']:
                pass
            else:
                pass

    return rightmost_position


def count_reads_per_feature(dgff3, dsam, featurebuffer=0):  #  target_feature='gene'
    print('Counting reads for all feature and Count reads per feature type (ex: exon, gene, CDS, etc.)')
    print(f'Modifying start and end feature coordinates by {featurebuffer}')
    dallfeatures = {}  # number of reads that overlap each feature(key == ID) (could even be one base)
    allnames = []  # every ID= for each feature with reads overlapping
    dfeaturetypes = {}  # number of reads that overlap each feature type. gene, mRNA, exon, CDS
    # overlap represent the percent of each sequence covered by each feature
    overlaptypes = {}  # for each feature type, record percent of overlap of each aligned seq to each feature
    # allfeaturelengths = []  # not used
    for scaffold in list(dgff3.keys()):
        for gff3line in dgff3[scaffold]:
            start_gff = int(gff3line[3]) + featurebuffer
            end_gff = int(gff3line[4]) - featurebuffer
            try:
                dsam[scaffold]
            except:
                pass  # nothing aligned to scaffold in dsam
            else:
                for aligned_seq in dsam[scaffold]:
                    ncRNAtype = ''  # will be equal to tRNA, rRNA, or snoRNA if present in ncRNA line
                    start_sam = int(aligned_seq[3])
                    end_sam = int(aligned_seq[3]) + len(aligned_seq[9]) - 1  # get_rightmost_reference_based_alignment_coordinate(aligned_seq[5], start_sam) # CIGAR = aligned_seq[5] # len(aligned_seq[9])
                    # if we assume that  ranges are well-formed (so that x1 <= x2 and y1 <= y2) then
                    # x1 <= y2 && y1 <= x2
                    # meaning the aligned sequence overlaps in some way to the gff3 feature
                    if (start_gff <= end_sam) and (start_sam <= end_gff):
                        # gff3ID = gff3line[8].split(';')[0]
                        allnames.append(gff3line[8].split(';')[0])  # [3:] would cut out 'ID=...'
                        dallfeatures[gff3line[8].split(';')[0]] = dallfeatures.setdefault(gff3line[8].split(';')[0], 0) + 1
                        dfeaturetypes[gff3line[2]] = dfeaturetypes.setdefault(gff3line[2], 0) + 1
                        if gff3line[2] == 'ncRNA':
                            # then iterate through metadata and check for type of ncRNA: rRNA, tRNA, snoRNA
                            for info in gff3line[8].split(';'):
                                if info[:5] == 'type=':
                                    ncRNAtype = info[5:]
                                    dfeaturetypes[ncRNAtype] = dfeaturetypes.setdefault(ncRNAtype, 0) + 1
                        # if gff3line[1] == '+':  # feature is 5' -> 3'
                        # elif gff3line[1] == '-':  # feature is 3' -> 5'
                        if (end_sam - start_sam) == 0:  # if (end_gff - start_gff) == 0:  # add one?
                            # avoid division by zero
                            print('Potential division by zero: %s, %s, %d, %d' % (scaffold, gff3line[2], start_gff, end_gff))
                        elif (start_sam <= start_gff) and (end_sam <= end_gff):
                            # 'right' side of aligned read overlaps feature but 'left' does not
                            #  if gff3line[2] == target_feature:
                            overlaptypes.setdefault(gff3line[2], []).append(((end_sam - start_gff) / (end_sam - start_sam)) * 100)
                            if ncRNAtype != '':
                                overlaptypes.setdefault(ncRNAtype, []).append(((end_sam - start_gff) / (end_sam - start_sam)) * 100)
                        elif (start_sam >= start_gff) and (end_sam >= end_gff):
                            # 'left' side of aligned read overlaps feature but 'right' does not
                            overlaptypes.setdefault(gff3line[2], []).append(((end_gff - start_sam) / (end_sam - start_sam)) * 100)
                            if ncRNAtype != '':
                                overlaptypes.setdefault(ncRNAtype, []).append(((end_gff - start_sam) / (end_sam - start_sam)) * 100)
                        elif (start_gff <= start_sam) and (end_sam <= end_gff):
                            # if read aligned is inside feature and smaller than feature
                            # (end_sam - start_sam) / (end_sam - start_sam)
                            overlaptypes.setdefault(gff3line[2], []).append(100)
                            if ncRNAtype != '':
                                overlaptypes.setdefault(ncRNAtype, []).append(100)
                        elif (start_sam <= start_gff) and (end_sam >= end_gff):
                            # if read aligned is larger than feature aligned to...
                            #     although read completely covers feature, the feature does not completely cover read
                            overlaptypes.setdefault(gff3line[2], []).append(100)
                            if ncRNAtype != '':
                                overlaptypes.setdefault(ncRNAtype, []).append(100)
                        else:
                            print('some form of overlap was not accounted for')
                            print('%s, SS: %d, SG: %d, ES: %d, EG: %d' % (scaffold, start_sam, start_gff, end_sam, end_gff))

        print(scaffold, end=', ')

    return dallfeatures, dfeaturetypes, allnames, overlaptypes


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


def read_gff3(GFF3file, features=['all']):
    print('Reading Gff3 file: %s' % GFF3file)
    d = {}
    dhead = {}
    dfeaturetypelengths = {}
    with open(GFF3file, 'r') as FILE:
        for line in FILE:
            if line[:5] == '##gff':
                pass
            elif line[:17] == '##sequence-region':
                dhead[line.strip().split()[1]] = int(line.strip().split()[3])
            elif line[0] == '#':
                pass
            elif line[2] == '*':
                pass  # skips reads unaligned to any scaffold
            elif features == ['all']:  # keep all lines
                # key is scaffold name (zeroth column)
                key = line.strip().split('\t')[0]
                d.setdefault(key, []).append(line.strip().split('\t'))
                start = int(line.strip().split('\t')[3])
                end = int(line.strip().split('\t')[4])
                dfeaturetypelengths[line.strip().split('\t')[2]] = dfeaturetypelengths.setdefault(line.strip().split('\t')[2], 0) + (end - start + 1)
                if line.strip().split('\t')[2] == 'ncRNA':
                    # then iterate through metadata and check for type of ncRNA: rRNA, tRNA, snoRNA
                    for info in line.strip().split('\t')[8].split(';'):
                        if info[:5] == 'type=':
                            ncRNAtype = info[5:]
                            dfeaturetypelengths[ncRNAtype] = dfeaturetypelengths.setdefault(ncRNAtype, 0) + (end - start + 1)
            elif line.strip().split('\t')[2] in features:  # keep lines only if in features
                # key is scaffold name (zeroth column)
                key = line.strip().split('\t')[0]
                d.setdefault(key, []).append(line.strip().split('\t'))
                start = int(line.strip().split('\t')[3])
                end = int(line.strip().split('\t')[4])
                dfeaturetypelengths[line.strip().split('\t')[2]] = dfeaturetypelengths.setdefault(line.strip().split('\t')[2], 0) + (end - start + 1)
                if line.strip().split('\t')[2] == 'ncRNA':
                    # then iterate through metadata and check for type of ncRNA: rRNA, tRNA, snoRNA
                    for info in line.strip().split('\t')[8].split(';'):
                        if info[:5] == 'type=':
                            ncRNAtype = info[5:]
                            dfeaturetypelengths[ncRNAtype] = dfeaturetypelengths.setdefault(ncRNAtype, 0) + (end - start + 1)
    print('Number of scaffolds: %d' % len(list(d.keys())))  # dhead.keys()
    # print(list(d.values())[:2])
    return d, dhead, dfeaturetypelengths


def mean(lst):
    return sum(lst) / len(lst)


def main(GFF3file, samfiles, featurebuffer=0):
    import os

    # featurebuffer reduces the size of features, changes feature borders (start+featurebuffer, end-featurebuffer)
    # dTE, dTEhead = read_sam(TEfile, positions=[0])  # sam file of aligned TEs
    # keys of dgff3 and dsRNA is scaffold name
    dgff3, dgff3head, dfeaturetypelengths = read_gff3(GFF3file)  # , features=['mRNA']
    # print(len(dmRNA['scaffold51_74']))
    # print(dmRNA['scaffold51_74'])

    for sam in samfiles:
        # sam files of aligned sRNAs
        dsRNA, dsRNAhead, countreads = read_sam(sam, positions=[0])  # dsRNA stands for dictionary small RNA
        if countreads == 0:
            # pass because there are no mapped reads, only reads available would be unaligned reads
            pass
        else:
            print('Mapped reads, so calculating seqs per feature and feature type')
            dmax = {**dsRNAhead, **dgff3}  # update dsRNAhead with values from dTEhead  # changed to dgff3
            dallfeatures, dfeaturetypes, allnames, overlaptypes = count_reads_per_feature(dgff3, dsRNA, featurebuffer)
            path, file = os.path.split(sam)

            # output = ['\t'.join(n.split('_') + [str(dallfeatures[n])]) for n in allnames]
            output = ['\t'.join([n, str(dallfeatures[n])]) for n in allnames]
            outpath = os.path.join(path, '.'.join(file.split('.')[:-1] + ['allfeatures', 'fixed', 'tsv']))
            write_out(outpath, output)

            sortvalues, sortkeys = sort_dictionary_by_value(dallfeatures)
            header = 'Gene\tRead Count\tPercent of Total Reads\tRPKM\tTotal Reads'
            # output = [header] + ['\t'.join([k, str(v), str((v / countreads) * 100), str(countreads)]) for k, v in zip(sortkeys, sortvalues)]
            # gene length = (int(k.split(':')[-1].split('..')[1]) - int(k.split(':')[-1].split('..')[0])) + 1
            output = [header] + ['\t'.join([k, str(v), str((v / countreads) * 100),
                                            str( (v / (((int(k.split(':')[-1].split('..')[1]) - int(k.split(':')[-1].split('..')[0]) + 1) / 1000) * (countreads / 1000000))) ),
                                            str(countreads)]) for k, v in zip(sortkeys, sortvalues)]
            outpath = os.path.join(path, '.'.join(file.split('.')[:-1] + ['allfeatures', 'sorted', 'fixed', 'tsv']))
            write_out(outpath, output)

            sortvalues, sortkeys = sort_dictionary_by_value(dfeaturetypes)
            header = 'Feature Type\tRead Count\tPercent of Total Reads\tRPKM\tTotal Reads\tAverage Percent Length of Overlap of Reads Overlapping Feature'
            # output = [header] + ['\t'.join([k, str(v), str((v / countreads) * 100), str(countreads), str(mean(overlaptypes[k]))]) for k, v in zip(sortkeys, sortvalues)]
            output = [header] + ['\t'.join([k, str(v), str((v / countreads) * 100), str( (v / ((dfeaturetypelengths[k] / 1000) * (countreads / 1000000))) ), str(countreads), str(mean(overlaptypes[k]))]) for k, v in zip(sortkeys, sortvalues)]
            outpath = os.path.join(path, '.'.join(file.split('.')[:-1] + ['featuretypes', 'fixed', 'tsv']))
            write_out(outpath, output)
            #print(output)


def write_out(outpath, output):
    print('Writing out file: %s' % outpath)
    with open(outpath, 'w') as OUT:
        OUT.write('\n'.join(output))


def rpkm_dictionary(d, dlengths):
    # How to calculate RPKM:

    # Count up the total reads in a sample and divide that number by 1,000,000 – this is our “per million” scaling factor.
    # Divide the read counts by the “per million” scaling factor. This normalizes for sequencing depth, giving you reads per million (RPM)
    # Divide the RPM values by the length of the gene, in kilobases. This gives you RPKM.

    totalreads = sum(d.values())
    scalefactor = totalreads / 1000000
    drpm = {k: v/scalefactor for k, v in d.items()}
    drpkm = {k: v/dlengths[k] for k, v in drpm.items()}

    return drpkm

def calculate_gene_length_gff3(d):
    ''' d is dictionary of gff3 file, scaffold as key, all entries to scaffold as values in list '''
    # Ex line= scaffold51_1	EuGene	gene	2175	2778	.	-	.	ID=PTET.51.1.G0010001;Name=PTET.51.1.G0010001
    dgenelengths = {}
    for k, v in d.items():
        for splitline in v:
            dgenelengths[splitline[8].split(';')[0][3:]] = int(splitline[4]) - int(splitline[3]) + 1
    return dgenelengths


def count_reads_per_GOterm_using_feature_counts(dGO, dallfeatures):
    dGOcounts = {}
    dGOdefineCounts = {}
    totalGOcounts = 0
    for gene in dGO.keys():
        try:
            dallfeatures[gene]
        except:
            pass
        else:
            if dGO[gene][8] != '':
                for GO in dGO[gene][8].split('; '):
                    # add the number of reads mapping to a gene feature to each GO term associated with that gene
                    dGOcounts[GO] = dGOcounts.setdefault(GO, 0) + dallfeatures[gene]
                for GO in dGO[gene][9].split('; '):
                    dGOdefineCounts[GO] = dGOdefineCounts.setdefault(GO, 0) + dallfeatures[gene]
                totalGOcounts += dallfeatures[gene]
    return dGOcounts, dGOdefineCounts, totalGOcounts


def count_reads_per_feature(dgff3, dsam):  #  target_feature='gene'
    print('Counting reads for all feature and Count reads per feature type (ex: exon, gene, CDS, etc.)')
    dallfeatures = {}  # number of reads that overlap each feature(key == ID) (could even be one base)
    allnames = []  # every ID= for each feature with reads overlapping
    dfeaturetypes = {}  # number of reads that overlap each feature type. gene, mRNA, exon, CDS
    # overlap represent the percent of each sequence covered by each feature
    overlaptypes = {}  # for each feature type, record percent of overlap of each aligned seq to each feature
    allfeaturelengths = []
    for scaffold in list(dgff3.keys()):
        for gff3line in dgff3[scaffold]:
            start_gff = int(gff3line[3])
            end_gff = int(gff3line[4])
            try:
                dsam[scaffold]
            except:
                pass  # nothing aligned to scaffold in dsam
            else:
                for aligned_seq in dsam[scaffold]:
                    start_sam = int(aligned_seq[3])
                    end_sam = int(aligned_seq[3]) + len(aligned_seq[9])
                    # if we assume that  ranges are well-formed (so that x1 <= x2 and y1 <= y2) then
                    # x1 <= y2 && y1 <= x2
                    # meaning the aligned sequence overlaps in some way to the gff3 feature
                    if (start_gff <= end_sam) and (start_sam <= end_gff):
                        # gff3ID = gff3line[8].split(';')[0]
                        allnames.append(gff3line[8].split(';')[0][3:])  # [3:] would cut out 'ID=...'
                        dallfeatures[gff3line[8].split(';')[0][3:]] = dallfeatures.setdefault(gff3line[8].split(';')[0][3:], 0) + 1
                        dfeaturetypes[gff3line[2]] = dfeaturetypes.setdefault(gff3line[2], 0) + 1
                        # if gff3line[1] == '+':  # feature is 5' -> 3'
                        # elif gff3line[1] == '-':  # feature is 3' -> 5'
                        if (end_sam - start_sam) == 0:  # if (end_gff - start_gff) == 0:  # add one?
                            # avoid division by zero
                            print('Potential division by zero: %s, %s, %d, %d' % (scaffold, gff3line[2], start_gff, end_gff))
                        elif (start_sam <= start_gff) and (end_sam <= end_gff):
                            # 'right' side of aligned read overlaps feature but 'left' does not
                            #  if gff3line[2] == target_feature:
                            overlaptypes.setdefault(gff3line[2], []).append(((end_sam - start_gff) / (end_sam - start_sam)) * 100)
                        elif (start_sam >= start_gff) and (end_sam >= end_gff):
                            # 'left' side of aligned read overlaps feature but 'right' does not
                            overlaptypes.setdefault(gff3line[2], []).append(((end_gff - start_sam) / (end_sam - start_sam)) * 100)
                        elif (start_gff <= start_sam) and (end_sam <= end_gff):
                            # if read aligned is inside feature and smaller than feature
                            # (end_sam - start_sam) / (end_sam - start_sam)
                            overlaptypes.setdefault(gff3line[2], []).append(100)
                        elif (start_sam <= start_gff) and (end_sam >= end_gff):
                            # if read aligned is larger than feature aligned to...
                            #     although read completely covers feature, the feature does not completely cover read
                            overlaptypes.setdefault(gff3line[2], []).append(100)
                        else:
                            print('some form of overlap was not accounted for')
                            print('%s, SS: %d, SG: %d, ES: %d, EG: %d' % (scaffold, start_sam, start_gff, end_sam, end_gff))

        print(scaffold)

    return dallfeatures, dfeaturetypes, allnames, overlaptypes


def sort_dictionary_by_value(d):
    ''' returns two sorted lists '''
    tkeys = list(d.keys())
    tvalues = [d[n] for n in tkeys]
    sortvalues, sortkeys = (list(t) for t in
                            zip(*sorted(zip(tvalues, tkeys), reverse=True)))  # *sorted would be increasing
    return sortvalues, sortkeys


def count_number_of_each_GO_term(dGO):
    d1, d2 = {}, {}  # d1 is number of times GO term is present in annotation file, d2 is the same for GO definition
    countall = 0
    for gene in dGO.keys():  # key is gene name
        if dGO[gene][8] != '':  # dGO[gene][8] is GO terms, dGO[gene][9] is GO term definition
            for GO in dGO[gene][8].split('; '):
                # add the number of reads mapping to a gene feature to each GO term associated with that gene
                d1[GO] = d1.setdefault(GO, 0) + 1
                countall += 1
            for GO in dGO[gene][9].split('; '):
                d2[GO] = d2.setdefault(GO, 0) + 1
    return(d1, d2, countall)


def read_GOterms(GOfile):
    d = {}
    count = 0
    with open(GOfile, 'r') as FILE:
        header = FILE.readline()
        for line in FILE:
            key = line.rstrip('\r\n').split('\t')[0]
            d[key] = line.rstrip('\r\n').split('\t')
            count += 1
    print('Number of GO term entries (lines in file) = %d' % (count))
    return d


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
    print('Number of scaffolds: %d' % len(list(dhead.keys())))
    print('Total number of lines in sam file (excluding header): %d' % countreads)
    return d, dhead, countreads


def read_gff3(GFF3file, features=['all']):
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
                key = line.strip().split('\t')[0]
                d.setdefault(key, []).append(line.strip().split('\t'))
            elif line.strip().split('\t')[2] in features:  # keep lines only if in features
                key = line.strip().split('\t')[0]
                d.setdefault(key, []).append(line.strip().split('\t'))
    print('Number of scaffolds: %d' % len(list(d.keys())))  # dhead.keys()
    # print(list(d.values())[:2])
    return d, dhead


def mean(lst):
    return sum(lst) / len(lst)


def main(GFF3file, samfiles, GOfile):
    import os

    # dTE, dTEhead = read_sam(TEfile, positions=[0])  # sam file of aligned TEs
    # keys of dgff3 and dsRNA is scaffold name
    dgff3, dgff3head = read_gff3(GFF3file)  # , features=['mRNA']
    # print(len(dmRNA['scaffold51_74']))
    # print(dmRNA['scaffold51_74'])

    for sam in samfiles:
        # sam files of aligned sRNAs
        dsRNA, dsRNAhead, countreads = read_sam(sam, positions=[0])  # dsRNA stands for dictionary small RNA
        dGO = read_GOterms(GOfile)

        # dNumGO = number of times GO term is present in annotation file, dNumGODefinition is the same for definition
        dNumGO, dNumGODefinition, AllGOCount = count_number_of_each_GO_term(dGO)

        dmax = {**dsRNAhead, **dgff3}  # update dsRNAhead with values from dTEhead  # changed to dgff3
        dallfeatures, dfeaturetypes, allnames, overlaptypes = count_reads_per_feature(dgff3, dsRNA)
        dGOcounts, dGOdefineCounts, totalGOcounts = count_reads_per_GOterm_using_feature_counts(dGO, dallfeatures)
        path, file = os.path.split(sam)

        dgenelengths = calculate_gene_length_gff3(dgff3)
        rpkm_dallfeatures = rpkm_dictionary(dallfeatures, dgenelengths)


        # output = ['\t'.join(n.split('_') + [str(dallfeatures[n])]) for n in allnames]
        output = ['\t'.join([n, str(dallfeatures[n])]) for n in allnames]
        outpath = os.path.join(path, '.'.join(file.split('.')[:-1] + ['allfeatures', 'fixed', 'tsv']))
        write_out(outpath, output)

        sortvalues, sortkeys = sort_dictionary_by_value(dallfeatures)
        header = 'Gene\tRead Count\tPercent of Total Reads\tTotal Reads'
        output = [header] + ['\t'.join([k, str(v), str((v / countreads) * 100), str(countreads)]) for k, v in zip(sortkeys, sortvalues)]
        outpath = os.path.join(path, '.'.join(file.split('.')[:-1] + ['allfeatures', 'sorted', 'fixed', 'tsv']))
        write_out(outpath, output)

        sortvalues, sortkeys = sort_dictionary_by_value(rpkm_dallfeatures)
        header = 'Gene\tRPKM'
        output = [header] + ['\t'.join([k, str(v)]) for k, v in zip(sortkeys, sortvalues)]
        outpath = os.path.join(path, '.'.join(file.split('.')[:-1] + ['allfeatures', 'sorted', 'RPKM', 'fixed', 'tsv']))
        write_out(outpath, output)

        sortvalues, sortkeys = sort_dictionary_by_value(dfeaturetypes)
        header = 'Feature Type\tRead Count\tPercent of Total Reads\tTotal Reads\tAverage Percent Length of Overlap of Reads Overlapping Feature'
        output = [header] + ['\t'.join([k, str(v), str((v / countreads) * 100), str(countreads), str(mean(overlaptypes[k]))]) for k, v in zip(sortkeys, sortvalues)]
        outpath = os.path.join(path, '.'.join(file.split('.')[:-1] + ['featuretypes', 'fixed', 'tsv']))
        write_out(outpath, output)

        # Output the number of times each go term is present in annotation file, decreasing
        sortvalues, sortkeys = sort_dictionary_by_value(dNumGO)
        header = 'GO Term\tNumber of Times GO Term is in File\tPercent of All GO Terms\tTotal Number of GO terms'
        output = [header] + ['\t'.join([k, str(v), str((v / AllGOCount) * 100), str(AllGOCount)]) for k, v in zip(sortkeys, sortvalues)]
        outpath = os.path.join(path, '.'.join(file.split('.')[:-1] + ['GOtermFromAnnotationFile', 'counts', 'tsv']))
        write_out(outpath, output)

        # Output the number of times each go term is present in annotation file, decreasing
        sortvalues, sortkeys = sort_dictionary_by_value(dNumGODefinition)
        header = 'GO Term Definition\tNumber of Times GO Term is in File\tPercent of All GO Terms\tTotal Number of GO terms'
        output = [header] + ['\t'.join([k, str(v), str((v / AllGOCount) * 100), str(AllGOCount)]) for k, v in zip(sortkeys, sortvalues)]
        outpath = os.path.join(path, '.'.join(file.split('.')[:-1] + ['GOtermDefinitionFromAnnotationFile', 'counts', 'tsv']))
        write_out(outpath, output)

        #output sorted GO terms, term with most reads first
        sortvalues, sortkeys = sort_dictionary_by_value(dGOcounts)
        header = 'GO Term\tRead Count\tPercent of Total Reads\tTotal Read count of Genes with GO term'
        output = [header] + ['\t'.join([k, str(v), str((v / totalGOcounts) * 100), str(totalGOcounts)]) for k, v in zip(sortkeys, sortvalues)]
        outpath = os.path.join(path, '.'.join(file.split('.')[:-1] + ['GOterm', 'counts', 'tsv']))
        write_out(outpath, output)

        #output sorted GO term descriptions, term with most reads first
        sortvalues, sortkeys = sort_dictionary_by_value(dGOdefineCounts)
        header = 'GO Term Definition\tRead Count\tPercent of Total Reads\tTotal Read count of Genes with GO term'
        output = [header] + ['\t'.join([k, str(v), str((v / totalGOcounts) * 100), str(totalGOcounts)]) for k, v in zip(sortkeys, sortvalues)]
        outpath = os.path.join(path, '.'.join(file.split('.')[:-1] + ['GOtermDefinition', 'counts', 'tsv']))
        write_out(outpath, output)

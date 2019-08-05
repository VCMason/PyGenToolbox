def plot_histogram(data, outpath):
    print('Plotting histogram to outfile: %s' % outpath)
    import matplotlib.pyplot as plt

    plt.hist(x=data, bins=50)  # bins='auto', color='#0504aa', alpha=0.7, rwidth=0.85
    plt.grid(axis='y')  # , alpha=0.75
    plt.xlabel('Relative Gene Position')
    plt.ylabel('Frequency')
    plt.title('Relative Seq Alignment For All Genes')
    #  maxfreq = n.max()
    # Set a clean upper y-axis limit.
    #  plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)
    plt.savefig(outpath)
    plt.show()
    plt.close()


def write_out(outpath, output):
    print('Writing out file: %s' % outpath)
    with open(outpath, 'w') as OUT:
        OUT.write('\n'.join(output))


def mean(lst):
    return sum(lst) / len(lst)


def count_reads_sense_antisense(dgff3, dsam):
    print('Counting reads per feature and feature type')
    dallfeatures = {}
    allnames = []
    dfeaturetypesensecounts = {}
    dfeaturetypeantisensecounts = {}
    dfeaturetypesenseproportions = {}
    dfeaturetypeantisenseproportions = {}
    dallsensecounts = {}
    dallsenseproportions = {}
    dallantisensecounts = {}
    dallantisenseproportions = {}
    for scaffold in list(dgff3.keys()):
        for gff3line in dgff3[scaffold]:
            countsensereads = 0
            countantisensereads = 0
            totalreads = 0
            start_gff = int(gff3line[3])
            end_gff = int(gff3line[4])
            feature = '_'.join([gff3line[8].split(';')[0][3:], gff3line[6]])  # [3:] cuts out 'ID=...' # adds orientation + or -
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
                        # feature types = exon, CDS, mRNA, 5' UTR, etc...
                        if gff3line[6] == '+':  # feature is 5' -> 3'
                            #####  This statement is not universal!!  #####
                            if aligned_seq[1] == 16:  # then the read is mapped in reverse
                                countantisensereads += 1
                            else:  # it is forward
                                countsensereads += 1
                        elif gff3line[6] == '-':  # feature is 3' -> 5'
                            #####  This statement is not universal!!  #####
                            if aligned_seq[1] == 16:  # then the read is mapped in reverse
                                countsensereads += 1
                            else:
                                countantisensereads += 1
                        totalreads += 1
            if totalreads != 0:
                allnames.append(feature)
                dallsensecounts.setdefault(feature, []).append(countsensereads)
                dallantisensecounts.setdefault(feature, []).append(countantisensereads)
                dallsenseproportions.setdefault(feature, []).append(countsensereads / totalreads)
                dallantisenseproportions.setdefault(feature, []).append(countantisensereads / totalreads)
                dfeaturetypesensecounts.setdefault(gff3line[2], []).append(countsensereads)
                dfeaturetypeantisensecounts.setdefault(gff3line[2], []).append(countantisensereads)
                dfeaturetypesenseproportions.setdefault(gff3line[2], []).append(countsensereads / totalreads)
                dfeaturetypeantisenseproportions.setdefault(gff3line[2], []).append(countantisensereads / totalreads)
        print(scaffold)

    return allnames, dallsensecounts, dallantisensecounts, dallsenseproportions, dallantisenseproportions, dfeaturetypesensecounts, dfeaturetypeantisensecounts, dfeaturetypesenseproportions, dfeaturetypeantisenseproportions


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
            elif features == ['all']:  # keep all lines
                key = line.strip().split('\t')[0]
                d.setdefault(key, []).append(line.strip().split('\t'))
            elif line.strip().split('\t')[2] in features:  # keep lines only if in features
                key = line.strip().split('\t')[0]
                d.setdefault(key, []).append(line.strip().split('\t'))
    print('Number of scaffolds: %d' % len(list(dhead.keys())))
    # print(list(d.values())[:2])
    return d, dhead


def main(GFF3file, samfiles):
    print('Start')
    import os

    # dTE, dTEhead = read_sam(TEfile, positions=[0])  # sam file of aligned TEs
    # keys of dgff3 and dsRNA is scaffold name
    dgff3, dgff3head = read_gff3(GFF3file)  # , features=['mRNA']

    for sam in samfiles:
        # sam files of aligned sRNAs
        dsRNA, dsRNAhead, countreads = read_sam(sam, positions=[0])  # dsRNA stands for dictionary small RNA

        dmax = {**dsRNAhead, **dgff3head}  # update dsRNAhead with values from dTEhead

        allnames, dallsensecounts, dallantisensecounts, dallsenseproportions, dallantisenseproportions, \
        dfeaturetypesensecounts, dfeaturetypeantisensecounts, dfeaturetypesenseproportions, \
        dfeaturetypeantisenseproportions = count_reads_sense_antisense(dgff3, dsRNA)

        path, file = os.path.split(sam)
        outpath = os.path.join(path, '.'.join(file.split('.')[:-1] + ['SenseVsAntiSense', 'allfeatures', 'tsv']))
        output = ['\t'.join([n, str(dallsensecounts[n]), str(dallantisensecounts[n]), str(dallsenseproportions[n]), str(dallantisenseproportions[n])]) for n in allnames]
        write_out(outpath, output)

        outpath = os.path.join(path, '.'.join(file.split('.')[:-1] + ['SenseVsAntiSense', 'featuretypes', 'tsv']))
        header = '\t'.join(['ID', 'Mean of Sense Counts', 'Mean of Antisense Counts', 'Mean of Sense Proportions', 'Mean of Antisense Proportions'])
        output = [header] + ['\t'.join([n, str(mean(dfeaturetypesensecounts[n])), str(mean(dfeaturetypeantisensecounts[n])), str(mean(dfeaturetypesenseproportions[n])), str(mean(dfeaturetypeantisenseproportions[n]))]) for n in list(dfeaturetypesensecounts.keys())]
        write_out(outpath, output)

    print('End')

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


def collect_all_covered_normalized_positions_per_feature(dgff3, dsam, target_feature='gene'):
    print('For each position in feature: %s, normalize by gene length, record value for histogram' % target_feature)

    allhistogram = []  # relative position for each gene (1st position = 1/len(gene), 5th = 5/len(gene))
    allfeaturelengths = []  # list of lengths of all features of type == target_feature
    for scaffold in list(dgff3.keys()):
        for gff3line in dgff3[scaffold]:
            start_gff = int(gff3line[3])
            end_gff = int(gff3line[4])
            featurelength = end_gff - start_gff
            if gff3line[2] == target_feature:
                allfeaturelengths.append(featurelength)
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
                        # if gff3line[1] == '+':  # feature is 5' -> 3'
                        # elif gff3line[1] == '-':  # feature is 3' -> 5'
                        if (end_sam - start_sam) == 0:  # if (end_gff - start_gff) == 0:  # add one?
                            # avoid division by zero
                            print('Potential division by zero: %s, %s, %d, %d' % (scaffold, gff3line[2], start_gff, end_gff))
                        elif (start_sam <= start_gff) and (end_sam <= end_gff):
                            # 'right' side of aligned read overlaps feature but 'left' does not
                            if gff3line[2] == target_feature:
                                # maybe round(i / featurelength, 2)
                                allhistogram.extend([i / featurelength for i in range(1, end_sam - start_gff)])
                        elif (start_sam >= start_gff) and (end_sam >= end_gff):
                            # 'left' side of aligned read overlaps feature but 'right' does not
                            if gff3line[2] == target_feature:
                                allhistogram.extend([i / featurelength for i in range(start_sam - start_gff, featurelength)])  # end_gff - start_gff
                        elif (start_gff <= start_sam) and (end_sam <= end_gff):
                            # if read aligned is inside feature and smaller than feature
                            if gff3line[2] == target_feature:
                                allhistogram.extend([i / featurelength for i in range(start_sam - start_gff, end_sam - start_gff)])
                        elif (start_sam <= start_gff) and (end_sam >= end_gff):
                            # if read aligned is larger than feature aligned to...
                            #     although read completely covers feature, the feature does not completely cover read
                            if gff3line[2] == target_feature:
                                allhistogram.extend([i / featurelength for i in range(1, featurelength)])
                        else:
                            print('some form of overlap was not accounted for')
                            print('%s, SS: %d, SG: %d, ES: %d, EG: %d' % (scaffold, start_sam, start_gff, end_sam, end_gff))

        # print(scaffold)
    print('Summary of feature: %s, Number of features: %d, Mean length: %.2f' % (target_feature, len(allfeaturelengths), mean(allfeaturelengths)))

    return allhistogram


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


def main(GFF3file, samfiles, target_feature='gene'):
    import os

    # dTE, dTEhead = read_sam(TEfile, positions=[0])  # sam file of aligned TEs
    # keys of dgff3 and dsRNA is scaffold name
    dgff3, dgff3head = read_gff3(GFF3file)  # , features=['mRNA']
    # print(len(dmRNA['scaffold51_74']))
    # print(dmRNA['scaffold51_74'])

    for sam in samfiles:
        # sam files of aligned sRNAs
        dsRNA, dsRNAhead, countreads = read_sam(sam, positions=[0])  # dsRNA stands for dictionary small RNA

        dmax = {**dsRNAhead, **dgff3head}  # update dsRNAhead with values from dTEhead
        allhistogram = collect_all_covered_normalized_positions_per_feature(dgff3, dsRNA, target_feature)

        path, file = os.path.split(sam)
        outpath = os.path.join(path, '.'.join(file.split('.')[:-1] + ['NormCov', target_feature, 'pdf']))
        plot_histogram(allhistogram, outpath)

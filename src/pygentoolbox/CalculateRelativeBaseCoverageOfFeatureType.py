def plot_sns_kde(df, outpath, target_feature='Gene'):
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt
    # allnames is 1d list of column names, data is 2d list: elements in container list are columns

    sns.kdeplot(data=df)  # bw_method=0.3
    plt.grid(axis='y')  # , alpha=0.75
    plt.xlabel('Relative Gene Position')
    plt.ylabel('Frequency')
    plt.title(f'Relative Coverage For All {target_feature}')
    #  maxfreq = n.max()
    # Set a clean upper y-axis limit.
    #  plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)
    plt.savefig(outpath)
    plt.show()
    plt.close()


def plot_multi_kde(df, outpath, target_feature='Gene'):
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt
    # allnames is 1d list of column names, data is 2d list: elements in container list are columns

    df.plot.kde()  # bw_method=0.3
    plt.grid(axis='y')  # , alpha=0.75
    plt.xlabel('Relative Gene Position')
    plt.ylabel('Frequency')
    plt.title(f'Relative Coverage For All {target_feature}')
    #  maxfreq = n.max()
    # Set a clean upper y-axis limit.
    #  plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)
    plt.savefig(outpath)
    plt.show()
    plt.close()


def plot_multihistogram(data, outpath, feature, numbins=10):  # names
    # data is 2d list
    # outpath is full output path
    # numbins = 10 or 50 etc...
    print('Plotting histogram to outfile: %s' % outpath)
    import matplotlib.pyplot as plt

    # bins='auto', color='#0504aa', alpha=0.7, rwidth=0.85  # , label=names
    plt.hist(x=data, bins=numbins, histtype='step', stacked=True, fill=False)
    #plt.legend(prop={'size': 10})
    # plt.hist(x, n_bins, density=True, histtype='bar', stacked=True)
    # or
    # colors = ['red', 'tan', 'lime']
    # plt.hist(x, n_bins, density=True, histtype='bar', color=colors, label=colors)
    # plt.legend(prop={'size': 10})
    # plt.set_title('bars with legend')
    # or of different sample sizes
    # x_multi = [np.random.randn(n) for n in [10000, 5000, 2000]]
    # ax3.hist(x_multi, n_bins, histtype='bar')
    # ax3.set_title('different sample sizes')
    plt.grid(axis='y')  # , alpha=0.75
    plt.xlabel('Relative Gene Position')
    plt.ylabel('Frequency')
    plt.title('Normalized %s Coverage' % feature)
    #  maxfreq = n.max()
    # Set a clean upper y-axis limit.
    #  plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)
    plt.savefig(outpath)
    plt.show()
    plt.close()


def plot_histogram(data, outpath, numbins=50):
    print('Plotting histogram to outfile: %s' % outpath)
    import matplotlib.pyplot as plt

    plt.hist(x=data, bins=numbins)  # bins='auto', color='#0504aa', alpha=0.7, rwidth=0.85
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


def collect_contacts_per_feature_normalized_by_feature_length(dgff3, dcontacts, target_feature='gene'):
    print('For each position in feature: %s, normalize by gene length, record value for histogram' % target_feature)

    allhistogram = []  # relative contact position inside a feature (scaled by strength of contact)
    countleft, counties, countright = 0, 0, 0
    posleft, posies, posright = [], [], []
    # contactpos - start_gff / len(feature) # one entry per 1 strength of contact (add 7 of these if strength = 7)
    # all histogram will be a 1d list of all contacts normalized by feautre length
    # allhistogram will contain all contacts 1 featurelength to the left of feature, within feature, and 1 featurelength to the right of the feature
    # <0 == left of feature. 0-1 == inside feature. >1 == right of feature
    allfeaturelengths = []  # list of lengths of all features of type == target_feature
    for scaffold in list(dgff3.keys()):
        for gff3line in dgff3[scaffold]:
            #print(gff3line)
            start_gff = int(gff3line[3])
            end_gff = int(gff3line[4])
            featurelength = end_gff - start_gff
            ncRNAtype = ''
            if gff3line[2] == 'ncRNA':
                # then iterate through metadata and check for type of ncRNA: rRNA, tRNA, snoRNA
                for info in gff3line[8].split(';'):
                    if info[:5] == 'type=':
                        ncRNAtype = info[5:]
            if (gff3line[2] == target_feature) or (ncRNAtype == target_feature):
                allfeaturelengths.append(featurelength)
            try:
                dcontacts[scaffold]
            except:
                pass  # nothing aligned to scaffold in dsam
            else:
                for contact in dcontacts[scaffold]:
                    #contact is list of DNA scaffold, position, strength of contact
                    contactpos = int(contact[1])
                    contactstrength = int(contact[2])
                    outsideleft = start_gff - featurelength - 1
                    outsideright = end_gff + featurelength
                    if (start_gff <= contactpos) and (contactpos <= end_gff):
                        # if contact inside feature
                        counties += 1
                        posies.append(f'{scaffold}_{contactpos}')
                        for i in range(contactstrength):
                            allhistogram.append((contactpos - start_gff) / featurelength)
                    elif (outsideleft <= contactpos) and (contactpos <= start_gff):
                        # if contact is left of feature
                        countleft += 1
                        posleft.append(f'{scaffold}_{contactpos}')
                        for i in range(contactstrength):
                            allhistogram.append(((contactpos - outsideleft) / featurelength) - 1)
                    elif (end_gff <= contactpos) and (contactpos <= outsideright):
                        # if contact is right of feature
                        countright += 1
                        posright.append(f'{scaffold}_{contactpos}')
                        for i in range(contactstrength):
                            allhistogram.append(((contactpos - end_gff) / featurelength) + 1)

    return allhistogram, counties, countleft, countright, posies, posleft, posright


def collect_all_covered_normalized_positions_per_feature(dgff3, dsam, target_feature='gene'):
    print('For each position in feature: %s, normalize by gene length, record value for histogram' % target_feature)

    allhistogram = []  # relative position for each gene (1st position = 1/len(gene), 5th = 5/len(gene))
    allfeaturelengths = []  # list of lengths of all features of type == target_feature
    for scaffold in list(dgff3.keys()):
        for gff3line in dgff3[scaffold]:
            start_gff = int(gff3line[3])
            end_gff = int(gff3line[4])
            featurelength = end_gff - start_gff
            ncRNAtype = ''
            if gff3line[2] == 'ncRNA':
                # then iterate through metadata and check for type of ncRNA: rRNA, tRNA, snoRNA
                for info in gff3line[8].split(';'):
                    if info[:5] == 'type=':
                        ncRNAtype = info[5:]
            if (gff3line[2] == target_feature) or (ncRNAtype == target_feature):
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
                            if (gff3line[2] == target_feature) or (ncRNAtype == target_feature):
                                # maybe round(i / featurelength, 2)
                                if gff3line[6] == '+':
                                    allhistogram.extend([float(i) / featurelength for i in range(1, end_sam - start_gff)])
                                elif gff3line[6] == '-':
                                    allhistogram.extend([float(i) / featurelength for i in range(end_gff - end_sam, featurelength)])
                        elif (start_sam >= start_gff) and (end_sam >= end_gff):
                            # 'left' side of aligned read overlaps feature but 'right' does not
                            if (gff3line[2] == target_feature) or (ncRNAtype == target_feature):
                                if gff3line[6] == '+':
                                    allhistogram.extend([float(i) / featurelength for i in range(start_sam - start_gff, featurelength)])  # end_gff - start_gff
                                elif gff3line[6] == '-':
                                    allhistogram.extend([float(i) / featurelength for i in range(1, end_gff - start_sam)])  # end_gff - start_gff
                        elif (start_gff <= start_sam) and (end_sam <= end_gff):
                            # if read aligned is inside feature and smaller than feature
                            if (gff3line[2] == target_feature) or (ncRNAtype == target_feature):
                                if gff3line[6] == '+':
                                    allhistogram.extend([float(i) / featurelength for i in range(start_sam - start_gff, end_sam - start_gff)])
                                elif gff3line[6] == '-':
                                    allhistogram.extend([float(i) / featurelength for i in range(end_gff - end_sam, end_gff - start_sam)])
                        elif (start_sam <= start_gff) and (end_sam >= end_gff):
                            # if read aligned is larger than feature aligned to...
                            #     although read completely covers feature, the feature does not completely cover read
                            if (gff3line[2] == target_feature) or (ncRNAtype == target_feature):
                                if gff3line[6] == '+':
                                    allhistogram.extend([float(i) / featurelength for i in range(1, featurelength)])
                                elif gff3line[6] == '-':  # it's the same
                                    allhistogram.extend([float(i) / featurelength for i in range(1, featurelength)])
                        else:
                            print('some form of overlap was not accounted for')
                            print('%s, SS: %d, SG: %d, ES: %d, EG: %d' % (scaffold, start_sam, start_gff, end_sam, end_gff))

        # print(scaffold, end=', ')
    print('Summary of feature: %s, Number of features: %d, Mean length: %.2f' % (target_feature, len(allfeaturelengths), mean(allfeaturelengths)))

    return allhistogram


def read_contacts(filename):
    import math

    print(f'Reading contact file: {filename}')
    print('Gathering DNA contact position and strength')
    # key is scaffold name
    d = {}
    count, countunique = 0, 0
    with open(filename, 'r') as FILE:
        for line in FILE:
            dnascaff = 'scaffold51_' + line.strip().split()[3][2:] + '_with_IES'
            dnapos = int(line.strip().split()[4])
            dnascore = int(math.ceil(float(line.strip().split()[6])))  # round up to nearest integer
            d.setdefault(dnascaff, []).append([dnascaff, dnapos, dnascore])
            count += dnascore
            countunique += 1
    print(f'Number of total contacts: {count}')
    print(f'Number of unique contacts: {countunique}')
    return d, count, countunique


def read_sam(sam, positions=[]):
    ''' key is scaffold, value is line.strip() '''
    import os
    
    print('Reading sam file: %s' % sam)
    p, f = os.path.split(sam)
    print('File: %s' % f)
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
    print('Number of scaffolds: %d' % len(list(d.keys())))
    # print(list(d.values())[:2])
    return d, dhead


def main(GFF3file, samfiles, sizefactors, datatype='sam', target_feature='gene', numbins=50, normalizationmultiple=1):
    # datatype is string. dattype can be sam, contacts, or bed
    # size factors are the value you multiple the cases by to change numbers to be comparable with control
    # sizes factors can be [1,1,1,1] if no changes are applies for 4 sam files
    import os
    import pandas as pd
    import random

    # dTE, dTEhead = read_sam(TEfile, positions=[0])  # sam file of aligned TEs
    # keys of dgff3 and dsRNA is scaffold name
    dgff3, dgff3head = read_gff3(GFF3file)  # , features=['mRNA']
    # print(len(dmRNA['scaffold51_74']))
    # print(dmRNA['scaffold51_74'])

    legendallfiles = []  # will be unique legend identifier for multi histogram # 1d list
    normcovallfiles = []  # will be a 2d array containing normalized feature coverage. One list per file [[], [], []...]
    allnames = []
    for sam, sizefactor in zip(samfiles, sizefactors):
        # sam files of aligned sRNAs
        if datatype == 'sam':
            dsRNA, dsRNAhead, countreads = read_sam(sam, positions=[0])  # dsRNA stands for dictionary small RNA

            if countreads == 0:
                # pass because there are no mapped reads, only reads available would be unaligned reads
                pass
            else:
                print('Mapped reads, so calculating feature type coverage normalized by feature length')
                dmax = {**dsRNAhead, **dgff3head}  # update dsRNAhead with values from dTEhead
                normcovdata = collect_all_covered_normalized_positions_per_feature(dgff3, dsRNA, target_feature)
                normcovallfiles.append(normcovdata)

                path, file = os.path.split(sam)
                # legendallfiles.append(sam.split('.')[0][-4:])  # extracting 15bp, or 27bp etc. NOT UNIVERSAL!
                outpath = os.path.join(path, '.'.join(file.split('.')[:-1] + ['NormCov', target_feature, 'pdf']))
                plot_histogram(normcovdata, outpath, numbins)
        elif datatype == 'contacts':
            dcontacts, allcountacts, alluniquecontacts = read_contacts(sam)
            normcovdata, counties, countleft, countright, posies, posleft, posright = collect_contacts_per_feature_normalized_by_feature_length(dgff3, dcontacts, target_feature)
            print(f'Number of unique contacts left of IESs: {countleft}')
            print(f'Number of unique contacts within IESs: {counties}')
            print(f'Number of unique contacts right of IESs: {countright}')

            print(f'Frequency of unique contacts left of IESs {countleft / alluniquecontacts}')
            print(f'Frequency of unique contacts within IESs {counties / alluniquecontacts}')
            print(f'Frequency of unique contacts right of IESs {countright / alluniquecontacts}')

            path, file = os.path.split(sam)
            # legendallfiles.append(sam.split('.')[0][-4:])  # extracting 15bp, or 27bp etc. NOT UNIVERSAL!
            outpath = os.path.join(path, '.'.join(file.split('.')[:-1] + ['NormCov', target_feature, 'pdf']))
            plot_histogram(normcovdata, outpath, numbins)

            print('Normalize central feature by normalization multiple:')
            outpath = os.path.join(path, '.'.join(file.split('.')[:-1] + ['NormCov', 'NormGATC', target_feature, 'pdf']))
            additionalentries = []
            for elem in normcovdata:
                if (elem >= 0.0001) and (elem <= 1.0):
                    additionalentries += [elem]*(normalizationmultiple - 1)  # subtract 1 b/c 1 already in normcovdata
            plot_histogram(normcovdata + additionalentries, outpath, numbins)

            print(f'Normalize by scaling factors: {sizefactors}')
            outpath = os.path.join(path, '.'.join(file.split('.')[:-1] + ['NormCov', 'NormGATC', 'NormEVKD', target_feature, 'pdf']))
            normevkd = normcovdata + additionalentries
            if sizefactor == 1:
                print(f'Size Factor = {sizefactor}: no adjustment needed')
            elif sizefactor > 1:  # increase entry number
                print(f'Size Factor = {sizefactor}: increasing number of results')
                keep = []  # these are element indeces that i will add to normevkd
                for count, elem in enumerate(normevkd):
                    if random.uniform(0, 1) >= (1 / sizefactor):
                        keep.append(count)
                normevkd = normevkd + [normevkd[i] for i in keep]
                plot_histogram(normevkd, outpath, numbins)
            elif sizefactor < 1:  # decrease entry number
                print(f'Size Factor = {sizefactor}: decreasing number of results')
                keep = []  # these are element indeces that i will keep from normevkd (excluding the remainder)
                for count, elem in enumerate(normevkd):
                    if random.uniform(0, 1) <= sizefactor:
                        keep.append(count)
                normevkd = [normevkd[i] for i in keep]
                plot_histogram(normevkd, outpath, numbins)

            # strposies = ', '.join(posies)
            # strposleft = ', '.join(posleft)
            # strposright = ', '.join(posright)
            # print(f'Contacts inside IESs: {strposies}\n')
            # print(f'Contacts inside IESs: {strposleft}\n')
            # print(f'Contacts inside IESs: {strposright}\n')
            normcovallfiles.append(normevkd)

            pathminusonedir, onedir = os.path.split(path)
            pathminustwodir, twodir = os.path.split(pathminusonedir)
            allnames.append(twodir)

        elif datatype == 'bed':
            pass

    if datatype == 'sam':
        outpath = os.path.join(path, '.'.join(['AllNormalizedCoverageHistograms', target_feature, 'pdf']))
        plot_multihistogram(normcovallfiles, outpath, target_feature, numbins)  # legendallfiles
    elif datatype == 'contacts':
        outpath = os.path.join(path, '.'.join(['AllNormalizedCoverageHistograms', target_feature, 'pdf']))
        # allnames is 1d list of column names, data is 2d list elements in container list are columns
        df = pd.DataFrame(list(map(list, zip(*normcovallfiles))), columns=allnames)
        print(df.head())
        plot_multi_kde(df, outpath, target_feature)




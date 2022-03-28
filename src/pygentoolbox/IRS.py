############################################################################################
##########        # v02 finds closest mRNA to each sRNA on each scaffold          ##########
############################################################################################


def plot_hexbins(x, y, outpath, nbins=20, xlab='mRNA position', ylab='sRNA position'):
    # not too slow
    print('Plotting scatter, hexbins, 2d histogram')
    # Libraries
    import numpy as np
    import matplotlib.pyplot as plt

    x, y = np.array(x), np.array(y)

    # Create a figure with 6 plot areas
    fig, axes = plt.subplots(ncols=3, nrows=1, figsize=(21, 5))

    # Everything sarts with a Scatterplot
    axes[0].set_title('Scatterplot')
    axes[0].plot(x, y, 'ko')
    # As you can see there is a lot of overplottin here!

    # Thus we can cut the plotting window in several hexbins
    # nbins = 20
    axes[1].set_title('Hexbin')
    axes[1].hexbin(x, y, gridsize=nbins, cmap=plt.cm.BuGn_r)

    # 2D Histogram
    axes[2].set_title('2D Histogram')
    axes[2].hist2d(x, y, bins=nbins, cmap=plt.cm.BuGn_r)

    for ax in axes.flat:
        ax.set(xlabel=xlab, ylabel=ylab)

    # Hide x labels and tick labels for top plots and y ticks for right plots.
    for ax in axes.flat:
        ax.label_outer()

    plt.show()
    plt.savefig(outpath)
    plt.close()


def plot_density(x, y, outpath, nbins=20):
    print('Plotting density')
    import matplotlib.pyplot as plt
    import numpy as np
    from scipy.stats import kde

    # Evaluate a gaussian kde on a regular grid of nbins x nbins over data extents
    x = np.array(x)
    y = np.array(y)
    # nbins = 300
    k = kde.gaussian_kde([x, y])
    xi, yi = np.mgrid[x.min():x.max():nbins * 1j, y.min():y.max():nbins * 1j]
    zi = k(np.vstack([xi.flatten(), yi.flatten()]))

    # Make the plot
    plt.pcolormesh(xi, yi, zi.reshape(xi.shape))

    # Change color palette
    plt.pcolormesh(xi, yi, zi.reshape(xi.shape), cmap=plt.cm.Greens_r)
    plt.show()
    plt.savefig(outpath)
    plt.close()


def mean(lst):
    return sum(lst) / len(lst)


def IRS(dgff3, dsam_mac, dsam_mac_ies):
    '''dgff3 and dsam have all lines saved as a list for values and the scaffold where they mapped as key'''
    '''for each ies feature, determine the IRS'''
    import re

    d_mac, d_left_mac_ies, d_right_ies_mac, d_both_mac_ies_mac = {}, {}, {}, {}
    for scaffold in list(dgff3.keys()):
        if (scaffold in list(dsam_mac.keys())) and (scaffold in list(dsam_mac_ies.keys())):
            print('%s is present in d_mac, dsam_mac, and dsam_mac_ies' % scaffold)
            for gff3line in dgff3[scaffold]:
                # Isolate name of ies and removes 'ID=' from front
                ies_name = gff3line[-1].split(';')[0][3:]
                # initiate values for each IES dictionary
                d_mac[ies_name] = 0
                # initiate values for each IES dictionary
                d_left_mac_ies[ies_name] = 0
                d_right_ies_mac[ies_name] = 0
                d_both_mac_ies_mac[ies_name] = 0
                for seq_mac in dsam_mac[scaffold]:
                    # Isolate coordinates for mac (somatic genome) with IES seqs eliminated
                    start_gff_mac, end_gff_mac = int(gff3line[-1].split(';')[0].split('.')[-1]), int(gff3line[-1].split(';')[0].split('.')[-1]) + 1
                    start_sam_mac, end_sam_mac, CIGAR_mac = seq_mac[0], seq_mac[1], seq_mac[2]

                    # print('%d, %d, %d, %d, %s' % (start_gff_mac, end_gff_mac, start_sam_mac, end_sam_mac, CIGAR_mac))
                    # if we assume that  ranges are well-formed (so that x1 <= x2 and y1 <= y2) then
                    # x1 <= y2 && y1 <= x2
                    # meaning the aligned sequence overlaps in some way to the gff3 feature
                    # overlaps the coordinates when ies's are eliminated
                    # print('### Count reads aligned to Mac (IESs eliminated) features ###')
                    if (start_gff_mac <= end_sam_mac) and (end_gff_mac <= end_gff_mac):
                        # now check CIGAR string if deletion
                        # re.findall() returns all the non-overlapping matches of patterns
                        # in a string as a list of strings
                        cigar = re.findall(r'\d+[A-Z]', CIGAR_mac)
                        # then the ends match and no hard or soft clipping
                        if cigar != []:  # if there was a match # sometimes CIGAR string == *
                            if (cigar[0][-1] == 'M') and (cigar[-1][-1] == 'M') and ('D' not in CIGAR_mac) and ('I' not in CIGAR_mac):
                                d_mac[ies_name] = d_mac.setdefault(ies_name, 0) + 1
                                #print(d_mac[ies_name])
                                #gate = 1
                                #for i in cigar:
                                #    if ('D' in i) or ('I' in i):  # for each deletion
                                #        gate = 0
                                #if gate == 1:
                                #    d_mac.setdefault(ies_name, 0) + 1

                # print('### Count reads aligned to Mac + IES features ###')
                for seq_mac_ies in dsam_mac_ies[scaffold]:
                    # Isolate coordinates for mac + IES ("germline" genome) when ies has been retained
                    start_gff, end_gff = int(gff3line[3]), int(gff3line[4])  # ies start coord, ies end coord
                    start_sam, end_sam, CIGAR = seq_mac_ies[0], seq_mac_ies[1], seq_mac_ies[2]
                    # meaning the aligned sequence overlaps in some way to the gff3 feature
                    # overlaps the ies feature when ies's are retained
                    if (start_gff <= end_sam) and (start_sam <= end_gff):
                        # now check CIGAR string if deletion
                        # re.findall() returns all the non-overlapping matches of patterns
                        # in a string as a list of strings
                        cigar = re.findall(r'\d+[A-Z]', CIGAR)
                        if cigar != []:  # if there was a match  # sometimes CIGAR string == *
                            # then the ends match and no hard or soft clipping
                            if (cigar[0][-1] == 'M') and (cigar[-1][-1] == 'M') and ('D' not in CIGAR) and ('I' not in CIGAR):
                                if (end_sam - start_sam) == 0:  # if (end_gff - start_gff) == 0:  # add one?
                                    # avoid division by zero
                                    print('Division by zero: %s, %s, %d, %d' % (scaffold, gff3line[2], start_gff, end_gff))
                                elif (start_sam <= start_gff) and (end_sam <= end_gff):
                                    # 'right' side of aligned read overlaps feature but 'left' does not
                                    #  if gff3line[2] == target_feature:
                                    d_left_mac_ies[ies_name] = d_left_mac_ies.setdefault(ies_name, 0) + 1
                                elif (start_sam >= start_gff) and (end_sam >= end_gff):
                                    # 'left' side of aligned read overlaps feature but 'right' does not
                                    d_right_ies_mac[ies_name] = d_right_ies_mac.setdefault(ies_name, 0) + 1
                                elif (start_gff <= start_sam) and (end_sam <= end_gff):
                                    # if read aligned is inside feature and smaller than feature
                                    # pass because i only want reads that overlap the left or right boundary (or both)
                                    pass
                                elif (start_sam <= start_gff) and (end_sam >= end_gff):
                                    # if read aligned is larger than feature aligned to...
                                    # although read completely covers feature, the feature does not completely cover read
                                    d_both_mac_ies_mac[ies_name] = d_both_mac_ies_mac.setdefault(ies_name, 0) + 1
                                else:
                                    print('some form of overlap was not accounted for')
                                    print('%s, SS: %d, SG: %d, ES: %d, EG: %d' % (scaffold, start_sam, start_gff, end_sam, end_gff))
                                    
    d_IRS, d_IRS_Alt = {}, {}  # classical way of calculating IRS, and alternative method
    print('Calculating IRS')
    countIRS = 0
    for k in list(d_mac.keys()):
        # classical
        if d_left_mac_ies[k] + d_mac[k] == 0:
            left_score = 0
        else:
            left_score = d_left_mac_ies[k] / (d_left_mac_ies[k] + d_mac[k])
        if d_right_ies_mac[k] + d_mac[k] == 0:
            right_score = 0
        else:
            right_score = d_right_ies_mac[k] / (d_right_ies_mac[k] + d_mac[k])

        d_IRS[k] = (left_score + right_score) / 2

        # alternative
        if d_left_mac_ies[k] + d_right_ies_mac[k] + d_both_mac_ies_mac[k] + d_mac[k] == 0:
            # If denominator is zero..
            d_IRS_Alt[k] = 0
        else:
            d_IRS_Alt[k] = (d_left_mac_ies[k] + d_right_ies_mac[k] + d_both_mac_ies_mac[k]) / (
                        d_left_mac_ies[k] + d_right_ies_mac[k] + d_both_mac_ies_mac[k] + d_mac[k])
        countIRS += 1
        if countIRS % 10000 == 0:
            print('IRS %d' % countIRS)

    return d_IRS, d_IRS_Alt


def subsample_dictionary(d, number=15000000):
    import random

    # print("choosing 2 random items from a dictionary using sample method ", random.sample(d.items(), k=2))
    # {k: v for (k, v) in random.sample(d.items(), k=3)}
    # random.sample returns list of tuples, so i need to convert back to dictionary
    new_d = {k: v for (k, v) in random.sample(d.items(), k=number)}
    return new_d


def read_sam_subsample(sam, freq=0.2, positions=[0]):  # , number=15000000
    ''' key is scaffold, value is line.strip() '''
    from random import random

    print('Reading sam file: %s' % sam)
    d = {}
    dhead = {}
    #countlines = 0
    #with open(sam, 'r') as FILE:
    #    for line in FILE:
    #        if line[0] != '@':
    #            countlines += 1
    #print('Number of lines in Sam file: %d' % countlines)

    with open(sam, 'r') as FILE:
        # lines = [line for line in FILE if random() <= .25]
        # print('Number of subsampled lines: %d' % len(lines))
        countlines = 0
        for line in FILE:
            if (line[0] != '@') and (random() <= freq):
                countlines += 1
                if int(line.strip().split('\t')[3]) not in positions:
                    key = line.strip().split('\t')[2]
                    d.setdefault(key, []).append(line.strip().split('\t'))
            elif line[:3] == '@SQ':
                dhead[line.strip().split('\t')[1][3:]] = int(line.strip().split('\t')[2].split(':')[1])
    print('Number of scaffolds: %d' % len(list(dhead.keys())))
    print('Number of sub-sampled reads (%.2f percent of reads) = %d' % (freq*100, countlines))
    return d, dhead


def read_sam(sam, positions=[]):
    ''' key is scaffold, value is line.strip() '''
    print('Reading sam file: %s' % sam)
    d = {}
    dhead = {}
    with open(sam, 'r') as FILE:
        for line in FILE:
            if line[0] != '@':
                if int(line.strip().split('\t')[3]) not in positions:
                    key = line.strip().split('\t')[2]
                    d.setdefault(key, []).append(line.strip().split('\t'))
            elif line[:3] == '@SQ':
                dhead[line.strip().split('\t')[1][3:]] = int(line.strip().split('\t')[2].split(':')[1])
    print('Number of scaffolds: %d' % len(list(dhead.keys())))
    return d, dhead


def read_gff3(GFF3file, features=['all']):
    print('Reading Gff3 file: %s' % GFF3file)
    d = {}
    dhead = {}
    with open(GFF3file, 'r') as FILE:
        for line in FILE:
            if line[0] == '#':
                pass
            # if line[:5] == '##gff':
            #     pass
            # elif line[:17] == '##sequence-region':
            #     dhead[line.strip().split()[1]] = int(line.strip().split()[3])
            elif features == ['all']:  # keep all lines
                key = line.strip().split('\t')[0][:-len('_with_IES')]  # remove the extra text from the 'with ies' .gff3 file
                d.setdefault(key, []).append(line.strip().split('\t'))
            elif line.strip().split('\t')[2] in features:  # keep lines only if in features
                key = line.strip().split('\t')[0][:-len('_with_IES')]  # remove the extra text from the 'with ies' .gff3 file
                d.setdefault(key, []).append(line.strip().split('\t'))
    print('Number of scaffolds: %d' % len(list(d.keys())))
    print(list(d.keys())[:10])
    return d, dhead


def read_sam_subsample_coordinates(sam, freq):
    ''' key is scaffold, value is line.strip() '''
    from random import random

    print('Reading, sub-sampling, and only coordinates for sam file: %s' % sam)
    d = {}
    dhead = {}
    with open(sam, 'r') as FILE:
        # lines = [line for line in FILE if random() <= .25]
        # print('Number of subsampled lines: %d' % len(lines))
        countlines = 0
        for line in FILE:
            if (line[0] != '@') and (random() <= freq):
                countlines += 1
                if line.strip().split('\t')[2][-len('_with_ies'):] == '_with_IES':
                    key = line.strip().split('\t')[2][:-len('_with_IES')]
                else:
                    key = line.strip().split('\t')[2]
                # key = str(scaffold), value = list of all [int(start coord), int(end coord), str(CIGAR)]
                d.setdefault(key, []).append([int(line.strip().split('\t')[3]), int(line.strip().split('\t')[3]) + len(line.strip().split('\t')[9]), line.strip().split('\t')[5]])
    print('Number of scaffolds: %d' % len(list(d.keys())))
    print(list(d.keys())[:10])
    print('Number of sub-sampled reads (%.2f percent of reads) = %d' % (freq*100, countlines))
    return d


def main(GFF3file, samfiles_mac, samfiles_mac_ies, freq=0.2, features=['all']):
    import os

    # dTE, dTEhead = read_sam(TEfile, positions=[0])  # sam file of aligned TEs

    dgff3, dgff3head = read_gff3(GFF3file, features)
    # print(len(dmRNA['scaffold51_74']))
    # print(dmRNA['scaffold51_74'])

    for sam_mac, sam_mac_ies in zip(samfiles_mac, samfiles_mac_ies):
        # dsam, dsamhead = read_sam(sam,
        #                           positions=[0])  # sam files of aligned sRNAs # dsRNA stands for dictionary small RNA
        # dsam, dsamhead = read_sam_subsample(sam, freq, positions=[0])
        dsam_mac = read_sam_subsample_coordinates(sam_mac, freq)
        dsam_mac_ies = read_sam_subsample_coordinates(sam_mac_ies, freq)

        # dmax = {**dsamhead, **dgff3head}  # update dsRNAhead with values from dTEhead

        print('##### Starting IRS Calculations #####')
        d_IRS, d_IRS_Alt = IRS(dgff3, dsam_mac, dsam_mac_ies)  # or dsam for all reads
        print('Number of IRS values: %d' % len(d_IRS.keys()))
        print(list(d_IRS.items())[:20])

        path, file = os.path.split(sam_mac)

        outpath = os.path.join(path, '.'.join(file.split('.')[:-1] + ['IRS', 'tsv']))
        with open(outpath, 'w') as OUT:
            OUT.write('\n'.join(['\t'.join([k, str(v)]) for (k, v) in d_IRS.items()]))
        outpath = os.path.join(path, '.'.join(file.split('.')[:-1] + ['IRS', 'Alternative', 'tsv']))
        with open(outpath, 'w') as OUT:
            OUT.write('\n'.join(['\t'.join([k, str(v)]) for (k, v) in d_IRS_Alt.items()]))

        # plot_hexbins(x, y, outpath, 20)

        print('###########################################')
        print('###########################################')
        print('###########################################')

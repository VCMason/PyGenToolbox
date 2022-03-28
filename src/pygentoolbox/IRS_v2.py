############################################################################################
##########        # v02                                                           ##########
############################################################################################


def mean(lst):
    return sum(lst) / len(lst)


def IRS(dgff3, dsam_mac, dsam_mac_ies):
    '''dgff3 and dsam have all lines saved as a list for values and the scaffold where they mapped as key'''
    '''for each ies feature, determine the IRS'''
    import re

    count_scaffolds = 0
    d_mac, d_left_mac_ies, d_right_ies_mac, d_both_mac_ies_mac = {}, {}, {}, {}
    for scaffold in list(dgff3.keys()):
        if (scaffold in list(dsam_mac.keys())) and (scaffold in list(dsam_mac_ies.keys())):
            # print('%s is present in d_mac, dsam_mac, and dsam_mac_ies' % scaffold)
            for gff3line in dgff3[scaffold]:
                # Isolate name of ies and removes 'ID=' from front
                ies_name = gff3line[-1].split(';')[0][3:]
                # # ies start coord, ies end coord in mac when IES is deleted
                start_gff_mac, end_gff_mac = int(gff3line[-1].split(';')[0].split('.')[-1]), int(
                    gff3line[-1].split(';')[0].split('.')[-1]) + 1
                # ies start coord, ies end coord when IES is retained
                start_gff, end_gff = int(gff3line[3])+5, int(gff3line[4])-5  # 5 bases excluded b/c alignment errors
                # initiate values for each IES dictionary
                d_mac[ies_name] = 0
                # initiate values for each IES dictionary
                d_left_mac_ies[ies_name] = 0
                d_right_ies_mac[ies_name] = 0
                d_both_mac_ies_mac[ies_name] = 0
                for seq_mac in dsam_mac[scaffold]:
                    # Isolate coordinates for mac (somatic genome) with IES seqs eliminated
                    start_sam_mac, end_sam_mac = seq_mac[1], seq_mac[2]
                    # print('%d, %d, %d, %d, %s' % (start_gff_mac, end_gff_mac, start_sam_mac, end_sam_mac, CIGAR_mac))
                    # if we assume that  ranges are well-formed (so that x1 <= x2 and y1 <= y2) then
                    # x1 <= y2 && y1 <= x2
                    # meaning the aligned sequence overlaps in some way to the gff3 feature
                    # overlaps the coordinates when ies's are eliminated
                    # print('### Count reads aligned to Mac (IESs eliminated) features ###')
                    if (start_gff_mac <= end_sam_mac) and (end_gff_mac <= end_gff_mac):
                        if (start_sam_mac <= start_gff_mac) and (end_sam_mac >= end_gff_mac):
                            # although read completely covers feature, the feature does not completely cover read
                            d_mac[ies_name] = d_mac[ies_name] + 1

                # print('### Count reads aligned to Mac + IES features ###')
                for seq_mac_ies in dsam_mac_ies[scaffold]:
                    # Isolate coordinates for mac + IES ("germline" genome) when ies has been retained
                    start_sam, end_sam = seq_mac_ies[1], seq_mac_ies[2]
                    # meaning the aligned sequence overlaps in some way to the gff3 feature
                    # overlaps the ies feature when ies's are retained
                    if (start_gff <= end_sam) and (start_sam <= end_gff):
                        if (end_sam - start_sam) == 0:  # if (end_gff - start_gff) == 0:  # add one?
                            # avoid division by zero
                            print('Division by zero: %s, %s, %d, %d' % (scaffold, gff3line[2], start_gff, end_gff))
                        elif (start_sam <= start_gff) and (end_sam <= end_gff):
                            # 'right' side of aligned read overlaps feature but 'left' does not
                            #  if gff3line[2] == target_feature:
                            d_left_mac_ies[ies_name] = d_left_mac_ies[ies_name] + 1
                        elif (start_sam >= start_gff) and (end_sam >= end_gff):
                            # 'left' side of aligned read overlaps feature but 'right' does not
                            d_right_ies_mac[ies_name] = d_right_ies_mac[ies_name] + 1
                        elif (start_gff <= start_sam) and (end_sam <= end_gff):
                            # if read aligned is inside feature and smaller than feature
                            # pass because i only want reads that overlap the left or right boundary (or both)
                            pass
                        elif (start_sam <= start_gff) and (end_sam >= end_gff):
                            # if read aligned is larger than feature aligned to...
                            # although read completely covers feature, the feature does not completely cover read
                            d_both_mac_ies_mac[ies_name] = d_both_mac_ies_mac[ies_name] + 1
                        else:
                            print('some form of overlap was not accounted for')
                            print('%s, SS: %d, SG: %d, ES: %d, EG: %d' % (scaffold, start_sam, start_gff, end_sam, end_gff))
        count_scaffolds += 1
        if count_scaffolds % 50 == 0:
            print('Processed %d scaffolds\n%s is present in d_mac, dsam_mac, and dsam_mac_ies' % (count_scaffolds, scaffold))

    d_IRS, d_IRS_Alt, d_counts = {}, {}, {}  # classical way of calculating IRS, and alternative method
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
        d_counts[k] = [str(d_left_mac_ies[k]), str(d_right_ies_mac[k]), str(d_both_mac_ies_mac[k]), str(d_mac[k])]
        countIRS += 1
        if countIRS % 5000 == 0:
            print('IRS %d' % countIRS)
            print('%d, %d, %d, %d' % (d_left_mac_ies[k], d_right_ies_mac[k], d_both_mac_ies_mac[k], d_mac[k]))
            print('%s %.2f, %.2f' % (k, d_IRS[k], d_IRS_Alt[k]))

    return d_IRS, d_IRS_Alt, d_counts


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

def read_sam_subsample_coordinates(sam, numreads):
    ''' key is scaffold, value is line.strip() '''
    from random import sample

    print('Reading, sub-sampling, trim by CIGAR, only coordinates,  for sam file: %s' % sam)
    samlist = []
    d = {}
    with open(sam, 'r') as FILE:
        # lines = [line for line in FILE if random() <= .25]
        # print('Number of subsampled lines: %d' % len(lines))
        subsampled_list = {}
        for line in FILE:
            if line[0] != '@':
                samlist.append([line.strip().split('\t')[2], int(line.strip().split('\t')[3]), int(line.strip().split('\t')[3]) + len(line.strip().split('\t')[9]) - 1, line.strip().split('\t')[5]])

        if len(samlist) <= numreads:
            # the total number of reads is less than the desired
            # just pass the dictionary so the return statement can return one
            gate = 0
            subsampled_list = samlist
        else:
            # here we need to randomly select reads from the total reads to reduce to desired number
            # print("choosing 2 random items from a dictionary using sample method ", random.sample(d.items(), k=2))
            # {k: v for (k, v) in random.sample(d.items(), k=3)}
            # random.sample returns list of tuples, so i need to convert back to dictionary
            gate = 1
            subsampled_list = sample(samlist, numreads)

        for line in subsampled_list:  # line is a list of int(start coord), int(end coord), and str(CIGAR)
            if line[0][-len('_with_ies'):] == '_with_IES':
                key = line[0][:-len('_with_IES')]
            else:
                key = line[0]
            # key = str(scaffold), value = list of all [str(scaffold), int(start coord), int(end coord), str(CIGAR)]
            d.setdefault(key, []).append(line)

    print('Number of scaffolds: %d' % len(list(d.keys())))
    print(list(d.keys())[:10])
    print('Desired Number of reads: %d' % numreads)
    if gate == 0:
        print('Number of sub-sampled reads = %d' % len(samlist))
    elif gate == 1:
        print('Number of sub-sampled reads = %d' % numreads)
    samlist = None  # free up memory
    subsampled_list = None  # free up memory

    return d


def blocks(files, size=65536):
    while True:
        b = files.read(size)
        if not b:
            break
        yield b


def read_sam_subsample_coordinates_cigar(sam, numreads):
    ''' key is scaffold, value is line.strip() '''
    from random import sample
    from random import random
    import re

    # count number of lines in sam file
    with open(sam, "r", encoding="utf-8", errors='ignore') as f:
        linecount = sum(bl.count("\n") for bl in blocks(f))

    if linecount > numreads:
        freq = float(numreads) / float(linecount)
    else:
        freq = 1.0

    print('Reading, sub-sampling, trim by CIGAR, only coordinates,  for sam file: %s' % sam)
    print('Total number of lines = %d. Sampling frequency: %f' % (linecount, freq))
    samlist = []
    d = {}
    with open(sam, 'r') as FILE:
        # lines = [line for line in FILE if random() <= .25]
        # print('Number of subsampled lines: %d' % len(lines))
        subsampled_list = {}
        for line in FILE:
            if (line[0] != '@') and (random() <= freq):  # this will only sample numlines = linecount * freq
                CIGAR = line.strip().split('\t')[5]
                # now check CIGAR string if deletion
                # re.findall() returns all the non-overlapping matches of patterns
                # in a string as a list of strings
                cigar = re.findall(r'\d+[A-Z]', CIGAR)
                # then the ends match and no hard or soft clipping
                if cigar != []:  # if there was a match # sometimes CIGAR string == * # this skips unmapped reads
                    if (cigar[0][-1] == 'M') and (cigar[-1][-1] == 'M') and ('D' not in CIGAR) and ('I' not in CIGAR) and ('N' not in CIGAR):
                        # I need to collect all relevant reads to subsample them randomly
                        samlist.append([line.strip().split('\t')[2], int(line.strip().split('\t')[3]), int(line.strip().split('\t')[3]) + len(line.strip().split('\t')[9]) - 1, CIGAR])

        #if len(samlist) <= numreads:
        #    # the total number of reads is less than the desired
        #    # just pass the dictionary so the return statement can return one
        #    gate = 0
        #    subsampled_list = samlist
        #else:
        #    # here we need to randomly select reads from the total reads to reduce to desired number
        #    # print("choosing 2 random items from a dictionary using sample method ", random.sample(d.items(), k=2))
        #    # {k: v for (k, v) in random.sample(d.items(), k=3)}
        #    # random.sample returns list of tuples, so i need to convert back to dictionary
        #    gate = 1
        #    subsampled_list = sample(samlist, numreads)

        #for line in subsampled_list:  # line is a list of int(start coord), int(end coord), and str(CIGAR)
        for line in samlist:  # line is a list of int(start coord), int(end coord), and str(CIGAR)
            if line[0][-len('_with_ies'):] == '_with_IES':
                key = line[0][:-len('_with_IES')]
            else:
                key = line[0]
            # key = str(scaffold), value = list of all [str(scaffold), int(start coord), int(end coord), str(CIGAR)]
            d.setdefault(key, []).append(line)

    print('Number of scaffolds: %d' % len(list(d.keys())))
    print(list(d.keys())[:10])
    print('Desired Number of reads: %d' % numreads)
    print('Number of sub-sampled reads = %d' % len(samlist))
    samlist = None  # free up memory

    return d


def main(GFF3file, samfiles_mac, samfiles_mac_ies, numreads=15000000, features=['all'], CIGARstring=True):
    import os

    # dTE, dTEhead = read_sam(TEfile, positions=[0])  # sam file of aligned TEs

    dgff3, dgff3head = read_gff3(GFF3file, features)
    # print(len(dmRNA['scaffold51_74']))
    # print(dmRNA['scaffold51_74'])

    for sam_mac, sam_mac_ies in zip(samfiles_mac, samfiles_mac_ies):
        # dsam, dsamhead = read_sam(sam,
        #                           positions=[0])  # sam files of aligned sRNAs # dsRNA stands for dictionary small RNA
        # dsam, dsamhead = read_sam_subsample(sam, freq, positions=[0])

        if CIGARstring == True:
            dsam_mac = read_sam_subsample_coordinates_cigar(sam_mac, numreads)
            dsam_mac_ies = read_sam_subsample_coordinates_cigar(sam_mac_ies, numreads)
        elif CIGARstring == False:
            dsam_mac = read_sam_subsample_coordinates(sam_mac, numreads)
            dsam_mac_ies = read_sam_subsample_coordinates(sam_mac_ies, numreads)

        # dmax = {**dsamhead, **dgff3head}  # update dsRNAhead with values from dTEhead

        print('##### Starting IRS Calculations #####')
        d_IRS, d_IRS_Alt, d_counts = IRS(dgff3, dsam_mac, dsam_mac_ies)  # or dsam for all reads
        print('Number of IRS values: %d' % len(d_IRS.keys()))
        print(list(d_IRS.items())[:20])

        path, file = os.path.split(sam_mac)

        outpath = os.path.join(path, '.'.join(file.split('.')[:-1] + ['IRS', 'tsv']))
        with open(outpath, 'w') as OUT:
            OUT.write('\n'.join(['\t'.join([k, str(round(v, 2)), '\t'.join(d_counts[k])]) for (k, v) in d_IRS.items()]))
        outpath = os.path.join(path, '.'.join(file.split('.')[:-1] + ['IRS', 'Alternative', 'tsv']))
        with open(outpath, 'w') as OUT:
            OUT.write('\n'.join(['\t'.join([k, str(round(v, 2)), '\t'.join(d_counts[k])]) for (k, v) in d_IRS_Alt.items()]))

        print('###########################################')
        print('###########################################')
        print('###########################################')

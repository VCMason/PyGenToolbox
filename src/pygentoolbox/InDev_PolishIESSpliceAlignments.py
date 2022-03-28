def read_sam_cigar(sam, numreads):
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


def polish_ies_splice_alignments(samfile, dgff3, dref):
    import re

    outsam = '.'.join(samfile.split('.')[:-1] + ['IESssPolished', 'sam'])
    with open(outsam, 'w') as OUT:
        with open(samfile, 'r') as FILE:
            for line in FILE:
                if line[0] == '@':
                    OUT.write(line.strip() + '\n')
                else:
                    CIGAR = line.strip().split('\t')[5]
                    if line.strip().split('\t')[2] == '*':
                        # skip the unmapped reads
                        OUT.write(line.strip() + '\n')
                    elif 'N' in CIGAR:
                        # if N is in CIGAR then there is some splicing
                        # re.findall() returns all the non-overlapping matches of patterns
                        # in a string as a list of strings
                        NM = int(line.strip().split('\t')[16].split(':')[-1]) # meta data NM:i:2 # number of mismatches
                        if NM > 0:
                            # then there are mismatches in the alignment
                            # does the alignment overlap IES feature? If so check and fix alignments at IES boundaries
                            cigar = re.findall(r'\d+[A-Z]', CIGAR)
                            MD = line.strip().split('\t')[17].split(':')[-1] # metadata MD:Z:7T0G139
                            md = re.findall(r'\d+[A-Z]', MD)
                            # then the ends match and no hard or soft clipping
                            #if cigar != []:  # if there was a match # sometimes CIGAR string == * # this skips unmapped reads
                            samstart = int(line.strip().split('\t')[3])
                            samend = int(line.strip().split('\t')[3]) + sum([int(c[:-1]) for c in cigar if (c[-1] != 'I') and (c[-1 != 'H'])])

                            if (cigar[0][-1] == 'M') and (cigar[-1][-1] == 'M') and ('D' not in CIGAR) and ('I' not in CIGAR):
                        else:
                            # then there are no mismatches
                            OUT.write(line.strip() + '\n')
                    else:
                        # if no 'splicing' (denoted as #N in CIGAR)
                        OUT.write(line.strip() + '\n')


def read_fasta(f):
    # f is full path to fasta file
    d = {}
    names = []
    with open(f, 'r') as FILE:
        for line in FILE:
            if line[0] == '>':
                scaffold = line.strip().split()[0][1:]
                d[scaffold] = []
                names.append(scaffold)
            else:
                d.setdefault(scaffold, []).append(line.strip())
    dnew = {k: ''.join(v) for k, v in d.items()}
    d = None
    return dnew, names


def read_gff3(GFF3file, features=['all']):
    # sam file aligned to MacAndIES genome
    # arguements suggested with Hisat2: --no-softclip --pen-noncansplice 0
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
                key = line.strip().split('\t')[0]  # [:-len('_with_IES')]  # remove the extra text from the 'with ies' .gff3 file
                d.setdefault(key, []).append(line.strip().split('\t'))
            elif line.strip().split('\t')[2] in features:  # keep lines only if in features
                key = line.strip().split('\t')[0]  # [:-len('_with_IES')]  # remove the extra text from the 'with ies' .gff3 file
                d.setdefault(key, []).append(line.strip().split('\t'))
    print('Number of scaffolds: %d' % len(list(d.keys())))
    print(list(d.keys())[:10])
    return d, dhead


def main(gff3file, samfilelist, refgenomefile, features=['all']):

    dgff3, dgff3head = read_gff3(GFF3file, features)

    dref, refscaffnames = read_fasta(refgenomefile)

    for samfile in samfilelist:
        polish_ies_splice_alignments(samfile, dgff3, dref)
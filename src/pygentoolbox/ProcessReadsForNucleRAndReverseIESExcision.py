############################################################################################
##########        Input are reads aligned to Mac+IES genome sam file              ##########
############################################################################################


def read_sam_add_n_write_out(sam):
    import re

    '''When read is spliced with 'N' in the CIGAR string I add N's to sequence in sam file'''
    '''Assuming splicing (not deletions) are results of reads spanning IES'''
    outfile = '.'.join(sam.split('.')[:-1] + ['ReversedExcision', 'sam'])
    with open(outfile, 'w') as OUT:
        OUT.write('')
    with open(outfile, 'a') as OUT:
        print('Reading sam file: %s' % sam)
        with open(sam, 'r') as FILE:
            for line in FILE:
                if line[0] == '@':
                    OUT.write(line.strip() + '\n')
                else:
                    #print(line)
                    CIGAR = line.strip().split('\t')[5]
                    cigar = re.findall(r'\d+[A-Z]', CIGAR)
                    if cigar == []:  # if there was a match # sometimes CIGAR string == * # this skips unmapped reads
                        OUT.write(line.strip() + '\n')
                    else:   # then read should be mapped
                        if 'N' in CIGAR:  # then the alignment of read was spliced.
                            seq = line.strip().split('\t')[9]
                            qual = line.strip().split('\t')[10]
                            newseq = []
                            newqual = []
                            newcigar = []
                            position = 0
                            for i in cigar:
                                if i[-1] in ['M', 'I', 'S', 'H']:
                                    x = int(i[:-1])
                                    newseq.append(seq[position:position + x])
                                    newqual.append(qual[position:position + x])
                                    newcigar.append(i)
                                    position += x
                                elif i[-1] in ['D', 'P']:
                                    # skip deletions and 'padding' because bases don't exist in sequence
                                    newcigar.append(i)
                                elif i[-1] == 'N':
                                    x = int(i[:-1])
                                    newseq.append('N'*x)
                                    newqual.append('F'*x)
                                    newcigar.append('%sM' % i[:-1])
                            newline = '\t'.join(line.strip().split('\t')[:5] + [''.join(newcigar)] + line.strip().split('\t')[6:9] + [''.join(newseq), ''.join(newqual)] + line.strip().split('\t')[11:])
                            OUT.write(newline + '\n')
                        else:
                            OUT.write(line.strip() + '\n')

    return


def fragment_midpoints_to_bed(samfile, dfive, dss, midwidth=50):
    outbedfile = '.'.join(samfile.split('.')[:-1] + [identifier, 'sam'])

    with open(outbedfile, 'w') as OUT:
        OUT.write('')
    with open(outbedfile, 'a') as OUT:
        with open(samfile, 'r') as SAM:
            for line in SAM:
                if line[0] != '@':
                    readname = line.strip().split('\t')[0]
                    try:
                        dbed[readname]
                    except:
                        # then no entry yet
                        if dfive[readname][0] < dfive[readname][1]:
                            start = 0
                            end = 1
                        elif dfive[readname][0] > dfive[readname][1]:
                            start = 1
                            end = 0
                        if sum(dss[readname]) == 0:
                            # then no splicing, so trim read to midpoint
                            reftlen = abs(dfive[end] - dfive[start]) + 1
                            trimlength = (reftlen - 50) / 2
                            startcoord = int(dfive[start] + trimlength)
                            endcoord = int(dfive[end] + trimlength)
                            OUT.write('\t'.join([line.strip().split('\t')[2], startcoord, endcoord, line.strip().split('\t')[0], '*']) + '\n')
                        elif sum(dss[readname]) > 0:
                            # then fragment was spliced, so do not trim the ends of the fragment to the midpoint
                            OUT.write('\t'.join([line.strip().split('\t')[2], dfive[start], dfive[end], line.strip().split('\t')[0], '*']) + '\n')
                    else:
                        # already done for fragment
                        pass


def filter_sam_by_read_name(samfile, dreadnames, identifier):
    # samfile is samfile name
    # readnames is list of readnames to keep in sam file output
    # identifier is unique identifier for filtration method/file
    outsamfile = '.'.join(samfile.split('.')[:-1] + [identifier, 'sam'])
    with open(outsamfile, 'w') as OUTSAM:
        OUTSAM.write('')
    with open(outsamfile, 'a') as OUTSAM:
        with open(samfile, 'r') as SAM:
            for line in SAM:
                if line[0] == '@':
                    OUTSAM.write(line.strip() + '\n')
                else:
                    try:
                        dreadnames[line.strip().split('\t')[0]]
                    except:
                        # there is not entry with that readname
                        pass
                    else:
                        # then readname present in dictionary
                        OUTSAM.write(line.strip() + '\n')
    print('Filtered sam file output to: %s' % outsamfile)

    return outsamfile


def get_lists_of_deletion_insertion_splice_lengths_cigar(CIGAR):
    import re

    splices = []
    insertions = []
    deletions = []
    cigar = re.findall(r'\d+[A-Z]', CIGAR)

    if cigar == []:  # if there was a match # sometimes CIGAR string == * # this skips unmapped reads
        print(f'Provided CIGAR string: {CIGAR} does not match CIGAR pattern \\d+[A-Z]')
    else:  # then read should be mapped
        for i in cigar:
            if i[-1] == 'N':
                splices.append(int(i[:-1]))
            elif i[-1] == 'I':
                insertions.append(int(i[:-1]))
            elif i[-i] == 'D':
                deletions.append(int(i[:-1]))
            else:
                pass
    return splices, insertions, deletions


def get_rightmost_reference_based_alignment_coordinate(CIGAR, leftmost_coordinate):
    import re
    cigar = re.findall(r'\d+[A-Z]', CIGAR)
    if cigar == []:  # if there was a match # sometimes CIGAR string == * # this skips unmapped reads
        print(f'Provided CIGAR string: {CIGAR} does not match CIGAR pattern \\d+[A-Z]')
        rightmost_position = 0  # assumes unmapped read
    else:  # then read should be mapped
        rightmost_position = leftmost_coordinate - 1  # subtract 1 because leftmost base is 1-based
        for i in cigar:
            if i[-1] in ['M', 'N', 'D']:
                rightmost_position += int(i[:-1])
            elif i[-1] in ['I', 'S', 'H', 'P']:
                pass
            else:
                pass
    return rightmost_position


def reference_based_tlen(sam):
    ''' assumes all reads are paired and mapped in proper pair -f 3 '''
    ''' assumes paired reads are read mapped, mate mapped, and only mapped to primary alignment -F 268 '''
    ''' assumes all reads are mapped in proper pair '''
    ''' input file is .sam file, output is dictionary: key is read name, value is ref based tlen value'''
    ''' calculates the fragment length between paired reads according to reference positions'''
    ''' so if the read is spliced this would be included in the fragment length (tlen) calculation'''
    ''' always rightmost 5' minus leftmost 5' reference coordinate'''
    ''' example: '''
    ''' if an IES was excised then the tlen value reported would include the IES length (as if it wasn't excised)'''
    d = {}
    # d is dictionary, key is read name, value is list of length 2
    # value = [5' ref-based coordinate R1, 5' ref-based coordinate R2]
    with open(sam, 'r') as FILE:
        for line in FILE:
            flag = int(line.strip().split('\t')[2])
            binaryflag = f'{flag:012b}'  # ex: flag = 128 would return: 000010000000
            try:
                d[line.strip().split('\t')[0]]
            except:
                # then dictionary entry for read name does not exist
                if binaryflag[-7] == 1:
                    # then it is first read in pair (R1) and on
                    if binaryflag[-5] == 0:
                        # read is on forward strand
                        # so add the leftmost position (5') end to first elem of list
                        d[line.strip().split('\t')[0]] = [int(line.strip().split('\t')[3]), 0]
                    else:
                        # if 1 then it is on the reverse strand
                        # since it is on reverse strand we need the rightmost coordinate to get the 5' end
                        # ïnt(line.strip().split('\t')[3]) is leftmost aligned position (excludes soft/hard clipping)
                        CIGAR = line.strip().split('\t')[5]
                        rightmost_coord = get_rightmost_reference_based_alignment_coordinate(CIGAR, int(line.strip().split('\t')[3]))
                        d[line.strip().split('\t')[0]] = [rightmost_coord, 0]
                elif binaryflag[-8] == 1:
                    # then it is second read in pair (R2)
                    if binaryflag[-5] == 0:
                        # read is on forward strand
                        # so add the leftmost position (5') end to first elem of list
                        d[line.strip().split('\t')[0]] = [0, int(line.strip().split('\t')[3])]
                    else:
                        # if 1 then it is on the reverse strand
                        # since it is on reverse strand we need the rightmost coordinate to get the 5' end
                        # ïnt(line.strip().split('\t')[3]) is leftmost aligned position (excludes soft/hard clipping)
                        CIGAR = line.strip().split('\t')[5]
                        rightmost_coord = get_rightmost_reference_based_alignment_coordinate(CIGAR, int(line.strip().split('\t')[3]))
                        d[line.strip().split('\t')[0]] = [0, rightmost_coord]
                else:
                    print('read is not first or second in pair.. sketchy')
            else:
                # then dictionary entry for read name exists
                if binaryflag[-7] == 1:
                    # then it is first read in pair (R1) and on
                    if binaryflag[-5] == 0:
                        # read is on forward strand
                        # so add the leftmost position (5') end to first elem of list
                        d[line.strip().split('\t')[0]][0] = int(line.strip().split('\t')[3])
                    else:
                        # if 1 then it is on the reverse strand
                        # since it is on reverse strand we need the rightmost coordinate to get the 5' end
                        # ïnt(line.strip().split('\t')[3]) is leftmost aligned position (excludes soft/hard clipping)
                        CIGAR = line.strip().split('\t')[5]
                        rightmost_coord = get_rightmost_reference_based_alignment_coordinate(CIGAR, int(line.strip().split('\t')[3]))
                        d[line.strip().split('\t')[0]][0] = rightmost_coord
                elif binaryflag[-8] == 1:
                    # then it is second read in pair (R2)
                    if binaryflag[-5] == 0:
                        # read is on forward strand
                        # so add the leftmost position (5') end to first elem of list
                        d[line.strip().split('\t')[0]][1] = int(line.strip().split('\t')[3])
                    else:
                        # if 1 then it is on the reverse strand
                        # since it is on reverse strand we need the rightmost coordinate to get the 5' end
                        # ïnt(line.strip().split('\t')[3]) is leftmost aligned position (excludes soft/hard clipping)
                        CIGAR = line.strip().split('\t')[5]
                        rightmost_coord = get_rightmost_reference_based_alignment_coordinate(CIGAR, int(line.strip().split('\t')[3]))
                        d[line.strip().split('\t')[0]][1] = rightmost_coord
                else:
                    print('read is not first or second in pair.. sketchy')
    # dreftlen abs() of reference based fragment length: rightmost 5' minus leftmost 5' reference coordinate
    dreftlen = {k: abs(v[1]-v[2]) for k, v in d.items()}

    return d, dreftlen


def tlen(sam):
    ''' assumes all reads are paired and mapped in proper pair -f 3 '''
    ''' assumes paired reads are read mapped, mate mapped, and only mapped to primary alignment -F 268 '''
    ''' input file is .sam file, output is dictionary: key is read name, value is ref based tlen value'''
    ''' calculates the actual fragment length between paired reads (subtracts spliced bases etc.)'''
    ''' so if the read is spliced this would be removed in the fragment length (tlen) calculation'''
    ''' always rightmost 5' minus leftmost 5' reference coordinate, but also subtracting splice sites (when appropriate)'''
    ''' all splice sites are subtracted if the reads do not overlap'''
    ''' if paired reads overlap and splice sites are in one read only it subtracts all splice sites'''
    ''' if paired reads overlap and splice sites are in both reads it subtracts:'''
    ''' 1.) the read with largest number of splicing events (when comparing R1 and R2)'''
    ''' 2.) if the read with fewer splice events is of unique length(s) compared to read with more splice events'''
    ''' (assumes reads are spliced properly, assumes reads with splicing of unique length cover different features'''
    ''' ex: if an IES was excised then the tlen value reported would notinclude the length of the IES (as if it was excised)'''

    import re
    import collections

    # read names of all sequences
    names = []
    d = {}
    # d is dictionary, key is read name, value is list of length 2
    # value = [5' ref-based coordinate R1, 5' ref-based coordinate R2]
    d3prime = {}
    # d is dictionary, key is read name, value is list of length 2
    # value = [3' ref-based coordinate R1, 3' ref-based coordinate R2]
    dss = {}
    # dss is dictionary, key is read name, value is list of splice site length to be subtracted from d[key]
    # keep all unique length splices, inserts, dels
    # if multiple splices, inserts, or dels of same length are present in the same read, i.e...
    ## two 27bp splices make sure there are two splices of the same length in the list
    ### if overlapping reads cover two DIFFERNT IESs of same length one will NOT be accounted for.
    di = {}
    # di is dictionary, key is read name, value is list of insertion lengths to be added to d[key]
    dd = {}
    # dd is dictionary, key is read name, value is list of deletion lengths to be subtracted from d[key]
    with open(sam, 'r') as FILE:
        for line in FILE:
            if line[0] != '@':
                names.append(readname)
                readname = line.strip().split('\t')[0]
                flag = int(line.strip().split('\t')[2])
                binaryflag = f'{flag:012b}'  # ex: flag = 128 would return: 000010000000
                CIGAR = line.strip().split('\t')[5]
                splices, insertions, deletions = get_lists_of_deletion_insertion_splice_lengths_cigar(CIGAR)
                #if ('N' in CIGAR) or ('I' in CIGAR) or ('D' in CIGAR):
                    # cigar = re.findall(r'\d+[A-Z]', CIGAR)
                try:
                    d[readname]
                except:
                    # then dictionary entry for read name does not exist
                    dss[readname] = splices
                    di[readname] = insertions
                    dd[readname] = deletions
                    if binaryflag[-7] == 1:
                        # then it is first read in pair (R1) and on
                        if binaryflag[-5] == 0:
                            # read is on forward strand
                            # so add the leftmost position (5') end to first elem of list
                            d[readname] = [int(line.strip().split('\t')[3]), 0]
                            rightmost_coord = get_rightmost_reference_based_alignment_coordinate(CIGAR, int(line.strip().split('\t')[3]))
                            d3prime[readname] = [rightmost_coord, 0]
                        else:
                            # if 1 then it is on the reverse strand
                            # since it is on reverse strand we need the rightmost coordinate to get the 5' end
                            # ïnt(line.strip().split('\t')[3]) is leftmost aligned position (excludes soft/hard clipping)
                            rightmost_coord = get_rightmost_reference_based_alignment_coordinate(CIGAR, int(line.strip().split('\t')[3]))
                            d[readname] = [rightmost_coord, 0]
                            d3prime[readname] = [int(line.strip().split('\t')[3]), 0]
                    elif binaryflag[-8] == 1:
                        # then it is second read in pair (R2)
                        if binaryflag[-5] == 0:
                            # read is on forward strand
                            # so add the leftmost position (5') end to first elem of list
                            d[readname] = [0, int(line.strip().split('\t')[3])]
                            rightmost_coord = get_rightmost_reference_based_alignment_coordinate(CIGAR, int(line.strip().split('\t')[3]))
                            d3prime[readname] = [0, rightmost_coord]
                        else:
                            # if 1 then it is on the reverse strand
                            # since it is on reverse strand we need the rightmost coordinate to get the 5' end
                            # ïnt(line.strip().split('\t')[3]) is leftmost aligned position (excludes soft/hard clipping)
                            rightmost_coord = get_rightmost_reference_based_alignment_coordinate(CIGAR, int(line.strip().split('\t')[3]))
                            d[readname] = [0, rightmost_coord]
                            d3prime[readname] = [0, int(line.strip().split('\t')[3])]
                    else:
                        print('read is not first or second in pair.. sketchy')
                else:
                    # then dictionary entry for read name exists
                    if binaryflag[-7] == 1:
                        # then it is first read in pair (R1) and on
                        if binaryflag[-5] == 0:
                            # read is on forward strand
                            # so add the leftmost position (5') end to first elem of list
                            d[readname][0] = int(line.strip().split('\t')[3])
                            rightmost_coord = get_rightmost_reference_based_alignment_coordinate(CIGAR, int(line.strip().split('\t')[3]))
                            d3prime[readname][0] = rightmost_coord
                        else:
                            # if 1 then it is on the reverse strand
                            # since it is on reverse strand we need the rightmost coordinate to get the 5' end
                            # ïnt(line.strip().split('\t')[3]) is leftmost aligned position (excludes soft/hard clipping)
                            CIGAR = line.strip().split('\t')[5]
                            rightmost_coord = get_rightmost_reference_based_alignment_coordinate(CIGAR, int(line.strip().split('\t')[3]))
                            d[readname][0] = rightmost_coord
                            d3prime[readname][0] = int(line.strip().split('\t')[3])
                    elif binaryflag[-8] == 1:
                        # then it is second read in pair (R2)
                        if binaryflag[-5] == 0:
                            # read is on forward strand
                            # so add the leftmost position (5') end to first elem of list
                            d[readname][1] = int(line.strip().split('\t')[3])
                            rightmost_coord = get_rightmost_reference_based_alignment_coordinate(CIGAR, int(line.strip().split('\t')[3]))
                            d3prime[readname][1] = rightmost_coord
                        else:
                            # if 1 then it is on the reverse strand
                            # since it is on reverse strand we need the rightmost coordinate to get the 5' end
                            # ïnt(line.strip().split('\t')[3]) is leftmost aligned position (excludes soft/hard clipping)
                            CIGAR = line.strip().split('\t')[5]
                            rightmost_coord = get_rightmost_reference_based_alignment_coordinate(CIGAR, int(line.strip().split('\t')[3]))
                            d[readname][1] = rightmost_coord
                            d3prime[readname][1] = int(line.strip().split('\t')[3])
                    else:
                        print('read is not first or second in pair.. sketchy')

                    # we have 5' coords to R1 and R2, and 3' coord for R1 and R2
                    # check if reads overlap, if they do only keep splice, insert, del if lengths are unique
                    if d[readname][0] < d3prime[readname][0]:
                        # then read R1 is mapped to forward strand
                        # if (5'_R1 <= 5'_R2) and (3'_R2 <= 3'_R1):
                        # R1: --------->
                        # R2:     <---------
                        if (d[readname][0] <= d[readname][1]) and (d3prime[readname][1] <= d3prime[readname][0]):
                            # then reads overlap to some degree and could share same splice sites/deletions/insertions
                            try:
                                dss[readname]
                            except:
                                # then this is the first read of the pair to be processed (NOT binaryflag[-7] == 1...)
                                dss[readname] = dss.setdefault(readname, []) + splices
                                di[readname] = di.setdefault(readname, []) + insertions
                                dd[readname] = dd.setdefault(readname, []) + deletions
                            else:
                                # then this is the second read of the pair to be processed
                                # then add all unique length splices, inserts, dels
                                # if multiple splices, inserts, or dels of same length are present in the same read
                                ## i.e. two 27bp splices make sure there are two splices of the same length in the list
                                # list example: [27, 101, 27, 52]
                                # counter object is dictionary Counter({27: 2, 101: 1, 52: 1})
                                # splice length is key and value is how often length is present in list
                                splicecounts = collections.Counter(splices)
                                for splicelength in splices:
                                    if splicelength not in dss[readname]:
                                        dss[readname] = dss[readname] + [splicelength]
                                    elif len(set(dss[readname])) < splicecounts[splicelength]:
                                        dss[readname] = dss[readname] + [splicelength]
                                insertcounts = collections.Counter(insertions)
                                for insertlength in insertions:
                                    if insertlength not in di[readname]:
                                        di[readname] = di[readname] + [insertlength]
                                    elif len(set(di[readname])) < insertcounts[insertlength]:
                                        di[readname] = di[readname] + [insertlength]
                                delcounts = collections.Counter(deletions)
                                for deletionlength in deletions:
                                    if deletionlength not in dd[readname]:
                                        dd[readname] = dd[readname] + [deletionlength]
                                    elif len(set(dd[readname])) < delcounts[deletionlength]:
                                        dd[readname] = dd[readname] + [deletionlength]
                        else:
                            # then reads do not overlap
                            dss[readname] = dss.setdefault(readname, []) + splices
                            di[readname] = di.setdefault(readname, []) + insertions
                            dd[readname] = dd.setdefault(readname, []) + deletions
                    else:
                        # then read R1 must be mapped to reverse strand (assumed proper paired reads)
                        # if (5'_R2 <= 5'_R1) and (3'_R1 <= 3'_R2):
                        # R2: --------->
                        # R1:     <---------
                        if (d[readname][1] <= d[readname][0]) and (d3prime[readname][0] <= d3prime[readname][1]):
                            # then reads overlap to some degree and could share same splice sites/deletions/insertions
                            try:
                                dss[readname]
                            except:
                                # then this is the first read of the pair to be processed (NOT binaryflag[-7] == 1...)
                                dss[readname] = dss.setdefault(readname, []) + splices
                                di[readname] = di.setdefault(readname, []) + insertions
                                dd[readname] = dd.setdefault(readname, []) + deletions
                            else:
                                # then this is the second read of the pair to be processed
                                # then add all unique length splices, inserts, dels
                                # if multiple splices, inserts, or dels of same length are present in the same read
                                ## i.e. two 27bp splices make sure there are two splices of the same length in the list
                                # list example: [27, 101, 27, 52]
                                # counter object is dictionary Counter({27: 2, 101: 1, 52: 1})
                                # splice length is key and value is how often length is present in list
                                splicecounts = collections.Counter(splices)
                                for splicelength in splices:
                                    if splicelength not in dss[readname]:
                                        dss[readname] = dss[readname] + [splicelength]
                                    elif len(set(dss[readname])) < splicecounts[splicelength]:
                                        dss[readname] = dss[readname] + [splicelength]
                                insertcounts = collections.Counter(insertions)
                                for insertlength in insertions:
                                    if insertlength not in di[readname]:
                                        di[readname] = di[readname] + [insertlength]
                                    elif len(set(di[readname])) < insertcounts[insertlength]:
                                        di[readname] = di[readname] + [insertlength]
                                delcounts = collections.Counter(deletions)
                                for deletionlength in deletions:
                                    if deletionlength not in dd[readname]:
                                        dd[readname] = dd[readname] + [deletionlength]
                                    elif len(set(dd[readname])) < delcounts[deletionlength]:
                                        dd[readname] = dd[readname] + [deletionlength]
                        else:
                            # then reads do not overlap
                            dss[readname] = dss.setdefault(readname, []) + splices
                            di[readname] = di.setdefault(readname, []) + insertions
                            dd[readname] = dd.setdefault(readname, []) + deletions

    # dreftlen abs() of reference based fragment length: rightmost 5' minus leftmost 5' reference coordinate plus one
    # if (v[0] != 0) and (v[1] != 0) if either value is zero then assume one of the reads did not map
    ## and if so do not calculate fragment length
    dreftlen = {k: abs(v[0]-v[1]) + 1 for k, v in d.items() if (v[0] != 0) and (v[1] != 0)}
    # dtlen is actual fragment length (not reference based)
    # rightmost 5' coordinate - leftmost 5' reference coordinate - sum(splices) - sum(deletions) + sum(insertions)
    # for overlapping reads, only unique length deletions, insertions, splices are accounted for
    dtlen = {k: abs(v[0]-v[1]) + 1 - sum(dss[k]) - sum(dd[k]) + sum(di[k]) for k, v in d.items() if (v[0] != 0) and (v[1] != 0)}

    return d, d3prime, dreftlen, dtlen, dss, dd, di


def process_reads_nucler_pear_assembled_se(samfile, lowlimit=100, uplimit=200, midwidth=50):
    ''' assume fragment length = length of sequence b/c all reads were assembled together with PEAR'''
    ''' read sam, remove reads > fragment length, ... '''
    ''' trim reads <= frag len to length=trim centered at midpoint of fragment IF read is NOT spliced over IES'''
    ''' if read is spliced over IES do not remove reads based upon fragment length and do not trim '''
    ''' if over IES 'reverse' the excision even by adding NNN's to sequence and CIGAR string'''
    ''' reads are not trimmed over IESs'''
    ''' this function should be run on sam file before ReadBAM.R and skip the processReads function in NucleR.R'''

    print('Starting fragment length calculations')
    # assumes all reads are paired and mapped in proper pair -f 3
    # assumes paired reads are read mapped, mate mapped, and only mapped to primary alignment -F 268
    # input file is .sam file, output is dictionary: key is read name, value is ref based tlen value
    # dfiveprime, d3prime, dreftlen, dtlen, dsplicesites, ddeletions, dinsertions
    dtlen, dss, dd, di = tlen_as_seq_len(samfile)
    print('Finished fragment length calculations')

    # determine which reads have splice sites (assuming b/c of IES)
    # if the length of splice sites sum to zero.. then no splice sites, else, it was spliced
    splicereadnames = [k for k, v in dss.items() if sum(v) > 0]

    print('Trimming reads based upon fragment length')
    # trim the reads to those with fragment length >=100 and <=200
    dfragtrim = {k: v for k, v in dtlen.items() if (v >= lowlimit) and (v <= uplimit)}
    print('Finished trimming reads by fragment length (not reference based)')

    print('Writing out sam file with pairs trimmed by fragment length')
    # output
    outsamfile = filter_sam_by_read_name(samfile, dfragtrim, '%dto%d' % (lowlimit, uplimit))
    print('Finished writing sam file, keeping only fragments between %d and %d in length' % (lowlimit, uplimit))

    #  trim to middle of each fragment (if not spliced), do not trim spliced reads
    print('Trimming non-spliced reads to midpoint of fragment with width = %d' % midwidth)
    fragment_midpoints_to_bed(outsamfile, dfive, dss, midwidth)

    # convert bedfile to .RData file (formatted as Granges object)


    print('###########################################')
    print('###################FIN#####################')
    print('###########################################')


def process_reads_nucler(samfile, lowlimit=100, uplimit=200, midwidth=50):
    ''' read sam, remove reads > fragment length, ... '''
    ''' trim reads <= frag len to length=trim centered at midpoint of fragment IF read is NOT spliced over IES'''
    ''' if read is spliced over IES do not remove reads based upon fragment length and do not trim '''
    ''' if over IES 'reverse' the excision even by adding NNN's to sequence and CIGAR string'''
    ''' reads are not trimmed over IESs'''
    ''' this function should be run on sam file before ReadBAM.R and skip the processReads function in NucleR.R'''

    print('Starting fragment length calculations')
    # assumes all reads are paired and mapped in proper pair -f 3
    # assumes paired reads are read mapped, mate mapped, and only mapped to primary alignment -F 268
    # input file is .sam file, output is dictionary: key is read name, value is ref based tlen value
    # dfiveprime, d3prime, dreftlen, dtlen, dsplicesites, ddeletions, dinsertions
    dfive, dthree, dreftlen, dtlen, dss, dd, di = tlen(samfile)
    print('Finished fragment length calculations')

    # determine which reads have splice sites (assuming b/c of IES)
    # if the length of splice sites sum to zero.. then no splice sites, else, it was spliced
    splicereadnames = [k for k, v in dss.items() if sum(v) > 0]

    print('Trimming reads based upon fragment length')
    # trim the reads to those with fragment length >=100 and <=200
    dfragtrim = {k: v for k, v in dtlen.items() if (v >= lowlimit) and (v <= uplimit)}
    print('Finished trimming reads by fragment length (not reference based)')

    print('Writing out sam file with pairs trimmed by fragment length')
    # output
    outsamfile = filter_sam_by_read_name(samfile, dfragtrim, '%dto%d' % (lowlimit, uplimit))
    print('Finished writing sam file, keeping only fragments between %d and %d in length' % (lowlimit, uplimit))

    #  trim to middle of each fragment (if not spliced), do not trim spliced reads
    print('Trimming non-spliced reads to midpoint of fragment with width = %d' % midwidth)
    fragment_midpoints_to_bed(outsamfile, dfive, dss, midwidth)

    # convert bedfile to .RData file (formatted as Granges object)


    print('###########################################')
    print('###################FIN#####################')
    print('###########################################')

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


def record_read_positions(bamfile, MQ, dreads={}):
    #bamfile is full path to bam file
    #MQ = int() = minimum mapping quality threshold
    #dreads is class object utilizing
    # if f'{flag:012b}'[-5] == '1':  # if 1 then it is reverse # '000000010000' # 16 # then the read is mapped in reverse
    import pysam

    # print(f'{bamfile}')
    with pysam.AlignmentFile(bamfile, 'rb') as FILE:
        for line in FILE:
            line = line.tostring()  # changes the pysam formatting stuff to the actual format of the sam file
            if line[0] == '@':
                pass
            elif int(line.strip().split('\t')[4]) >= MQ:
                flag = int(line.strip().split('\t')[1])  # sam flag
                # if () and (f'{flag:012b}'[-9] == '0')
                if f'{flag:012b}'[-3] != '1':  # if the read is not unmapped # '000000000100' == 4
                    if f'{flag:012b}'[-5] == '1':  # if 1 then it is reverse # '000000010000' == 16 read mapped reverse
                        # key is read name, value is scaffold name and then 5' position of read mapped
                        # read name is reduced to final three fields which should be a unique identifier
                        # scaffold name is reduced to only scaffold number
                        # '_'.join(line.strip().split('\t')[0].split(':')[-3:])
                        fiveprimecoord = get_rightmost_reference_based_alignment_coordinate(
                            line.strip().split('\t')[5],
                            int(line.strip().split('\t')[3]))
                        dreads.setdefault(line.strip().split('\t')[0], []).append(
                            [line.strip().split('\t')[2], str(fiveprimecoord)])
                        # str(int(line.strip().split('\t')[3]) + len(line.strip().split('\t')[9]))
                    else:
                        # key is read name, value is scaffold name and then 5' position of read mapped
                        dreads.setdefault(line.strip().split('\t')[0], []).append(
                            [line.strip().split('\t')[2], line.strip().split('\t')[3]])

    return dreads
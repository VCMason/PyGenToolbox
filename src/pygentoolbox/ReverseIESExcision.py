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


def main(samfile):

    read_sam_add_n_write_out(samfile)

    print('###########################################')
    print('###################FIN#####################')
    print('###########################################')

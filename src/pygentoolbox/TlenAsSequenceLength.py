
def add_tlen_as_seq_length_and_filter(samfile, outsam, low=0, high=0):
    # no filter applied if low and high both == 0
    tlens = []
    print('reading sam file: %s' % samfile)
    with open(outsam, 'w') as OUT:
        OUT.write('')
    with open(outsam, 'a') as OUT:
        with open(samfile, 'r') as FILE:
            for line in FILE:
                if line[0] == '@':
                    OUT.write(line.strip() + '\n')
                elif line[0] != '@':
                    seqlen = len(line.strip().split('\t')[9])
                    if (low == 0) and (high == 0):
                        # replace tlen value with sequence length value
                        # assuming PEAR assembled reads so fragment length == sequence length
                        OUT.write('\t'.join(line.strip().split('\t')[:8] + [str(seqlen)] + line.strip().split('\t')[9:]) + '\n')
                        tlens.append(seqlen)
                    elif (seqlen >= low) and (seqlen <= high):
                        # replace tlen value with sequence length value
                        # assuming PEAR assembled reads so fragment length == sequence length
                        OUT.write('\t'.join(line.strip().split('\t')[:8] + [str(seqlen)] + line.strip().split('\t')[9:]) + '\n')
                        tlens.append(seqlen)
                    else:
                        pass

    print('Finished replacing tlen value with sequence length')
    print('output file: %s' % outsam)
    avgseqlen = sum(tlens)/len(tlens)
    print('Mean sequence length: %f' % avgseqlen)
    print('Min sequence length: %d\nMax sequence length: %d' % (min(tlens), max(tlens)))

def main(samfile, outsam, low=0, high=0):
    add_tlen_as_seq_length_and_filter(samfile, outsam, low, high)
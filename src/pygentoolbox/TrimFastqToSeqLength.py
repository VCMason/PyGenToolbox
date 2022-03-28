
def output_file(outfile, outlist):
    with open(outfile, 'w') as OUT:
        OUT.write('\n'.join(outlist))
    print(f'Output file: {outfile}')


def read_fastq_trim_to_sequence_length(f, seqlength=25):
    # d = {}  # fastq names are key, values are sequences as string
    outlist = []
    with open(f, 'r') as FILE:
        for count, line in enumerate(FILE):
            if count % 4 == 0:
                name = line.strip()[1:]
            elif count % 4 == 1:
                seq = line.strip()
            elif (count % 4 == 3) and (len(seq) == seqlength):
                qual = line.strip()
                outseq = f'@{name}\n{seq}\n+\n{qual}'
                outlist.append(outseq)

    return outlist


def main(fastqfile, outfile, seqlength=25):
    outlist = read_fastq_trim_to_sequence_length(fastqfile, seqlength)
    output_file(outfile, outlist)
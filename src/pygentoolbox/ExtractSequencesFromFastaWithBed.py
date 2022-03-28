
def write_out(output, outfilename):
    with open(outfilename, 'w') as OUT:
        OUT.write('\n'.join(output))
    print(f'Output file: {outfilename}')


def extract_seqs_from_fasta_with_bed(dfasta, fastanames, bedlist):
    # bed file coordinates are assumed to be 1 based
    output = []
    for bed in bedlist:
        for name in fastanames:
            if bed[0] == name:
                start = int(bed[1]) - 1
                end = int(bed[2])
                id = bed[3]
                seq = dfasta[name][start:end]
                output.append(f'>{name}:{start+1}:{end}:{id}\n{seq}')
                break

    return output

def read_bedfile_as_list(f, delim='\t'):
    with open(f, 'r') as FILE:
        bedlist = [line.strip().split(delim) for line in FILE]

    return bedlist


def read_fasta_as_dict(f):
    d = {}  # fasta names are key, values are sequences as string
    namelist = []
    with open(f, 'r') as FILE:
        for line in FILE:
            if line[0] == '>':
                if ' ' in line:
                    name = line.strip().split()[0][1:]  # line.strip().split()[0][1:-len('_with_IES')]
                    namelist.append(name)
                    d[name] = []
                else:
                    name = line.strip()[1:]
                    namelist.append(name)
                    d[name] = []
            elif line.strip() != '':  # else: # trying to prevent some crap happening on the last line
                d[name].append(line.strip())
    for name in namelist:
        d[name] = ''.join(d[name])  # join list of partial sequences. Useful if interleaved fasta

    return d, namelist


def main(fastafile, bedfile, outfilename):
    dfasta, fastanames = read_fasta_as_dict(fastafile)
    bedlist = read_bedfile_as_list(bedfile)

    output = extract_seqs_from_fasta_with_bed(dfasta, fastanames, bedlist)
    write_out(output, outfilename)
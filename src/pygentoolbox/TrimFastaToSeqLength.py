
def output_file(outfile, d, names):
    output = [f'>{n}\n{d[n]}' for n in names]
    with open(outfile, 'w') as OUT:
        OUT.write('\n'.join(output))
    print(f'Output file: {outfile}')


def trim_fasta_to_sequence_length(dfasta, names, seqlength):
    dtrim = {}
    trimnames = []
    for n in names:
        if len(dfasta[n]) == seqlength:
            dtrim[n] = dfasta[n]
            trimnames.append(n)
    return dtrim, trimnames


def read_fasta_as_dict(f):
    d = {}  # fasta names are key, values are sequences as string
    namelist = []
    with open(f, 'r') as FILE:
        for line in FILE:
            if line[0] == '>':
                name = line.strip()[1:]
                namelist.append(name)
                d[name] = []
            elif line.strip() != '':  # else: # trying to prevent some crap happening on the last line
                d[name].append(line.strip())
    for name in namelist:
        d[name] = ''.join(d[name])  # join list of partial sequences. Useful if interleaved fasta

    return d, namelist


def main(fastafile, outfile, seqlength=25):
    dfasta, names = read_fasta_as_dict(fastafile)
    dtrim, trimnames = trim_fasta_to_sequence_length(dfasta, names, seqlength)
    output_file(outfile, dtrim, trimnames)
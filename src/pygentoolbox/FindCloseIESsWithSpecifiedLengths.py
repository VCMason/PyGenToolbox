class IES:
    def __init__(self, line):
        # line is line form annotation file
        self.id = line.strip().split('\t')[8].split(';')[0][len('ID='):]
        self.scaffold = line.strip().split('\t')[0]
        self.feature = line.strip().split('\t')[2]
        self.start = int(line.strip().split('\t')[3])
        self.end = int(line.strip().split('\t')[4])
        self.length = len(line.strip().split('\t')[8].split(';')[-1][len('sequence='):])
        self.sequence = line.strip().split('\t')[8].split(';')[-1][len('sequence='):]


def design_primers(seq, start, minlength, maxlength, maximum=25, minimum=18, optimal=20, mingcpercent=20.0):
    # seq is string sequence
    # start is start position of desired amplicon that can not have a primer within it
    # length is length of desired amplicon that can not have a primer within it
    # max is max length of primer to be designed
    # min is minimum length of primer to be designed
    # opt is optimal length of primer to be designed
    print('Starting primer design')

    import primer3plus
    import json
    from primer3plus.utils import reverse_complement

    design = primer3plus.Design()
    design.settings.template(seq)
    # <start>,<length> pairs where <start> is the index of the first base of a target,and <length> is its length.
    design.settings.target((start, minlength))
    # you can set primer3 parameters directly with .set()
    design.set('PRIMER_MAX_SIZE', maximum)
    design.set('PRIMER_MIN_SIZE', minimum)
    design.set('PRIMER_OPT_SIZE', optimal)
    design.set('PRIMER_MIN_GC', mingcpercent)
    design.set('PRIMER_PRODUCT_SIZE_RANGE', [minlength, maxlength])
    results, explain = design.run()
    if len(results) == 0:
        lloc = (0, 0)
        rloc = (0, 0)
        lseq = 'NOPRIMERSFOUND'
        rseq = 'NOPRIMERSFOUND'
        productsize = 0
    else:
        result = results[0]

        # ['location'] returns tuple (2,30) i.e. start position (zero based numbering) and length of primer
        # 'TCATGAA...' start = 2 is 'TC ^ ATGAA' so it starts at ATGAA
        lloc = result['LEFT']['location']
        lseq = result['LEFT']['SEQUENCE']
        rloc = result['RIGHT']['location']
        rseq = result['RIGHT']['SEQUENCE']
        productsize = result['PAIR']['PRODUCT_SIZE']

        ### NOTE indexing differences for reverse sequences using rloc
        # print('LEFT')
        # print(lseq)
        # print(design.SEQUENCE_TEMPLATE.value[lloc[0]:lloc[0] + lloc[1]])
        # print()
        # print('RIGHT')
        # print(rseq)
        # print(reverse_complement(design.SEQUENCE_TEMPLATE.value[rloc[0] + 1 - rloc[1]:rloc[0] + 1]))

        # to print all results
        # print(json.dumps(results, indent=1))
        # print(json.dumps(explain, indent=1))

    return lseq, rseq, lloc, rloc, productsize


def parse_fasta(f):
    print(f)
    seqs = {}
    FILE = open(f, 'r')
    l = FILE.readline()
    while l:
        if l[0] == '>':
            n = l[1:].strip().split()[0]
            seqs[n] = ''
        else:
            seqs[n] += l.strip().upper()
        l = FILE.readline()
    FILE.close()

    return seqs


def read_gff3(gff3file, scaffname, features=['all']):
    # print('Reading Gff3 file: %s' % gff3file)
    scaffoldlines = []
    with open(gff3file, 'r') as FILE:
        for line in FILE:
            if (line[0] == '#') or (line.strip() == ''):
                pass
            elif features == ['all']:  # keep all lines
                # key = line.strip().split('\t')[8].split(';')[0][len('ID='):]
                # d[key] = line.strip().split('\t')
                scaffoldlines.append(line)
            elif (line.strip().split('\t')[0] == scaffname) and (line.strip().split('\t')[2] in features):  # keep lines only if in features
                scaffoldlines.append(line.strip())
    # print(f'Number of features in scaffold {scaffname}: {len(scaffoldlines)}')
    # print(list(d.keys())[:10])
    return scaffoldlines


def main(gff3file='/media/sf_LinuxShare/Ciliates/Genomes/Annotations/internal_eliminated_sequence_PGM_ParTIES.pt_51_with_ies.gff3', seqfile='/media/sf_LinuxShare/Ciliates/Genomes/Seqs/ptetraurelia_mac_51_with_ies.fa', repattern='GATC', numscaff=696, window=3000, outeriesminlen=100, innerieslen=27):
    # annotfile is string as full path to annotation file
    # window is int as maximum allowed length between outer most IESs
    # outeriesminlen is int as minimum length of two outer IESs
    # innerieslen is int as length of inner IES

    import os

    p, f = os.path.split(gff3file)
    outfile = os.path.join(p, '.'.join(f.split('.')[:-1] + [f'{repattern}', f'win{window}', f'outmin{outeriesminlen}', f'in{innerieslen}', 'txt']))
    dseqs = parse_fasta(seqfile)

    outputiess = []
    for n in range(1, numscaff+1):
        scaffname = '_'.join(['scaffold51', str(n), 'with', 'IES'])
        scafflines = read_gff3(gff3file, scaffname, features=['internal_eliminated_sequence'])
        for count, line in enumerate(scafflines):
            if count <= len(scafflines) - 3:
                # iess will be list containing classes of IES
                iess = [IES(line)]
                iess.append(IES(scafflines[count + 1]))
                iess.append(IES(scafflines[count + 2]))
                # if
                # 1) distance between start and end of two outer IESs is <= window size
                # 2) length of both outer IESs are >= outeriesminlen
                # 3) length of one inner IES must be less than innerieslen
                # 4) restriction site is inside sequence of both outer IESs
                # 5) the innermost sequence of outer IESs (after splitting by restriction site) is >= 100 bps in length
                # 6) the restriction pattern does not split inner IES sequence
                if (iess[-1].end - iess[0].start <= window) and (iess[0].length >= outeriesminlen) \
                    and (iess[-1].length >= outeriesminlen) and (iess[1].length < innerieslen) \
                    and (len(iess[0].sequence.split(repattern)) > 1) and (len(iess[-1].sequence.split(repattern)) > 1) \
                    and (len(iess[0].sequence.split(repattern)[-1]) >= 100) \
                    and (len(iess[-1].sequence.split(repattern)[0]) >= 100) \
                    and (len(iess[1].sequence.split(repattern)) == 1):
                    # and (iess[1].length == innerieslen)
                    # if re site does not cut within the outermost re cut sites (which are inside outermost IESs)
                    # start and end are 1 based numbers
                    leftiesREposition = iess[0].end - len(iess[0].sequence.split(repattern)[-1])
                    rightiesREposition = iess[-1].start + len(iess[-1].sequence.split(repattern)[0]) - 1
                    seq = dseqs[scaffname][leftiesREposition: rightiesREposition]
                    print('checking if inner sequence does not have RE site')
                    if len(seq.split(repattern)) == 1:  # if it does not split the sequence
                        pcrproductminlength = iess[-1].start - iess[0].end
                        pcrproductmaxlength = len(seq)
                        print(f'Length of sequence after being cut by RE: {len(seq)}')
                        print(f'Desired: PCR product length: {pcrproductminlength}')
                        lseq, rseq, lloc, rloc, size = design_primers(seq, len(iess[0].sequence.split(repattern)[-1]),
                                                                pcrproductminlength, pcrproductmaxlength,
                                                                maximum=35, minimum=25, optimal=30, mingcpercent=15.0)

                        outputiess.append(scafflines[count:count + 3] + [
                                                                         f'Middle IES length = {iess[1].length}',
                                                                         iess[0].sequence.split(repattern)[-1], seq,
                                                                         iess[-1].sequence.split(repattern)[0],
                                                                         f'LEFT_PRIMER: {lseq}',
                                                                         f'LEFT_PRIMER_START_LENGTH: {lloc}',
                                                                         f'RIGHT_PRIMER: {rseq}',
                                                                         f'RIGHT_PRIMER_START_LENGTH: {rloc}',
                                                                         f'PCR Product size = {size}'
                                                                         ])
                        print(f'Left Primer: {lseq}')
                        print(f'Right Primer: {rseq}')
                        print('RESULT FOUND')
                        print(f'Inner IES length = {iess[1].length}')
                    else:
                        print('FAILED')
                        print(f'it had {len(seq.split(repattern)) - 1} RE sites present')

    print(f'Number of IES groups that match criterion = {len(outputiess)}')
    with open(outfile, 'w') as OUT:
        OUT.write('\n'.join(['\n'.join(iess) + '\n#####\n' for iess in outputiess]))
    print(f'output written to file: {outfile}')




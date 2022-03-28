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


def main(gff3file='/media/sf_LinuxShare/Ciliates/Genomes/Annotations/internal_eliminated_sequence_PGM_ParTIES.pt_51_with_ies.gff3', seqfile='/media/sf_LinuxShare/Ciliates/Genomes/Seqs/ptetraurelia_mac_51_with_ies.fa', repattern='GATC', numscaff=696, iesminlen=3000,  minrelen=500, minfraglen=300):
    # annotfile is string as full path to annotation file
    # window is int as maximum allowed length between outer most IESs
    # outeriesminlen is int as minimum length of two outer IESs
    # innerieslen is int as length of inner IES

    import os

    # minrelen = 500  # int(0.1 * iesminlen)  # minimum RE fragment length for left and right RE DNA fragments

    p, f = os.path.split(gff3file)
    outfile = os.path.join(p, '.'.join(f.split('.')[:-1] + [f'{repattern}', f'minieslen{iesminlen}', f'outmin{minrelen}', f'in{minrelen}', 'txt']))
    dseqs = parse_fasta(seqfile)

    outputiess = []
    for n in range(1, numscaff+1):
        scaffname = '_'.join(['scaffold51', str(n), 'with', 'IES'])
        scafflines = read_gff3(gff3file, scaffname, features=['internal_eliminated_sequence'])
        for count, line in enumerate(scafflines):
            # iess will be a class of an IES
            iess = IES(line)
            # if
            # 2) length of both outer IESs are >= iesminlen
            # 4) restriction site is inside the IES and present 2 times
            #   4) len(iess.sequence.split(repattern)) == 3 means three fragments from two cuts
            # 5) the leftmost re fragment and rightmost re fragment >= minrelen bps in length

            if (iess.length >= iesminlen) and (len(iess.sequence.split(repattern)) >= 3):
                # and (len(iess.sequence.split(repattern)[0]) >= minrelen) \
                # and (len(iess.sequence.split(repattern)[-1]) >= minrelen):
                # start and end are 1 based numbers
                numfrags = len(iess.sequence.split(repattern))
                # isolate the middle-most RE fragment inside the IES
                mid = round(numfrags / 2)  # if even (4) then value is 2 if odd (7) value is 4
                if (len(repattern.join(iess.sequence.split(repattern)[:mid-1])) >= minrelen) \
                        and (len(repattern.join(iess.sequence.split(repattern)[mid:])) >= minrelen):
                    leftREposition = iess.start + len(repattern.join(iess.sequence.split(repattern)[:mid-1])) - 1
                    rightREposition = iess.end - len(repattern.join(iess.sequence.split(repattern)[mid:]))
                    #leftREposition = iess.start + len(iess.sequence.split(repattern)[0]) - 1
                    seq = dseqs[scaffname][leftREposition: rightREposition]
                    #if len(seq.split(repattern)) == 3:  # the ies must have 3 fragments (two RE sites)
                    if len(seq) > minfraglen:
                        pcrproductminlength = len(seq)-100
                        pcrproductmaxlength = len(seq)
                        print('\nSUCCESS')
                        print(f'it had {numfrags - 1} RE sites present')
                        print(f'Length of sequence after being cut by RE: {len(seq)}')
                        print(f'Desired: PCR product length: {pcrproductminlength}')
                        lseq, rseq, lloc, rloc, size = design_primers(seq, 0, pcrproductminlength, pcrproductmaxlength,
                                                                maximum=35, minimum=25, optimal=30, mingcpercent=5.0)
                         # start... len(iess[0].sequence.split(repattern)[-1]
                        outputiess.append([scafflines[count]] + [
                                                               f'IES length = {iess.length}',
                                                               f'length of target RE framgent {len(seq)}',
                                                               f'Number RE sites: {numfrags - 1}',
                                                               f'min_RE_fragment_length {minrelen}',
                                                               f'left_RE_position: {leftREposition}',
                                                               f'right_RE_position: {rightREposition}',
                                                               seq,
                                                               ])
                                                                         # f'LEFT_PRIMER: {lseq}',
                                                                         # f'LEFT_PRIMER_START_LENGTH: {lloc}',
                                                                         # f'RIGHT_PRIMER: {rseq}',
                                                                         # f'RIGHT_PRIMER_START_LENGTH: {rloc}',
                                                                         # f'PCR Product size = {size}'
                        print(f'Left Primer: {lseq}')
                        print(f'Right Primer: {rseq}')
                        print('RESULT FOUND')
                        print(f'IES length = {iess.length}')
                        print(f'min_RE_fragment_length: {minrelen}')
                        print(f'left_RE_position: {leftREposition}')
                        print(f'right_RE_position: {rightREposition}')
                        print(f'length of target RE fragment {len(seq)}')
                        print(f'Number RE sites: {numfrags - 1}')
                    else:
                        print(f'\nFAILED, target fragment is to short < {minfraglen} bps')
                        print(f'it had {numfrags - 1} RE sites present')
                else:
                    print('FAILED, RE fragment ends are too close to ends of IES')

    print(f'Number of IESs that match criterion = {len(outputiess)}')
    with open(outfile, 'w') as OUT:
        OUT.write('\n'.join(['\n'.join(iess) + '\n#####\n' for iess in outputiess]))
    print(f'output written to file: {outfile}')


print('DONE')

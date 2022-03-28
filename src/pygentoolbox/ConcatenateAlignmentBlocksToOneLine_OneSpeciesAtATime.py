

def main(alignmentfile, speciesfile, referencespecies):
    ''' filename (full path) '''
    with open(alignmentfile + '.reformat.step.fa', 'w') as OUT:  # clear output file for appending later
        OUT.write('')

    print(alignmentfile)
    print(speciesfile)
    print(referencespecies)
    # read in all species from species file, one species name per line
    with open(speciesfile, 'r') as FILE:
        specieslist = [line.strip() for line in FILE]

    # d is a dictionary number of keys = number of species and values (will) = all concatenated seq blocks
    for speciesname in specieslist:
        speciesseq = []
        speciescoord = 1
        refcoord = 1
        with open(alignmentfile, 'r') as FILE:
            line = FILE.readline()
            while line:
                if (line[0] == '>') and (line.split('-')[0] == '>' + referencespecies):
                    # >Homo_sapiens-chr22(+)/10509923-10509970
                    # for humans only extract the start and end coordinates
                    start, end = [int(i) for i in line.strip().split('/')[1].split('-')]
                    # extend all sequences by some number of N's if the human start coordinate for the present block is
                    # some distance from the previous block (start - end) 1001 - 1001 = 0, 1005 - 1001 = 4
                    # the end coordinates are +1 from the actual end base for the block (seq is 1000 bases, end is 1001)
                    species = line.split('-')[0][1:]
                    line = FILE.readline()  # read next line
                elif line[0] == '>':
                    # for all other species
                    species = line.split('-')[0][1:]
                    line = FILE.readline()  # read next line
                # here you add the N extension to the beginning of the seq block (if needed as a spacer)
                refcoord = end
                if speciesname == species:
                    seq = 'N'*(start - speciescoord)  # also corrects for if species was missing from all prev blocks
                    speciescoord = end
                    if line[0] != '>':
                        # should only enter on line immediately after name lines
                        # gather the sequence (which could be interleaved) to non-interleaved
                        while line[0] != '>':
                            seq += line.strip()
                            line = FILE.readline()  # it will break when it is at a name line where line[0] == '>'
                        # add the sequence for each species
                        # creates keys = species name (if not already made)
                        # concatenates sequences
                        # if concatenating strings (slow) d[species] = d.setdefault(species, '') + seq
                        speciesseq.append(seq)
                    # chrom 22 is about 51 million bps long which is about 71 Gb in size
                # if end > 10509971:
                #    break

        # add N's to the end of each sequence if needed
        if speciescoord < refcoord:
            speciesseq.append('N' * (refcoord - speciescoord))

        with open(alignmentfile + '.reformat.step.fa', 'a') as OUT:
            # + '\n' just add a new line to the end of the file
            OUT.write('>%s\n%s\n' % (speciesname, ''.join(speciesseq)))
        print('Species finished: %s' % speciesname)

main('chr22_version1.aln', 'Species251.txt', 'Homo_sapiens')




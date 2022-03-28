

def main(alignmentfile, speciesfile, referencespecies):
    ''' filename (full path) '''
    print(alignmentfile)
    print(speciesfile)
    print(referencespecies)
    # read in all species from species file, one species name per line
    with open(speciesfile, 'r') as FILE:
        specieslist = [line.strip() for line in FILE]

    # d is a dictionary number of keys = number of species and values (will) = all concatenated seq blocks
    d = {species: [] for species in specieslist}  # initiate keys for all species in species list file
    d_coords = {species: 1 for species in specieslist}  # tracks the last long N extension when missing from prev block
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
            seq = 'N'*(start - d_coords[species])  # also corrects for if species was missing from all prev blocks
            d_coords[species] = end
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
                d.setdefault(species, []).append(seq)
            # chrom 22 is about 51 million bps long which is about 71 Gb in size
            if end > 10509971:
                break

    for species in specieslist:
        if d_coords[species] < d_coords[referencespecies]:
            d.setdefault(species, []).append('N' * (d_coords[referencespecies] - d_coords[species]))

    with open(alignmentfile + '.reformat.fa', 'w') as OUT:
        # + '\n' just add a new line to the end of the file
        OUT.write('\n'.join(['>%s\n%s' % (species, ''.join(d[species])) for species in specieslist]) + '\n')


main('chr22_version1.aln', 'Species251.txt', 'Homo_sapiens')




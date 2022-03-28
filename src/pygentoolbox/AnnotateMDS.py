
def write_out(outlines, outfile):
    with open(outfile, 'w') as OUT:
        OUT.write('\n'.join(outlines))
    print(f'Output file written to: {outfile}')


def create_metadata_mds(line, start, end):
    # create metadata for MDS
    metadata = line[8].split(';')
    for data in metadata:
        if data[:len('ID=')] == 'ID=':
            mdsid = 'ID=MDS' + '.'.join(data[len('ID=IES'):].split('.')[:-1] + [f'{start}', f'{end}'])
        if data[:len('Name=')] == 'Name=':
            mdsname = 'Name=MDS' + '.'.join(
                data[len('Name=IES'):].split('.')[:-1] + [f'{start}', f'{end}'])
    return mdsid, mdsname


def create_mds_annotations(dies, ids, dseqlens, scaffsingenomenotingff3):
    # last mds (betweeen ies and end of chromosome is excluded. I don't have chromosome length in this script
    # input is dictionary, key is ies ID
    # value is a list, representing the ies line from gff3 file
    prevscaff = ''
    mds = []  # list of output lines for mds only gff3
    mdsandies = []  # list of outputlines for mds and ies features for gff3
    for count, id in enumerate(ids):
        ies = dies[id]  # ies line of .gff3 file
        scaffold = ies[0]
        if count == 0: # if it is the first ies we need to make mds annotation before and after ies annotation
            # mds before IES
            mdsstart = 1
            mdsend = int(ies[3]) - 1
            mdsid, mdsname = create_metadata_mds(ies, mdsstart, mdsend)
            mdslist = ies[:2] + ['mds', str(mdsstart), str(mdsend)] + ies[5:8] + [';'.join([mdsid, mdsname])]
            mds.append('\t'.join(mdslist))
            # first ies
            mdsandies.append('\t'.join(mdslist))
            mdsandies.append('\t'.join(ies))
            # mds after IES
            mdsstart = int(ies[4]) + 1  # start of mds is one base passed the end of previous ies
            prevscaff = scaffold
        elif scaffold != prevscaff:  # if it is not the first ies but the first one of a new scaffold
            # make last MDS annotation for previous scaffold.
            # use previous mdsstart coordinate (last time scaff == prevscaff)  # end is length of previous scaffold
            mdsend = dseqlens[prevscaff]
            mdsid, mdsname = create_metadata_mds(previes, mdsstart, mdsend)
            mdslist = previes[:2] + ['mds', str(mdsstart), str(mdsend)] + previes[5:8] + [';'.join([mdsid, mdsname])]
            mds.append('\t'.join(mdslist))
            mdsandies.append('\t'.join(mdslist))
            # we need to make first mds annotation of scaffold before first ies
            # mds before IES
            mdsstart = 1
            mdsend = int(ies[3]) - 1
            mdsid, mdsname = create_metadata_mds(ies, mdsstart, mdsend)
            mdslist = ies[:2] + ['mds', str(mdsstart), str(mdsend)] + ies[5:8] + [';'.join([mdsid, mdsname])]
            mds.append('\t'.join(mdslist))
            # now make first ies annotation
            mdsandies.append('\t'.join(mdslist))
            mdsandies.append('\t'.join(ies))
            # mds after IES
            mdsstart = int(ies[4]) + 1  # start of mds is one base passed the end of previous ies
            prevscaff = scaffold
        elif scaffold == prevscaff:
            # mdsstart is already specified
            mdsend = int(ies[3]) - 1  # end of mds is start of this ies - 1
            mdsid, mdsname = create_metadata_mds(ies, mdsstart, mdsend)
            mdslist = ies[:2] + ['mds', str(mdsstart), str(mdsend)] + ies[5:8] + [';'.join([mdsid, mdsname])]
            mds.append('\t'.join(mdslist))
            mdsandies.append('\t'.join(mdslist))
            mdsandies.append('\t'.join(ies))
            mdsstart = int(ies[4]) + 1
            prevscaff = scaffold
            previes = ies

    # create last mds annotation for last scaffold that has ies annotation
    mdsend = dseqlens[prevscaff]
    mdsid, mdsname = create_metadata_mds(previes, mdsstart, mdsend)
    mdslist = previes[:2] + ['mds', str(mdsstart), str(mdsend)] + previes[5:8] + [';'.join([mdsid, mdsname])]
    mds.append('\t'.join(mdslist))
    mdsandies.append('\t'.join(mdslist))

    # create MDS annotation for each scaffold that doesn't have ies annotation
    for scaffold in scaffsingenomenotingff3:
        scaffnumber = scaffold[len('scaffold51_'):-len('_with_IES')]
        end = dseqlens[scaffold]
        mdslist = [scaffold, 'MICA', 'mds', str(1), str(end), '.', '.', '.', f'ID=MDSPGM.PTET51.1.{scaffnumber}.{1}.{end};Name=MDSPGM.PTET51.1.{scaffnumber}.{1}.{end}']
        mds.append('\t'.join(mdslist))
        mdsandies.append('\t'.join(mdslist))

    print(f'Number of mds features: {len(mds)}')
    return mds, mdsandies


def read_fasta_as_dict(f):
    d = {}  # fasta names are key, values are sequences as string
    namelist = []
    with open(f, 'r') as FILE:
        for line in FILE:
            if line[0] == '>':
                if ' ' in line:
                    name = line.strip().split()[0][1:]  # [1:-len('_with_IES')]
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


def length_of_fasta_sequences(genomefile):
    import os

    print('Counting lengths of all scaffolds')
    path, f = os.path.split(genomefile)
    dgenome, names = read_fasta_as_dict(genomefile)
    d = {k: len(v) for k, v in dgenome.items()}

    return d, names


def read_ies_gff3(gff3file):
    print('Reading Gff3 file: %s' % gff3file)
    d = {}
    iesids = []  # list of transcript IDs of mRNA features
    gff3scaffolds = []
    with open(gff3file, 'r') as FILE:
        for line in FILE:
            if (line[0] == '#') or (line.strip() == ''):
                pass
            elif line.strip().split('\t')[2] == 'internal_eliminated_sequence':  # keep all lines
                key = line.strip().split('\t')[8].split(';')[0][len('ID='):]
                iesids.append(key)
                d[key] = line.strip().split('\t')
                gff3scaffolds.append(line.strip().split('\t')[0])
    gff3scaffolds = set(gff3scaffolds)
    print('Number of ies features: %d' % len(list(d.keys())))
    print('Number of scaffolds in gff3: %d' % len(list(gff3scaffolds)))
    return d, iesids, gff3scaffolds


def main(gff3file, genomefile):
    # Annotate all regions in between IESs as MDS
    dies, ids, gff3scaffolds = read_ies_gff3(gff3file)
    dseqlens, seqnames = length_of_fasta_sequences(genomefile)
    setseqnames = set(seqnames)
    # these scaffolds do not have an ies in them
    scaffsingenomenotingff3 = {name for name in setseqnames if name not in gff3scaffolds}
    mds, mdsandies = create_mds_annotations(dies, ids, dseqlens, scaffsingenomenotingff3)

    outfile = '.'.join(gff3file.split('.')[:-1] + ['MDS', 'gff3'])
    write_out(mds, outfile)
    outfile = '.'.join(gff3file.split('.')[:-1] + ['MDS', 'IES', 'gff3'])
    write_out(mdsandies, outfile)
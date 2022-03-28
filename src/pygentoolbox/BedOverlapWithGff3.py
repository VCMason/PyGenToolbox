def bed_overlaping_gff3(dbed, bednames, dgff3, gff3names):
    overlapoutput = []
    for bedkey in bednames:
        bedscaff = dbed[bedkey][0]
        bedstart = int(dbed[bedkey][1])
        bedend = int(dbed[bedkey][2])
        bedid = dbed[bedkey][3]

        for gff3key in gff3names:
            gff3scaff = dgff3[gff3key][0]
            gff3start = int(dgff3[gff3key][3])
            gff3end = int(dgff3[gff3key][4])

            if (bedscaff == gff3scaff) and (gff3start <= bedend) and (bedstart <= gff3end):
                # if true then the features overlap
                overlapoutput.append('\t'.join([bedid, bedkey, gff3key, dgff3[gff3key][8].split(';')[0][3:]])) #  + dgff3[gff3key]

    return(overlapoutput)


def readgff3(f):
    dgff3 = {}
    gff3names = []
    with open(f, 'r') as FILE:
        for line in FILE:
            if line[0] != '#':
                s = line.strip().split('\t')
                dgff3['_'.join([s[0], s[3], s[4]])] = line.strip().split('\t')
                gff3names.append('_'.join([s[0], s[3], s[4]]))

    return(dgff3, gff3names)


def readbed(f):
    dbed = {}
    bednames = []
    with open(f, 'r') as FILE:
        for line in FILE:
            dbed['_'.join(line.strip().split('\t')[:3])] = line.strip().split('\t')
            bednames.append('_'.join(line.strip().split('\t')[:3]))

    return(dbed, bednames)


def main(fbed, fgff3, outfile):
    from natsort import natsorted
    print('Input files:\n%s\n%s' % (fbed, fgff3))
    dbed, bednames = readbed(fbed)
    dgff3, gff3names = readgff3(fgff3)
    overlapoutput = bed_overlaping_gff3(dbed, bednames, dgff3, gff3names)

    sortoutput = natsorted(overlapoutput)
    with open(outfile, 'w') as OUT:
        OUT.write('\n'.join(sortoutput))
    print('Output file:\n%s' % outfile)


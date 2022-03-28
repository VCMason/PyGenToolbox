def main(sam1, sam2):
    # FilterSam1ByReadNamesFromSam2
    outsamfile = '.'.join(sam1.split('.')[:-1] + ['RNAReadsConnectedDNAOver300winExtIES', 'sam'])

    readnames = []
    with open(sam2, 'r') as FILE:
        for line in FILE:
            if line[0] != '@':
                readnames.append(line.strip().split('\t')[0])

    print('Number of reads in sam2: %d' % len(readnames))

    d = {}
    with open(sam1, 'r') as FILE:
        for line in FILE:
            if line[0] != '@':
                d[line.strip().split('\t')[0]] = line.strip()

    count = 0
    outsamlist = []
    for name in readnames:
        try:
            d[name]
        except:
            pass
        else:
            outsamlist.append(d[name])
            count += 1
            if count % 10000 == 0:
                print(count)

    with open(outsamfile, 'w') as OUT:
        OUT.write('\n'.join(outsamlist))

    print('Number of reads in sam2: %d' % len(readnames))
    print('Number of reads from sam2 matching in sam1 (could include repeats) : %d' % count)
    print('output file: %s' % outsamfile)
    print('Finished')



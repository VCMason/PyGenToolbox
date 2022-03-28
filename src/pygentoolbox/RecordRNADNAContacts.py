
def main(fRNA, fDNA, MQ, sn):
    #fRNA full path to RNA sam file
    #fDNA full path to DNa sam file
    #MQ = int() = minimum mapping quality threshold
    #sn = str() = Species name
    from natsort import natsorted, ns
    import os
    import math

    sntrim = sn[0].lower() + sn.split(' ')[1][0]
    drna = {}
    with open(fRNA, 'r') as FILE:
        for line in FILE:
            if line[0] == '@':
                pass
            elif int(line.strip().split('\t')[4]) >= MQ:
                drna.setdefault(line.strip().split('\t')[0], []).append([line.strip().split('\t')[2], line.strip().split('\t')[3]])
    ddna = {}
    with open(fDNA, 'r') as FILE:
        for line in FILE:
            if line[0] == '@':
                pass
            elif int(line.strip().split('\t')[4]) >= MQ:
                ddna.setdefault(line.strip().split('\t')[0], []).append([line.strip().split('\t')[2], line.strip().split('\t')[3]])

    #rnaorderedkeys = natsorted(drna.keys(), alg=ns.IGNORECASE)
    #print(rnaorderedkeys[:10])
    output = []
    for k in drna.keys():
        try:
            ddna[k]
        except:
            pass
        else:
            # if the read maps in the RNA and DNA sam files
            for rnapos in drna[k]:  # the read could map to multiple locations 'equally' well (according to MQ cutoff)
                # print(rnapos)
                for dnapos in ddna[k]:  # the read could map to multiple locations 'equally' well (according to MQ)
                    # print(dnapos)
                    output.append([sntrim + rnapos[0].split('_')[1], rnapos[1], str(int(rnapos[1])+1), sntrim + dnapos[0].split('_')[1], dnapos[1], str(int(dnapos[1])+1)])  # , 'thickness=' + str(len(ddna[k])) # str(math.log(len(drna[k]) + len(ddna[k]), 10))
    sortoutput = natsorted(output, key=lambda y: (y[0], y[1]))
    path, f = os.path.split(fDNA)
    outpath = os.path.join(path, 'RNA.DNA.Contacts.%s.raw.txt' % sntrim)
    with open(outpath, 'w') as OUT:
        OUT.write('\n'.join([' '.join(l) for l in sortoutput]))
    print('Finished')

    d = {}
    for l in sortoutput:
        d[' '.join(l)] = d.setdefault(' '.join(l), 0) + 1  # count how many times each RNA - DNA contact appears
    # the keys of d are now unique and while the thickness value will represent how often the contacts appeared
    # need to split the keys again so we can sort them for output
    newsortoutput = natsorted([k.split(' ') for k in d.keys()], key=lambda y: (y[0], y[1]))
    newsortoutput2 = [l + ['thickness=' + str(d[' '.join(l)]) + 'p'] for l in newsortoutput]
    outpath = os.path.join(path, 'RNA.DNA.Contacts.%s.wthickness.txt' % sntrim)
    with open(outpath, 'w') as OUT:
        OUT.write('\n'.join([' '.join(l) for l in newsortoutput2]))
    print('Finished2')

    # take the log of each number of connections to reduce high values to <100
    # remove low counts less than 2
    d = {k: math.log(v, 2) for k, v in d.items() if v >= 2}
    # the keys of d are now unique and while the thickness value will represent how often the contacts appeared
    # need to split the keys again so we can sort them for output
    newsortoutput = natsorted([k.split(' ') for k in d.keys()], key=lambda y: (y[0], y[1]))
    newsortoutput2 = [l + ['thickness=' + str(d[' '.join(l)]) + 'p'] for l in newsortoutput]
    outpath = os.path.join(path, 'RNA.DNA.Contacts.%s.wthickness.log2.d2.txt' % sntrim)
    with open(outpath, 'w') as OUT:
        OUT.write('\n'.join([' '.join(l) for l in newsortoutput2]))
    print('Finished3')


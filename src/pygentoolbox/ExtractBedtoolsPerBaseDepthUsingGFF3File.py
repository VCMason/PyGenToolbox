
def ks_test():
    from scipy.stats import ks_2samp
    val = ks_2samp(x, y)  # val is D-statistic, and pvalue, val[0] = dstat


def main(beddepthfile, gff3file):
    iesdepths = {}
    with open(beddepthfile, 'r') as FILE:
        ddepth = {'_'.join(line.strip().split('\t')[:2]): int(line.strip().split('\t')[2]) for line in FILE}
    iesnames = []
    with open(gff3file, 'r') as FILE:
        for line in FILE:
            start = int(line.strip().split('\t')[3])
            end = int(line.strip().split('\t')[4])
            scaff = line.strip().split('\t')[0]
            iesnames.append('%s_%d_%d' % (scaff, start, end))
            for i in range(start, end):
                iesdepths.setdefault('%s_%d_%d' % (scaff, start, end), []) + [ddepth['%s_%d' % (scaff, i)]]
    for ies in iesnames:
        # sum all values for same element numbers and calculate "average" CPD
        dprob = {}
        for depth in iesdepths[ies]:
            dprob[ies] = depth / sum(iesdepths[ies])  # make each depth value for each position of each ies a probabilit
    # now i need to find the average probability per position, OORR







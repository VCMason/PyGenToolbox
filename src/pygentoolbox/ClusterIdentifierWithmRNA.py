
def plot_one_normalized_per_base_strandedness_y1_minus_y2(dnorm, halfbins, prefix, color1='blue', label1='F - R'):
    import matplotlib.pyplot as plt
    import os

    path, plottitle = os.path.split(prefix)

    # dnorm is a dictionary
    # dnorm key is a float representig midpoint of range ex: 0.025 for 0.00-0.05
    # dnorm value is list [sumf, sumr], sumf is sum of all forward reads for this bin, and same for sumr for reverse

    x = halfbins
    y1 = [dnorm[k][0] + dnorm[k][1] for k in halfbins]  # dnorm[k][1] # reverse reads are already negative
    print(x)
    print(y1)
    plt.plot(x, y1, color=color1, label=label1)  # 'F'- 'R' or 'S' - 'A'
    plt.grid()
    plt.legend(loc='upper left')
    plt.title(plottitle)

    outpath = prefix + '.perbase.strandedness.AllClusters.Normalized.plots.pdf'
    print(outpath)
    plt.savefig(outpath)
    plt.show()
    plt.close()

    return


def plot_one_normalized_per_base_strandedness(dnorm, halfbins, prefix, color1='blue', color2='orange', label1='F', label2='R'):
    import matplotlib.pyplot as plt
    import os

    path, plottitle = os.path.split(prefix)

    # dnorm is a dictionary
    # dnorm key is a float representig midpoint of range ex: 0.025 for 0.00-0.05
    # dnorm value is list [sumf, sumr], sumf is sum of all forward reads for this bin, and same for sumr for reverse

    x = halfbins
    y1 = [dnorm[k][0] for k in halfbins]
    y2 = [dnorm[k][1] for k in halfbins]  # reverse reads are already negative
    print(x)
    print(y1)
    print(y2)
    plt.plot(x, y1, color=color1, label=label1)  # 'F' or 'S'
    plt.plot(x, y2, color=color2, label=label2)  # # 'R' or 'A'
    plt.grid()
    plt.legend(loc='upper left')
    plt.title(plottitle)

    outpath = prefix + '.perbase.strandedness.AllClusters.Normalized.plots.pdf'
    print(outpath)
    plt.savefig(outpath)
    plt.show()
    plt.close()

    return


def make_sns_pair_plot(matrix, columnlabels, pcaindex, prefix):
    import pandas as pd
    import matplotlib.pyplot as plt
    # Seaborn visualization library
    import seaborn as sns
    from sklearn.preprocessing import StandardScaler

    scaler = StandardScaler()
    scaler.fit(matrix)
    scaledmatrix = scaler.transform(matrix)

    dfx = pd.DataFrame(data=scaledmatrix, columns=columnlabels, index=pcaindex)

    #pcaindex = [int(i) for i in pcaindex]
    # pcaindex is list of '0' and '1's to specify hue
    pairdf = dfx.assign(mRNA=pcaindex)
    print(pairdf.head())
    # Create the default pairplot
    fig = sns.pairplot(pairdf, hue='mRNA')  # , hue='mRNA' # vars=pairdf.columns[:-1],
    outpath = prefix + '.pairplot.pdf'
    print(outpath)
    fig.savefig(outpath)
    plt.show()
    plt.close()

    return


def make_pca(matrix, columnlabels, pcaindex, prefix):
    # matrix is list of lists, first column is rowIDs, first row is column names
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from sklearn import datasets
    from pca import pca
    from sklearn.preprocessing import StandardScaler

    ## iris = datasets.load_iris()
    ## label = iris.feature_names
    ## y = iris.target
    ## x = iris.data
    ## x = matrix
    #scaler = StandardScaler()
    #scaler.fit(matrix)
    #scaledmatrix = scaler.transform(matrix)
    ## matrix is list of lists
    ## columnlabels are the column names for the matrix
    ## pcaindex will be color of dots, i.e. if cluster is covered by >= number of mRNA
    #dfx = pd.DataFrame(data=scaledmatrix, columns=columnlabels, index=pcaindex)
    dfx = pd.DataFrame(data=matrix, columns=columnlabels, index=pcaindex)

    # model = pca()
    # The number of components are extracted that cover at least 95% of the explained variance.
    model = pca(n_components=0.95, normalize=True)  # , normalize=True
    results = model.fit_transform(dfx)

    print(model.results['topfeat'])
    outpath = prefix + '.pca.TopFeatures.tsv'
    results['topfeat'].to_csv(outpath, sep='\t')  # pandas dataframe

    # Cumulative explained variance
    print(model.results['explained_var'])
    outpath = prefix + '.pca.ExplainedVariance.tsv'
    results['explained_var'].tofile(outpath)  #numpy.array
    # [0.92461872 0.97768521 0.99478782]

    # Explained variance per PC
    print(model.results['variance_ratio'])
    outpath = prefix + '.pca.VarianceRatio.tsv'
    results['variance_ratio'].tofile(outpath)
    # [0.92461872, 0.05306648, 0.01710261]

    print(model.results['loadings'])
    outpath = prefix + '.pca.Loadings.tsv'
    results['loadings'].to_csv(outpath, sep='\t')

    print(model.results['PC'])
    outpath = prefix + '.pca.PC.tsv'
    results['PC'].to_csv(outpath, sep='\t')

    print(model.results['outliers'])
    outpath = prefix + '.pca.Outliers.tsv'
    results['outliers'].to_csv(outpath, sep='\t')

    # Make plot
    fig, ax = model.plot()
    outpath = prefix + '.pca.VarianceExplained.pdf'
    print(outpath)
    fig.savefig(outpath)  ##### use FIG.savfig, not PLT.savefig... #####
    plt.show()
    plt.close()

    # 2D plot
    fig, ax = model.scatter()
    outpath = prefix + '.pca.scatter.pdf'
    print(outpath)
    fig.savefig(outpath)
    plt.show()
    plt.close()
    # 3d Plot
    fig, ax = model.scatter3d()
    outpath = prefix + '.pca.scatter.3D.pdf'
    print(outpath)
    fig.savefig(outpath)
    plt.show()
    plt.close()

    # 2D plot
    fig, ax = model.biplot(n_feat=8, PC=[0, 1])  # n_feat controls number of vectors
    outpath = prefix + '.pca.biplot.pdf'
    print(outpath)
    fig.savefig(outpath)
    plt.show()
    plt.close()
    # 3d Plot
    fig, ax = model.biplot3d(n_feat=8, PC=[0, 1, 2])
    outpath = prefix + '.pca.biplot.3D.pdf'
    print(outpath)
    fig.savefig(outpath)
    plt.show()
    plt.close()

    return


def normalize_clusters_for_histogram_sense_antisense(dclusters, clusternames, dgff3, gff3names, nbins=20):
    import matplotlib.pyplot as plt
    import numpy as np

    # dnorm is a dictionary
    # dnorm key is a float (midpoint of range 0.00-0.05, ...)
    # dnorm value is list [sumA, sumS], sumA is sum of all sense reads for this bin, and same for sumS for antisense
    bins, halfbins = make_histogram_bins(nbins)
    dnorm = {midpoint: [0, 0] for midpoint in halfbins}  # initiate one key per bin for histogram with frequency value 0
    for name in clusternames:
        scaff = name.split('|')[0]
        start = int(name.split('|')[1])
        end = int(name.split('|')[2])
        geneids, geneorientations, genebinary = gff3_gene_overlap(dgff3, gff3names, scaff, start, end, coveragecutoff=50.0)
        # only gather information for clusters overlapping genes # gid is [] if no gene overlap
        for orientation in geneorientations:
            coordinates = dclusters[name][0]  # list of coordinates
            flst = dclusters[name][1]  # list of number of forward reads per coordinate
            rlst = dclusters[name][2]  # list of number of reverse reads per coordinate
            clusterlength = len(coordinates)
            for i in range(clusterlength):  # coordinates are actual scaffold positions, we need 0 ... len(cluster)
                # long clusters will have more entries / bin than small clusters so divide by clusterlength
                f = flst[i] / clusterlength
                r = rlst[i] / clusterlength
                # make cluster positions relative 0-1
                normpos = i / clusterlength
                for count, b in enumerate(bins):  # for each histogram bin, add the f and r values per normposition
                    left, right = b[0], b[1]
                    k = halfbins[count]  # k is name and midpoint of bin (0.050 if bin is 0.0-0.1)
                    if (normpos >= left) and (normpos < right):
                        # f is positive, r is negative. I want sense positive and antisense negative
                        if orientation == '+':  # sense reads in 0 position, antisense in 1 position
                            dnorm[k] = [dnorm[k][0] + f, dnorm[k][1] + r]
                        elif orientation == '-':
                            dnorm[k] = [dnorm[k][0] - r, dnorm[k][1] - f]  # -r bc -(-number) == + and -f bc -(+num) = -
                    # only need to add f and r values to one histogram bin

    return dnorm, halfbins


def make_histogram_bins(nbins, numdecimals=4):
    bins = []
    halfbins = []
    step = 1 / nbins
    halfstep = step / 2
    count = 0
    for i in range(nbins):
        bins.append([round(count, numdecimals), round(count + step, numdecimals)])
        halfbins.append(round(count + halfstep, numdecimals))
        count += step

    return bins, halfbins


def normalize_clusters_for_histogram(dclusters, clusternames, nbins=20):
    import matplotlib.pyplot as plt
    import numpy as np

    # dnorm is a dictionary
    # dnorm key is a float representing the midpoint of a range, ex: 0.025 for 0.00-0.05, ...)
    # dnorm value is list [sumf, sumr], sumf is sum of all forward reads for this bin, and same for sumr for reverse
    bins, halfbins = make_histogram_bins(nbins)
    # initiate one key per bin for histogram with frequency value 0 and 0 for f and r
    dnorm = {midpoint: [0, 0] for midpoint in halfbins}
    for name in clusternames:
        coordinates = dclusters[name][0]  # list of coordinates
        flst = dclusters[name][1]  # list of number of forward reads per coordinate
        rlst = dclusters[name][2]  # list of number of reverse reads per coordinate
        clusterlength = len(coordinates)
        for i in range(clusterlength):  # coordinates are actual scaffold positions, we need 0 ... len(cluster)
            # long clusters will have more entries / bin than small clusters so divide by clusterlength
            f = flst[i] / clusterlength
            r = rlst[i] / clusterlength
            # make cluster positions relative 0-1
            normpos = i / clusterlength
            for count, b in enumerate(bins):  # for each histogram bin, add the f and r values per normposition
                left, right = b[0], b[1]
                k = halfbins[count]  # round(start + halfstep, 4) # k is name and midpoint of bin (0.05 if bin 0.0-0.1)
                if (normpos >= left) and (normpos < right):
                    dnorm[k] = [dnorm[k][0] + f, dnorm[k][1] + r]
                # only need to add f and r values to one histogram bin

    return dnorm, halfbins


def plot_per_base_strandedness(dclusters, clusternames, prefix, dbed, bednames, numplots=4):
    import matplotlib.pyplot as plt

    for i in range(0, len(clusternames), numplots):
        numsubplots = int(numplots / 2)
        fig, ax = plt.subplots(numsubplots, numsubplots, figsize=(10, 10))  # sharex=False
        a, b = 0, 0
        for j in range(i, i + numplots):
            if j == len(clusternames):
                break  # break to avoid out of range clusternames[j]
            if j == i:
                pass
            elif j == i + numplots:
                b += 1
            elif j % 2 == 1:
                b += 1
            elif j % 2 == 0:
                a += 1
                b -= 1
            name = clusternames[j]
            scaff = name.split('|')[0]
            start = name.split('|')[1]
            end = name.split('|')[2]
            # srcids is list of all src cluster ids that overlap window
            srcids = src_bed_overlap(dbed, bednames, scaff, int(start), int(end))
            ids = ','.join(srcids)
            if ids == '':
                ids = 'NewCluster'
            x = dclusters[name][0]  # genomic positions, x axis coordinates
            y = dclusters[name][1]  # number of forward reads per x axis coordinate
            z = dclusters[name][2]  # number of reverse reads per x axis coordinate
            ax[a, b].plot(x, y, color='blue', label='F')
            ax[a, b].plot(x, z, color='orange', label='R')
            ax[a, b].set_title(ids + ': ' + name[len('scaffold51_'):])
            ax[a, b].set_xlim([int(start), int(end)])
            ax[a, b].legend(loc='upper left')
            for k, label in enumerate(ax[a, b].get_xticklabels()):
                if k % 2 == 1:
                    label.set_visible(True)
                elif k % 2 == 0:
                    label.set_visible(False)
        outpath = prefix + '.perbase.strandedness.plots' + str(i) + '.pdf'
        print(outpath)
        plt.savefig(outpath)
        plt.show()
        plt.close()

    return


def per_base_strandedness(bamfile, windowfile):
    import pysam

    with open(windowfile) as FILE:
        windows = [line.strip().split('\t') for line in FILE]

    samfile = pysam.AlignmentFile(bamfile, "rb")

    dclusters = {}
    clusternames = []
    for window in windows:
        scaff = window[0]
        start = int(window[1])
        end = int(window[2])
        name = '|'.join([scaff, window[1], window[2]])
        clusternames.append(name)
        coordinates, f, r = [], [], []  # f = sum of all forward reads per base, r = sum of all reverse reads per base
        # positions = [i for i in range(start, end)]
        # for pos in range(start, end):
        for pileupcolumn in samfile.pileup(scaff, start, end, stepper='nofilter'):  # pileup gives depth values outside region
            # if pileupcolumn.pos == pos:  # it's a bitch but need to check that the position is correct
            if (pileupcolumn.pos <= end) and (pileupcolumn.pos >= start):  # it's a bitch but need to check that the position is correct
                forward, reverse = 0, 0
                for pileupread in pileupcolumn.pileups:
                    if not pileupread.is_del and not pileupread.is_refskip:
                        if pileupread.alignment.is_reverse == True:
                            reverse -= 1
                        elif pileupread.alignment.is_reverse == False:
                            forward += 1
                coordinates.append(pileupcolumn.pos)
                f.append(forward)
                r.append(reverse)
        dclusters[name] = [coordinates, f, r]

    return dclusters, clusternames


def window_read_count(samfile, scaff, start, end):
    # samfile was already opened with (samfile = pysam.AlignmentFile(bamfile, "rb")
    count = 0
    for read in samfile.fetch(scaff, start, end):
        count += 1
    norm = (count / (end - start + 1)) * 100  # 100 is a scaling factor to make values more readable
    normbysize = round(norm, 2)  # normalizes mRNA count by cluster size
    if normbysize >= 10:
        coveredbymRNA = 1
    else:
        coveredbymRNA = 0

    return count, normbysize, coveredbymRNA


def gff3_gene_overlap(dgff3, gff3names, scaff, start, end, coveragecutoff=25.0):
    geneids = []
    orientations = []
    genebinary = 0  # 0 if overlaps no genes, 1 if overlaps at least one gene
    for k in gff3names:
        genescaff = dgff3[k][0]
        genestart = int(dgff3[k][3])
        geneend = int(dgff3[k][4])
        geneorientation = dgff3[k][6]
        id = dgff3[k][8].split(';')[0][3:]  # isolate gene id, from metadata column
        if (scaff == genescaff) and (start <= geneend) and (end >= genestart):
            # below calculates amount of overlap, if they don't overlap, returns 0
            # len(range(max(x[0], y[0]), min(x[-1], y[-1]) + 1))
            overlap = len(range(max(start, genestart), min(end, geneend) + 1))
            percentoverlap = overlap / (end - start + 1) * 100
            if percentoverlap >= coveragecutoff:  # if gene covers more than 25% of cluster record it
                geneids.append(id)
                orientations.append(geneorientation)
                genebinary = 1

    return geneids, orientations, genebinary


def src_bed_overlap(dbed, bednames, scaff, start, end):
    srcids = []
    for k in bednames:
        bedscaff = dbed[k][0]
        bedstart = int(dbed[k][1])
        bedend = int(dbed[k][2])
        id = dbed[k][3]
        if (scaff == bedscaff) and (start <= bedend) and (end >= bedstart):
            srcids.append(id)

    return srcids


def fasta_seq_gc_percent(d, scaff, start, end):
    GCcount = 0
    seq = d[scaff][start:end+1]
    for base in seq:
        if (base == 'G') or (base == 'C'):
            GCcount += 1
    refGC = round((GCcount / len(seq)) * 100, 2)

    return refGC, seq


def strandedness(bamfile, mRNAfile, windowfile, dbed, bednames, dgff3, gff3names, dgenome):
    import pysam

    with open(windowfile) as FILE:
        windows = [line.strip().split('\t') for line in FILE]

    newwindows = []
    pcatable = []
    pcarownames = []
    pcaindex = []
    countnewcluster = 0
    samfile = pysam.AlignmentFile(bamfile, "rb")  # pysam.AlignmentFile("ex1.bam", "rb")
    rnafile = pysam.AlignmentFile(mRNAfile, "rb")
    for window in windows:
        forward, reverse = 0, 0
        scaff = window[0]
        start = int(window[1])
        end = int(window[2])
        averagedepth = round(float(window[3]), 2)
        readsgc = round(float(window[4]), 2)

        length = end - start + 1
        refseqGCpercent, seq = fasta_seq_gc_percent(dgenome, scaff, start, end)

        for read in samfile.fetch(scaff, start, end):
            if read.is_reverse == True:
                reverse += 1
            elif read.is_reverse == False:
                forward += 1
        totalreads = forward + reverse
        strandfor = round((forward / totalreads) * 100, 2)
        strandrev = round((reverse / totalreads) * 100, 2)
        # srcids is list of all src cluster ids that overlap window
        srcids = src_bed_overlap(dbed, bednames, scaff, start, end)
        ids = ','.join(srcids)
        if ids == '':
            countnewcluster += 1
            ids = f'n{countnewcluster}'

        geneids, geneorientations, genebinary = gff3_gene_overlap(dgff3, gff3names, scaff, start, end)
        polarity = []
        if (strandfor >= 75.0) or (strandrev >= 75.0):
            strandbinary = 1  # 0 if cluster is not stranded (<75% reads of one orientation) 1 if stranded
        else:
            strandbinary = 0
        # 0 if cluster reads are sense to all coevered genes, or cluster reads not stranded
        # 1 if cluster reads are antisense to at least one covered gene
        # if multiple genes cover cluster, value is biased to 1, if one gene is antisense return 1
        antisensebinary = 0
        gids = ','.join(geneids)
        orientations = ','.join(geneorientations)
        if gids == '':
            gids = 'NA'  # no gene overlap
            orientations = 'NA'
            polarity.append('X')  # polarity will be S, A, N, -, or X
        else:
            for gorientation in geneorientations:
                if strandfor >= 75.0:
                    if gorientation == '+':
                        polarity.append('S')  # S == sense
                    elif gorientation == '-':
                        polarity.append('A')  # A == antisense
                        antisensebinary = 1
                    elif gorientation == '.':
                        polarity.append('N')  # N == no gene orientation information
                elif strandrev >= 75.0:
                    if gorientation == '+':
                        polarity.append('A')  # S == sense
                        antisensebinary = 1
                    elif gorientation == '-':
                        polarity.append('S')  # A == antisense
                    elif gorientation == '.':
                        polarity.append('N')  # N == no gene orientation information
                else:  # if if strandfor and strandrev are both < 75.0 (no majority of reads pointing in same direction)
                    polarity.append('-')  # - == no strong strand bias for cluster (can't say sense or antisense)
        pols = ','.join(polarity)

        # get number of mRNA reads overlapping window
        mRNAcount, countnormbysize, coveredbymRNA = window_read_count(rnafile, scaff, start, end)  # raw count

        newwindows.append('\t'.join([ids] + window + [str(refseqGCpercent), str(length), str(forward), str(reverse), str(strandfor), str(strandrev), str(strandbinary), gids, orientations, pols, str(antisensebinary), str(mRNAcount), str(countnormbysize), seq]))
        pcatable.append([str(averagedepth), str(readsgc), str(refseqGCpercent), str(length), str(strandfor), str(genebinary), str(antisensebinary), str(countnormbysize)])
        pcarownames.append(ids)
        pcaindex.append(str(coveredbymRNA))
    header = '\t'.join(['srcID', 'scaffold', 'start', 'end', 'average_depth', 'GCPercent_reads', 'GCPercent_reference', 'length', 'forward_reads', 'reverse_reads', 'forward_percent', 'stranded', 'reverse_percent', 'geneIDs', 'gene_orientations', 'clusterVSgene_polarity', 'antisense', 'mRNA_count', 'mRNA_count_normbyclusterlength', 'sequence'])
    pcacolumnnames = ['average_depth', 'ReadGC', 'ReferenceGC', 'length', 'forward_percent', 'covered_by_gene', 'antisense', 'mRNA_NormCount']
    samfile.close()

    outwindows = [header] + newwindows
    # pcaoutput = [pcaheader] + pcatable

    # x = [1, 2, 3]
    # y = [[4, 5, 6], [7, 8, 9]]
    # x + y = [1, 2, 3, [4, 5, 6], [7, 8, 9]]
    # [x] + y = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]

    print(f'Number of new clusters identified: {countnewcluster}')
    print(f'Number of total clusters identified: {len(newwindows)}')

    return outwindows, pcatable, pcacolumnnames, pcarownames, pcaindex


def refine_window_coordinates(sortwindows, bamfile):
    import pysam
    import numpy as np

    refinedwindows = []
    positions = set()
    # d = {}  # i don't keep all read depth values
    samfile = pysam.AlignmentFile(bamfile, "rb")  # pysam.AlignmentFile("ex1.bam", "rb")
    for win in sortwindows:
        scaff = win.split('\t')[0]
        winstart = int(win.split('\t')[1])
        winend = int(win.split('\t')[2])
        # array.arrays A, C, G, T. bases, per coordinate
        coverage = samfile.count_coverage(scaff, winstart, winend, read_callback='nofilter')
        numbases = np.sum(coverage)
        depth = list(np.sum(coverage, axis=0))  # sum four arrays vertically with np.sum axis=0
        # numATs = np.sum(coverage[0] + coverage[3])  # number of AT bases in window
        numGCs = np.sum(coverage[1] + coverage[2])  # total number of GC bases for all aligned reads in window
        countstart, countend = 0, 0
        while depth[0] == 0:
            countstart += 1
            del depth[0]
        while depth[-1] == 0:
            countend += 1
            del depth[-1]
        refinedstart = winstart + countstart
        refinedend = winend - countend
        refinedwinlength = refinedend - refinedstart + 1
        average_depth = mean_list(depth, denominator=refinedwinlength)
        GCpercent = round((numGCs / numbases) * 100, 2)  # GC percentage of all aligned reads
        refinedwindows.append(f'{scaff}\t{refinedstart}\t{refinedend}\t{average_depth}\t{GCpercent}')

    samfile.close()

    #refinedwindows = ['\t'.join(lst) for lst in refinedwindows]
    return refinedwindows


def write_table(table, outpath):
    # tale is list of lists
    print('Writing table to: %s' % outpath)
    output = ['\t'.join(l) for l in table]
    with open(outpath, 'w') as OUT:
        OUT.write('\n'.join(output))

def write_to_bed(windows, outpath):
    print('Writing windows to: %s' % outpath)
    with open(outpath, 'w') as OUT:
        OUT.write('\n'.join(windows))


def collapse_windows(sortwindows):
    # this function is janky
    # sortwindows is a naturally sorted list of strings
    # the strings inside the list are formatted as scaffold\tstartposition\tendposition
    collapsedsortwindows = []
    collapsewindepth = []
    collapseGC = []
    overlapgate = 1  # record start position of first window in overlap
    for count, win in enumerate(sortwindows):
        scaff = win.split('\t')[0]
        winstart = int(win.split('\t')[1])
        winend = int(win.split('\t')[2])
        windepth = float(win.split('\t')[3])
        GCpercent = float(win.split('\t')[4])
        if count == 0:
            pass
        elif (scaff == prevscaff) and (winstart <= prevwinend) and (prevwinstart <= winend):
            # some type of overlap # also equivalent to (i think) a[1] > b[0] and a[0] < b[1]
            # determine if there is overlap between previous window and current window, then combine them
            if overlapgate == 1:
                collapsewindow = '%s\t%s' % (scaff, prevwinstart)  # will add end coordinate when i find end of overlap
                overlapgate = 0
                collapsewindepth.append(prevwindepth)
                collapsewindepth.append(windepth)
                collapseGC.append(prevwinGC)
                collapseGC.append(GCpercent)
        elif overlapgate == 0:
            # if it makes it to this point, there is no overlap with previous window
            # but if overlapgate == 0 then there was overlap between at least two windows and we combine coordinates
            collapsedsortwindows.append('%s\t%s\t%.2f\t%.2f' % (collapsewindow, prevwinend, sum(collapsewindepth) / len(collapsewindepth), sum(collapseGC) / len(collapseGC)))
            collapsewindepth = []  # reset for next series of overlapping windows
            collapseGC = []
            overlapgate = 1  # reset overlapgate to allow for first start coordinates to be recorded
        elif overlapgate == 1:
            # if no overlap with previous, and overlapgate == 1 then no windows were overlapping
            # check if this window WILL overlap with the NEXT window and add it to collapsedsortwindows if NO
            if win != sortwindows[-1]:  # need to have this to avoid list indexing error
                if (scaff == sortwindows[count+1].split('\t')[0]) and (winstart <= int(sortwindows[count+1].split('\t')[2])) and (int(sortwindows[count+1].split('\t')[1]) <= winend):
                    pass
                else:
                    # no overlap with NEXT window
                    collapsedsortwindows.append(win)
                    collapsewindepth = []
                    collapseGC = []
            else:  # if it is the last window here output the window
                collapsedsortwindows.append(win)
                collapsewindepth = []
                collapseGC = []

        prevscaff = win.split('\t')[0]
        prevwinstart = int(win.split('\t')[1])
        prevwinend = int(win.split('\t')[2])
        prevwindepth = float(win.split('\t')[3])
        prevwinGC = float(win.split('\t')[4])

    return collapsedsortwindows


def mean_list(lst, denominator):
    mean = sum(lst) / denominator
    return mean


def mean_with_denominator(value, denominator):
    # mean = sum(lst) / denominator
    mean = value / denominator
    return mean


def calculate_depth_with_pysam_sliding_windows(bamfile, dscaff_max, step=25, winsize=100, cutoff=100):
    import pysam
    import numpy as np

    windows = []
    positions = set()
    # d = {}  # i don't keep all read depth values
    samfile = pysam.AlignmentFile(bamfile, "rb")  #pysam.AlignmentFile("ex1.bam", "rb")
    for scaff, maxcoord in dscaff_max.items():
        print(scaff, end=', ')
        for i in range(0, maxcoord, step):
            ## tempdepth = []  # temppos = []
            ## .pileup skips bases when depth = 0  # need to figure out mpileup -a option in pysam
            ## tempdepth = [pileupcolumn.nsegments for pileupcolumn in samfile.pileup(scaff, i, i + winsize)]
            #tempdepth = []
            #for pos in range(i, i + winsize):
            #    depth = 0
            #    for pileupcolumn in samfile.pileup(scaff, pos, pos+1):  # pileup gives depth values outside region
            #        if pileupcolumn.pos == pos:  # it's a bitch but need to check that the position is correct
            #            #for pileupread in pileupcolumn.pileups:
            #            #    if not pileupread.is_del and not pileupread.is_refskip:
            #            #        depth += 1
            #            #tempdepth.append(depth)
            #            ##tempdepth.append(pileupcolumn.nsegments)  # pileupcolumn.nsegments # n depricated
            coverage = samfile.count_coverage(scaff, i, i + winsize, read_callback='nofilter')  # array.arrays A, C, G, T. bases, per coordinate
            numbases = np.sum(coverage)
            # numATs = np.sum(coverage[0] + coverage[3])  # number of AT bases in window
            numGCs = np.sum(coverage[1] + coverage[2])  # number of GC bases in window

            if numbases >= 10:
                average_depth = mean_with_denominator(numbases, denominator=winsize)
                GCpercent = round((numGCs / numbases) * 100, 2)
            else:
                average_depth = 0
                GCpercent = 0
            if average_depth >= cutoff:
                # print(f'{average_depth}, {scaff}, {i+1}, {i+winsize}, {tempdepth}')
                # print(str(average_depth), end=', ')
                # scaff, start, end, averagedepth
                windows.append([scaff, str(i + 1), str(i + winsize), str(average_depth), str(GCpercent)])
                for j in range(i + 1, i + winsize):
                    if j <= maxcoord:
                        positions.add('_'.join([scaff, str(j)]))
    samfile.close()

    return windows, positions


def read_gff3(f):
    # assumes gff3 is trimmed to 'gene' lines only (no exon/CDS/whatever else), one line per gene
    dgff3 = {}
    gff3names = []
    with open(f, 'r') as FILE:
        for line in FILE:
            if line[0] != '#':
                s = line.strip().split('\t')
                dgff3['|'.join([s[0], s[3], s[4]])] = line.strip().split('\t')
                gff3names.append('|'.join([s[0], s[3], s[4]]))

    return dgff3, gff3names


def read_bed(f):
    dbed = {}
    bednames = []
    with open(f, 'r') as FILE:
        header = FILE.readline()
        for line in FILE:
            dbed['|'.join(line.strip().split('\t')[:3])] = line.strip().split('\t')
            bednames.append('|'.join(line.strip().split('\t')[:3]))

    return dbed, bednames


def read_fasta_as_dict(f):
    d = {}  # fasta names are key, values are sequences as string
    namelist = []
    with open(f, 'r') as FILE:
        for line in FILE:
            if line[0] == '>':
                if ' ' in line:
                    name = line.strip().split()[0][1:-len('_with_IES')]
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

    return d, dgenome, names


def max_coords_per_scaffold(d):
    dmax_coords = {}
    dscaff_max = {}
    for n in list(d.keys()):
        # collect all coordinates in list, for each scaffold
        dmax_coords.setdefault('_'.join(n.split('_')[:-1]), []).append(int(n.split('_')[-1]))
    for scaff in list(dmax_coords.keys()):
        # get maximum coordinate per scaffold
        dscaff_max[scaff] = max(dmax_coords[scaff])
    return dscaff_max


def main(bamfiles, mRNAbamfiles, genomefile, srcbedfile, gff3file, step=25, winsize=100, cutoff=100, numplots=4):
    # input is sam file # output does not include positions of 0 coverage
    # key is scafffold_position # value is depth
    import os
    from natsort import natsorted, ns
    import pandas as pd

    for bam, mRNA in zip(bamfiles, mRNAbamfiles):
        path, file = os.path.split(bam)

        prefix = os.path.join(path, '.'.join(file.split('.')[:-1]))
        print('step: %d, window size: %d, cutoff: %d' % (step, winsize, cutoff))

        dscafflengths, dgenome, names = length_of_fasta_sequences(genomefile)
        windows, positions = calculate_depth_with_pysam_sliding_windows(bam, dscafflengths, step, winsize, cutoff)

        sortwindows = natsorted(['\t'.join(w) for w in windows], key=lambda y: y.lower())
        write_to_bed(sortwindows, '%s.windows.s%d.w%d.d%d.pileup.bed' % (prefix, step, winsize, cutoff))
        collapsedsortwindows = collapse_windows(sortwindows)
        #collapsedwindowfilename = '%s.collapsedwindows.s%d.w%d.d%d.pileup.bed' % (prefix, step, winsize, cutoff)
        print(f'Starting to refine window coordinates.')
        refinedsortwindows = refine_window_coordinates(collapsedsortwindows, bam)
        refinedwindowfilename = '%s.refinedwindows.s%d.w%d.d%d.pileup.bed' % (prefix, step, winsize, cutoff)
        write_to_bed(refinedsortwindows, refinedwindowfilename)

        dbed, bednames = read_bed(srcbedfile)
        dgff3, gff3names = read_gff3(gff3file)

        strandwindows, pcatable, pcacolumnnames, pcarownames, pcaindex = strandedness(bam, mRNA, refinedwindowfilename, dbed, bednames, dgff3, gff3names, dgenome)
        strandedfilename = '%s.refinedwindows.stranded.s%d.w%d.d%d.pileup.tsv' % (prefix, step, winsize, cutoff)
        write_to_bed(strandwindows, strandedfilename)
        pcatablefilename = '%s.refinedwindows.stranded.s%d.w%d.d%d.pileup.pcatable.tsv' % (prefix, step, winsize, cutoff)
        write_table(pcatable, pcatablefilename)

        dclusters, clusternames = per_base_strandedness(bam, refinedwindowfilename)
        plot_per_base_strandedness(dclusters, clusternames, prefix, dbed, bednames, numplots)

        # dfx = pd.DataFrame(data=pcatable, columns=pcacolumnnames, index=pcaindex)
        # print(dfx.head())
        make_pca(pcatable, pcacolumnnames, pcaindex, prefix)
        # pairdf = pd.DataFrame(data=pcatable, columns=pcacolumnnames, index=pcarownames)
        make_sns_pair_plot(pcatable, pcacolumnnames, pcaindex, prefix)

        dnorm, halfbins = normalize_clusters_for_histogram(dclusters, clusternames, nbins=50)
        plot_one_normalized_per_base_strandedness(dnorm, halfbins, prefix)
        plot_one_normalized_per_base_strandedness_y1_minus_y2(dnorm, halfbins, prefix, color1='blue', label1='F - R')

        prefix = prefix + '.sense.antisense'
        dnorm, halfbins = normalize_clusters_for_histogram_sense_antisense(dclusters, clusternames, dgff3, gff3names, nbins=50)
        plot_one_normalized_per_base_strandedness(dnorm, halfbins, prefix, color1='cyan', color2='red', label1='Sense', label2='Antisense')
        prefix = prefix + '.sense.antisense.difference'
        plot_one_normalized_per_base_strandedness_y1_minus_y2(dnorm, halfbins, prefix, color1='teal', label1='Sense - Antisense')


    print('###FIN###')

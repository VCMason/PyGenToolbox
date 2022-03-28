############################################################################################
##########        # finds closest sRNA to each mRNA on each scaffold          ##########
############################################################################################


def plot_hexbins(x, y, outpath, nbins=20, xlab='mRNA position', ylab='sRNA position'):
    # not too slow
    print('Plotting scatter, hexbins, 2d histogram')
    # Libraries
    import numpy as np
    import matplotlib.pyplot as plt

    x, y = np.array(x), np.array(y)

    # Create a figure with 6 plot areas
    fig, axes = plt.subplots(ncols=3, nrows=1, figsize=(21, 5))

    # Everything sarts with a Scatterplot
    axes[0].set_title('Scatterplot')
    axes[0].plot(x, y, 'ko')
    # As you can see there is a lot of overplottin here!

    # Thus we can cut the plotting window in several hexbins
    # nbins = 20
    axes[1].set_title('Hexbin')
    axes[1].hexbin(x, y, gridsize=nbins, cmap=plt.cm.BuGn_r)

    # 2D Histogram
    axes[2].set_title('2D Histogram')
    axes[2].hist2d(x, y, bins=nbins, cmap=plt.cm.BuGn_r)

    for ax in axes.flat:
        ax.set(xlabel=xlab, ylabel=ylab)

    # Hide x labels and tick labels for top plots and y ticks for right plots.
    for ax in axes.flat:
        ax.label_outer()
        
    plt.show()
    plt.savefig(outpath)
    plt.close()


def plot_density(x, y, outpath, nbins=20):
    print('Plotting density')
    import matplotlib.pyplot as plt
    import numpy as np
    from scipy.stats import kde

    # Evaluate a gaussian kde on a regular grid of nbins x nbins over data extents
    x = np.array(x)
    y = np.array(y)
    # nbins = 300
    k = kde.gaussian_kde([x, y])
    xi, yi = np.mgrid[x.min():x.max():nbins * 1j, y.min():y.max():nbins * 1j]
    zi = k(np.vstack([xi.flatten(), yi.flatten()]))

    # Make the plot
    plt.pcolormesh(xi, yi, zi.reshape(xi.shape))

    # Change color palette
    plt.pcolormesh(xi, yi, zi.reshape(xi.shape), cmap=plt.cm.Greens_r)
    plt.show()
    plt.savefig(outpath)
    plt.close()


def plot_scatter(x, y, outpath):
    print('Plotting')

    import matplotlib.pyplot as plt
    plt.scatter(x, y)
    plt.show()
    plt.savefig(outpath)
    plt.close()


def mean(lst):
    return sum(lst) / len(lst)


def compare_positions_gff3_sam(dgff3, dsRNA, dmax={}, analysis='allpairs', type='genomecoordinates'):  # dsRNA stands for dictionary small RNA
    # make list with format [[x,y], [x,y], [x,y]...]
    # xy = [[int(TEline.split('\t')[3]), int(sRNAline.split('\t')[3])] for TEline in dgff3[scaffold] for sRNAline in dsRNA[scaffold]]
    from scipy.stats.stats import pearsonr
    print('Collecting x and y data')
    print('Performing analysis: %s and of type: %s' % (analysis, type))

    xy = []
    for scaffold in list(dgff3.keys()):
        for TEline in dgff3[scaffold]:
            if type == 'proportions':
                TEcoord = float(TEline[3]) / dmax[scaffold]
            else:
                TEcoord = int(TEline[3])  # position is 5' of forward reads, 3' of reverse reads
            try:
                dsRNA[scaffold]
            except:
                pass  # no sRNA mapped to the same scaffold as the TE
            else:  # some sRNA mapped to the same scaffold as the TE
                if analysis == 'closestpairs':
                    for i in range(len(dsRNA[scaffold])):  # find closest mRNA to each sRNA
                        sRNAline = dsRNA[scaffold][i]
                        if i == 0:  # for first line only
                            sRNAcoord = mean([int(sRNAline[3]), int(sRNAline[4])])
                        else:  # for all other lines
                            currentdist = abs(TEcoord - mean([int(sRNAline[3]), int(sRNAline[4])]))
                            previousdist = abs(TEcoord - mean([int(dsRNA[scaffold][i - 1][3]), int(dsRNA[scaffold][i - 1][4])]))
                            if currentdist < previousdist:  # could exclude mRNAs of equal distance after first entry
                                sRNAcoord = mean([int(sRNAline[3]), int(sRNAline[4])])
                            else:
                                pass
                    if type == 'proportions':
                        sRNAcoord = sRNAcoord / dmax[scaffold]
                    else:
                        #  gff3coord = gff3coord   # so don't need to rename it
                        pass
                    xy.append([TEcoord, sRNAcoord])  # record only one mrna coordinate for each sRNA
                else:
                    for gff3line in dgff3[scaffold]:
                        if type == 'proportions':
                            sRNAcoord = mean([int(sRNAline[3]), int(sRNAline[4])]) / dmax[scaffold]
                        else:
                            sRNAcoord = mean([int(sRNAline[3]), int(sRNAline[4])])
                        xy.append([TEcoord, sRNAcoord])


    print('Number of plot points = %d' % len(xy))
    x = [plotpoint[0] for plotpoint in xy]
    y = [plotpoint[1] for plotpoint in xy]
    print('Pearson Correlation:')
    print(pearsonr(x, y))
    return x, y


def read_sam(sam, positions=[]):
    ''' key is scaffold, value is line.strip() '''
    print('Reading sam file: %s' % sam)
    d = {}
    dhead = {}
    with open(sam, 'r') as FILE:
        for line in FILE:
            if line[0] != '@':
                if int(line.strip().split('\t')[3]) not in positions:
                    key = line.strip().split('\t')[2]
                    d.setdefault(key, []).append(line.strip().split('\t'))
            elif line[:3] == '@SQ':
                dhead[line.strip().split('\t')[1][3:]] = int(line.strip().split('\t')[2].split(':')[1])
    print('Number of scaffolds: %d' % len(list(dhead.keys())))
    return d, dhead


def read_gff3(GFF3file, features=['all']):
    print('Reading Gff3 file: %s' % GFF3file)
    d = {}
    dhead = {}
    with open(GFF3file, 'r') as FILE:
        for line in FILE:
            if line[:5] == '##gff':
                pass
            elif line[:17] == '##sequence-region':
                dhead[line.strip().split()[1]] = int(line.strip().split()[3])
            elif features == ['all']:  # keep all lines
                key = line.strip().split('\t')[0]
                d.setdefault(key, []).append(line.strip().split('\t'))
            elif line.strip().split('\t')[2] in features:  # keep lines only if in features
                key = line.strip().split('\t')[0]
                d.setdefault(key, []).append(line.strip().split('\t'))
    print('Number of scaffolds: %d' % len(list(dhead.keys())))
    # print(list(d.values())[:2])
    return d, dhead


def main(GFF3file, samfiles):
    import os

    # dTE, dTEhead = read_sam(TEfile, positions=[0])  # sam file of aligned TEs

    dmRNA, dmRNAhead = read_gff3(GFF3file, features=['mRNA'])

    for sam in samfiles:
        dsRNA, dsRNAhead = read_sam(sam, positions=[0])  # sam files of aligned sRNAs # dsRNA stands for dictionary small RNA

        dmax = {**dsRNAhead, **dmRNAhead}  # update dsRNAhead with values from dTEhead

        # x, y = compare_positions_gff3_sam(dmRNA, dsRNA, analysis='allpairs')
        # path, file = os.path.split(sam)
        # outpath = os.path.join(path, '.'.join(file.split('.')[:-1] + ['hexbins', 'png']))
        # plot_hexbins(x, y, outpath, 50)
        # 
        # x, y = compare_positions_gff3_sam(dmRNA, dsRNA, dmax, analysis='allpairs', type='proportions')
        # outpath = os.path.join(path, '.'.join(file.split('.')[:-1] + ['prop', 'hexbins', 'png']))
        # plot_hexbins(x, y, outpath, 50)

        print('##### closest pairs #####')

        x, y = compare_positions_gff3_sam(dmRNA, dsRNA, analysis='closestpairs')
        path, file = os.path.split(sam)
        outpath = os.path.join(path, '.'.join(file.split('.')[:-1] + ['hexbins', 'png']))
        plot_hexbins(x, y, outpath, 20)

        x, y = compare_positions_gff3_sam(dmRNA, dsRNA, dmax, analysis='closestpairs', type='proportions')
        outpath = os.path.join(path, '.'.join(file.split('.')[:-1] + ['prop', 'hexbins', 'png']))
        plot_hexbins(x, y, outpath, 20)

        print('###########################################')
        print('###########################################')
        print('###########################################')

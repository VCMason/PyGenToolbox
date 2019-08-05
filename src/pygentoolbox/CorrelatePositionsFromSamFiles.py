
def plot_hexbins(x, y, outpath, nbins=20):
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


def compare_positions(dTE, dsRNA, dmax={}, analysis='allpairs', type='genomecoordinates'):  # dsRNA stands for dictionary small RNA
    # make list with format [[x,y], [x,y], [x,y]...]
    # xy = [[int(TEline.split('\t')[3]), int(sRNAline.split('\t')[3])] for TEline in dTE[scaffold] for sRNAline in dsRNA[scaffold]]
    from scipy.stats.stats import pearsonr
    print('Collecting x and y data')
    print('Performing analysis: %s and of type: %s' % (analysis, type))

    xy = []
    for scaffold in list(dTE.keys()):
        for TEline in dTE[scaffold]:
            if type == 'proportions':
                TEcoord = float(TEline.split('\t')[3]) / dmax[scaffold]
            else:
                TEcoord = int(TEline.split('\t')[3])  # position is 5' of forward reads, 3' of reverse reads
            try:
                dsRNA[scaffold]
            except:
                pass  # no sRNA mapped to the same scaffold as the TE
            else:  # some sRNA mapped to the same scaffold as the TE
                if analysis == 'closestpairs':
                    distances = [abs(TEcoord - int(sRNAline.split('\t')[3])) for sRNAline in dsRNA[scaffold]]
                    if type == 'proportions':
                        sRNAcoord = min(distances) / dmax[scaffold]
                    else:
                        sRNAcoord = min(distances)
                    xy.append([TEcoord, sRNAcoord])
                else:
                    for sRNAline in dsRNA[scaffold]:
                        if type == 'proportions':
                            sRNAcoord = float(sRNAline.split('\t')[3]) / dmax[scaffold]
                        else:
                            sRNAcoord = int(sRNAline.split('\t')[3])
                        xy.append([TEcoord, sRNAcoord])

    print('Number of plot points = %d' % len(xy))
    x = [plotpoint[0] for plotpoint in xy]
    y = [plotpoint[1] for plotpoint in xy]
    print('Pearson Correlation:')
    print(pearsonr(x, y))
    return x, y


def max_of_maximums(dTEmax, dsRNAmax):
    print('Finding max of maximums')
    dmax = {}
    for i in list(dTEmax.keys()):
        for j in list(dsRNAmax.keys()):
            try:
                dTEmax[j]
            except:
                dmax[j] = dsRNAmax[j]
            else:
                dmax[j] = max(dTEmax[j], dsRNAmax[j])

            try:
                dsRNAmax[i]
            except:
                dmax[i] = dTEmax[i]  # if no i in dsRNAmax
            else:  # if scaffold is present in both dictionaries, find the max between them, i and j should be same
                dmax[i] = max(dTEmax[i], dsRNAmax[i])
    for k, v in list(dmax.items()):
        if v == 0:
            print(k)
    return dmax


def max_coords_per_scaffold(d):
    print('Calcualting maximum coordinates per scaffold')
    dmax_coords = {}
    dscaff_max = {}
    for scaff in list(d.keys()):
        # collect all coordinates in list, for each scaffold
        for line in d[scaff]:
            dmax_coords.setdefault(scaff, []).append(int(line.split('\t')[3]))
    for scaff in list(dmax_coords.keys()):
        # get maximum coordinate per scaffold
        dscaff_max[scaff] = max(dmax_coords[scaff])
    return dscaff_max


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
                    d.setdefault(key, []).append(line.strip())
            elif line[:3] == '@SQ':
                    dhead[line.strip().split('\t')[1][3:]] = int(line.strip().split('\t')[2].split(':')[1])
    print('Number of scaffolds: %d' % len(list(dhead.keys())))
    return d, dhead


def main(TEfile, samfiles):
    import os

    dTE, dTEhead  = read_sam(TEfile, positions=[0])  # sam file of aligned TEs

    for sam in samfiles:
        dsRNA, dsRNAhead = read_sam(sam, positions=[0])  # sam files of aligned sRNAs # dsRNA stands for dictionary small RNA

        dmax = {**dsRNAhead, **dTEhead}  # update dsRNAhead with values from dTEhead

        x, y = compare_positions(dTE, dsRNA, analysis='allpairs')
        path, file = os.path.split(sam)
        outpath = os.path.join(path, '.'.join(file.split('.')[:-1] + ['hexbins', 'png']))
        plot_hexbins(x, y, outpath, 50)
        
        x, y = compare_positions(dTE, dsRNA, dmax, analysis='allpairs', type='proportions')
        outpath = os.path.join(path, '.'.join(file.split('.')[:-1] + ['prop', 'hexbins', 'png']))
        plot_hexbins(x, y, outpath, 50)

        print('##### closest pairs #####')

        x, y = compare_positions(dTE, dsRNA, analysis='closestpairs')
        path, file = os.path.split(sam)
        outpath = os.path.join(path, '.'.join(file.split('.')[:-1] + ['hexbins', 'png']))
        plot_hexbins(x, y, outpath, 50)

        x, y = compare_positions(dTE, dsRNA, dmax, analysis='closestpairs', type='proportions')
        outpath = os.path.join(path, '.'.join(file.split('.')[:-1] + ['prop', 'hexbins', 'png']))
        plot_hexbins(x, y, outpath, 50)


        print('###########################################')
        print('###########################################')
        print('###########################################')
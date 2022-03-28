# use ConverBedToRangeCoordinates.py to make bedrangecoord file
# use DoMyRNATasteLikePi to get forandrevcountfile


def seaborn_lineplot_with_barh(x, y, z, outpath, names, ys, starts, lengths):
    # x is x axis position, # y is number of forward read values per position, # z is reverse read values (negative)
    import numpy as np
    import seaborn as sns
    import matplotlib.pyplot as plt

    print('Plotting heatmap with horizontal bar plot')

    fig, axes = plt.subplots(2, 1, figsize=(8, 8), sharex=True, gridspec_kw={'height_ratios': [1, 32]})
    sns.lineplot(x, y, ax=axes[1])
    sns.lineplot(x, z, ax=axes[1])
    # axes[1].set(ylabel=None)
    axes[0].get_yaxis().set_visible(False)
    axes[0].barh(y=ys, width=lengths, left=starts)
    ends = [starts[i] + lengths[i] for i in range(len(starts))]
    for s, e in zip(starts, ends):
        axes[1].axvline(s, 0, len(x), linestyle='--', linewidth=0.25)
        axes[1].axvline(e, 0, len(x), linestyle='--', linewidth=0.25)
    if outpath != '':
        fig.savefig(outpath)
        print('Finished plotting line and features: %s' % outpath)
    else:
        print('Finished plotting line and features')

    plt.show()
    plt.close()

    return


def plot_scatter(x, y, outpath, title_analysis, exitplotcount, countplot, prefix='NaN', xlab='Distance From Read End', ylab='Depth', figlength=10, figheight=5):
    ''' x is list of x values, y is list of y values, path is full path to output file '''
    import matplotlib.pyplot as plt
    import os
    # font = {'family': 'normal',
    #         'weight': 'normal',
    #         'size': 20}
    if countplot == 1:
        plt.rcParams.update({'font.size': 16})
        plt.figure(figsize=(figlength, figheight))

    plt.plot(x, y, label=prefix)
    plt.grid(True)
    plt.title(title_analysis)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.legend(loc='upper right')
    plt.savefig(outpath)
    if countplot == exitplotcount:  # 2 = len(infiles), infiles passed to phasing function
        path, outname = os.path.split(outpath)
        print(outname)
        plt.show()
        plt.close()


def main(bedrangecoordfile, forandrevcountfile):
    # from natsort import natsorted, ns
    import os

    d = {}
    d_for = {}
    d_rev = {}
    with open(forandrevcountfile, 'r') as FILE:
        for line in FILE:
            d[line.strip().split('\t')[0]] = line.strip().split('\t')
            d_for[line.strip().split('\t')[0]] = int(line.strip().split('\t')[1])
            d_rev[line.strip().split('\t')[0]] = int(line.strip().split('\t')[2])
    output = []
    sortnames = []
    # each bed range file represetns one feature
    feature_for_depth = []  # list of forward read depth at each position
    feature_rev_depth = []  # list of reverse read depth at each position
    with open(bedrangecoordfile, 'r') as FILE:
        #feature_key = '.'.join(bedrangecoordfile.split('.')[4:7])
        for line in FILE:
            sortnames.append(line.strip().split('\t')[0])
            try:
                # limit coverage values to only coordinates present in bed range file
                d[line.strip().split('\t')[0]]
            except:
                output.append(line.strip())
                # d[line.strip().split('\t')[0]] = line.strip().split('\t')
                # if there is no coverage add zero coverage entries from the bedrange file to d_for
                d_for[line.strip().split('\t')[0]] = int(line.strip().split('\t')[1])
                d_rev[line.strip().split('\t')[0]] = int(line.strip().split('\t')[2])
                feature_for_depth.append(0)
                feature_rev_depth.append(0)
            else:
                # if there is coverage information for this postion then output it
                output.append('\t'.join(d[line.strip().split('\t')[0]]))
                feature_for_depth.append(int(d[line.strip().split('\t')[0]][1]))
                feature_rev_depth.append(int(d[line.strip().split('\t')[0]][2]))

    print('###################################')
    print('###################################')
    print('###################################')

    total_per_base_depth = sum(feature_for_depth) + sum(feature_rev_depth)
    if total_per_base_depth > 50:
        print(f'sum all base depths = {total_per_base_depth}')
        avg_for = sum(feature_for_depth) / len(feature_for_depth)
        avg_rev = sum(feature_rev_depth) / len(feature_rev_depth)
        avg_strandedness = avg_for / (avg_for + avg_rev)
        print('Average strandedness per feature, for features with > 25  = sum(for) + sum(rev)')
        print('0 = all forward reads, 1 = all reverse reads, 0.5 = half forward half reverse')
        print(f'Average strandedness of this feature = {avg_strandedness}')
    else:
        print('depth was too low to calculate average strandedness')
        print(f'sum all base depths = {total_per_base_depth}')


    print('###################################')
    print('###################################')
    print('###################################')

    outfile = forandrevcountfile + '.out'
    with open(outfile, 'w') as OUT:
        OUT.write('\n'.join(output))

    outfileprefix = 'SingleStrandRNA.FvsR.Depth.CompleteXaxis'
    analysis = 'Single Stranded RNA Analysis'

    # names should already be sorted, but in case you need to do it, you can use this below
    # sortnames = natsorted(names, key=lambda y: y.lower())  # or # natsorted(x, alg=ns.IGNORECASE)  # or alg=ns.IC
    # plot forward reads with positive value
    path, file = os.path.split(bedrangecoordfile)
    plot_scatter(list(range(len(sortnames))), [d_for[n] for n in sortnames],
                     os.path.join(path, '.'.join(file.split('.')[:-1] + [outfileprefix, 'pdf'])), analysis,
                     2, 1, '_'.join(file.split('.')[:4]), 'Genomic Position', figlength=20, figheight=10)

    path, file = os.path.split(bedrangecoordfile)
    # plot reverse reads with negative value
    plot_scatter(list(range(len(sortnames))), [-1 * d_rev[n] for n in sortnames],
                     os.path.join(path, '.'.join(file.split('.')[:-1] + [outfileprefix, 'pdf'])), analysis,
                     2, 2, '_'.join(file.split('.')[:4]), 'Genomic Position', figlength=20, figheight=10)
    print('done plotting')

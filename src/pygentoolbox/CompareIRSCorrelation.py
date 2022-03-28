

def plot_sns_grouped_histogram(x, y, xlab='', ylab='', figtitle='', outpath='', labels=[]):
    import os
    import seaborn as sns
    import matplotlib.pyplot as plt
    import numpy as np

    # make grouped histogram chart
    fig, axes = plt.subplots(figsize=(8.27, 11.7))
    for a in zip([x, y]):
        sns.distplot(a, ax=axes)  # , label=labels, kde=False # bins=range(1, 110, 10)
    # ax.set_xlim([0, 100])
    if figtitle != '':
        plt.suptitle(figtitle)
    if (xlab != '') and (ylab != ''):
        axes.set(xlabel=xlab, ylabel=ylab)
    if outpath != '':
        fig.savefig(outpath)
    if labels != []:
        plt.legend(labels)
    plt.show()
    plt.close()


def plot_sns_scatter_and_kde_heatmap(x, y, xlab='', ylab='', figtitle='', outpath='', colorcategories=[], nbins=100):
    import seaborn as sns
    import matplotlib.pyplot as plt
    from scipy.stats.stats import pearsonr
    from scipy.stats import linregress

    # x = np.array(l1)
    # y = np.array(l2)
    a4size = (11.7,8.27)
    # set fig and axes
    fig, axes = plt.subplots(1, 2, sharex=True, figsize=(8, 4))  # figsize(x,x)

    # plot first subplot scatter to axis 0
    if figtitle != '':
        plt.suptitle(figtitle)
    if (xlab != '') and (ylab != ''):
        axes[0].set(xlabel=xlab, ylabel=ylab)
    if colorcategories != []:
        sns.scatterplot(x, y, hue=colorcategories, legend="brief", ax=axes[0])
    else:
        sns.scatterplot(x, y, ax=axes[0])

    # plot second subplot density to axis 1
    # df = pd.DataFrame(list(zip(x, y, colorcategories)), columns=[xlab, ylab, 'color'])
    sns.kdeplot(x, y, levels=nbins, cmap="mako", ax=axes[1])  # cbar=True, # fill=True, thresh=0, (these arg should work, but they are not utilized by contour)
    if (xlab != '') and (ylab != ''):
        axes[1].set(xlabel=xlab, ylabel=ylab)
    # now both subplots finished so show them and save them
    if outpath != '':
        fig.savefig(outpath)
    plt.show()
    plt.close()

    print('Simple Linear Regression:')
    print(linregress(x, y))
    # print('Pearson Correlation:')
    # print(pearsonr(x, y))

def plot_sns_scatter(x, y, xlab='', ylab='', figtitle='', outpath='', colorcategories=[], nbins=100):
    import seaborn as sns
    import matplotlib.pyplot as plt
    from scipy.stats import linregress

    print('Simple Linear Regression:')
    print(linregress(x, y))

    a4size = (11.7,8.27)
    # set fig and axes
    fig, axes = plt.subplots(figsize=(8, 4))  # figsize(x,x)

    # plot first subplot scatter to axis 0
    if figtitle != '':
        plt.suptitle(figtitle)
    if (xlab != '') and (ylab != ''):
        axes.set(xlabel=xlab, ylabel=ylab)

    if colorcategories != []:
        sns.scatterplot(x, y, hue=colorcategories, legend="brief", ax=axes)
    else:
        sns.scatterplot(x, y, ax=axes)

    if outpath != '':
        fig.savefig(outpath)
    plt.show()
    plt.close()


def main(multiirstablefile, singleirsfile, colname='CAF1'):
    import pandas as pd

    # df = pd.read_csv(multiirstablefile, delimiter='\t')
    # mydf = pd.read_csv(singleirsfile, delimiter='\t')
    # x = df[colname].tolist()
    # y = mydf.iloc[:, 1].tolist()  # first column should be IRS name, second should be IRS

    with open(multiirstablefile, 'r') as FILE:
        header = FILE.readline()
        for count, elem in enumerate(header.strip().split('\t')):
            if elem == colname:
                colnumber = count
        d = {line.strip().split('\t')[0]: line.strip().split('\t')[colnumber] for line in FILE}

    with open(singleirsfile, 'r') as FILE:
        names = [line.strip().split('\t')[0] for line in FILE]
    with open(singleirsfile, 'r') as FILE:
        d2 = {line.strip().split('\t')[0]: line.strip().split('\t')[1] for line in FILE}

    print('Number IRS in my file: %d' % len([names]))
    x = []
    y = []
    for k in names:
        try:
            d[k]
        except:
            pass
        else:
            y.append(float(d[k]))
            x.append(float(d2[k]))
    # y = [d[k] for k in names]
    # x = [d2[k] for k in names]
    print('number of IRS in both files: %d' % len(y))

    plot_sns_scatter(x, y, 'Victor_Ptcaf1_IRS', 'Xyrus_Ptcaf1_IRS', 'Ptcaf1 IRS Values')
    # plot_sns_scatter_and_kde_heatmap(x, y, 'Xyrus_Ptcaf1_IRS', 'Victor_Ptcaf1_IRS', 'Ptcaf1 IRS Values')
    plot_sns_grouped_histogram(x, y, 'Victor_Ptcaf1_IRS', 'Xyrus_Ptcaf1_IRS', 'Ptcaf1 IRS Values', labels=['Xyrus_IRS', 'Victor_IRS'])

    print('#### FIN ####')
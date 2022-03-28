''' idea of the program '''
''' mac prep DNA sequenced, and IRS are calculated '''
''' nucleosome prep DNA sequenced (MNase treated Mac prep DNA), and IRS are calculated '''
''' assumption, we can only loose fragments during MNase digestions '''
# if log2(IRS.mp / IRS.np) is positive, then loss of retained IES fragments during MNase digestion
# if log2(IRS.mp / IRS.np) is negative, then loss of excised IES fragments during MNase digestion
#
# IRS.mp / IRS.np = ###
# if < 1.0 , then loss of retained IES fragments during MNase digestion
# if > 1.0, then loss of excised IES fragments during MNase digestion
# if = 1.0, then no change
#
# if IRS.mp = 0.6
# if IRS.np = 0.4
# IRS.mp / IRS.np = 0.6 / 0.4 = 1.5
# log2 fold change = log2(IRS.mp / IRS.np) or log2(1.5) = log(1.5) / log(2) = 0.58 = positive fold change ==
### then loss of retained IES fragments during MNase digestion
# inverse of log2(IRS.mp / IRS.np) == 1 / log2(IRS.mp / IRS.np) == 2^log2(IRS.mp / IRS.np)
#
# if IRS.mp = 0.6
# if IRS.np = 0.8
# IRS.mp / IRS.np = 0.6 / 0.8 = 0.75
# log2 fold change = log2(IRS.mp / IRS.np) or log2(0.75) = log(0.75) / log(2) = -0.42 = negative fold change ==
### then loss of excised IES fragments during MNase digestion

### log2(ΔIRS) or log(#)/log(2) == log.base.2(#)


def ks_test(x, y):
    # x is list
    # y is list

    from scipy.stats import ks_2samp

    ksresult = ks_2samp(x, y)
    # ksresult[0] is d statistic
    # ksresult[1] is p-value

    return ksresult


def organize_data_for_two_sample_ks_test(alllog2fc, orderedcommonies, dlencategorical):
    print('\n####\n####\nCalculate Two-sample Kolmogorov-Smirnov Test (KS Test)')
    # for each KD compare log2FC(irs.mp/irs.np) from KD_PGM to EV_PGM
    # EV_log2fc is list of log2FCs, assumes EV_PGM is first file specified in input
    evlog2fc = [alllog2fc[0][k] for k in orderedcommonies]
    print('KS-test of all IESs that pass filtration')
    for kdfc in alllog2fc[1:]:
        kdlog2fc = [kdfc[k] for k in orderedcommonies]
        ksresult = ks_test(evlog2fc, kdlog2fc)
        print(ksresult)

    print('KS-test of SHORT IESs that pass filtration')
    shorties = [k for k in orderedcommonies if dlencategorical[k] == 'short']
    evlog2fc = [alllog2fc[0][k] for k in shorties]
    for kdfc in alllog2fc[1:]:
        kdlog2fc = [kdfc[k] for k in shorties]
        ksresult = ks_test(evlog2fc, kdlog2fc)
        print(ksresult)

    mediumies = [k for k in orderedcommonies if dlencategorical[k] == 'medium']
    print('KS-test of MEDIUM IESs that pass filtration')
    evlog2fc = [alllog2fc[0][k] for k in mediumies]
    for kdfc in alllog2fc[1:]:
        kdlog2fc = [kdfc[k] for k in mediumies]
        ksresult = ks_test(evlog2fc, kdlog2fc)
        print(ksresult)

    longies = [k for k in orderedcommonies if dlencategorical[k] == 'long']
    print('KS-test of LONG IESs that pass filtration')
    evlog2fc = [alllog2fc[0][k] for k in longies]
    for kdfc in alllog2fc[1:]:
        kdlog2fc = [kdfc[k] for k in longies]
        ksresult = ks_test(evlog2fc, kdlog2fc)
        print(ksresult)


def seaborn_pairplot(filename):
    '''  give tab delimited file, 1st column is ies names, all other columns have header and values '''
    import pandas as pd
    # Seaborn visualization library
    import seaborn as sns

    df = pd.read_csv(filename, sep="\t")

    # Create the default pairplot
    sns_plt = sns.pairplot(df, hue='Length_Category', plot_kws={'edgecolor': 'k'})  #
    sns_plt.savefig(filename[:-len('.tsv')] + '.pairplot.pdf')


def log2_irs(dirsmac, dirsnuc):
    # irs_mac and irs_nuc are dictionaries of IRS values, keys are ies names, same IESs are present in each dictionary
    import math

    log2fc = {k: math.log2( dirsnuc[k] / dirsmac[k]) for k in list(dirsmac.keys())}

    return log2fc


def simulate_noise(irsmaclist=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9], depthlist=[20]*9, irsnuclist=[0.5]*9, noiselevel=0.5, iterations=10):
    ### not finished ###
    # irsmaclist is list containing irs values from mac prep (assumed the 'real' starting IRS)
    # depthlist contains the total depth values from nucprep for each IES
    # irsnuclist is list containing irs values from nuc prep
    # I assume number of fragments must have been higher in Mac prep before digestion
    # i increase the depth in the nuc prep by nucdepth*(1/noiselevel) (this is theoretical number of frags in mac)
    # i simulate random loss of ret and exc frags from theoretical mac depth (using mac IRS) to nucdepth values
    ## to see what IRS values are obtained simply due to noise
    # does the nucprep irs fall within or outside the 95% confidence interval (if outside then signal might not be due to noise)
    # noiselevel is a float (frequency). Will cause loss of this freq of fragments from total depth of IES.
    # noiselevel is low @ 0 (no noise), high at 1 (complete loss of all fragments due to noise).
    import random
    import seaborn as sns
    import matplotlib.pyplot as plt
    import math

    fclist = []
    allnoiseirs, allnoisefc = [], []
    x1, y, = [], []
    for irsmac, depthnuc, irsnuc in zip(irsmaclist, depthlist, irsnuclist):
        depth = depthnuc * (1 / noiselevel)
        ret = int(depth * irsmac)
        exc = int(depth - ret)
        lfc = math.log2(ret / exc)
        fclist.append(lfc)
        irsnoise, fcnoise = [], []
        for i in range(iterations):
            ies = [0] * exc + [1] * ret  # 1 is retained fragment 0 is excised
            for j in range(int(depth*noiselevel)):
                ies.pop(random.randint(0, len(ies) - 1))
            # collect noisy scores
            newret = sum(ies)
            newexc = len(ies) - newret
            irsnoise.append(newret / len(ies))
            x1.append(irsmac)
            y.append(newret / len(ies))
            if newexc != 0:
                newfc = newret / newexc
                if newfc != 0:
                    fcnoise.append(math.log2(newfc))
                else:
                    fcnoise.append(lfc)
            else:
                # avoid div by zero
                if newfc != 0:
                    newfc = newret / 1
                else:
                    fcnoise.append(lfc)
        allnoiseirs.append(irsnoise)
        allnoisefc.append(fcnoise)

    names = [str(irs) for irs in irsmaclist]
    namesfc = [str(round(fc, 1)) for fc in fclist]

    plot_sns_kde_heatmap_scatter_overlay(x1, y, irsmaclist, irsnuclist, 'IRS.MP', 'Noisy IRS.NP & Actual IRS.NP', 'Noise Effect On IRS.MP, Noiselevel=%.2f' % noiselevel)


def plot_sns_kde_heatmap_scatter_overlay(x1, y, x2, z, xlab='', ylab='', figtitle='', outpath='', colorcategories=[], nbins=100):
    import seaborn as sns
    import matplotlib.pyplot as plt
    # from scipy.stats.stats import pearsonr
    # from scipy.stats import linregress

    a4size = (11.7,8.27)
    # set fig and axes
    fig, axes = plt.subplots(sharex=True, figsize=(8, 8))  # 1, 2, # figsize(x,x)

    if figtitle != '':
        plt.suptitle(figtitle)
    if (xlab != '') and (ylab != ''):
        axes.set(xlabel=xlab, ylabel=ylab)

    # plot density plot
    # df = pd.DataFrame(list(zip(x, y, colorcategories)), columns=[xlab, ylab, 'color'])
    sns.kdeplot(x1, y, levels=nbins, cmap="mako", ax=axes)  # cbar=True, # fill=True, thresh=0, (these arg should work, but they are not utilized by contour)

    if colorcategories != []:
        sns.scatterplot(x2, z, hue=colorcategories, legend="brief", ax=axes)
    else:
        sns.scatterplot(x2, z, ax=axes)

    # now both plots finished so show them and save them
    if outpath != '':
        fig.savefig(outpath)
    plt.show()
    plt.close()

    # print('Simple Linear Regression:')
    # print(linregress(x, y))
    # print('Pearson Correlation:')
    # print(pearsonr(x, y))


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


def plot_kde_heatmap(l1, l2, nbins=300):
    import statistics
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.stats import kde

    x = np.array(l1)
    y = np.array(l2)
    # Evaluate a gaussian kde on a regular grid of nbins x nbins over data extents
    k = kde.gaussian_kde([x, y])
    xi, yi = np.mgrid[x.min():x.max():nbins * 1j, y.min():y.max():nbins * 1j]
    zi = k(np.vstack([xi.flatten(), yi.flatten()]))

    # Make the plot
    plt.pcolormesh(xi, yi, zi.reshape(xi.shape))
    # Add color bar
    # plt.pcolormesh(xi, yi, zi.reshape(xi.shape), cmap=plt.cm.Greens_r)
    plt.colorbar()
    plt.show()
    plt.close()


def plot_sns_scatter_hue(x, y, colorcategories, figtitle):
    import seaborn as sns
    import matplotlib.pyplot as plt
    from scipy.stats.stats import pearsonr

    sns.scatterplot(x, y, hue=colorcategories, legend="brief")
    plt.title(figtitle)
    plt.show()
    plt.close()
    print('Pearson Correlation:')
    print(pearsonr(x, y))


def show_bias_by_excision(allirsmacfilter, allirsnucfilter, alldfcmacfilter, alldfcnucfilter, alllog2fc, alldmacret, alldmacexc, alldnucret, alldnucexc, files, orderednames, dcolor, allresiduals, nbins):
    # plot correlation between KD_PGM.mp and KD_PGM.np
    # if very strong correlation, then the IRS in np is determined by mp IRS not by MNase digestion
    import os
    import math

    print('\n#####\nRepresenting the IES excision bias with non-normalized values, and the normalized values')
    colorcategories = [dcolor[k] for k in orderednames]

    print('\n#####\nPlotting correlations between EV_PGM and first KD_PGM')

    path, file = os.path.split(files[0])
    pathkd, filekd = os.path.split(files[1])
    prefix = '_'.join(file.split('_')[:2] + ['IRS.MP', 'vs'] + filekd.split('_')[:2] + ['IRS.MP', 'NotNorm'])
    outpath = os.path.join(path, prefix + '.pdf')
    print(prefix)
    xmacev = [allirsmacfilter[0][k] for k in orderednames]
    ymackd = [allirsmacfilter[1][k] for k in orderednames]
    plot_sns_scatter_and_kde_heatmap(xmacev, ymackd, 'IRS.MP.EV', 'IRS.MP.KD', prefix, outpath, colorcategories, nbins)
    # plot_sns_scatter_hue(xmacev, ymackd, colorcategories, prefix)
    # plot_kde_heatmap(xmacev, ymackd)

    path, file = os.path.split(files[0])
    pathkd, filekd = os.path.split(files[1])
    prefix = '_'.join(file.split('_')[:2] + ['IRS.NP', 'vs'] + filekd.split('_')[:2] + ['IRS.NP', 'NotNorm'])
    outpath = os.path.join(path, prefix + '.pdf')
    print(prefix)
    xnucev = [allirsnucfilter[0][k] for k in orderednames]
    ynuckd = [allirsnucfilter[1][k] for k in orderednames]
    plot_sns_scatter_and_kde_heatmap(xnucev, ynuckd, 'IRS.NP.EV', 'IRS.NP.KD', prefix, outpath, colorcategories, nbins)

    path, file = os.path.split(files[0])
    pathkd, filekd = os.path.split(files[1])
    prefix = '_'.join(file.split('_')[:2] + ['IRS.MP', 'vs'] + filekd.split('_')[:2] + ['IRS.NP', 'NotNorm'])
    outpath = os.path.join(path, prefix + '.pdf')
    print(prefix)
    plot_sns_scatter_and_kde_heatmap(xmacev, ynuckd, 'IRS.MP.EV', 'IRS.NP.KD', prefix, outpath, colorcategories, nbins)

    # path, file = os.path.split(files[0])
    # pathkd, filekd = os.path.split(files[1])
    # prefix = '_'.join(file.split('_')[:2] + ['FC.MP', 'vs'] + filekd.split('_')[:2] + ['FC.MP', 'NotNorm'])
    # print(prefix)
    # outpath = os.path.join(path, prefix + '.pdf')
    # xmacevfc = [alldfcmacfilter[0][k] for k in orderednames]
    # ymackdfc = [alldfcmacfilter[1][k] for k in orderednames]
    # plot_sns_scatter_and_kde_heatmap(xmacevfc, ymackdfc, 'FC.MP.EV', 'FC.NP.KD', prefix, outpath, colorcategories, nbins=200)
    # # plot_sns_scatter_hue(xmacevfc, ymackdfc, colorcategories, prefix)
    # # plot_kde_heatmap(xmacevfc, ymackdfc)

    # path, file = os.path.split(files[0])
    # pathkd, filekd = os.path.split(files[1])
    # prefix = '_'.join(file.split('_')[:2] + ['Log2FC.MP', 'vs'] + filekd.split('_')[:2] + ['Log2FC.MP', 'NotNorm'])
    # print(prefix)
    # outpath = os.path.join(path, prefix + '.pdf')
    # xmacevlfc = [math.log2(alldfcmacfilter[0][k]) for k in orderednames]
    # ymackdlfc = [math.log2(alldfcmacfilter[1][k]) for k in orderednames]
    # plot_sns_scatter_and_kde_heatmap(xmacevlfc, ymackdlfc, 'Log2(FC.MP.EV)', 'Log2(FC.MP.KD)', prefix, outpath, colorcategories, nbins)
    #
    # path, file = os.path.split(files[0])
    # pathkd, filekd = os.path.split(files[1])
    # prefix = '_'.join(file.split('_')[:2] + ['Log2FC.MP', 'vs'] + filekd.split('_')[:2] + ['Log2FC.NP', 'NotNorm'])
    # print(prefix)
    # outpath = os.path.join(path, prefix + '.pdf')
    # ynuckdlfc = [math.log2(alldfcnucfilter[1][k]) for k in orderednames]
    # plot_sns_scatter_and_kde_heatmap(xmacevlfc, ynuckdlfc, 'Log2(FC.MP.EV)', 'Log2(FC.NP.KD)', prefix, outpath, colorcategories, nbins)
    print('\n#####\nPlotting correlations within each knockdown (EV_PGM, first KD_PGM, second KD_PGM...)')

    for dmacirs, dnucirs, dfcmac, dfcnuc, dlog2fc, dmacret, dmacexc, dnucret, dnucexc, dres, nucfilename in zip(allirsmacfilter, allirsnucfilter, alldfcmacfilter, alldfcnucfilter, alllog2fc, alldmacret, alldmacexc, alldnucret, alldnucexc, allresiduals, files):

        # plot_sns_scatter_hue(xmac, ynuc, colorcategories, prefix)
        # plot_kde_heatmap(xmac, ynuc)
        mac = [dmacirs[k] for k in orderednames]
        nuc = [dnucirs[k] for k in orderednames]
        fcmac = [math.log2(dfcmac[k]) for k in orderednames]
        fcnuc = [math.log2(dfcnuc[k]) for k in orderednames]
        lfc = [dlog2fc[k] for k in orderednames]
        rawfcmac = [dfcmac[k] for k in orderednames]
        rawfcnuc = [dfcnuc[k] for k in orderednames]
        macret = [dmacret[k] for k in orderednames]
        macexc = [dmacexc[k] for k in orderednames]
        nucret = [dnucret[k] for k in orderednames]
        nucexc = [dnucexc[k] for k in orderednames]
        lmactotaldepth = [dmacret[k] + dmacexc[k] for k in orderednames]  # math.log2()
        lnuctotaldepth = [dnucret[k] + dnucexc[k] for k in orderednames]  # math.log2()
        res = [dres[k] for k in orderednames]

        path, nucfile = os.path.split(nucfilename)
        prefix = '_'.join(nucfile.split('_')[:2] + ['IRS.MP', 'vs', 'IRS.NP', 'NotNorm'])
        print(prefix)
        outpath = os.path.join(path, prefix + '.pdf')
        irsnucminusmac = [dnucirs[k] - dmacirs[k] for k in orderednames]
        plot_sns_scatter_and_kde_heatmap(mac, nuc, 'IRS.MP', 'IRS.NP', prefix, outpath, colorcategories, nbins)

        # this is still the sum of the MP IRS and the MNase digestion, not MNase digestion only. so not normalized.
        # first two plots show dispraportionate effect of noise (equal loss of fragments) on FC for different IRS values

        # prefix = '_'.join(nucfile.split('_')[:2] + ['FC.MP', 'vs', 'Log2FC.MP', 'Transformed'])
        # print(prefix)
        # outpath = os.path.join(path, prefix + '.pdf')
        # plot_sns_scatter_and_kde_heatmap(rawfcmac, fcmac, 'FC.MP', 'Log2(FC.MP)', prefix, outpath, colorcategories, nbins)
        # 
        # prefix = '_'.join(nucfile.split('_')[:2] + ['FC.NP', 'vs', 'Log2FC.NP', 'Transformed'])
        # print(prefix)
        # outpath = os.path.join(path, prefix + '.pdf')
        # plot_sns_scatter_and_kde_heatmap(rawfcnuc, fcnuc, 'FC.NP', 'Log2(FC.NP)', prefix, outpath, colorcategories, nbins)

        # prefix = '_'.join(nucfile.split('_')[:2] + ['IRS.MP', 'vs', 'FC.MP', 'NotNorm'])
        # print(prefix)
        # outpath = os.path.join(path, prefix + '.pdf')
        # plot_sns_scatter_and_kde_heatmap(mac, rawfcmac, 'IRS.MP', 'FC.MP', prefix, outpath, colorcategories, nbins)
        #
        # prefix = '_'.join(nucfile.split('_')[:2] + ['IRS.NP', 'vs', 'FC.NP', 'NotNorm'])
        # print(prefix)
        # outpath = os.path.join(path, prefix + '.pdf')
        # plot_sns_scatter_and_kde_heatmap(nuc, rawfcnuc, 'IRS.NP', 'FC.NP', prefix, outpath, colorcategories, nbins)

        prefix = '_'.join(nucfile.split('_')[:2] + ['IRS.MP', 'vs', 'MP.Depth', 'NotNorm'])
        print(prefix)
        outpath = os.path.join(path, prefix + '.pdf')
        plot_sns_scatter_and_kde_heatmap(mac, lmactotaldepth, 'IRS.MP', 'MP.Depth', prefix, outpath, colorcategories, nbins)

        prefix = '_'.join(nucfile.split('_')[:2] + ['IRS.NP', 'vs', 'NP.Depth', 'NotNorm'])
        print(prefix)
        outpath = os.path.join(path, prefix + '.pdf')
        plot_sns_scatter_and_kde_heatmap(nuc, lnuctotaldepth, 'IRS.NP', 'NP.Depth', prefix, outpath, colorcategories, nbins)

        # prefix = '_'.join(nucfile.split('_')[:2] + ['IRS.MP', 'vs', 'MP.Depth.Retained', 'NotNorm'])
        # print(prefix)
        # outpath = os.path.join(path, prefix + '.pdf')
        # plot_sns_scatter_and_kde_heatmap(mac, macret, 'IRS.MP', 'MP.Depth.Retained', prefix, outpath, colorcategories, nbins)
        # 
        # prefix = '_'.join(nucfile.split('_')[:2] + ['IRS.NP', 'vs', 'NP.Depth.Retained', 'NotNorm'])
        # print(prefix)
        # outpath = os.path.join(path, prefix + '.pdf')
        # plot_sns_scatter_and_kde_heatmap(nuc, nucret, 'IRS.NP', 'NP.Depth.Retained', prefix, outpath, colorcategories, nbins)
        # 
        # prefix = '_'.join(nucfile.split('_')[:2] + ['IRS.MP', 'vs', 'MP.Depth.Excised', 'NotNorm'])
        # print(prefix)
        # outpath = os.path.join(path, prefix + '.pdf')
        # plot_sns_scatter_and_kde_heatmap(mac, macexc, 'IRS.MP', 'MP.Depth.Excised', prefix, outpath, colorcategories, nbins)
        # 
        # prefix = '_'.join(nucfile.split('_')[:2] + ['IRS.NP', 'vs', 'NP.Depth.Excised', 'NotNorm'])
        # print(prefix)
        # outpath = os.path.join(path, prefix + '.pdf')
        # plot_sns_scatter_and_kde_heatmap(nuc, nucexc, 'IRS.NP', 'NP.Depth.Excised', prefix, outpath, colorcategories, nbins)

        # prefix = '_'.join(nucfile.split('_')[:2] + ['Log2FC.MP', 'vs', 'MP.Depth', 'NotNorm'])
        # print(prefix)
        # outpath = os.path.join(path, prefix + '.pdf')
        # plot_sns_scatter_and_kde_heatmap(fcmac, lmactotaldepth,  'Log2(FC.MP)', 'MP.Depth', prefix, outpath, colorcategories, nbins)
        # 
        # prefix = '_'.join(nucfile.split('_')[:2] + ['Log2FC.NP', 'vs', 'NP.Depth', 'NotNorm'])
        # print(prefix)
        # outpath = os.path.join(path, prefix + '.pdf')
        # plot_sns_scatter_and_kde_heatmap(fcnuc, lnuctotaldepth,  'Log2(FC.NP)', 'NP.Depth', prefix, outpath, colorcategories, nbins)

        prefix = '_'.join(nucfile.split('_')[:2] + ['IRS.MP', 'vs', 'Log2FC.MP', 'NotNorm'])
        print(prefix)
        outpath = os.path.join(path, prefix + '.pdf')
        plot_sns_scatter_and_kde_heatmap(mac, fcmac, 'IRS.MP', 'Log2(FC.MP)', prefix, outpath, colorcategories, nbins)

        prefix = '_'.join(nucfile.split('_')[:2] + ['IRS.NP', 'vs', 'Log2FC.NP', 'NotNorm'])
        print(prefix)
        outpath = os.path.join(path, prefix + '.pdf')
        plot_sns_scatter_and_kde_heatmap(nuc, fcnuc, 'IRS.NP', 'Log2(FC.NP)', prefix, outpath, colorcategories, nbins)

        # prefix = '_'.join(nucfile.split('_')[:2] + ['IRS.MP', 'vs', 'Log2FC.NP', 'NotNorm'])
        # print(prefix)
        # outpath = os.path.join(path, prefix + '.pdf')
        # plot_sns_scatter_and_kde_heatmap(mac, fcnuc, 'IRS.MP', 'Log2(FC.NP)', prefix, outpath, colorcategories, nbins)
        # 
        # prefix = '_'.join(nucfile.split('_')[:2] + ['IRS.NP', 'vs', 'Log2FC.MP', 'NotNorm'])
        # print(prefix)
        # outpath = os.path.join(path, prefix + '.pdf')
        # plot_sns_scatter_and_kde_heatmap(nuc, fcmac, 'IRS.NP', 'Log2(FC.MP)', prefix, outpath, colorcategories, nbins)

        prefix = '_'.join(nucfile.split('_')[:2] + ['IRS.MP', 'vs', 'ResidualsOfIRS.NP-IRS.MP' 'Log2FC.MP', 'Norm'])
        print(prefix)
        outpath = os.path.join(path, prefix + '.pdf')
        plot_sns_scatter_and_kde_heatmap(mac, res, 'IRS.MP', 'ResidualsOfIRS.NP-IRS.MP', prefix, outpath, colorcategories, nbins)

        prefix = '_'.join(nucfile.split('_')[:2] + ['IRS.NP', 'vs', 'ResidualsOfIRS.NP-IRS.MP' 'Log2FC.MP', 'Norm'])
        print(prefix)
        outpath = os.path.join(path, prefix + '.pdf')
        plot_sns_scatter_and_kde_heatmap(nuc, res, 'IRS.NP', 'ResidualsOfIRS.NP-IRS.MP', prefix, outpath, colorcategories, nbins)

        prefix = '_'.join(nucfile.split('_')[:2] + ['IRS.MP', 'vs', 'IRS.NP-IRS.MP' 'Log2FC.MP', 'Norm'])
        print(prefix)
        outpath = os.path.join(path, prefix + '.pdf')
        plot_sns_scatter_and_kde_heatmap(mac, irsnucminusmac, 'IRS.MP', 'IRS.NP-IRS.MP', prefix, outpath, colorcategories, nbins)

        prefix = '_'.join(nucfile.split('_')[:2] + ['IRS.NP', 'vs', 'IRS.NP-IRS.MP' 'Log2FC.MP', 'Norm'])
        print(prefix)
        outpath = os.path.join(path, prefix + '.pdf')
        plot_sns_scatter_and_kde_heatmap(nuc, irsnucminusmac, 'IRS.NP', 'IRS.NP-IRS.MP', prefix, outpath, colorcategories, nbins)

        prefix = '_'.join(nucfile.split('_')[:2] + ['IRS.MP', 'vs', 'Log2FC.NPvsMP', 'Norm'])
        print(prefix)
        outpath = os.path.join(path, prefix + '.pdf')
        plot_sns_scatter_and_kde_heatmap(mac, lfc, 'IRS.MP', 'Log2(FC.NP/FC.MP)', prefix, outpath, colorcategories, nbins)

        prefix = '_'.join(nucfile.split('_')[:2] + ['IRS.NP', 'vs', 'Log2FC.NPvsMP', 'Norm'])
        print(prefix)
        outpath = os.path.join(path, prefix + '.pdf')
        plot_sns_scatter_and_kde_heatmap(nuc, lfc, 'IRS.NP', 'Log2(FC.NP/FC.MP)', prefix, outpath, colorcategories, nbins)

        prefix = '_'.join(nucfile.split('_')[:2] + ['Log2FC.NPvsMP', 'vs', 'MP.Depth', 'Norm'])
        print(prefix)
        outpath = os.path.join(path, prefix + '.pdf')
        plot_sns_scatter_and_kde_heatmap(lfc, lmactotaldepth, 'Log2(FC.NP/FC.MP)', 'MP.Depth', prefix, outpath, colorcategories, nbins)

        prefix = '_'.join(nucfile.split('_')[:2] + ['Log2FC.NPvsMP', 'vs', 'NP.Depth', 'Norm'])
        print(prefix)
        outpath = os.path.join(path, prefix + '.pdf')
        plot_sns_scatter_and_kde_heatmap(lfc, lnuctotaldepth, 'Log2(FC.NP/FC.MP)', 'NP.Depth', prefix, outpath, colorcategories, nbins)

        # prefix = '_'.join(nucfile.split('_')[:2] + ['IRS.NP-IRS.MP', 'vs', 'Log2FC.MP', 'Norm'])
        # print(prefix)
        # outpath = os.path.join(path, prefix + '.pdf')
        # plot_sns_scatter_and_kde_heatmap(irsnucminusmac, fcmac, 'IRS.NP-IRS.MP', 'Log2(FC.MP)', prefix, outpath, colorcategories, nbins)
        #
        # prefix = '_'.join(nucfile.split('_')[:2] + ['IRS.NP-IRS.MP', 'vs', 'Log2FC.NP', 'Norm'])
        # print(prefix)
        # outpath = os.path.join(path, prefix + '.pdf')
        # plot_sns_scatter_and_kde_heatmap(irsnucminusmac, fcnuc, 'IRS.NP-IRS.MP', 'Log2(FC.NP)', prefix, outpath, colorcategories, nbins)

        prefix = '_'.join(nucfile.split('_')[:2] + ['IRS.NP-IRS.MP', 'vs', 'Log2FC.NPvsMP', 'Norm'])
        print(prefix)
        outpath = os.path.join(path, prefix + '.pdf')
        plot_sns_scatter_and_kde_heatmap(irsnucminusmac, lfc, 'IRS.NP-IRS.MP', 'Log2(FC.NP/FC.MP)', prefix, outpath, colorcategories, nbins)

        # prefix = '_'.join(nucfile.split('_')[:2] + ['Log2FC.MP', 'vs', 'Log2FC.NP', 'NotNorm'])
        # print(prefix)
        # outpath = os.path.join(path, prefix + '.pdf')
        # plot_sns_scatter_and_kde_heatmap(fcmac, fcnuc, 'Log2(FC.MP)', 'Log2(FC.NP)', prefix, outpath, colorcategories, nbins)

        # prefix = '_'.join(nucfile.split('_')[:2] + ['Log2FC.MP', 'vs', 'Log2FC.NPvsMP', 'Norm'])
        # print(prefix)
        # outpath = os.path.join(path, prefix + '.pdf')
        # plot_sns_scatter_and_kde_heatmap(fcmac, lfc, 'Log2(FC.MP)', 'Log2(FC.NP/FC.MP)', prefix, outpath, colorcategories, nbins)
        #
        # prefix = '_'.join(nucfile.split('_')[:2] + ['Log2FC.NP', 'vs', 'Log2FC.NPvsMP', 'Norm'])
        # print(prefix)
        # outpath = os.path.join(path, prefix + '.pdf')
        # plot_sns_scatter_and_kde_heatmap(fcnuc, lfc, 'Log2(FC.NP)', 'Log2(FC.NP/FC.MP)', prefix, outpath, colorcategories, nbins)
        # print('Done with within knockdown plots for this knockdown\n#####')


def plot_sns_grouped_histogram(l, xlab='', ylab='', figtitle='', outpath='', labels=[]):
    # l is list of lists, l = zip([x, y, z, ...])

    import os
    import seaborn as sns
    import matplotlib.pyplot as plt
    import numpy as np

    # make grouped histogram chart
    fig, axes = plt.subplots(figsize=(8.27, 11.7))
    for a in l:
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


def plot_irs_histogram(dirs, irsfile, coveragelimit, irslimit, binnumber=50):
    import os
    import matplotlib.pyplot as plt

    path, file = os.path.split(irsfile)
    fileprefix = file.split('.')[0]
    outpath = os.path.join(path, '%s.Alt.hist.CovLim%d.IrsLim%.2f.pdf' % (fileprefix, coveragelimit, irslimit))
    plt.figure(figsize=(11.7, 8.27))  # figsize=(3,4))
    plt.hist(list(dirs.values()), bins=binnumber, edgecolor='black', linewidth=1.2)
    plt.title('%s IRS Alternative Histogram' % fileprefix)
    plt.xlabel('IRS Alternative')
    plt.ylabel('Frequency')
    plt.savefig(outpath)
    plt.show()
    plt.close()


def limit_ies_by_dkd(allirsmac, limitiesbydkdinverse=False, dkddiffthreshold=0.2):
    # allirsmac is list of dictionaries. Element is dictionary. Each dict key is IES name, value is IRS
    # if dkddiffthreshold is positive (or zero) filters for KD IRS.mp >= EV IRS.mp
    # if dkddiffthreshold is negative filters for KD IRS.mp <= EV IRS.mp
    # function will also only keep IESs that passed IRS filters in EV_PGM and each KD_PGM (different # for diff KD possible)
    # assumes first element in list is dictionary for mac prep of EV_PGM
    evmac = allirsmac[0]
    # alldkddep is list of dictionaries. each dict only has IESs also present in evmac.
    ## alldkddep also requires KD_PGM IRS must be >= EV_PGM(evmac) IRS by dkddiffthreshold
    alldkddep = [evmac]
    # skip EV_PGM (element zero)
    if limitiesbydkdinverse is False:
        if dkddiffthreshold >= 0:
            print('dkddiff is positive or zero so filtering for KD IRS.mp >= EV IRS.mp')
        elif dkddiffthreshold < 0:
            print('dkddiff is negative so filtering for KD IRS.mp <= EV IRS.mp')
    elif limitiesbydkdinverse is True:
        print('Finding KD independent IESs')

    for dirs in allirsmac[1:]:
        print('Number of IESs for KD_PGM before DKD dependent filter')
        print(len(list(dirs.keys())))
        d = {}
        for k, v in evmac.items():
            try:
                dirs[k]
            except:
                pass
            else:
                diff = dirs[k] - v
                if limitiesbydkdinverse is False:
                    if dkddiffthreshold >= 0:
                        if diff >= dkddiffthreshold:
                            d[k] = dirs[k]
                    elif dkddiffthreshold < 0:
                        if diff <= dkddiffthreshold:
                            d[k] = dirs[k]
                elif limitiesbydkdinverse is True:
                    absdiffthreshold = abs(dkddiffthreshold)
                    if (diff < absdiffthreshold) and (diff > -absdiffthreshold):
                        d[k] = dirs[k]

        print('Number of IESs for KD_PGM after DKD dependent filter')
        print(len(list(d.keys())))
        alldkddep.append(d)

    return alldkddep


def sturges_rule(n):
    import math
    # n is number of observations
    numbins = round(1 + 3.322 * math.log(n, 10))

    print('Number if bins claculated by Sturge\'s Rule : %d' % numbins)

    return numbins


def filter_irs(irsfile, retexclimit=1, coveragelimit=1, irslimit=0.01):
    # coveragelimit has default value of 1 so that denominator if IRS must be > 0
    dirs = {}
    dfc = {}
    dret, dexc = {}, {}
    with open(irsfile, 'r') as FILE:
        for line in FILE:
            # if sum of leftmost, rightmost, and both (left&right) molecules supporting retention plus the number of ...
            # molecules support excision are >= limit (and integer value)
            numreadsretained = int(line.strip().split('\t')[2]) + int(line.strip().split('\t')[3]) + int(line.strip().split('\t')[4])
            numreadsexcised = int(line.strip().split('\t')[5])
            if (numreadsexcised >= retexclimit) and (numreadsretained >= retexclimit):
                irs = float(line.strip().split('\t')[1])
                iesname = line.strip().split('\t')[0]
                if (numreadsretained + numreadsexcised >= coveragelimit) and (irs >= irslimit):
                    dirs[iesname] = irs
                    dfc[iesname] = numreadsretained / numreadsexcised
                    dret[iesname] = numreadsretained
                    dexc[iesname] = numreadsexcised

    print('#####  For file: %s #####' % irsfile)
    print('%d IES elements with counts >= than %d and IRS >=%.2f' % (len(list(dirs.keys())), coveragelimit, irslimit))
    print('%d IES elements with IRS less than 0.1' % len([k for k, v in dirs.items() if v < 0.1]))
    print('%d IES elements with IRS greater than 0.1 and less than 0.25' % len([k for k, v in dirs.items() if (v > 0.1) and (v < 0.25)]))
    print('%d IES elements with IRS above 0.25 and less than 0.75' % len([k for k, v in dirs.items() if (v < 0.75) and (v > 0.25)]))
    print('%d IES elements with IRS greater than 0.75' % len([k for k, v in dirs.items() if v > 0.75]))
    print('%d IES elements with IRS greater than 0.95' % len([k for k, v in dirs.items() if v > 0.95]))

    binnumber = sturges_rule(len(list(dirs.keys())))
    plot_irs_histogram(dirs, irsfile, coveragelimit, irslimit, binnumber)

    return dirs, dfc, dret, dexc


def read_gff3(GFF3file, features=['all']):
    print('Reading Gff3 file: %s' % GFF3file)
    d = {}
    with open(GFF3file, 'r') as FILE:
        for line in FILE:
            if (line[0] == '#') or (line.strip() == ''):
                pass
            elif features == ['all']:  # keep all lines
                key = line.strip().split('\t')[8].split(';')[0][len('ID='):]
                d[key] = line.strip().split('\t')
            elif line.strip().split('\t')[2] in features:  # keep lines only if in features
                key = line.strip().split('\t')[8].split(';')[0][len('ID='):]
                d[key] = line.strip().split('\t')
    print('Number of features: %d' % len(list(d.keys())))
    # print(list(d.keys())[:10])
    return d


def log2fc_mp_vs_np(MacPrepIRSFiles, NucPrepIRSFiles, iesannotationfile, limitiesfiles=[], limitiesbydkd=False, limitiesbydkdinverse=False, nbins=100, dkddiff=0.2, skdlimit=0.3, retexclimit=1, coveragelimit=1, irslimit=0.01):
    # output files written to directory of first nucleosome prep file in NucPrepIRSFiles
    # limitiesbydkdinverse # True or False # requires limitiesbydkd=True # True if you want to test IESs that are not dependent on KD: (IRS.MP in KD_PGM - IRS.MP in EV_PGM <= +dkddiff) and (IRS.MP in KD_PGM - IRS.MP in EV_PGM >= -dkddiff)
    # dkddiff  # float value # if dkddiff is positive, then: IRS.MP in KD_PGM - IRS.MP in EV_PGM >= dkddiff # if dkddiff is negative, then: IRS.MP in KD_PGM - IRS.MP in EV_PGM <= dkddiff
    import os
    import pandas as pd
    import math

    dgff3 = read_gff3(iesannotationfile, features=['internal_eliminated_sequence'])
    # dictionary with length of each IES, k is ies name, v is list of a line (for k) split by tabs from ies .gff3 file
    dlen = {k: len(v[8].split(';')[-1][len('sequence='):]) for k, v in dgff3.items()}
    print('Average IES length = %.2f' % (sum(list(dlen.values()))/len(list(dlen.keys()))))
    dlencategorical = {k: 'short' if v < 100 else 'medium' if (v >= 100) and (v <= 200) else 'long' if v > 200 else 'NA' for k, v in dlen.items()}
    #dlencategorical = {k: 0 if v < 100 else 1 if (v >= 100) and (v <= 200) else 2 if v > 200 else 9999 for k, v in dlen.items()}

    allirsmac = []
    allirsnuc = []
    alldfcmac, alldfcnuc = [], []
    alldmacret, alldmacexc, alldnucret, alldnucexc = [], [], [], []
    for macfile, nucfile in zip(MacPrepIRSFiles, NucPrepIRSFiles):
        dirsmac, dfcmac, dmacret, dmacexc = filter_irs(macfile, retexclimit, coveragelimit, irslimit)
        dirsnuc, dfcnuc, dnucret, dnucexc = filter_irs(nucfile, retexclimit, coveragelimit, irslimit)
        allirsmac.append(dirsmac)
        allirsnuc.append(dirsnuc)
        alldfcmac.append(dfcmac)
        alldfcnuc.append(dfcnuc)
        alldmacret.append(dmacret)
        alldmacexc.append(dmacexc)
        alldnucret.append(dnucret)
        alldnucexc.append(dnucexc)

    if limitiesbydkd == True:
        # send list of dictionaries allirsmac
        alldkddepirs = limit_ies_by_dkd(allirsmac, limitiesbydkdinverse, dkddiff)

    # this would use an IRS file from a single KD (like Ptcaf1) to further filter the IESs dependent on a single KD
    if len(limitiesfiles) > 0:
        alllimits = []
        for limitfile in limitiesfiles:
            dlimit, dlimitfc, dretlimit, dexclimit = filter_irs(limitfile, retexclimit, coveragelimit, skdlimit)
            alllimits.append(dlimit)
            print('%d = Number of dependent IESs from limitfile: %s' % (len(list(dlimit.keys())), limitfile))

    # filter IESs to those that pass the filter in all files (so we have the same IES in each file and the same number)
    # dcommonies, key is ies name, value is number of times that key is present after filtration for all input files
    dcommonies = {}
    if (len(limitiesfiles) == 0) and (limitiesbydkd is False):
        for dirs in allirsmac + allirsnuc:
            for k, v in dirs.items():
                dcommonies[k] = dcommonies.setdefault(k, 0) + 1
        # list listing names of IES that passed filtration in all input files
        commonies = [k for k, v in dcommonies.items() if v == len(MacPrepIRSFiles + NucPrepIRSFiles)]
    elif (len(limitiesfiles) == 0) and (limitiesbydkd is True):
        for dirs in alldkddepirs + allirsnuc:
            for k, v in dirs.items():
                dcommonies[k] = dcommonies.setdefault(k, 0) + 1
        # list listing names of IES that passed filtration in all input files
        commonies = [k for k, v in dcommonies.items() if v == len(MacPrepIRSFiles + NucPrepIRSFiles)]
    elif (len(limitiesfiles) > 0) and (limitiesbydkd is False):
        for dirs in allirsmac + allirsnuc + alllimits:
            for k, v in dirs.items():
                dcommonies[k] = dcommonies.setdefault(k, 0) + 1
        # list listing names of IES that passed filtration in all input files
        commonies = [k for k, v in dcommonies.items() if v == len(MacPrepIRSFiles + NucPrepIRSFiles + alllimits)]
    elif (len(limitiesfiles) > 0) and (limitiesbydkd is True):
        for dirs in alldkddepirs + allirsnuc + alllimits:
            for k, v in dirs.items():
                dcommonies[k] = dcommonies.setdefault(k, 0) + 1
        # list listing names of IES that passed filtration in all input files
        commonies = [k for k, v in dcommonies.items() if v == len(MacPrepIRSFiles + NucPrepIRSFiles + alllimits)]

    # order commonies list to the order of the first macprep input file
    with open(MacPrepIRSFiles[0], 'r') as FILE:
        orderedcommonies = [line.strip().split('\t')[0] for line in FILE if line.strip().split('\t')[0] in commonies]

    if len(limitiesfiles) == 0:
        print('%d = number of IES that passed filtration in all %d input files' % (len(orderedcommonies),
                                                                               len(MacPrepIRSFiles + MacPrepIRSFiles)))
    elif len(limitiesfiles) > 0:
        print('%d = number of IES that passed filtration in all %d input files' % (len(orderedcommonies),
                                                            len(MacPrepIRSFiles + MacPrepIRSFiles + limitiesfiles)))

    allirsmacfilter = []
    allirsnucfilter = []
    alldfcmacfilter, alldfcnucfilter = [], []
    # these are list of dictionaries that only contain IRS values for IES that pass filtration in all input files
    for dirsmac, dirsnuc, dfcmac, dfcnuc in zip(allirsmac, allirsnuc, alldfcmac, alldfcnuc):
        dirsmacfilter = {k: dirsmac[k] for k in orderedcommonies}
        dirsnucfilter = {k: dirsnuc[k] for k in orderedcommonies}
        dfcmacfilter = {k: dfcmac[k] for k in orderedcommonies}
        dfcnucfilter = {k: dfcnuc[k] for k in orderedcommonies}
        allirsmacfilter.append(dirsmacfilter)
        allirsnucfilter.append(dirsnucfilter)
        alldfcmacfilter.append(dfcmacfilter)
        alldfcnucfilter.append(dfcnucfilter)

    for dirsmac, dirsnuc, macfile, nucfile in zip(allirsmac, allirsnuc, MacPrepIRSFiles, NucPrepIRSFiles):
        x = [dirsmac[k] for k in orderedcommonies]
        y = [dirsnuc[k] for k in orderedcommonies]
        peeath1, macfilename = os.path.split(macfile)
        peeath2, nucfilename = os.path.split(nucfile)
        label = ['%s_IRS.MP' % '.'.join(macfilename.split('_')[:2]), '%s_IRS.NP' % '.'.join(nucfilename.split('_')[:2])]
        l = zip([x, y])
        plot_sns_grouped_histogram(l, xlab='', ylab='', figtitle='', outpath='', labels=label)

    allirslog2fc = []
    alllog2fc = []
    for dirsmacfilter, dirsnucfilter, dfcmacfilter, dfcnucfilter in zip(allirsmacfilter, allirsnucfilter, alldfcmacfilter, alldfcnucfilter):
        # calculates log2(nuc / mac)
        d_irs_log2fc = log2_irs(dirsmacfilter, dirsnucfilter)
        dlog2fc = log2_irs(dfcmacfilter, dfcnucfilter)
        allirslog2fc.append(d_irs_log2fc)
        alllog2fc.append(dlog2fc)

    ### important calculation! ###
    # calculate differences between EV_PGM nucleosome content and each KD_PGM log2FC.KD(FC.np/FC.mp) - log2FC.EV(FC.np/FC.mp)
    alllog2fckdev = []
    evlog2fc = alllog2fc[0]
    for dlog2fc in alllog2fc[1:]:
        dlog2kdev = {k: dlog2fc[k] - evlog2fc[k] for k in orderedcommonies}
        alllog2fckdev.append(dlog2kdev)

    ##### fit 'independent' IRS.mp to IRS.np-IRS.mp #####
    allirsnucminusmac = []
    for dirsmac, dirsnuc in zip(allirsmacfilter, allirsnucfilter):
        dirsnucminusmac = {k: dirsnuc[k] - dirsmac[k] for k in orderedcommonies}
        allirsnucminusmac.append(dirsnucminusmac)

    print('Histograms of IRS.mp EV vs IRS.mp KD')
    x1 = [allirsmacfilter[0][k] for k in orderedcommonies]
    y1 = [allirsmacfilter[1][k] for k in orderedcommonies]
    label = ['EV_PGM', 'KD_PGM']
    l = zip([x1, y1])
    plot_sns_grouped_histogram(l, xlab='', ylab='', figtitle='', outpath='', labels=label)

    import statsmodels.formula.api as smf
    allresiduals = []
    for dirsmac, dirsnucminusmac in zip(allirsmacfilter, allirsnucminusmac):
        mac = [dirsmac[k] for k in orderedcommonies]
        nucminusmac = [dirsnucminusmac[k] for k in orderedcommonies]
        df = pd.DataFrame({'irsmac': mac, 'irsnucminusmac': nucminusmac})
        reg = smf.ols('irsnucminusmac ~ irsmac', data=df).fit()
        print(reg.summary())
        pred_val = reg.fittedvalues.copy()
        true_val = df['irsnucminusmac'].values.copy()
        residual = true_val - pred_val
        dres = {k: res for k, res in zip(orderedcommonies, residual)}  # assumed IESs maintained order
        allresiduals.append(dres)
    #####  #####

    ##### fit 'independent' IRS.mp to log2FC(irs.np/irs.mp) #####
    allfcresiduals = []
    for dirsmac, dlog2fc in zip(allirsmacfilter, alllog2fc):
        mac = [dirsmac[k] for k in orderedcommonies]
        lfc = [dlog2fc[k] for k in orderedcommonies]
        df = pd.DataFrame({'irsmac': mac, 'log2fcNPvsMP': lfc})
        reg = smf.ols('log2fcNPvsMP ~ irsmac', data=df).fit()
        print(reg.summary())
        pred_val = reg.fittedvalues.copy()
        true_val = df['log2fcNPvsMP'].values.copy()
        residual = true_val - pred_val
        dres = {k: res for k, res in zip(orderedcommonies, residual)}  # assumed IESs maintained order
        allfcresiduals.append(dres)
    #####  #####

    # calculate total depth then SIMULATE NOISE for each IES
    alldmactotal, alldnuctotal = [], []
    for dmacret, dmacexc, dnucret, dnucexc in zip(alldmacret, alldmacexc, alldnucret, alldnucexc):
        dmactotal = {k: dmacret[k] + dmacexc[k] for k in orderedcommonies}
        dnuctotal = {k: dnucret[k] + dnucexc[k] for k in orderedcommonies}
        alldmactotal.append(dmactotal)
        alldnuctotal.append(dnuctotal)
    for dirsmacfilter, dirsnucfilter, dmactotal, dnuctotal in zip(allirsmacfilter, allirsnucfilter, alldmactotal, alldnuctotal):
        macirs = [dirsmacfilter[k] for k in orderedcommonies]
        nucirs = [dirsnucfilter[k] for k in orderedcommonies]
        # mactotal = [dmactotal[k] for k in orderedcommonies]
        nuctotal = [dnuctotal[k] for k in orderedcommonies]
        simulate_noise(irsmaclist=macirs, depthlist=nuctotal, irsnuclist=nucirs, noiselevel=0.5, iterations=10)

    # not really desired output, but proof that there is some bias by excision that needs to be accounted for
    show_bias_by_excision(allirsmacfilter, allirsnucfilter, alldfcmacfilter, alldfcnucfilter, alllog2fc, alldmacret, alldmacexc, alldnucret, alldnucexc, NucPrepIRSFiles, orderedcommonies, dlencategorical, allresiduals, nbins)

    # plot histograms of log2fc(IRS_NP/IRS_MP)
    print('###\n###\nNow plotting log2FC(FC_NP/FC_MP)\n###\n###')
    for d_log2fc, irsfile in zip(alllog2fc, NucPrepIRSFiles):
        binnumber = sturges_rule(len(list(d_log2fc.keys())))
        plot_irs_histogram(d_log2fc, irsfile, coveragelimit, irslimit, binnumber)

    print('Histograms of log2(FC.np/FC.mp) for EV and KD')
    x1 = [alllog2fc[0][k] for k in orderedcommonies]
    y1 = [alllog2fc[1][k] for k in orderedcommonies]
    label = ['EV_PGM', 'KD_PGM']
    l = zip([x1, y1])
    plot_sns_grouped_histogram(l, xlab='', ylab='', figtitle='', outpath='', labels=label)

    print('Histograms of resiudals log2(FC.NP/FC.MP) for EV and KD')
    x2 = [allfcresiduals[0][k] for k in orderedcommonies]
    y2 = [allfcresiduals[1][k] for k in orderedcommonies]
    label = ['EV_PGM', 'KD_PGM']
    l = zip([x2, y2])
    plot_sns_grouped_histogram(l, xlab='', ylab='', figtitle='', outpath='', labels=label)

    plot_sns_scatter_and_kde_heatmap(x1, x2, 'Log2FCNPvsMP', 'ResidualsOfXaxis', 'EV_PGM')
    plot_sns_scatter_and_kde_heatmap(y1, y2, 'Log2FCNPvsMP', 'ResidualsOfXaxis', 'KD_PGM')

    print('Histograms of IRS.np - IRS.mp for EV and KD')
    x = [allirsnucminusmac[0][k] for k in orderedcommonies]
    y = [allirsnucminusmac[1][k] for k in orderedcommonies]
    label = ['EV_PGM', 'KD_PGM']
    l = zip([x, y])
    plot_sns_grouped_histogram(l, xlab='', ylab='', figtitle='', outpath='', labels=label)

    print('Histograms of resiudals Delta IRS for EV and KD')
    x = [allresiduals[0][k] for k in orderedcommonies]
    y = [allresiduals[1][k] for k in orderedcommonies]
    label = ['EV_PGM', 'KD_PGM']
    l = zip([x, y])
    plot_sns_grouped_histogram(l, xlab='', ylab='', figtitle='', outpath='', labels=label)

    # get file prefixes from nuc prep files
    nucfilenames = ['Name']
    for nucfile in NucPrepIRSFiles:
        path, file = os.path.split(nucfile)
        nucfilenames.append(file)
    header = '\t'.join(['_'.join(nucfile.split('_')[:2]) for nucfile in nucfilenames] + ['Length', 'Length_Category'])

    # make output file names
    path, file = os.path.split(NucPrepIRSFiles[0])
    if len(limitiesfiles) == 0:
        outfile = os.path.join(path, 'Log2FC_NPvsMP_%dFiles_%dIES_%dCovLim_%.2fIRSLim.LimDKD%s.LimDKDInv%s.DKDdiff%.2f.RetExcLim%d.tsv' % (
            len(MacPrepIRSFiles + MacPrepIRSFiles), len(orderedcommonies), coveragelimit, irslimit, limitiesbydkd, limitiesbydkdinverse, dkddiff, retexclimit))
    elif len(limitiesfiles) > 0:
        outfile = os.path.join(path, 'Log2FC_NPvsMP_%dFiles_%dIES_%dCovLim_%.2fIRSLim.LimDKD%s.LimDKDInv%s.DKDdiff%.2f.SingleKDFilter.RetExcLim%d.tsv' % (
            len(MacPrepIRSFiles + MacPrepIRSFiles + limitiesfiles), len(orderedcommonies), coveragelimit, irslimit, limitiesbydkd, limitiesbydkdinverse, dkddiff, retexclimit))

    # write out data
    with open(outfile, 'w') as OUT:
        OUT.write(header + '\n')
    with open(outfile, 'a') as OUT:
        for k in orderedcommonies:
            outline = '\t'.join([k] + [str(round(d_log2fc[k], 2)) if isinstance(d_log2fc[k], (float, int)) else d_log2fc[k] if isinstance(d_log2fc[k], str) else 'NA' for d_log2fc in alllog2fc + [dlen, dlencategorical]]) + '\n'
            OUT.write(outline)

    print('output file: %s' % outfile)

    # give tab delimited file, 1st column is ies names, all other columns have header and values
    seaborn_pairplot(outfile)

    print('\n####\n####\nCalculate Two-sample Kolmogorov-Smirnov Test (KS Test)')
    # for each KD compare log2FC(irs.np/irs.mp) from KD_PGM to EV_PGM
    # EV_log2fc is list of log2FCs, assumes EV_PGM is first file specified in input
    organize_data_for_two_sample_ks_test(alllog2fc, orderedcommonies, dlencategorical)

    # summarize
    print('\n####\n####\nSummary of output Log2FCs:')
    df = pd.read_csv(outfile, sep="\t")
    print(df.describe())

    tempfiles = []
    for name in NucPrepIRSFiles:
        path, file = os.path.split(name)
        tempfiles.append(os.path.join(path, '_'.join(file.split('_')[:2] + ['Log2FC.KD.Minus.Log2FC.EV'])))

    # scatter plots of comparing results between files,
    xev_fc_log2macnuc = [alllog2fc[0][k] for k in orderedcommonies]
    colorcategories = [dlencategorical[k] for k in orderedcommonies]
    ieslens = [math.log2(dlen[k]) for k in orderedcommonies]
    for dlog2fc, tempfile in zip(alllog2fc[1:], tempfiles[1:]):
        path, macfile = os.path.split(tempfiles[0])
        pathkd, macfilekd = os.path.split(tempfile)
        prefix = '_'.join(macfile.split('_')[:2] + ['Log2FCNPvsMP.EV', 'vs'] + macfilekd.split('_')[:2] + ['Log2FCNPvsMP.KD', 'Normalized'])
        print(prefix)
        outpath = os.path.join(path, prefix + '.pdf')
        outpathevlength = os.path.join(path, 'Log2FCNPvsMP.EV.IESLength.pdf')
        outpathkdlength = os.path.join(path, 'Log2FCNPvsMP.KD.IESLength.pdf')
        ykd_fc_log2macnuc = [dlog2fc[k] for k in orderedcommonies]
        plot_sns_scatter_and_kde_heatmap(xev_fc_log2macnuc, ykd_fc_log2macnuc, 'Log2(FC.NP/FC.MP) EV', 'Log2(FC.NP/FC.MP) KD', prefix, outpath, colorcategories, nbins)
        plot_sns_scatter_and_kde_heatmap(xev_fc_log2macnuc, ieslens, 'Log2(FC.NP/FC.MP) EV', 'log2(IES.Length)', 'IES length vs EV', outpathevlength, colorcategories, nbins)
        plot_sns_scatter_and_kde_heatmap(ykd_fc_log2macnuc, ieslens, 'Log2(FC.NP/FC.MP) KD', 'log2(IES.Length)', 'IES length vs KD', outpathkdlength, colorcategories, nbins)

        # plot_sns_scatter_hue(xev_fc_log2macnuc, ykd_fc_log2macnuc, colorcategories, prefix)
        # plot_kde_heatmap(xev_fc_log2macnuc, ykd_fc_log2macnuc)

    print('\n###\n###\n Comparison between EV_PGM and KD_PGM: log2FC.KD(FC.np/FC.mp) - log2FC.EV(FC.np/FC.mp)')
    print('Histograms represents change in nucleosome content for every IES between EV and KD')
    print('Zero == no change')
    print('Negative: proportionally more nucleosomes bound to excised IES DNA fragments in KD_PGM')
    print('Positive: proportionally more nucleosomes bound to retained IES DNA fragments in KD_PGM')
    for dlog2kdev, tempfile in zip(alllog2fckdev, tempfiles[1:]):
        binnumber = sturges_rule(len(list(dlog2kdev.keys())))
        plot_irs_histogram(dlog2kdev, tempfile, coveragelimit, irslimit, binnumber)
        outpath = tempfile + 'IESLen.pdf'
        log2kdev = [dlog2kdev[k] for k in orderedcommonies]
        plot_sns_scatter_and_kde_heatmap(log2kdev, ieslens, 'Log2(FC.NP/FC.MP) KD - Log2(FC.NP/FC.MP) EV', 'log2(IES.Length)', 'IES length vs KD - EV', outpath, colorcategories, nbins)

    print('For Log2FC residuals')
    for dlog2kdev, tempfile in zip(alllog2fckdev, tempfiles[1:]):
        binnumber = sturges_rule(len(list(dlog2kdev.keys())))
        plot_irs_histogram(allfcresiduals, tempfile, coveragelimit, irslimit, binnumber)
        outpath = tempfile + 'IESLen.pdf'
        fcresiduals = [allfcresiduals[k] for k in orderedcommonies]
        plot_sns_scatter_and_kde_heatmap(fcresiduals, ieslens, 'Log2(FC.NP/FC.MP) KD - Log2(FC.NP/FC.MP) EV', 'log2(IES.Length)', 'IES length vs KD - EV', outpath, colorcategories, nbins)

    # make output file names
    path, file = os.path.split(NucPrepIRSFiles[0])
    if len(limitiesfiles) == 0:
        outfile = os.path.join(path, 'Log2FCNPvsMP.KD.Minus.Log2FCNPvsMP.EV_%dFiles_%dIES_%dCovLim_%.2fIRSLim.LimDKD%s.LimDKDInv%s.DKDdiff%.2f.RetExcLim%d.tsv' % (
            len(MacPrepIRSFiles + MacPrepIRSFiles), len(orderedcommonies), coveragelimit, irslimit, limitiesbydkd, limitiesbydkdinverse, dkddiff, retexclimit))
    elif len(limitiesfiles) > 0:
        outfile = os.path.join(path, 'Log2FCNPvsMP.KD.Minus.Log2FCNPvsMP.EV_%dFiles_%dIES_%dCovLim_%.2fIRSLim.LimDKD%s.LimDKDInv%s.DKDdiff%.2f.SingleKDFilter.RetExcLim%d.tsv' % (
            len(MacPrepIRSFiles + MacPrepIRSFiles + limitiesfiles), len(orderedcommonies), coveragelimit, irslimit, limitiesbydkd, limitiesbydkdinverse, dkddiff, retexclimit))

    nucfilenames = []
    for nucfile in NucPrepIRSFiles:
        path, file = os.path.split(nucfile)
        nucfilenames.append(file)
    header = '\t'.join(['Name'] + ['_'.join(nucfile.split('_')[:2] + ['Minus'] + nucfilenames[0].split('_')[:2]) for nucfile in nucfilenames[1:]] + ['Length', 'Length_Category'])

    # write out data
    with open(outfile, 'w') as OUT:
        OUT.write(header + '\n')
    with open(outfile, 'a') as OUT:
        for k in orderedcommonies:
            outline = '\t'.join([k] + [str(round(d_log2fckdev[k], 2)) if isinstance(d_log2fckdev[k], (float, int)) else d_log2fckdev[k] if isinstance(d_log2fckdev[k], str) else 'NA' for d_log2fckdev in alllog2fckdev + [dlen, dlencategorical]]) + '\n'
            OUT.write(outline)

    print('output file: %s' % outfile)

    # summarize
    print('\n####\n####\nSummary of output Log2FCs:')
    df = pd.read_csv(outfile, sep="\t")
    print(df.describe())



    print('### Finished ###')


def main(MacPrepIRSFiles, NucPrepIRSFiles, iesannotationfile, limitiesfiles=[], limitiesbydkd=False, limitiesbydkdinverse=False, nbins=100, dkddiff=0.2, skdlimit=0.3, retexclimit=1, coveragelimit=1, irslimit=0.01):

    log2fc_mp_vs_np(MacPrepIRSFiles, NucPrepIRSFiles, iesannotationfile, limitiesfiles, limitiesbydkd, limitiesbydkdinverse, nbins, dkddiff, skdlimit, retexclimit, coveragelimit, irslimit)








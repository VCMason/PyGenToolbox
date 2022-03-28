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


def organize_data_for_two_sample_ks_test(allirslog2fc, orderedcommonies, dlencategorical):
    print('\n####\n####\nCalculate Two-sample Kolmogorov-Smirnov Test (KS Test)')
    # for each KD compare log2FC(irs.mp/irs.np) from KD_PGM to EV_PGM
    # EV_log2fc is list of log2FCs, assumes EV_PGM is first file specified in input
    evlog2fc = [allirslog2fc[0][k] for k in orderedcommonies]
    print('KS-test of all IESs that pass filtration')
    for kdirs in allirslog2fc[1:]:
        kdlog2fc = [kdirs[k] for k in orderedcommonies]
        ksresult = ks_test(evlog2fc, kdlog2fc)
        print(ksresult)

    print('KS-test of SHORT IESs that pass filtration')
    shorties = [k for k in orderedcommonies if dlencategorical[k] == 'short']
    evlog2fc = [allirslog2fc[0][k] for k in shorties]
    for kdirs in allirslog2fc[1:]:
        kdlog2fc = [kdirs[k] for k in shorties]
        ksresult = ks_test(evlog2fc, kdlog2fc)
        print(ksresult)

    mediumies = [k for k in orderedcommonies if dlencategorical[k] == 'medium']
    print('KS-test of MEDIUM IESs that pass filtration')
    evlog2fc = [allirslog2fc[0][k] for k in mediumies]
    for kdirs in allirslog2fc[1:]:
        kdlog2fc = [kdirs[k] for k in mediumies]
        ksresult = ks_test(evlog2fc, kdlog2fc)
        print(ksresult)

    longies = [k for k in orderedcommonies if dlencategorical[k] == 'long']
    print('KS-test of LONG IESs that pass filtration')
    evlog2fc = [allirslog2fc[0][k] for k in longies]
    for kdirs in allirslog2fc[1:]:
        kdlog2fc = [kdirs[k] for k in longies]
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

    irs_log2fc = {k: math.log2(dirsmac[k] / dirsnuc[k]) for k in list(dirsmac.keys())}

    return irs_log2fc

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


def show_bias_by_excision(allirsmacfilter, allirsnucfilter, allirslog2fc, MacPrepIRSFiles, orderednames, dcolor):
    # plot correlation between KD_PGM.mp and KD_PGM.np
    # if very strong correlation, then the IRS in np is determined by mp IRS not by MNase digestion
    import os
    import seaborn as sns
    import matplotlib.pyplot as plt
    from scipy.stats.stats import pearsonr

    colorcategories = [dcolor[k] for k in orderednames]

    for dmacirs, dnucirs, dlog2fc, macfilename in zip(allirsmacfilter, allirsnucfilter, allirslog2fc, MacPrepIRSFiles):
        path, macfile = os.path.split(macfilename)
        prefix = '_'.join(macfile.split('_')[:2] + ['MP', 'vs', 'NP', 'NoExcisionNormalization'])
        print(prefix)
        xmac = [dmacirs[k] for k in orderednames]
        ynuc = [dnucirs[k] for k in orderednames]
        sns.scatterplot(xmac, ynuc, hue=colorcategories, legend="brief")
        plt.title(prefix)
        plt.show()
        plt.close()
        print('Pearson Correlation:')
        print(pearsonr(xmac, ynuc))
        plot_kde_heatmap(xmac, ynuc)

        prefix = '_'.join(macfile.split('_')[:2] + ['MP', 'vs', 'Log2FC', 'ExcisionNormalized'])
        print(prefix)
        ylfc = [dlog2fc[k] for k in orderednames]
        sns.scatterplot(xmac, ylfc, hue=colorcategories, legend="brief")
        plt.title(prefix)
        plt.show()
        plt.close()
        print('Pearson Correlation:')
        print(pearsonr(xmac, ylfc))
        plot_kde_heatmap(xmac, ylfc)


def plot_irs_grouped_histogram():

    ## not working yet
    import os
    import matplotlib.pyplot as plt
    import numpy as np

    # make grouped histogram chart
    outpath = os.path.join(path, 'IRS.MultiHist.Alt.pdf')
    plt.figure(figsize=(11.7, 8.27))
    colors = ['red', 'black']
    # colors = ['blue', 'blue', 'red', 'red', 'black']
    labels = ['EV_PGM', 'CAF_PGM']
    t_allIRS2 = np.array(
        list(map(list, zip(*allIRS2))))  # transposing list of lists, lists are IRS values for each file
    plt.hist(t_allIRS2, binnumbermulti, density=True, histtype='bar', color=colors, label=labels)
    plt.legend(prop={'size': 10})
    plt.title('IRSAlt Multi-Histogram')
    plt.xlabel('IRSAlt')
    plt.ylabel('Frequency')
    plt.gca().set_ylim([0, 4])  # ymin,ymax
    plt.savefig(outpath)
    plt.show()
    plt.close()


def plot_irs_histogram(dirs, irsfile, coveragelimit, irslimit, binnumber=50):
    import os
    import matplotlib.pyplot as plt

    path, file = os.path.split(irsfile)
    fileprefix = file.split('.')[0]
    outpath = os.path.join(path, '%s.Alt.hist.CovLim%d.IrsLim%.2f.pdf' % (fileprefix, coveragelimit, irslimit))
    plt.figure(figsize=(11.7,8.27)) # figsize=(3,4))
    plt.hist(list(dirs.values()), bins=binnumber)
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
        print('Finding non KD dependent IESs')

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


def filter_irs(irsfile, coveragelimit=1, irslimit=0.01):
    # coveragelimit has default value of 1 so that denominator if IRS must be > 0
    dirs = {}
    with open(irsfile, 'r') as FILE:
        for line in FILE:
            # if sum of leftmost, rightmost, and both (left&right) molecules supporting retention plus the number of ...
            # molecules support excision are >= limit (and integer value)
            numreadsretained = int(line.strip().split('\t')[2]) + int(line.strip().split('\t')[3]) + int(line.strip().split('\t')[4])
            numreadsexcised = int(line.strip().split('\t')[5])
            irs = float(line.strip().split('\t')[1])
            iesname = line.strip().split('\t')[0]
            if (numreadsretained + numreadsexcised >= coveragelimit) and (irs >= irslimit):
                dirs[iesname] = irs

    print('#####  For file: %s #####' % irsfile)
    print('%d IES elements with counts >= than %d and IRS >=%.2f' % (len(list(dirs.keys())), coveragelimit, irslimit))
    print('%d IES elements with IRS less than 0.1' % len([k for k, v in dirs.items() if v < 0.1]))
    print('%d IES elements with IRS greater than 0.1 and less than 0.25' % len([k for k, v in dirs.items() if (v > 0.1) and (v < 0.25)]))
    print('%d IES elements with IRS above 0.25 and less than 0.75' % len([k for k, v in dirs.items() if (v < 0.75) and (v > 0.25)]))
    print('%d IES elements with IRS greater than 0.75' % len([k for k, v in dirs.items() if v > 0.75]))
    print('%d IES elements with IRS greater than 0.95' % len([k for k, v in dirs.items() if v > 0.95]))

    binnumber = sturges_rule(len(list(dirs.keys())))
    plot_irs_histogram(dirs, irsfile, coveragelimit, irslimit, binnumber)

    return dirs


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


def filter_retention_excision_MP_vs_NP(macfile, nucfile, coveragelimit):
    ''' I CAN'T DO THIS WITHOUT FIRST NORMALIZING COUNTS or DEPTH'''
    ''' calcuates  frequency of fragments kept, and frequency of excised fragments kept from MP to NP '''

    dret, dexc = {}, {}
    with open(macfile, 'r') as FILE:
        for line in FILE:
            # if sum of leftmost, rightmost, and both (left&right) molecules supporting retention plus the number of ...
            # molecules support excision are >= limit (and integer value)
            numreadsretained = int(line.strip().split('\t')[2]) + int(line.strip().split('\t')[3]) + int(line.strip().split('\t')[4])
            numreadsexcised = int(line.strip().split('\t')[5])
            irs = float(line.strip().split('\t')[1])
            iesname = line.strip().split('\t')[0]
            if (numreadsretained + numreadsexcised >= coveragelimit) and (irs >= irslimit):
                dirs[iesname] = irs

    return dret, dexc


def log2fc_mp_vs_np(MacPrepIRSFiles, NucPrepIRSFiles, iesannotationfile, limitiesfiles=[], limitiesbydkd=True, limitiesbydkdinverse=False, dkddiff=0.2, coveragelimit=1, irslimit=0.01):
    # output files written to directory of first nucleosome prep file in NucPrepIRSFiles
    # limitiesbydkdinverse # True or False # requires limitiesbydkd=True # True if you want to test IESs that are not dependent on KD: (IRS.MP in KD_PGM - IRS.MP in EV_PGM <= +dkddiff) and (IRS.MP in KD_PGM - IRS.MP in EV_PGM >= -dkddiff)
    # dkddiff  # float value # if dkddiff is positive, then: IRS.MP in KD_PGM - IRS.MP in EV_PGM >= dkddiff # if dkddiff is negative, then: IRS.MP in KD_PGM - IRS.MP in EV_PGM <= dkddiff
    import os
    import pandas as pd


    dgff3 = read_gff3(iesannotationfile, features=['internal_eliminated_sequence'])
    # dictionary with length of each IES, k is ies name, v is list of a line (for k) split by tabs from ies .gff3 file
    dlen = {k: len(v[8].split(';')[-1][len('sequence='):]) for k, v in dgff3.items()}
    print('Average IES length = %.2f' % (sum(list(dlen.values()))/len(list(dlen.keys()))))
    dlencategorical = {k: 'short' if v < 100 else 'medium' if (v >= 100) and (v <= 200) else 'long' if v > 200 else 'NA' for k, v in dlen.items()}
    #dlencategorical = {k: 0 if v < 100 else 1 if (v >= 100) and (v <= 200) else 2 if v > 200 else 9999 for k, v in dlen.items()}

    allirsmac = []
    allirsnuc = []
    allret, allexc = [], []
    for macfile, nucfile in zip(MacPrepIRSFiles, NucPrepIRSFiles):
        dirsmac = filter_irs(macfile, coveragelimit, irslimit)
        dirsnuc = filter_irs(nucfile, coveragelimit, irslimit)
        dret, dexc = filter_retention_excision_MP_vs_NP(macfile, nucfile, coveragelimit)
        allirsmac.append(dirsmac)
        allirsnuc.append(dirsnuc)
        allret.append(dret)
        allex.append(dexc)

    if limitiesbydkd == True:
        # send list of dictionaries allirsmac
        alldkddepirs = limit_ies_by_dkd(allirsmac, limitiesbydkdinverse, dkddiff)

    # this would use an IRS file from a single KD (like Ptcaf1) to further filter the IESs dependent on a single KD
    if len(limitiesfiles) > 0:
        alllimits = []
        for limitfile in limitiesfiles:
            dlimit = filter_irs(limitfile, coveragelimit, irslimit)
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
        commonies = [k for k, v in dcommonies.items() if v == len(MacPrepIRSFiles + MacPrepIRSFiles)]
    elif (len(limitiesfiles) == 0) and (limitiesbydkd is True):
        for dirs in alldkddepirs + allirsnuc:
            for k, v in dirs.items():
                dcommonies[k] = dcommonies.setdefault(k, 0) + 1
        # list listing names of IES that passed filtration in all input files
        commonies = [k for k, v in dcommonies.items() if v == len(MacPrepIRSFiles + MacPrepIRSFiles)]
    elif (len(limitiesfiles) > 0) and (limitiesbydkd is False):
        for dirs in allirsmac + allirsnuc + alllimits:
            for k, v in dirs.items():
                dcommonies[k] = dcommonies.setdefault(k, 0) + 1
        # list listing names of IES that passed filtration in all input files
        commonies = [k for k, v in dcommonies.items() if v == len(MacPrepIRSFiles + MacPrepIRSFiles + alllimits)]
    elif (len(limitiesfiles) > 0) and (limitiesbydkd is True):
        for dirs in alldkddepirs + allirsnuc + alllimits:
            for k, v in dirs.items():
                dcommonies[k] = dcommonies.setdefault(k, 0) + 1
        # list listing names of IES that passed filtration in all input files
        commonies = [k for k, v in dcommonies.items() if v == len(MacPrepIRSFiles + MacPrepIRSFiles + alllimits)]

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
    # these are list of dictionaries that only contain IRS values for IES that pass filtration in all input files
    for dirsmac, dirsnuc in zip(allirsmac, allirsnuc):
        dirsmacfilter = {k: dirsmac[k] for k in orderedcommonies}
        dirsnucfilter = {k: dirsnuc[k] for k in orderedcommonies}
        allirsmacfilter.append(dirsmacfilter)
        allirsnucfilter.append(dirsnucfilter)

    allirslog2fc = []
    for dirsmacfilter, dirsnucfilter in zip(allirsmacfilter, allirsnucfilter):
        d_irs_log2fc = log2_irs(dirsmacfilter, dirsnucfilter)
        allirslog2fc.append(d_irs_log2fc)

    # not really desired output, but proof that there is some bias by excision that needs to be accounted for
    show_bias_by_excision(allirsmacfilter, allirsnucfilter, allirslog2fc, MacPrepIRSFiles, orderedcommonies, dlencategorical)

    # plot histograms of log2fc(IRS_MP/IRS_NP)
    print('###\n###\nNow plotting log2FC(IRS_MP/IRS_NP)\n###\n###')
    for d_irs_log2fc, irsfile in zip(allirslog2fc, NucPrepIRSFiles):
        binnumber = sturges_rule(len(list(d_irs_log2fc.keys())))
        plot_irs_histogram(d_irs_log2fc, irsfile, coveragelimit, irslimit, binnumber)

    # get file prefixes from nuc prep files
    nucfilenames = ['Name']
    for nucfile in NucPrepIRSFiles:
        path, file = os.path.split(nucfile)
        nucfilenames.append(file)
    header = '\t'.join(['_'.join(nucfile.split('_')[:2]) for nucfile in nucfilenames] + ['Length', 'Length_Category'])

    # make output file names
    path, file = os.path.split(NucPrepIRSFiles[0])
    if len(limitiesfiles) == 0:
        outfile = os.path.join(path, 'Log2FC_IRS_MPvsNP_%dFiles_%dIES_%dCovLim_%.2fIRSLim.tsv' % (
            len(MacPrepIRSFiles + MacPrepIRSFiles), len(orderedcommonies), coveragelimit, irslimit))
    elif len(limitiesfiles) > 0:
        outfile = os.path.join(path, 'Log2FC_IRS_MPvsNP_%dFiles_%dIES_%dCovLim_%.2fIRSLim.tsv' % (
        len(MacPrepIRSFiles + MacPrepIRSFiles + limitiesfiles), len(orderedcommonies), coveragelimit, irslimit))

    # write out data
    with open(outfile, 'w') as OUT:
        OUT.write(header + '\n')
    with open(outfile, 'a') as OUT:
        for k in orderedcommonies:
            outline = '\t'.join([k] + [str(round(d_irs_log2fc[k], 2)) if isinstance(d_irs_log2fc[k], (float, int)) else d_irs_log2fc[k] if isinstance(d_irs_log2fc[k], str) else 'NA' for d_irs_log2fc in allirslog2fc + [dlen, dlencategorical]]) + '\n'
            OUT.write(outline)

    print('output file: %s' % outfile)

    # give tab delimited file, 1st column is ies names, all other columns have header and values
    seaborn_pairplot(outfile)

    print('\n####\n####\nCalculate Two-sample Kolmogorov-Smirnov Test (KS Test)')
    # for each KD compare log2FC(irs.mp/irs.np) from KD_PGM to EV_PGM
    # EV_log2fc is list of log2FCs, assumes EV_PGM is first file specified in input
    organize_data_for_two_sample_ks_test(allirslog2fc, orderedcommonies, dlencategorical)

    # summarize
    print('\n####\n####\nSummary of output Log2FCs:')
    df = pd.read_csv(outfile, sep="\t")
    print(df.describe())

    print('### Finished ###')


def main(MacPrepIRSFiles, NucPrepIRSFiles, iesannotationfile, limitiesfiles=[], limitiesbydkd=True, limitiesbydkdinverse=False, dkddiff=0.2, coveragelimit=1, irslimit=0.01):

    log2fc_mp_vs_np(MacPrepIRSFiles, NucPrepIRSFiles, iesannotationfile, limitiesfiles, limitiesbydkd, limitiesbydkdinverse, dkddiff, coveragelimit, irslimit)








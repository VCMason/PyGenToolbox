
def barplot(clusters, strengths, outpath, title='', xlab='Ptiwi08 Bound sRNA Cluster', ylab='PingPong Signal Strength', figlength=10, figheight=5):
    # strengths is list of integers representing height of bar plot
    # clusters is list of cluster names (strings), representing bar plot categories
    import matplotlib.pyplot as plt
    import os

    p, f = os.path.split(outpath)
    title = f[:-len('.F.sort.53.PingPong.FvsR.55.DicerSignal.SortPosition.pdf')]
    #general figure setup
    plt.rcParams.update({'font.size': 12})
    plt.rcParams.update({'figure.autolayout': True})  # automatically adjust fig so x axis labels are visible
    plt.figure(figsize=(figlength, figheight))
    plt.xticks(rotation=45)  # rotation='vertical' # , fontsize=12 # rotation=45, ha='right'
    # plot data and axis labels
    plt.bar(clusters, strengths)
    plt.grid(True)
    plt.title(title)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.savefig(outpath)
    plt.show()
    plt.close()
    return


def summarize_ping_ping_signal_per_cluster(bed, dpp):
    dclusterstrength = {}
    names = []
    for i in bed:
        scaffold, start, end = i[0], int(i[1]), int(i[2])
        key = i[3]  # cluster id
        names.append(key)
        dclusterstrength[key] = 0
        for position, strength in dpp.items():
            ppscaff = '_'.join(position.split('_')[:-1])
            ppcoord = int(position.split('_')[-1])
            if (ppscaff == scaffold) and (ppcoord >= start) and (ppcoord <= end):
                dclusterstrength[key] = dclusterstrength.setdefault(key, 0) + int(strength)

    return dclusterstrength, names


def collect_ping_pong_strength(pingpingfile):
    # scaffold51_8_23351	18
    with open(pingpingfile, 'r') as FILE:
        dpingpong = {line.strip().split('\t')[0]: int(line.strip().split('\t')[1]) for line in FILE if line.strip() != ''}

    return dpingpong


def read_bed(bedfile):
    # scaffold51_8	23301	23738   cluster8
    with open(bedfile, 'r') as FILE:
        # scaffold = line.strip().split('\t')[0]
        # start = int(line.strip().split('\t')[1])
        # end = int(line.strip().split('\t')[2])
        # id = line.strip().split('\t')[3]
        bed = [line.strip().split('\t') for line in FILE if line.strip() != '']

    return bed


def main(pingpingfiles, bedfile):
    import os

    # pingpingfiles is list of files.
    # bedfile has at least four tab delimited columns. Scaffold start   end id
    bed = read_bed(bedfile)
    for filename in pingpingfiles:
        dpingpong = collect_ping_pong_strength(filename)
        dclusterstrength, names = summarize_ping_ping_signal_per_cluster(bed, dpingpong)

        path, f = os.path.split(filename)
        outpath = os.path.join(path, '.'.join(filename.split('.')[:-1] + ['pdf']))
        strengths = [dclusterstrength[n] for n in names]
        barplot(names, strengths, outpath)
def plot_scatter(x, y, outpath):
    ''' x is list of x values, y is list of y values, path is full path to output file '''
    import matplotlib.pyplot as plt
    # import os
    plt.plot(x, y)
    plt.grid(True)
    # plt.title(os.path.split(outpath)[1])
    plt.xlabel('All Scaffold Positions')
    plt.ylabel('Depth')
    plt.savefig(outpath)
    plt.close()

def format_data_separate(d):
    perscaffd = {}
    for k, valuelist in d.items():
        x_values = []
        y_values = []
        # iterate for each value, each value is a list of lists
        #  v[1] == depth, v[0] relative scaffold position
        for v in valuelist:
            # iterate for each list within each within value
            # record all x and y coordinates in large lists
            x_values.append(v[0])
            y_values.append(v[1])
        perscaffd[k] = perscaffd.get(k, []) + [x_values, y_values]

    return perscaffd


def format_data_together(d):
    x_values = []
    y_values = []
    for valuelist in d.values():
        # iterate for each value, each value is a list of lists
        #  v[1] == depth, v[0] relative scaffold position
        for v in valuelist:
            # iterate for each list within each within the list of lists value
            # record all x and y coordinates in large lists
            x_values.append(v[0])
            y_values.append(v[1])
    return x_values, y_values


def parse_bed_coord_depth_to_dict(bedfile, cutoff):
    ''' key is unique first column entries. values are position and depth int values in a list of list of lists '''
    d = {}
    with open(bedfile) as BED:
        for line in BED:
            if int(line.strip().split('\t')[2]) >= cutoff:
                d[line.strip().split('\t')[0]] = d.get(line.strip().split('\t')[0], []) + [[int(line.strip().split('\t')[1]), int(line.strip().split('\t')[2])]]
    return d


def main(bedfile, outpath, cutoff=5, separate=False):
    ''' bedfile is full path to input bed file, cutoff is int value >= desired depth, path is full path to output file '''
    print('Reading in Bedfile')
    d = parse_bed_coord_depth_to_dict(bedfile, cutoff)

    if separate is False:
        print('Formatting Data Together')
        x_values, y_values = format_data_together(d)
        print('Plotting data')
        plot_scatter(x_values, y_values, outpath)

    if separate is True:
        print('Formatting Data Separate and Plotting')
        perscaffd = format_data_separate(d)
        for k in perscaffd.keys():
            outfilename = '.'.join(outpath.split('.')[:-1] + [k] + [str(cutoff)] + [outpath.split('.')[-1]])
            plot_scatter(perscaffd[k][0], perscaffd[k][1], outfilename)
    print('Finished')
    # xaxislength = max(map(len, d.values())) # length of longest dict value

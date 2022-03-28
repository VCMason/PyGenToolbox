def mean(lst):
    return sum(lst) / len(lst)


def plot_histogram(x, lowerlimit=1000):
    import matplotlib.pyplot as plt
    print('Mean of fragment size: %f' % (mean(x)))
    # x = [value1, value2, value3, ....]
    print('Largest fragment length (in bps): %d' % max(x))
    x = sorted(x)
    plt.hist(x, bins=100)
    plt.show()
    plt.close()

    x2 = [i for count, i in enumerate(x) if (count >= 0.01*len(x)) and (count <= 0.99*len(x))]
    print('Mean of fragment sizes, excluding lower 1 percent and upper 1 percent: %f' % (mean(x2)))
    plt.hist(x2, bins=100)
    plt.show()
    plt.close()

    x3 = [i for count, i in enumerate(x) if (count >= 0.05*len(x)) and (count <= 0.95*len(x))]
    print('Mean of fragment sizes, excluding lower 5 percent and upper 5 percent: %f' % (mean(x3)))
    plt.hist(x3, bins=100)
    plt.show()
    plt.close()

    x4 = [i for i in x if i >= lowerlimit]
    print('Mean of mean fragment size > %d: %f' % (lowerlimit, mean(x4)))
    plt.hist(x4, bins=100)
    plt.show()
    plt.close()

    x5 = [i for i in x if i >= lowerlimit * 2]
    print('Mean of mean fragment size > %d: %f' % (lowerlimit * 2, mean(x5)))
    plt.hist(x5, bins=100)
    plt.show()
    plt.close()


def calculate_distribution_of_length_of_genomic_seqs_after_cutting_by_feature(names, efficiency=0.9):
    # this function might fuck up if the features are overlapping (calculate a negative fragment size)
    # this would happen if the start of the "next" feature is before the end of the previous one
    distances = []
    gate = 0
    import random
    print('Cuttting efficiency = %f' % efficiency)
    prevscaffold = 'INITIATE'
    for n in names:
        scaffold = n.split('.')[0]
        # start = int(n.split('.')[1])
        # end = int(n.split('.')[2])
        if scaffold != prevscaffold:
            gate = 0
        if gate == 0:
            current_pos = int(n.split('.')[1])
            distances.append(current_pos - 0)
            previous_pos = int(n.split('.')[2])
            gate = 1
        elif (gate == 1) and (random.randint(1, 10) <= 10 * efficiency):
            current_pos = int(n.split('.')[1])
            distances.append(current_pos - previous_pos)
            previous_pos = int(n.split('.')[2])
        prevscaffold = n.split('.')[0]
    print(distances[:50])
    return distances


def ParseGff3SkipComments(f):
    print(f)
    dgff3 = {}
    names = []
    with open(f, 'r') as FILE:
        for l in FILE:
            if l[0] == '#':
                pass
            else:
                n = '.'.join([l.strip().split('\t')[0]] + l.strip().split('\t')[3:5])  # n == scaffold_start_end
                names.append(n)
                dgff3[n] = l.strip().split('\t')
    return dgff3, names


def main_MakeDistributionOfCutGenomeSeqsByFeature(f, feature):
    # import pkg_resources
    # DATA_PATH = pkg_resources.resource_filename('pygentoolbox', 'data/seqs/ptetraurelia_mac_51.fa')

    #f = input('Enter full path to file: ')  # f should be full path to annotation file .gff3
    #feature = input('Enter feature type to cut genome sequence, ex. internal_eliminated_sequence: ')
    # f = 'D:\\LinuxShare\\Ciliates\\Genomes\\Seqs\\ptetraurelia_mac_51.fa'
    # f = 'D:\\LinuxShare\\Ciliates\\Genomes\\Seqs\\ptetraurelia_mac_51_with_ies.fa'
    # f = 'D:\\LinuxShare\\Ciliates\\Genomes\\Seqs\\ptetraurelia_mic2.fa'
    # f = 'D:\\LinuxShare\\Ciliates\\Genomes\\Seqs\\contigs_ABK_COSP_best_k51_no_scaf.fa'

    dgff3, names = ParseGff3SkipComments(f)
    print('Total number of features: %d' % (len(dgff3.keys())))
    dredgff3 = {k: v for k, v in dgff3.items() if v[2] == feature}  # isolate only relevant features
    rednames = [n for n in names if n in dredgff3.keys()]
    print('Number of features matching %s: %d' % (feature, len(dredgff3.keys())))
    distances = calculate_distribution_of_length_of_genomic_seqs_after_cutting_by_feature(rednames, efficiency=0.9)
    plot_histogram(distances)

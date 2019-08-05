
def compare_ies_to_features_up_and_down(f, depthcutoff):
    ''' should run clone_features_up_down() in pygentoolbox.Tools before this code'''
    ''' f is fullpath to file, depthcutoff is int variable. Keep depth values >= depthcutoff'''
    import math
    output = []
    with open(f, 'r') as FILE:
        lines = FILE.readlines()
        for i in range(len(lines)):
            if lines[i][0] == '#':
                pass
            else:
                l = lines[i].strip().split('\t')
                if l[2] == 'internal_eliminated_sequence':
                    up = int(lines[i-1].strip().split('\t')[9])
                    ies = int(l[9])
                    down = int(lines[i+1].strip().split('\t')[9])
                    if (up >= depthcutoff) and (down >= depthcutoff):  # and (ies > 0)
                        diffup, diffdown = up-ies, down-ies
                        # log2up, log2down = math.log2(up/ies), math.log2(down/ies)
                        x = [str(j) for j in [up, ies, down, diffup, diffdown]]  # , log2up, log2down
                        output.append(l[:5] + x + l[5:-1])

    fout =  '.'.join(f.split('.')[:-1] + ['diff'] + [f.split('.')[-1]])
    with open(fout, 'w') as OUT:
        OUT.write('\n'.join(['\t'.join(y) for y in output]))
    print('Added fields: up, ies, down, diffup, diffdown')  # , log2up, log2down'
    # print('Requiring ies to have > 0 depth')
    print('Output file: %s' % (fout))
    print('Number of ies with up and down features with >= %d coverage: %d' % (depthcutoff, len(output)))

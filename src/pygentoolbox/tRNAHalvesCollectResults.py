''' This program is retarded '''

def collect_files_from_dir_with_extension(directory, fileextension):
    import os

    # fullpath = os.path.join(directory, '*%s' % fileextension)
    # glob.glob(fullpath)
    filelist = []
    for file in os.listdir(directory):
        if file.endswith(fileextension):
            filelist.append(os.path.join(directory, file))

    return filelist


def main(tRNAAnnotationfile, gff3file, directorylist, countcutoff=100):
    dtrnaannot = {}
    with open(tRNAAnnotationfile, 'r') as FILE:
        header = FILE.readline()
        # '.'.join([l.split('\t').strip()[0]] + l.split('\t').strip()[2:4]) == scaffold51_9.start.end
        # note, start could be bigger then end if '-' orientation
        # '_'.join(l.split('\t').strip()[4:6] == Trp_CCA for example
        #dtrnaannot = {'.'.join([l.split('\t').strip()[0]] + l.split('\t').strip()[2:4]): l for l in FILE}
        for line in FILE:
            s = line.strip().split('\t')
            sclean = [elem.strip() for elem in s]  # some columns have spaces after values....
            if int(sclean[2]) > int(sclean[3]):  # '-' orientation
                dtrnaannot['.'.join([sclean[0], sclean[3], sclean[2]])] = '_'.join(sclean[4:6])
            else:  # it is forward '+' orientation
                dtrnaannot['.'.join([sclean[0]] + sclean[2:4])] = '_'.join(sclean[4:6])

    dgff3annot = {}
    dgff3names = []
    dtrnafam = {}
    newannotfile = []
    with open(gff3file, 'r') as FILE:
        # line.strip().split('\t')[8].strip().split(';')[0] == 'ID=PTET... gene id ex: ID=PTET.51.1.G0390204
        for line in FILE:
            trnakey = '.'.join([line.strip().split('\t')[0]] + line.strip().split('\t')[3:5])
            if line[0] != '#':
                if line.strip().split('\t')[2] == 'ncRNA':
                    try:
                        dtrnaannot[trnakey]  # line.strip().split('\t')[8].strip().split(';')[0]
                    except:
                        newannotfile.append(line.strip())
                    else:
                        dgff3names.append(line.strip().split('\t')[8].split(';')[0])
                        dgff3annot[line.strip().split('\t')[8].split(';')[0]] = line.strip() + ';' + 'tRNAFamily=%s' % dtrnaannot[trnakey]
                        newannotfile.append(line.strip() + ';' + 'tRNAFamily=%s' % dtrnaannot[trnakey])
                        dtrnafam[line.strip().split('\t')[8].split(';')[0]] = dtrnaannot[trnakey]

    with open('.'.join(gff3file.split('.')[:-1] + ['tRNAFam', 'gff3']), 'w') as OUT:
        OUT.write('\n'.join(newannotfile))

    for directory in directorylist:
        flagfiles = collect_files_from_dir_with_extension(directory, '.sort.sam.flagstat')
        featurefiles = collect_files_from_dir_with_extension(directory, '.sort.allfeatures.sorted.fixed.tsv')

        for file in featurefiles:
            print(file)
            outfeatures = []
            # Keep only gene feature IDs, per gene feature add:
            # "tRNA family" as AA_anticodon, scaffold, start, end, orientation
            with open(file, 'r') as FILE:
                header = FILE.readline().strip().split('\t')
                header = '\t'.join(header + ['Scaffold', 'start', 'end', 'orientation', 'tRNAFamily'])
                for feature in FILE:
                    featsplit = feature.strip().split('\t')
                    # dgff3annot[featsplit[0]] == ID=PTET....
                    if int(featsplit[1]) > countcutoff:
                        if featsplit[0] in dgff3names:
                            outfeatures.append('\t'.join(featsplit + [dgff3annot[featsplit[0]].split('\t')[0], dgff3annot[featsplit[0]].split('\t')[3], dgff3annot[featsplit[0]].split('\t')[4], dgff3annot[featsplit[0]].split('\t')[6], dtrnafam[featsplit[0]]]))
                        else:
                            outfeatures.append('\t'.join(featsplit + [dgff3annot[featsplit[0]].split('\t')[0], dgff3annot[featsplit[0]].split('\t')[3], dgff3annot[featsplit[0]].split('\t')[4], dgff3annot[featsplit[0]].split('\t')[6], dtrnafam[featsplit[0]]]))
            print('writing to file: %s' % file + '.out')
            with open(file + '.out', 'w') as OUT:
                OUT.write('\n'.join([header] + outfeatures))






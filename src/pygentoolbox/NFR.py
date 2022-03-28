# run on linux
# need Nucleosome Dynamics installed
# specifically readBAM.R, nucleR.R, and NFR.R (+ dependencies...)


def nfr(nucleroutgff, nfroutgff, minwidth):
    import subprocess
    # cmd = "Rscript /media/sf_LinuxShare/Programs/nucleosome_dynamics/bin/NFR.R  --input $PREFIX.RData.gff --output $PREFIX.NFR.RData.gff"
    cmd = "Rscript /media/sf_LinuxShare/Programs/nucleosome_dynamics/bin/NFR.R  --input %s --output %s --minwidth %d" % (nucleroutgff, nfroutgff, minwidth)
    cmdlist = cmd.split()
    subprocess.call(cmdlist)


def read_scaffolds(scaffoldfile):
    with open(scaffoldfile, 'r') as FILE:
        scaffolds = [line.strip() for line in FILE]
    return scaffolds


def main(bamfile, scaffoldfile, minwidth=110):
    # bamfile is full path to file, bamfile was aligned with paired end reads
    # scaffoldfile is a file with one name of one scaffold per line

    # scaffolds is list of all scaffold names that we want to analyze in bamfile
    scaffolds = read_scaffolds(scaffoldfile)

    for scaff in scaffolds:
        print(scaff)
        nucleroutgff = '.'.join(bamfile.split('.')[:-1] + [scaff, 'nucleR', 'gff'])
        nfroutgff = '.'.join(bamfile.split('.')[:-1] + [scaff, 'NFR', str(minwidth), 'gff'])

        nfr(nucleroutgff, nfroutgff, minwidth)
    print('### FIN ###')
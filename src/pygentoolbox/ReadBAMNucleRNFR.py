# run on linux
# need Nucleosome Dynamics installed
# specifically readBAM.R, nucleR.R, and NFR.R (+ dependencies...)


def nfr(nucleroutgff, nfroutgff, minwidth):
    import subprocess
    # cmd = "Rscript /media/sf_LinuxShare/Programs/nucleosome_dynamics/bin/NFR.R  --input $PREFIX.RData.gff --output $PREFIX.NFR.RData.gff"
    cmd = "Rscript /media/sf_LinuxShare/Programs/nucleosome_dynamics/bin/NFR.R  --input %s --output %s --minwidth %d" % (nucleroutgff, nfroutgff, minwidth)
    cmdlist = cmd.split()
    subprocess.call(cmdlist)


def nucler(rdata, nucleroutgff, fraglen, scaff):
    import subprocess
     # cmd = "Rscript /media/sf_LinuxShare/Programs/nucleosome_dynamics/bin/nucleR.R --input $PREFIX.RData --output $PREFIX.RData.gff --type paired --fragmentLen 200 --chr scaffold51_1_with_IES"
    cmd = "Rscript /media/sf_LinuxShare/Programs/nucleosome_dynamics/bin/nucleR.R --input %s --output %s --type paired --fragmentLen %d --chr %s" % (rdata, nucleroutgff, fraglen, scaff)
    cmdlist = cmd.split()
    subprocess.call(cmdlist)


def readbam_r(bamfile, rdata):
    import subprocess
    # cmd = "Rscript /media/sf_LinuxShare/Programs/nucleosome_dynamics/bin/readBAM.R --input $PREFIX.bam --output $PREFIX.RData --type paired"
    cmd = "Rscript /media/sf_LinuxShare/Programs/nucleosome_dynamics/bin/readBAM.R --input %s --output %s --type paired" % (bamfile, rdata)
    cmdlist = cmd.split()
    subprocess.call(cmdlist)


def read_scaffolds(scaffoldfile):
    with open(scaffoldfile, 'r') as FILE:
        scaffolds = [line.strip() for line in FILE]
    return scaffolds


def main(bamfile, scaffoldfile, fraglen=200, minwidth=110):
    # bamfile is full path to file, bamfile was aligned with paired end reads
    # scaffoldfile is a file with one name of one scaffold per line

    # only need to convert bam file to Rdata file once
    print(bamfile)
    rdata = '.'.join(bamfile.split('.')[:-1] + ['RData'])
    readbam_r(bamfile, rdata)

    # scaffolds is list of all scaffold names that we want to analyze in bamfile
    scaffolds = read_scaffolds(scaffoldfile)

    for scaff in scaffolds:
        print(scaff, end=',')
        nucleroutgff = '.'.join(bamfile.split('.')[:-1] + [scaff, 'nucleR', 'gff'])
        nfroutgff = '.'.join(bamfile.split('.')[:-1] + [scaff, 'NFR', str(minwidth), 'gff'])

        nucler(rdata, nucleroutgff, fraglen, scaff)
        nfr(nucleroutgff, nfroutgff, minwidth)
    print('### FIN ###')
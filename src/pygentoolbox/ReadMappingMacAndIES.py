def run_samtools(samfile, cleanup=False):
    # samfile is full path to sam file
    print('Starting samtools converting sam to sorted indexed bam')
    import os
    import subprocess

    path, sam = os.path.split(samfile)
    bamfile = '.'.join(samfile.split('.')[:-1] + ['bam'])
    sortbamfile = '.'.join(samfile.split('.')[:-1] + ['sort', 'bam'])
    sortsamfile = '.'.join(samfile.split('.')[:-1] + ['sort', 'sam'])
    flagstatfile = '.'.join(samfile.split('.')[:-1] + ['sort', 'sam', 'flagstat'])

    with open(bamfile, 'w') as OUT:
        # cmd = 'samtools view -h -b %s > %s' % (samfile, bamfile)
        cmd = 'samtools view -h -b %s' % samfile
        ps = subprocess.Popen(cmd.split(), stdout=OUT)
        ps.wait()
    # cmd = 'rm %s' % samfile
    with open(sortbamfile, 'w') as OUT:
        # cmd = 'samtools sort %s > %s' % (bamfile, sortbamfile)
        cmd = 'samtools sort %s' % (bamfile)
        ps = subprocess.Popen(cmd.split(), stdout=OUT)
        ps.wait()
    if cleanup == True:
        os.remove(samfile)  # delete sam file
        os.remove(bamfile)  # delete bam file
    cmd = 'samtools index %s' % sortbamfile
    ps = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
    ps.wait()
    # with open(sortsamfile, 'w') as OUT:
    #     # cmd = 'samtools view -h %s > %s' % (sortbamfile, sortsamfile)
    #     cmd = 'samtools view -h %s' % sortbamfile
    #     ps = subprocess.Popen(cmd.split(), stdout=OUT)
    #     ps.wait()
    with open(flagstatfile, 'w') as OUT:
        # cmd = 'samtools flagstat %s > %s' % (sortbamfile, flagstatfile)
        cmd = 'samtools flagstat %s' % sortbamfile
        ps = subprocess.Popen(cmd.split(), stdout=OUT)
        ps.wait()

    print('Finished with Samtools, returning %s\n' % sortbamfile)
    return sortbamfile


def run_hista2_PE(aligndatabase, forwardfile, reversefile, threads, knownsplicesites, paired=True):
    import os
    import subprocess
    from pygentoolbox.Tools import make_directory

    print('Aligning reads to Reference:\n%s' % aligndatabase)
    refpath, referenceprefix = os.path.split(aligndatabase)
    path, f = os.path.split(forwardfile)
    # pathminusonedir, dir = os.path.split(path)
    make_directory(os.path.join(path, 'hisat2'))

    # assuming read files are .fastq.gz
    if paired == True:
        # then data is paired end
        print('Starting Hisat2: aligning\n%s, %s' % (forwardfile, reversefile))
        outsamfile = os.path.join(path, 'hisat2', '.'.join(f.split('.')[:-2] + [referenceprefix, 'sam']))
        cmd = 'hisat2 -q -p %d --no-softclip --pen-noncansplice 0 --known-splicesite-infile %s -x %s -1 %s -2 %s -S %s' % (threads, knownsplicesites, aligndatabase, forwardfile, reversefile, outsamfile)
    else:
        # then data is single end
        print('Starting Hisat2: aligning\n%s' % forwardfile)
        outsamfile = os.path.join(path, 'hisat2', '.'.join(f.split('.')[:-2] + [referenceprefix, 'sam']))
        cmd = 'hisat2 -q -p %d --no-softclip --pen-noncansplice 0 --known-splicesite-infile %s -x %s -U %s -S %s' % (threads, knownsplicesites, aligndatabase, forwardfile, outsamfile)

    subprocess.call(cmd.split())

    print('Finished with Hisat2\n')
    return outsamfile


def run_fastqc(fullpath):
    # fullpath is full path to input file
    # '/media/sf_LinuxShare/Programs/FastQC/fastqc -o /media/sf_LinuxShare/Projects/Lyna/DATA/fastqc -f fastq fastq 200107_NB501850_A_L1-4_ADPF-98_R1.fastq'
    print('Starting fastqc')
    import os
    import subprocess
    from pygentoolbox.Tools import make_directory

    path, f = os.path.split(fullpath)
    # pathminusonedir, dir = os.path.split(path)
    outpath = os.path.join(path, 'fastqc')
    print(outpath)
    make_directory(outpath)
    # /media/sf_LinuxShare/Programs/FastQC/fastqc is path to executable
    cmd = '/media/sf_LinuxShare/Programs/FastQC/fastqc -o %s -f fastq fastq %s' % (outpath, fullpath)
    print(cmd)
    subprocess.call(cmd.split())
    outputfile = os.path.join(outpath, f)
    print('Finished fastqc, output at directory:\n%s\n' % outpath)

    return


def run_fastp_single_end(forwardreadfile):
    import subprocess
    import os
    from pygentoolbox.Tools import make_directory

    ### start fastp ###
    path1, f1 = os.path.split(forwardreadfile)

    print('Start trimming of:\n%s\n' % forwardreadfile)
    make_directory(os.path.join(path1, 'fastp'))
    # example: cd /media/sf_LinuxShare/Projects/Lyna/DATA
    #cmd = "cd %s" % path1
    #print(cmd)
    #cmdlist = cmd.split()
    #p = subprocess.call(cmdlist)

    ## fastp -i 500_LK_L1_R1.fastq.gz -I 500_LK_L1_R2.fastq.gz -o 500_LK_R1.trim.fastq.gz -O 500_LK_R2.trim.fastq.gz
    ## fastp -i /media/sf_LinuxShare/Projects/Lyna/TestMyPipe/500_LK_L1_R1.fastq.gz -I /media/sf_LinuxShare/Projects/Lyna/TestMyPipe/500_LK_L1_R2.fastq.gz -o /media/sf_LinuxShare/Projects/Lyna/TestMyPipe/fastp/500_LK_L1_R1.fastq.gz -O /media/sf_LinuxShare/Projects/Lyna/TestMyPipe/fastp/500_LK_L1_R2.fastq.gz
    forwardout = os.path.join(path1, "fastp", '.'.join(f1.split('.')[:-2] + ['trim', 'fastq', 'gz']))  # output file for
    cmd = "fastp -i %s -o %s" % (forwardreadfile, forwardout)
    print(cmd)
    cmdlist = cmd.split()
    p = subprocess.Popen(cmdlist, stdout=subprocess.PIPE)
    cmdout, err = p.communicate()
    print(cmdout)
    print('Finished trimming, files output to:\n%s\n' % forwardout)
    ### end fastp ###

    return forwardout


def run_fastp_paired_end(forwardreadfile, reversereadfile):
    import subprocess
    import os
    from pygentoolbox.Tools import make_directory

    ### start fastp ###
    print('Start trimming of:\n%s\n%s' % (forwardreadfile, reversereadfile))

    path1, f1 = os.path.split(forwardreadfile)
    path2, f2 = os.path.split(reversereadfile)

    make_directory(os.path.join(path1, 'fastp'))
    make_directory(os.path.join(path2, 'fastp'))
    # example: cd /media/sf_LinuxShare/Projects/Lyna/DATA
    #cmd = "cd %s" % path1
    #print(cmd)
    #cmdlist = cmd.split()
    #p = subprocess.call(cmdlist)

    ## fastp -i 500_LK_L1_R1.fastq.gz -I 500_LK_L1_R2.fastq.gz -o 500_LK_R1.trim.fastq.gz -O 500_LK_R2.trim.fastq.gz
    ## fastp -i /media/sf_LinuxShare/Projects/Lyna/TestMyPipe/500_LK_L1_R1.fastq.gz -I /media/sf_LinuxShare/Projects/Lyna/TestMyPipe/500_LK_L1_R2.fastq.gz -o /media/sf_LinuxShare/Projects/Lyna/TestMyPipe/fastp/500_LK_L1_R1.fastq.gz -O /media/sf_LinuxShare/Projects/Lyna/TestMyPipe/fastp/500_LK_L1_R2.fastq.gz
    forwardout = os.path.join(path1, "fastp", '.'.join(f1.split('.')[:-2] + ['trim', 'fastq', 'gz']))  # output file for
    reverseout = os.path.join(path2, "fastp", '.'.join(f2.split('.')[:-2] + ['trim', 'fastq', 'gz']))  # output file rev
    cmd = "fastp -i %s -I %s -o %s -O %s" % (forwardreadfile, reversereadfile, forwardout, reverseout)
    print(cmd)
    cmdlist = cmd.split()
    p = subprocess.Popen(cmdlist, stdout=subprocess.PIPE)
    cmdout, err = p.communicate()
    print(cmdout)
    print('Finished trimming, files output to:\n%s\n%s\n' % (forwardout, reverseout))
    ### end fastp ###

    return forwardout, reverseout


def main(forward, reverse, directory='', threads=1, skipsteps=[], paired=True,
         reference='/media/sf_LinuxShare/Ciliates/Genomes/Hisat2_Indexes/Pt_51_Mac',
         knownsplicesites='/media/sf_LinuxShare/Ciliates/Genomes/Annotations/ss.IES.pt_51_MacAndIES.table'):
    ### input .fastq.gz files
    ### forwardreadfile='/media/sf_LinuxShare/Projects/Lyna/TestMyPipe/500_LK_L1_R1.fastq.gz',
    ### reversereadfile='/media/sf_LinuxShare/Projects/Lyna/TestMyPipe/500_LK_L1_R2.fastq.gz',

    # run script from directory with flypipe scripts and data
    # assumes data is paired end
    # forwardreadfile == full path to forward read file .fastq.gz
    # reversereadfile == full path to reverse read file .fastq.gz
    import os
    import glob
    import pandas as pd
    from pygentoolbox.Tools import run_gzip_decompress

    if directory == '':  # os if they do not specify an input directory... then use cwd
        directory = os.getcwd()
        forwardlist = [os.path.join(directory, f) for f in forward]
        reverselist = [os.path.join(directory, f) for f in reverse]
    else:
        forwardlist = [os.path.join(directory, f) for f in forward]
        reverselist = [os.path.join(directory, f) for f in reverse]

    d = {}  # will have referencename as key and list as value with counts
    dnames = []
    data = [[], [], []]  # data[0] == list of seq lengths, data[1] == list of reference genomes (molecule types),
                         # data[2] == list of number of reads that map to each reference for that seq length
    for count, forwardfile in enumerate(forwardlist):
        reversefile = reverselist[count]
        if 'fastp' not in skipsteps:
            ### start fastp ###
            forwardclean, reverseclean = run_fastp_paired_end(forwardfile, reversefile)
            #forwardclean = run_fastp_single_end(file)
        else:
            print('yep your reads must already be clean')
            forwardclean = forwardfile
            reverseclean = reversefile
        if 'gzip' not in skipsteps:
            ### start gzip -dk file.trim.fastq.gz
            forwardtrim = run_gzip_decompress(forwardclean)
            reversetrim = run_gzip_decompress(reverseclean)
        else:
            forwardtrim = forwardclean
            reversetrim = reverseclean

        if 'fastqc' not in skipsteps:
            ### start fastqc ###
            run_fastqc(forwardtrim)
            run_fastqc(reversetrim)

            ### delete read files, still have .tar.gz files ##
            os.remove(forwardtrim)
            os.remove(reversetrim)

        if 'hisat2' not in skipsteps:
            ### align with hisat2 ###
            #samfile = run_hisat2(reference, forwardtrim)
            samfile = run_hista2_PE(reference, forwardclean, reverseclean, threads, knownsplicesites, paired=True)

        if 'samtools' not in skipsteps:
            ### run samtools convert to bam, sort, index ###
            bamsortfile = run_samtools(samfile, cleanup=True)

    print('##########\n###FIN####\n##########')
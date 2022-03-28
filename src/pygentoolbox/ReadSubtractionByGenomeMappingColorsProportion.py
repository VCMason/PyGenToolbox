# order of read subtraction
# kleb pneu -> Ptet mtDNA -> Mac -> Mac+IES -> Mic2.0 or OES
# input files .fastq.gz

# samtools view -b -f 4 file.bam > unmapped.bam
# samtools view -b -F 4 file.bam > unmapped.bam
# number of reads in bam file
# samtools view -c SAMPLE.bam
# counting only mapped (primary aligned) reads
# samtools view -c -F 260 SAMPLE.bam
# convert bam to fastq
# samtools bam2fq SAMPLE.bam > SAMPLE.fastq


def make_pandas_stacked_barplot(pivot_df, order, outpath, colors):
    ### pivot_df is product of:
    # rows = zip(data[0], data[1], data[2])
    # headers = ['Sequence_Length', 'Reference', 'Mapped_Reads']
    # df = pd.DataFrame(rows, columns=headers)
    # pivot_df = df.pivot(index='Sequence_Length', columns='Reference', values='Mapped_Reads')
    ###
    # order is a list of column names to order 'stacking' in barplot
    # outpath should be full path with ending like '.pdf' or '.png
    # Note: .loc[:,['Jan','Feb', 'Mar']] is used here to rearrange the layer ordering
    import matplotlib.pyplot as plt
    # import pandas as pd

    pivot_df.loc[:, order].plot.bar(stacked=True, figsize=(10, 7), color=colors)  #  color=colors,
    plt.savefig(outpath)
    print('Output Stacked Barplot to:\n%s' % outpath)


def read_subtraction_by_mapping_single_end(data, d, dnames, aligndatabaselist, readfile, seqlength, totalreads):
    # data[0] == list of seq lengths, data[1] == list of reference genomes (molecule types),
    # data[2] == list of number of reads that map to each reference for that seq length
    # d is dictionary keys are reference sequence, values are list of numbers aligned to each reference
    # dnames is list of references in d
    # aligndatabaselist is list of full paths to indexed reference files
    # readfile is full path to .fastq seq file (assumed single end)
    # startseqlength == the smallest fragment size
    # hisat2(options)... | \
    # tee > (samtools flagstat - > hisat2_output.flagstat) | \
    # samtools sort - O BAM | \
    # tee hisat2_output.bam | \
    # samtools index - hisat2_output.bam.bai
    import os

    for count, reference in enumerate(aligndatabaselist):

        if count == 0:
            ### align with hisat2 ###
            samfile = run_hisat2(reference, readfile)
        else:
            ### align with hisat2 ###
            samfile = run_hisat2(reference, unmappedfastqfile)

        ### run samtools convert to bam, sort, index ###
        bamsortfile = run_samtools(samfile, cleanup=True)

        ### samtools, seperate unmapped reads from mapped reads ###
        unmappedbamfile, mappedbamfile = run_samtools_filer_by_flag(bamsortfile)

        ### samtools, count reads in file ###
        numunmappedreads = run_samtools_count_reads_in_file(unmappedbamfile)
        nummappedreads = run_samtools_count_reads_in_file(mappedbamfile)

        ### record number of mapped reads to reference ###
        path, ref = os.path.split(reference)
        # '_' + readfile.split('_')[-1].split('.')[0] should be _15bp or 23bp etc.
        d.setdefault(ref, []).append(nummappedreads/totalreads)
        dnames.append(ref)
        data[0].append(seqlength)  # count start at zero
        data[1].append(ref)
        data[2].append(nummappedreads/totalreads)

        ### samtools, convert bam to fastq ###
        unmappedfastqfile = run_samtools_bam2fastq(unmappedbamfile)

        if reference == aligndatabaselist[-1]:
            d.setdefault('Unmapped_Reads', []).append(numunmappedreads/totalreads)
            dnames.append('Unmapped_Reads')
            data[0].append(seqlength)  # count start at zero
            data[1].append('Unmapped_Reads')
            data[2].append(numunmappedreads/totalreads)

        ### cleanup ###
        os.remove(unmappedbamfile)
        os.remove(mappedbamfile)
        #os.remove(unmappedfastqfile)

    return data, d, dnames


def run_samtools_bam2fastq(bamfile, cleanup=False):
    import subprocess
    import os

    print('Start converting bam to fastq')
    fastqfile = '.'.join(bamfile.split('.')[:-1] + ['fastq'])
    with open(fastqfile, 'w') as OUT:
        cmd = 'samtools bam2fq %s' % (bamfile)  # samtools bam2fq %s > %s
        ps = subprocess.Popen(cmd.split(), stdout=OUT)
        ps.wait()
        print('Finished')
    print('Converted file:\n%s\nto\n%s' % (bamfile, fastqfile))
    if cleanup == True:
        os.remove(bamfile)  # delete bam file
    print('Finished')
    return fastqfile


def run_samtools_count_reads_in_file(bamfile, cleanup=False):
    import subprocess
    import os

    print('Start counting reads in file:\n%s' % bamfile)
    cmd = 'samtools view -c %s' % bamfile
    ps = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
    cmdout, err = ps.communicate()
    ps.wait()
    print('Finished')
    print('Number of reads mapped in file:\n%s = %d' % (bamfile, int(cmdout)))
    if cleanup == True:
        os.remove(bamfile)  # delete bam file
    print('Finished counting reads in file')
    return int(cmdout)


def run_samtools_filer_by_flag(bamsortfile, flag=4, cleanup=False):
    # bamsortfile is full path to file # assumes file ends with .sort.bam
    # flag 4 == read unmapped
    # samtools view -b -f 4 file.bam > unmapped.bam
    # samtools view -b -F 4 file.bam > mapped.bam # F == Not Matched
    # number of reads in bam file
    # samtools view -c SAMPLE.bam
    import subprocess
    import os

    print('Start Filtering file:\n%s\nwith flag: %d' % (bamsortfile, flag))
    matchflagfile = '.'.join(bamsortfile.split('.')[:-2] + ['Match%d' % flag, 'sort',  'bam'])  # unmapped
    notmatchflagfile = '.'.join(bamsortfile.split('.')[:-2] + ['NotMatch%d' % flag, 'sort',  'bam'])  # mapped
    with open(matchflagfile, 'w') as OUT:
        cmd = 'samtools view -h -b -f %d %s' % (flag, bamsortfile)
        ps = subprocess.Popen(cmd.split(), stdout=OUT)
        ps.wait()
    with open(notmatchflagfile, 'w') as OUT:
        cmd = 'samtools view -h -b -F %d %s' % (flag, bamsortfile)
        ps = subprocess.Popen(cmd.split(), stdout=OUT)
        ps.wait()

    if cleanup == True:
        os.remove(bamsortfile)  # delete bam file
    print('Finished')
    return matchflagfile, notmatchflagfile


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
    #     cmd = 'samtools view -h %s' % (sortbamfile)
    #     ps = subprocess.Popen(cmd.split(), stdout=OUT)
    #     ps.wait()
    # with open(flagstatfile, 'w') as OUT:
    #     # cmd = 'samtools flagstat %s > %s' % (sortsamfile, flagstatfile)
    #     cmd = 'samtools flagstat %s' % (sortsamfile)
    #     ps = subprocess.Popen(cmd.split(), stdout=OUT)
    #     ps.wait()

    print('Finished with Samtools\n')
    return sortbamfile


def run_hisat2(aligndatabase, fastqfile):
    print('Starting Hisat2: aligning\n%s' % fastqfile)
    print('Aligning reads to Reference:\n%s' % aligndatabase)
    import os
    import subprocess
    from pygentoolbox.Tools import make_directory

    refpath, referenceprefix = os.path.split(aligndatabase)
    path, f = os.path.split(fastqfile)
    # pathminusonedir, dir = os.path.split(path)
    make_directory(os.path.join(path, 'hisat2'))
    outsamfile = os.path.join(path, 'hisat2', '.'.join(f.split('.')[:-1] + [referenceprefix, 'sam']))
    cmd = 'hisat2 -q -x %s -U %s -S %s' % (aligndatabase, fastqfile, outsamfile)
    subprocess.call(cmd.split())

    print('Finished with Hisat2\n')
    return outsamfile


def run_blastn_match_db(fastafile, database, outformat=6, percentidentitiythreshold=98.00, bestbitscoreonly=True):
    # fasta file is full path to .fasta file
    # database is full path to blastn database
    # calculate % of human RNA spike-in by with BLASTn # this represents expected % free-floating RNA in sample
    print('Start BLASTn on file:\n%s\nTo Database:\n%s\n' % (fastafile, database))
    import os
    import subprocess
    from pygentoolbox.Tools import make_directory

    path, f = os.path.split(fastafile)
    pathminusonedir, dir = os.path.split(path)
    outpath = os.path.join(pathminusonedir, 'blastn')
    make_directory(outpath)
    if bestbitscoreonly == True:
        # takes the best BLAST result by bit score
        outfile = 'best_bit_score_per_query.blastn.RNA.tsv'
        fulloutpath = os.path.join(outpath, outfile)
        # full cmd = 'blastn -query %s -db %s -outfmt %d | sort -k1,1 -k12,12nr -k11,11n | sort -u -k1,1 --merge > %s' % (fastafile, database, outformat, fulloutpath)
        cmdpipe = ['blastn -query %s -db %s -outfmt %d' % (fastafile, database, outformat), 'sort -k1,1 -k12,12nr -k11,11n', 'sort -u -k1,1 --merge -']
        for count, cmd in enumerate(cmdpipe):
            if count == 0:
                ps = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
            elif count != len(cmdpipe)-1: # if it is not the last command
                ps = subprocess.Popen(cmd.split(), stdin=ps.stdout, stdout=subprocess.PIPE)
                ps.wait()
            else:  # it must be the last command
                with open(fulloutpath, 'w') as OUT:
                    ps = subprocess.Popen(cmd.split(), stdin=ps.stdout, stdout=OUT)
                    ps.wait()

    else:
        ### won't work right now
        outfile = 'All_hits_per_query.blastn.RNA.tsv'
        fulloutpath = os.path.join(outpath, outfile)
        cmd = 'blastn -query %s -db %s -outfmt %d > %s' % (fastafile, database, outformat, fulloutpath)
        cmdlist = cmd.split()
        subprocess.call(cmdlist)
    outdbmatching = os.path.join(outpath, 'human.rna.freefloating.tsv')
    # cmd = 'awk -F \"\t\" \'$3 > %f {print $1}\' %s > %s' % (percentidentitiythreshold, fulloutpath, outdbmatching)
    # subprocess.call(cmd.split())  # doesnt work with ' characters somehow
    with open(fulloutpath, 'r') as FILE:
        output = [line.strip() for line in FILE if float(line.split('\t')[2]) > percentidentitiythreshold]
    with open(outdbmatching, 'w') as OUT:
        OUT.write('\n'.join(output))
    print('Finished BLASTn, output results to:\n%s\n%s\n' % (fulloutpath, outdbmatching))

    return fulloutpath, outdbmatching


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


def main(totalreads, fileextension='.fastq.gz', directory='', outplotname='MyOutputPlot.pdf', startseqlength=15,
         cleanreads=True, blastdatabase='/media/sf_LinuxShare/Humans/Genome/Seqs/GRCh38_top_level.fa',
         aligndatabaselist=['/media/sf_LinuxShare/Ciliates/Genomes/Hisat2_Indexes/Kpne_CompleteGenome',
                            '/media/sf_LinuxShare/Ciliates/Genomes/Hisat2_Indexes/Pt_mtGenome',
                            '/media/sf_LinuxShare/Ciliates/Genomes/Hisat2_Indexes/Pt_51_Mac',
                            '/media/sf_LinuxShare/Ciliates/Genomes/Hisat2_Indexes/Pt_51_MacAndIES',
                            '/media/sf_LinuxShare/Ciliates/Genomes/Hisat2_Indexes/Pt_51_Mic2'],
         colors=['black', 'orange', 'forestgreen', 'red', 'dodgerblue', 'grey']):
    ### forwardreadfile='/media/sf_LinuxShare/Projects/Lyna/TestMyPipe/500_LK_L1_R1.fastq.gz',
    ### reversereadfile='/media/sf_LinuxShare/Projects/Lyna/TestMyPipe/500_LK_L1_R2.fastq.gz',

    # colors is length of align databases + 1 because we need a color for unmapped
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
        files = glob.glob('*%s' % fileextension)
        filelist = [os.path.join(directory, f) for f in files]
    else:
        # fullpath = os.path.join(directory, '*%s' % fileextension)
        # glob.glob(fullpath)
        filelist = []
        for file in os.listdir(directory):
            if file.endswith(fileextension):
                filelist.append(os.path.join(directory, file))

    d = {}  # will have referencename as key and list as value with counts
    dnames = []
    data = [[], [], []]  # data[0] == list of seq lengths, data[1] == list of reference genomes (molecule types),
                         # data[2] == list of number of reads that map to each reference for that seq length
    for count, file in enumerate(filelist):
        if cleanreads == False:
            ### start fastp ###
            #forwardout, reverseout = run_fastp_paired_end(forwardreadfile, reversereadfile)
            forwardclean = run_fastp_single_end(file)
        else:
            print('yep your reads must already be clean')
            forwardclean = file

        ### start gzip -dk file.trim.fastq.gz
        forwardtrim = run_gzip_decompress(forwardclean)
        # reversetrim = run_gzip_decompress(reverseout)

        ### start fastqc ###
        run_fastqc(forwardtrim)
        # run_fastqc(reversetrim)

        ### align reads to each genome sequentially, use only unmapped reads for next genome ###
        seqlength = startseqlength + count
        data, d, dnames = read_subtraction_by_mapping_single_end(data, d, dnames, aligndatabaselist, forwardtrim,
                                                                 seqlength, totalreads)
        # fullpath to fasta file, outformat == 6, percent identity threshold to determine if read is from a species,
        # BestBitScoreOnly?? == True
        # blastoutpath, blastoutdbmatching = run_blastn_match_db(rnafastafile, blastdatabase, 6, 98.00, True)

    df = pd.DataFrame(d)
    outpathpandas = os.path.join(directory, 'ReadsMappedToEachReference.tsv')
    df.to_csv(outpathpandas, sep='\t')
    print('wrote out table of counts aligned to each genome to:\n%s\n' % outpathpandas)

    # with pandas make stacked bar plot
    rows = zip(data[0], data[1], data[2])
    headers = ['sRNA_Size', 'Reference', 'Mapped_Reads']
    df = pd.DataFrame(rows, columns=headers)
    df.to_csv(outpathpandas + '.pivot', sep='\t')
    pivot_df = df.pivot(index='sRNA_Size', columns='Reference', values='Mapped_Reads')
    ref_prefixes = []
    # keeps only the ref prefix name 'Pt_51_Mac' from '/media/sf_LinuxShare/Ciliates/Genomes/Hisat2_Indexes/Pt_51_Mac'
    for fullpath in aligndatabaselist:
        path, prefix = os.path.split(fullpath)
        ref_prefixes.append(prefix)
    ref_prefixes = ref_prefixes + ['Unmapped_Reads']
    outpath = os.path.join(directory, outplotname)
    make_pandas_stacked_barplot(pivot_df, ref_prefixes, outpath, colors)





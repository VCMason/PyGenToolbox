
def make_directory(dirName):
    # dirName is directory in cwd or full path to directory
    import os

    if not os.path.exists(dirName):
        os.mkdir(dirName)
        print("Directory ", dirName, " Created ")
    else:
        print("Directory ", dirName, " already exists")


def run_blastn_match_db(fastafile, database, outformat=6):
    # fasta file is full path to .fasta file
    # database is full path to blastn database
    # minlen is minimum length allow for match
    # pid is percent identitiy threshold (match must be greater than this float value
    # bbs is best bit score only, True or False
    # calculate % of human RNA spike-in by with BLASTn # this represents expected % free-floating RNA in sample
    print('Start BLASTn on file:\n%s\nTo Database:\n%s\n' % (fastafile, database))
    import subprocess
    import os

    path, f = os.path.split(fastafile)
    pathminusonedir, dir = os.path.split(path)
    outpath = os.path.join(pathminusonedir, 'blastn')
    make_directory(outpath)
    if bbs == True:
        # takes the best BLAST result by bit score
        outfile = 'best_bit_score_per_query.blastn.RNA.tsv'
        fulloutpath = os.path.join(outpath, outfile)
        # full cmd = 'blastn -query %s -db %s -outfmt %d | sort -k1,1 -k12,12nr -k11,11n | sort -u -k1,1 --merge > %s' % (fastafile, database, outformat, fulloutpath)
        cmdpipe = ['blastn -query %s -db %s -outfmt %d' % (fastafile, database, outformat), 'sort -k1,1 -k12,12nr -k11,11n %s' % (fulloutpath + '.blastn'), 'sort -u -k1,1 --merge %s' % (fulloutpath + '.sort1')]
        for count, cmd in enumerate(cmdpipe):
            if count == 0:
                with open(fulloutpath + '.blastn', 'w') as OUT:
                    print('Pipe step 1.1')
                    print(cmd)
                    ps = subprocess.Popen(cmd.split(), bufsize=-1, stdout=OUT)
                    print('Pipe step 1.2')
                    ps.wait()
                    print('Pipe step 1.3')
            elif count != len(cmdpipe) - 1:  # if it is not the last command
                with open(fulloutpath + '.sort1', 'w') as OUT:
                    print('Pipe step 2.1')
                    ps = subprocess.Popen(cmd.split(), bufsize=-1, stdout=OUT)
                    print('Pipe step 2.2')
                    ps.wait()
                    print('Pipe step 2.3')
            else:  # it must be the last command
                with open(fulloutpath, 'w') as OUT:
                    print('Pipe step 3.1')
                    ps = subprocess.Popen(cmd.split(), bufsize=-1, stdout=OUT)
                    print('Pipe step 3.2')
                    ps.wait()
                    print('Pipe step 3.3')

    outdbmatching = os.path.join(outpath, 'human.rna.freefloating.tsv')
    cmd = 'awk -F \"\t\" \'$3 > %f {print $1}\' %s > %s' % (pid, fulloutpath, outdbmatching)
    subprocess.call(cmd.split())  # doesnt work with ' characters somehow
    with open(fulloutpath, 'r') as FILE:
        output = [line.strip() for line in FILE if (float(line.split('\t')[2]) > pid) and (float(line.split('\t')[3] > minlen))]
    with open(outdbmatching, 'w') as OUT:
        OUT.write('\n'.join(output))
    print('Finished BLASTn, output results to:\n%s\n%s\n' % (fulloutpath, outdbmatching))

    return fulloutpath, outdbmatching


def run_blast(fastafile, database, task='blastn', outformat=6, minlen=11, pid=90.00):
    # fasta file is full path to .fasta file
    # database is full path to blastn database
    # minlen is minimum length allow for match
    # pid is percent identitiy threshold (match must be greater than this float value
    # bbs is best bit score only, True or False

    print('Start BLASTn on file:\n%s\nTo Database:\n%s\n' % (fastafile, database))
    import subprocess
    import os

    path, f = os.path.split(fastafile)
    pathminusonedir, dir = os.path.split(path)
    outpath = os.path.join(pathminusonedir, 'blastn')
    make_directory(outpath)

    ### won't work right now
    outfile = 'All_hits_per_query.blastn.RNA.tsv'
    fulloutpath = os.path.join(outpath, outfile)
    cmd = 'blastn -task %s -query %s -db %s -word_size %d -perc_identity %f -outfmt %d > %s' % (task, fastafile, database, minlen, pid, outformat, fulloutpath)
    cmdlist = cmd.split()
    subprocess.call(cmdlist)

    return outpath


def main(fastafile, blastdatabase, task='blastn', outformat=6, minlen=11, pid=90.00, bbs=False):
    # makeblastdb -in IES.pt51.MacIES.DCL23.0.1.up.fa -dbtype nucl
    # fastafile is query # RNA from Dcr1 dependent IESs
    # blastdatabase is blast database # RNA from contacts near IESs
    # bbs is bestbitscoreonly limits the results to only the first match with the highest bit score

    if bbs == False:
        outpath = run_blast(fastafile, blastdatabase, task, outformat, minlen, pid)
    elif bbs == True:
        blastoutpath, blastoutdbmatching = run_blastn_match_db(fastafile, blastdatabase, 6)  # fullpath to fasta file, outformat == 6, percent identity threshold to determine if read is from a species, BestBitScoreOnly?? == True

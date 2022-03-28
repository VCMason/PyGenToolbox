

def read_fastq_and_split_by_length(fullpath, minlength=0, maxlength=0, onefile=False):  # fullpath = string, full path to file
    """ Reads .fastq file by every four lines, outputs each four lines to files by sequence length """
    import os

    path, f = os.path.split(fullpath)
    FILE = open(fullpath, 'r')
    line = FILE.readline()
    count = 0
    outtemp = []
    d = {}
    while line:
        if count % 4000000 == 0:
            print('On sequence: %d' % (count/4))
        if count % 4 == 0:  # if the entire line is only a +
            # if os.path.exists(fout):  # if the output file exists delete it before appending
            #    os.remove(fout)
            output = line.strip() + '\n'
        elif count % 4 == 1:
            length = len(line.strip())
            output += line.strip() + '\n'
        elif count % 4 == 2:
            output += line.strip() + '\n'
        elif count % 4 == 3:
            output += line.strip() + '\n'
            outtemp.append(output)
            if len(outtemp) == 1000:
                if (minlength == 0) and (maxlength == 0):  # then write out all sequences
                    if onefile is False:
                        fout = os.path.join(path, '.'.join(f.split('.')[:-1]) + '%dbp.fastq' % length)
                        with open(fout, 'a') as OUT:
                            OUT.write('\n'.join(outtemp) + '\n')
                    elif onefile is True:
                        fout = os.path.join(path, '.'.join(f.split('.')[:-1]) + '%dbp%dbp.fastq' % (minlength, maxlength))
                        with open(fout, 'a') as OUT:
                            OUT.write('\n'.join(outtemp) + '\n')
                elif (length >= minlength) and (length <= maxlength):  # then write out only sequences of length = desiredlength
                    if onefile is False:
                        fout = os.path.join(path, '.'.join(f.split('.')[:-1]) + '%dbp.fastq' % length)
                        with open(fout, 'a') as OUT:
                            OUT.write('\n'.join(outtemp) + '\n')
                    elif onefile is True:
                        fout = os.path.join(path, '.'.join(f.split('.')[:-1]) + '%dbp%dbp.fastq' % (minlength, maxlength))
                        with open(fout, 'a') as OUT:
                            OUT.write('\n'.join(outtemp) + '\n')
                outtemp = []  # reset temporary output chunk to zero
        line = FILE.readline()
        count += 1
    FILE.close()

    # output any straglers
    if len(outtemp) > 0:
        if (minlength == 0) and (maxlength == 0):  # then write out all sequences
            if onefile is False:
                fout = os.path.join(path, '.'.join(f.split('.')[:-1]) + '%dbp.fastq' % length)
                with open(fout, 'a') as OUT:
                    OUT.write('\n'.join(outtemp) + '\n')
            elif onefile is True:
                fout = os.path.join(path, '.'.join(f.split('.')[:-1]) + '%dbp%dbp.fastq' % (minlength, maxlength))
                with open(fout, 'a') as OUT:
                    OUT.write('\n'.join(outtemp) + '\n')
        elif (length >= minlength) and (length <= maxlength):  # then write out only sequences of length = desiredlength
            if onefile is False:
                fout = os.path.join(path, '.'.join(f.split('.')[:-1]) + '%dbp.fastq' % length)
                with open(fout, 'a') as OUT:
                    OUT.write('\n'.join(outtemp) + '\n')
            elif onefile is True:
                fout = os.path.join(path, '.'.join(f.split('.')[:-1]) + '%dbp%dbp.fastq' % (minlength, maxlength))
                with open(fout, 'a') as OUT:
                    OUT.write('\n'.join(outtemp) + '\n')

    print('Made one fastq file for each sequence length:\nInput file: %s\nExample output filename: %s' % (fullpath, fout))
    print('Number of lines in input file: %d' % count)


def main(fullpath, minlength=0, maxlength=0, onefile=False):
    # fullpath is full path to .fastq file
    # set min and max lengths to zero to output all seqs of all lengths, slow...
    read_fastq_and_split_by_length(fullpath, minlength, maxlength, onefile)

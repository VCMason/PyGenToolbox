

def read_fastq_and_split_by_length(fullpath):  # fullpath = string, full path to file
    """ Reads .fastq file by every four lines, outputs each four lines to files by sequence length """
    import os
    import re

    path, f = os.path.split(fullpath)
    FILE = open(fullpath, 'r')
    line = FILE.readline()
    count = 0
    countout = 0
    d = {}
    while line:
        if count % 400000 == 0:
            print('On sequence: %d' % (count/4))
        if count % 4 == 0:  # if the entire line is only a +
            # if os.path.exists(fout):  # if the output file exists delete it before appending
            #    os.remove(fout)
            output = line.strip() + '\n'
        elif count % 4 == 1:
            length = len(line.strip())
            if re.match(r"T[ATCG]G", line.strip()[:3]):
            #if re.match(r"A[ATCG]C", line.strip()[:3]):
                m = 'forward'
            elif re.match(r"C[ATCG]A", line.strip()[-3:]):
                m = 'reverse'
            else:
                m = 'none'
            output += line.strip() + '\n'
        elif count % 4 == 2:
            output += line.strip() + '\n'
        elif count % 4 == 3:
            output += line.strip() + '\n'
            if length == 25:
                # if (m == 'forward') or (m == 'reverse'):
                if m == 'forward':
                    fout = os.path.join(path, '.'.join(f.split('.')[:-1]) + '.scnRNA.fastq')
                    with open(fout, 'a') as OUT:
                        OUT.write(output)
                    countout += 1

        line = FILE.readline()
        count += 1

    FILE.close()

    if countout == 0:
        fout = 'No scnRNA found'

    print('Made one fastq file for each sequence length:\nInput file: %s\nExample output filename: %s' % (fullpath, fout))
    print('Number of sequences in input file: %d' % (count/4))
    print('Number of sequences in output file: %d' % countout)


def main(fullpath):
    # fullpath is full path to .fastq file
    read_fastq_and_split_by_length(fullpath)

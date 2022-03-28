def main(f, chunk=4000000, gz=False):
    # f = str, full path to file
    # chunk = int ## chunk is number of lines per split file, default 400000, multiples of 4 work well for fastq files
    # f is full path to file # limit is line number limit 1st line is 1 NOT 0 # still outputs the linecount <= limit

    if gz is False:
        with open(f, 'r') as FILE:
            splitcount = 0
            splitfilenames = []
            outlines = []
            for count, line in enumerate(FILE, start=1):
                if count % chunk < chunk:
                    outlines.append(line)
                elif count % chunk == 0:
                    # then append last line and split file
                    outlines.append(line)
                    outf = '.'.join(f.split('.')[:-1] + [f'split{splitcount}', f.split('.')[-1]])
                    splitfilenames.append(outf)
                    with open(outf, 'w') as OUT:
                        OUT.write(''.join(outlines))
                    outlines = []
                    splitcount += 1
        if len(outlines) > 0:  # if some lines were not written to file b/c len(outlines) < chunk value
            outf = '.'.join(f.split('.')[:-1] + [f'split{splitcount}', f.split('.')[-1]])
            splitfilenames.append(outf)
            with open(outf, 'w') as OUT:
                OUT.write(''.join(outlines))
            outlines = []

    elif gz is True:
        import gzip
        with gzip.open(f, 'rb') as FILE:
            splitcount = 0
            splitfilenames = []
            outlines = []
            for count, line in enumerate(FILE, start=1):
                if count % chunk < chunk:
                    outlines.append(line)
                elif count % chunk == 0:
                    # then append last line and split file
                    outlines.append(line)
                    outf = '.'.join(f.split('.')[:-1] + [f'split{splitcount}', f.split('.')[-1]])
                    splitfilenames.append(outf)
                    with gzip.open(outf, 'wb') as OUT:
                        OUT.write(b''.join(outlines))
                    outlines = []
                    splitcount += 1
        if len(outlines) > 0:  # if some lines were not written to file b/c len(outlines) < chunk value
            outf = '.'.join(f.split('.')[:-1] + [f'split{splitcount}', f.split('.')[-1]])
            splitfilenames.append(outf)
            with gzip.open(outf, 'wb') as OUT:
                OUT.write(b''.join(outlines))
            outlines = []

    else:
        print('gz parameter incorrectly specified, nothing happened')

    print(f)
    print(f'Split into {splitcount+1} files')
    print('Finished')

    return splitfilenames, splitcount+1

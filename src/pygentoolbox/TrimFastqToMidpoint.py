

def trim_fastq_to_midpoint_gz(fastqfile, outfile, n=50, seqlenmin=125, seqlenmax=175):
    # assume all sequences are complete fragments.
    # n is length of midpoint sequence
    print('Assuming all sequences are complete sequences fragments.')
    print(f'Trimming front and end of all sequences to midpoints. Midpoint length = {n}')
    print(f'Input file: {fastqfile}')

    import gzip
    import math

    with gzip.open(outfile, 'wb') as OUT:
        OUT.write(b'')

    with gzip.open(outfile, 'ab') as OUT:
        outlines = []
        count, total = 0, 0
        with gzip.open(fastqfile, 'rb') as GZ:
            for lnum, line in enumerate(GZ):
                # ord('i') converts ASCII to int value, chr(105) converts int to ASCII
                if (lnum % 4 == 0) or (lnum % 4 == 2):  # (line[0] == ord('@')) or (line[0] == ord('+')):
                    name = line.strip()
                elif (lnum % 4 == 1) or (lnum % 4 == 3):
                    seqlength = len(line.strip())
                    if (seqlength >= seqlenmin) and (seqlength <= seqlenmax):
                        outlines.append(name)
                        if seqlength / n <= 1:
                            outlines.append(line.strip())
                        else:
                            extrabases = seqlength - n
                            trimlen = extrabases / 2
                            # rounding ensures that 1 extra base will be removed from end of seq if trimlen is float ex:50.5
                            outlines.append(line.strip()[math.floor(trimlen):-1 * math.ceil(trimlen)])
                count += 1

                # output data every 10,000 lines
                if (count >= 1000000) and (count % 1000000 >= 1):
                    if total == 0:
                        OUT.write(b'\n'.join(outlines))
                    else:
                        OUT.write(b'\n' + b'\n'.join(outlines))
                    total += count
                    print(total)
                    outlines = []
                    count = 0
        # if data left over output last lines
        if len(outlines) > 0:
            OUT.write(b'\n' + b'\n'.join(outlines) + b'\n')

    print(f'Midpoint sequences output to file: {outfile}')

    return


def trim_fastq_to_midpoint(fastqfile, outfile, n=50):
    # assume all sequences are complete fragments.
    # n is length of midpoint sequence
    print('Assuming all sequences are complete sequences fragments.')
    print(f'Trimming front and end of all sequences to midpoints. Midpoint length = {n}')
    print(f'Input file: {fastqfile}')

    import math

    outlines = []
    with open(fastqfile, 'r') as FILE:
        for line in FILE:
            if line[0] == '>':
                outlines.append(line.strip())
            else:
                seqlength = len(line.strip())
                if seqlength / n <= 1:
                    outlines.append(line.strip())
                else:
                    extrabases = seqlength - n
                    trimlen = extrabases / 2
                    # rounding ensures that 1 extra base will be removed from end of seq if trimlen is float (ex: 50.5)
                    outlines.append(line.strip()[math.floor(trimlen):-1*math.ceil(trimlen)])

    with open(outfile, 'w') as OUT:
        OUT.write('\n'.join(outlines))

    print(f'Midpoint sequences output to file: {outfile}')

    return


def main(fastqfile, outfile, midpointlength=50, seqlenmin=125, seqlenmax=175):
    if fastqfile[-3:] == '.gz':
        trim_fastq_to_midpoint_gz(fastqfile, outfile, midpointlength, seqlenmin, seqlenmax)
    else:
        trim_fastq_to_midpoint(fastqfile, outfile, midpointlength)


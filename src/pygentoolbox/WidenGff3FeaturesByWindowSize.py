
def main(f, outfile, windowlen):
    # f == full path to annotation file
    # outfile is full path
    # window len == total length of window to widen feature, if 400 then the feature will be widened s-200 and e+200
    output = []
    with open(f, 'r') as FILE:
        for line in FILE:
            if line[0] == '#':
                output.append(line.strip())
            else:
                newstart = int(int(line.strip().split('\t')[3]) - (0.5 * windowlen))  # [3] is start coord
                if newstart < 1:  # gff3 files are 1 based
                    newstart = 1
                newend = int(int(line.strip().split('\t')[4]) + (0.5 * windowlen))  # [4] is end coord # theoretically this could be longer than the scaffold
                # new line with new start and end coordinates for the feature
                output.append('\t'.join(line.strip().split('\t')[:3] + [str(newstart), str(newend)] + line.strip().split('\t')[5:]))
    with open(outfile, 'w') as OUT:
        OUT.write('\n'.join(output))


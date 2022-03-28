
infile = 'D:/LinuxShare/Ciliates/Genomes/Annotations/ptetraurelia_mac_51_annotation_v2.0.gff3'
outfile = '.'.join(infile.split('.')[:-1] + ['geneid', 'gff3'])
with open(infile, 'r') as FILE:
    output = []
    for line in FILE:
        if line[0] == '#':
            output.append(line.strip())
        else:
            if line.strip().split('\t')[2] == 'gene':
                geneid = line.strip().split('\t')[8].split(';')[0][3:]
            outline = line.strip() + ';gene_id=%s' % geneid
            output.append(outline)

with open(outfile, 'w') as OUT:
    OUT.write('\n'.join(output))

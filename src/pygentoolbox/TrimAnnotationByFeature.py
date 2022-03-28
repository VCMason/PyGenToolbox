
def keep_one_mRNA_with_same_stop_codon_position(f):
    output = []
    nameset = set()
    with open(f, 'r') as FILE:
        for line in FILE:
            if line[0] == '#':
                output.append(line.strip())
            else:
                s = line.strip().split('\t')
                if s[6] == '+':  # then stop codons are on the right
                    name = s[8].split(';')[1] + '.' + s[4]
                    if name in nameset:
                        pass
                    else:
                        nameset.add(name)
                        output.append(line.strip())
                elif s[6] == '-':  # then stop codons are on the left
                    name = s[8].split(';')[1] + '.' + s[3]
                    if name in nameset:
                        pass
                    else:
                        nameset.add(name)
                        output.append(line.strip())
    outfile = '.'.join(f.split('.')[:-1] + ['onestop', 'gff3'])
    with open(outfile, 'w') as OUT:
        OUT.write('\n'.join(output))
    print('Output file: %s' % outfile)


def trim_to_mRNA_and_add_protein_id(f, feature='mRNA'):
    output = []
    with open(f, 'r') as FILE:
        for line in FILE:
            if line[0] == '#':
                output.append(line.strip())

            else:
                if line.strip().split('\t')[2] == feature:
                    outline = line.strip()
                    gate = 'open'  # allow it to find the next CDS line, which should have protein ID
                elif (line.strip().split('\t')[2] == 'CDS') and (gate == 'open'):
                    # extract: XP_011540840.1
                    # from: ID=cds4;Parent=rna43;Dbxref=GeneID:105378947,Genbank:XP_011540840.1;Name=XP_011540840.1;
                    outline += ';%s' % (line.strip().split('\t')[8].split(';')[3][5:])
                    gate = 'closed'
                    output.append(outline)

    outfile = '.'.join(f.split('.')[:-1] + [feature, 'wProtID', 'gff3'])
    with open(outfile, 'w') as OUT:
        OUT.write('\n'.join(output))
    print('Output file: %s' % outfile)


def trim_ensembl_gtf_to_3utrs_of_target_genes(f, feature='three_prime_utr', targetgenenames=[], minpeplength=5):
    # this is the worst piece of shit i have ever written <3
    # f is full path to ensembl .gtf file
    # example ensembl .gtf 3' UTR feature:
    # 1	havana	three_prime_utr	944154	944259	.	+	.	gene_id "ENSG00000187634"; gene_version "13"; transcript_id "ENST00000455979"; transcript_version "1"; gene_name "SAMD11"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "SAMD11-204"; transcript_source "havana"; transcript_biotype "protein_coding"; tag "cds_start_NF"; tag "mRNA_start_NF"; transcript_support_level "2";

    collectedfeatures = []
    stopcodons = []
    with open(f, 'r') as FILE:
        for line in FILE:
            if line[0] == '#':
                collectedfeatures.append(line.strip())
            else:
                if line.strip().split('\t')[2] == feature:
                    orientation = line.strip().split('\t')[6]
                    metadata = line.strip().split('\t')[8].split('"; ')
                    genename = ''
                    for meta in metadata:
                        if meta[:len('gene_name "')] == 'gene_name "':
                            genename = meta[len('gene_name "'):]
                    if (genename != '') and (genename in targetgenenames):
                        collectedfeatures.append(line.strip())

    output = []
    outputgenes = []
    # enename_startposition value is line of longest 3'UTR for each unique start
    for line in collectedfeatures:
        if line[0] == '#':
            output.append(line)
        else:
            gate = 'open'
            if line.split('\t')[6] == '+':
                start = line.split('\t')[3]
                end = line.split('\t')[4]
                orientation = line.split('\t')[6]
                difference = int(end) - int(start) + 1
            elif line.split('\t')[6] == '-':
                end = line.split('\t')[3]
                start = line.split('\t')[4]
                orientation = line.split('\t')[6]
                difference = int(start) - int(end) + 1
            if difference >= minpeplength * 3:  # it needs to be at least minpeplength AA long
                metadata = line.split('\t')[8].split('"; ')
                for meta in metadata:
                    if meta[:len('gene_name "')] == 'gene_name "':
                        genename = meta[len('gene_name "'):]
                for line2 in collectedfeatures:
                    if line2[0] == '#':
                        pass
                    else:
                        metadata2 = line2.split('\t')[8].split('"; ')
                        for meta2 in metadata2:
                            if meta2[:len('gene_name "')] == 'gene_name "':
                                genename2 = meta2[len('gene_name "'):]
                        if line2.split('\t')[6] == '+':
                            start2 = line2.split('\t')[3]
                            end2 = line2.split('\t')[4]
                            if (genename == genename2) and (start == start2) and (int(end) <= int(end2)):
                                if int(end) < int(end2):
                                    gate = 'closed'
                                elif int(end) == int(end2):
                                    for out in output:
                                        if out[0] == '#':
                                            pass
                                        else:
                                            metadata3 = out.split('\t')[8].split('"; ')
                                            for meta3 in metadata3:
                                                if meta3[:len('gene_name "')] == 'gene_name "':
                                                    genename3 = meta3[len('gene_name "'):]
                                            if line.split('\t')[6] == '+':
                                                start3 = out.split('\t')[3]
                                                end3 = out.split('\t')[4]
                                                # if the same gene with same 3'UTR start & end coordinates is already in output then
                                                # do not output a duplicate 3UTR feature with same start and end (close gate)
                                                # However if no gene with exact start and end is in output then output the first one
                                                if (genename == genename3) and (start == start3) and (end == end3):
                                                    gate = 'closed'
                        elif line2.split('\t')[6] == '-':
                            end2 = line2.split('\t')[3]
                            start2 = line2.split('\t')[4]
                            if (genename == genename2) and (start == start2) and (int(end) >= int(end2)):
                                if int(end) > int(end2):
                                    gate = 'closed'
                                elif int(end) == int(end2):
                                    for out in output:
                                        if out[0] == '#':
                                            pass
                                        else:
                                            metadata3 = out.split('\t')[8].split('"; ')
                                            for meta3 in metadata3:
                                                if meta3[:len('gene_name "')] == 'gene_name "':
                                                    genename3 = meta3[len('gene_name "'):]
                                            if line.split('\t')[6] == '-':
                                                end3 = out.split('\t')[3]
                                                start3 = out.split('\t')[4]
                                                # if the same gene with same 3'UTR start & end coordinates is already in output then
                                                # do not output a duplicate 3UTR feature with same start and end (close gate)
                                                # However if no gene with exact start and end is in output then output the first one
                                                if (genename == genename3) and (end == end3) and (start == start3):
                                                    gate = 'closed'

                        # if (genename == genename2) and (start == start2) and (x == ''):  # int(end) <= int(end2)
                        #     if int(end) < int(end2):
                        #         gate = 'closed'
                        #     elif int(end) == int(end2):
                        #         for out in output:
                        #             if out[0] == '#':
                        #                 pass
                        #             else:
                        #                 start3 = out.split('\t')[3]
                        #                 end3 = out.split('\t')[4]
                        #                 metadata3 = out.split('\t')[8].split('"; ')
                        #                 for meta3 in metadata3:
                        #                     if meta3[:len('gene_name "')] == 'gene_name "':
                        #                         genename3 = meta3[len('gene_name "'):]
                        #                 # if the same gene with same 3'UTR start & end coordinates is already in output then
                        #                 # do not output a duplicate 3UTR feature with same start and end (close gate)
                        #                 # However if no gene with exact start and end is in output then output the first one
                        #                 if (genename == genename3) and (start == start3) and (int(end) == int(end3)):
                        #                     gate = 'closed'
            else:
                gate = 'closed'

            if (gate == 'open') and (line[0] != '#'):
                output.append(line)
                outputgenes.append(genename)

    print(f'Number of gene names entered and searched for {len(targetgenenames)}')
    print(f'Number of genes with a transcript in output: {len(set(outputgenes))}')
    print(f'Genes with a transcript in output: {set(outputgenes)}')
    outfile = '.'.join(f.split('.')[:-1] + [feature, 'all', 'gtf'])
    with open(outfile, 'w') as OUT:
        OUT.write('\n'.join(collectedfeatures))
    print('Output file: %s' % outfile)

    outfile = '.'.join(f.split('.')[:-1] + [feature, 'longest', 'gtf'])
    with open(outfile, 'w') as OUT:
        OUT.write('\n'.join(output))
    print('Output file: %s' % outfile)


##### broken #####
class GENE:
    def __init__(self, line):
        # line is line form annotation file
        if line.strip().split('\t')[2] == 'gene':
            self.genename = line.strip().split('\t')[8].split(';')[3][len('Name='):].upper()
            self.scaffold = line.strip().split('\t')[0]
            self.orientation = line.strip().split('\t')[6]
            self.biotype = line.strip().split('\t')

        elif line.strip().split('\t')[2] == 'mRNA':
            self.mrnafeaturelength = int(line.strip().split('\t')[4]) - int(line.strip().split('\t')[3])
        self.start = int(line.strip().split('\t')[3])
        self.end = int(line.strip().split('\t')[4])
        self.length = len(line.strip().split('\t')[8].split(';')[-1][len('sequence='):])
        self.sequence = line.strip().split('\t')[8].split(';')[-1][len('sequence='):]


def main(f, feature='gene'):
    output = []
    with open(f, 'r') as FILE:
        for line in FILE:
            if line[0] == '#':
                output.append(line.strip())
            else:
                if line.strip().split('\t')[2] == feature:
                    output.append(line.strip())
    outfile = '.'.join(f.split('.')[:-1] + [feature, 'gff3'])
    with open(outfile, 'w') as OUT:
        OUT.write('\n'.join(output))
    print('Output file: %s' % outfile)


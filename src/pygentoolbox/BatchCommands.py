def call_consensus_sequences_mask_reference_bases(roots, ref, bams):
    ''' roots is a list of base names '52_Late', ref is full path to reference .fa (or just filename if in cwd) '''
    ''' bams is list of full paths to alignment files (or filenames if in cwd) '''
    ''' output files to cwd '''
    ''' requires bcftools, and cat command'''
    ''' unfinished function'''
    import os
    for root in roots:

        # call variants
        os.system('bcftools mpileup -Ou -f %s %s | bcftools call -c -Oz -o %s_calls.vcf.gz' % (ref, ' '.join(bams), root))
        # os.system('bcftools mpileup -Ou -f %s %s | bcftools call -mv -Oz -o %s_calls.vcf.gz' % (ref, ' '.join(bams), root))
        os.system('bcftools index %s_calls.vcf.gz' % (root))

        # normalize indels
        # os.system('bcftools norm - f %s %s_calls.vcf.gz - Ob - o %s_calls.norm.bcf' % (ref, root,  root))

        # filter adjacent indels within 5bp
        # os.system('bcftools filter --IndelGap 5 %s_calls.norm.bcf - Ob - o %s_calls.norm.flt-indels.bcf' % (root, root))

        #call consensus
        os.system('cat %s | bcftools consensus %s_calls.vcf.gz > %s_cns.fa' % (ref, root, root))

def align_star_se(nthreads, genome, fqs):
    '''  map reads to reference genome with STAR with single end '''
    ''' genome is full path to indexed genome file /path/to/genomeprefix '''
    ''' fqs is list of fastq files '''
    ''' assumes a unique prefix/root name is present right before .fastq in fq files blah.blah.UNINQUEPREFIX.fastq '''
    ''' get roots from fastq files '''

    roots = [fq.split('.')[-2:-1] for fq in fqs]
    print('Unique prefix sequences: %s' % ('\n'.join(roots)))
    print('Start STAR single end alignment.')
    
    import os
    for root in roots:
        os.system('STAR --runThreadN %d --genomeDir %s --readFilesIn %s --outSAMtype BAM SortedByCoordinate --outTmpDir ./STARTempFiles/%s --outFileNamePrefix %s' % (nthreads, genome, root, root, root))

    return roots
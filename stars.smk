###############################################################################
## Module: star_2nd -- Run STAR with 2-pass mode                             ##
## Version 0.1.0                                                             ##
## Written by Sen Li, sen.li_at_sund.ku.dk                                   ##
## Notes: alway use *.fq.gz to save space                                    ##
###############################################################################

rule star_2nd:
	input:
		read1="{WorkDir}/01_QualityFiltered/{sample}_filtered_1.fq.gz",
		read2="{WorkDir}/01_QualityFiltered/{sample}_filtered_2.fq.gz"
	params:
		sjtab=expand("{sjtab}", sjtab=config['sjtab']),
		prefix="{WorkDir}/03_StarMap2nd/{sample}",
		StarMinN=expand("{StarMinN}", StarMinN=config['StarMinN']),
		StarMinS=expand("{StarMinS}", StarMinS=config['StarMinS']),
		StarRef=expand("{StarRef}", StarRef=config['StarRef'])
	output:
		endfile="{WorkDir}/03_StarMap2nd/{sample}.done"
	threads: 36
	shell:
		"""
		module unload gcc/8.2.0 star/2.7.2b
		module load gcc/8.2.0 star/2.7.2b pigz/2.3.4
		STAR --genomeDir {params.StarRef} --readFilesIn {input.read1} {input.read2} --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread {params.StarMinN} --outFilterScoreMinOverLread {params.StarMinS} --twopassMode Basic --sjdbFileChrStartEnd {params.sjtab} --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 16000000000 --runThreadN {threads} --outFileNamePrefix {params.prefix}_
		mv {params.prefix}_Unmapped.out.mate1 {params.prefix}_unmapped_1.fq
		mv {params.prefix}_Unmapped.out.mate2 {params.prefix}_unmapped_2.fq
		pigz {params.prefix}_unmapped_1.fq
		pigz {params.prefix}_unmapped_2.fq
		touch {output.endfile}
		"""

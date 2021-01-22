###############################################################################
## Module: rm_adpt -- AdapterRemoval, star_1st -- STAR 1st pass              ##
## Version 0.1.0                                                             ##
## Written by Sen Li, sen.li_at_sund.ku.dk                                   ##
## Notes: alway use *.fq.gz to save space                                    ##
###############################################################################

###
rule rm_adpt:
	input:
		read1="{WorkDir}/00_RawData/{sample}_1.fq.gz",
		read2="{WorkDir}/00_RawData/{sample}_2.fq.gz"
	output:
		settings="{WorkDir}/01_QualityFiltered/{sample}.settings",
		singleton="{WorkDir}/01_QualityFiltered/{sample}.singleton.gz",
		discarded="{WorkDir}/01_QualityFiltered/{sample}.discarded.gz",
		read1="{WorkDir}/01_QualityFiltered/{sample}_filtered_1.fq.gz",
		read2="{WorkDir}/01_QualityFiltered/{sample}_filtered_2.fq.gz"
	threads: 4
	params:
		adapter1=expand("{adapter1}", adapter1=config['adapter1']),
		adapter2=expand("{adapter2}", adapter2=config['adapter2']),
		maxns=expand("{maxns}", maxns=config['maxns']),
		minquality=expand("{minquality}", minquality=config['minquality']),
		minlength=expand("{minlength}", minlength=config['minlength'])
	shell:
		"""
		module load gcc/8.2.0 AdapterRemoval/2.2.4
		AdapterRemoval --file1 {input.read1} --file2 {input.read2} --settings {output.settings} --singleton {output.singleton} --discarded {output.discarded} --output1 {output.read1} --output2 {output.read2} --trimqualities --trimns --maxns {params.maxns} --minquality {params.minquality} --minlength {params.minlength} --threads {threads} --adapter1 {params.adapter1} --adapter2 {params.adapter2} --gzip
		"""

###
rule star_1st:
	input:
		read1="{WorkDir}/01_QualityFiltered/{sample}_filtered_1.fq.gz",
		read2="{WorkDir}/01_QualityFiltered/{sample}_filtered_2.fq.gz"
	params:
		prefix="{WorkDir}/02_StarMap1st/{sample}",
		StarMinN=expand("{StarMinN}", StarMinN=config['StarMinN']),
		StarMinS=expand("{StarMinS}", StarMinS=config['StarMinS']),
		StarRef=expand("{StarRef}", StarRef=config['StarRef'])
	output:
		endfile="{WorkDir}/02_StarMap1st/{sample}.done"
	threads: 32
	shell:
		"""
		module unload gcc/8.2.0 star/2.7.2b
		module load gcc/8.2.0 star/2.7.2b
		STAR --genomeDir {params.StarRef} --readFilesIn {input.read1} {input.read2} --readFilesCommand zcat --outReadsUnmapped None --outFilterMatchNminOverLread {params.StarMinN} --outFilterScoreMinOverLread {params.StarMinS} --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 16000000000 --runThreadN {threads} --outFileNamePrefix {params.prefix}_
		touch {output.endfile}
		"""
###

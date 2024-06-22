rule all:
	input:
		"results/variantes_anotadas.csv"

rule anotar_variantes:
	input:
		"data/NIST.vcf.gz"
	output:
		"results/variantes_anotadas.csv"
	script:
		"scripts/anotacao_de_variantes2.py"
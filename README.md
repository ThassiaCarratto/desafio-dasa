# desafio-dasa

Esse é um workflow com o snakemake para anotação de variantes a partir de um VCF.

A partir dos rs id (ou cromossomo e posição, caso o id esteja faltando), é feita uma requisição ao REST API do Ensembl.

Com os dados retornados, são buscadas informações do dbSNP id, gene e frequência, e os resultados são anotados em um arquivo CSV. 
A frequência apresentada é referente a populações AMR do 1K Genomes (preferencialmente), ou GnomeAD (caso o primeiro esteja ausente).

Basicamente, o workflow:

Separa as variantes entre aquelas que possuem ou não um rsId;
Para as variantes com rsId, faz uma requisição POST para obter dados de VEP Ensembl;
Para as variantes sem rsId, é realizada uma requisição GET para obter dados do VEP Ensembl a partir do cromossomo e posição;
Confere se o rsId da variante está de acordo com o que consta no dbsnp, e, se necessário, faz a atualização;
Se for uma variante sem rsId, procura se há algum identificador no registro obtido;
Anota o gene em que a variante está localizada;
Verifica se há dados de frequência;
Compara se o alelo para qual a frequência é dada está presente no vcf para aquela variante;
Se sim, anota a frequencia.

Para utilizar esse workflow, recomenda-se criar um novo ambiente conda:
mamba env create --name anotacao_de_variantes --file environment.yaml

seguido de sua ativação:
conda activate anotacao_de_variantes

para rodar o pipeline, tenha certeza de que o vcf input está na pasta data e atualize o nome do arquivo input no Snakefile.

snakemake --cores N

em que N é o número de cores a ser utilizado.

Obs: o arquivo api não está funcionando.

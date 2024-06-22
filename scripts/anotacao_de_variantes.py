import gzip, requests, sys, json, re
import pandas as pd

def solicitar_dados_ensemble(snps_dict):
    # Definindo o servidor
    server = "http://grch37.rest.ensembl.org"
    
    # Construindo a extensão da URL
    ext = "/vep/human/id"
    headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
    
    # Fazendo a solicitação POST para os snps contidos em snps_dict
    r = requests.post(server+ext, headers=headers, json=snps_dict)
    
    # Verificando a resposta
    if r.ok:
        response_data = r.json()
    else:
        response_data = False
    
    return response_data

def elementos_sao_iguais(lista):
    # confere se todos os elementos de uma lista são iguais
    return len(set(lista)) == 1

def atualizar_rs(lista, id):
    # remove o elemento id da lista, deixando só um rs, que é o atualizado do dbsnp
    match_set = set(lista)
    match_set.remove(id)
    match_list = list(match_set)
    return match_list[0]

def encontrar_gene_snps_desconhecidos(chr, pos, i):
    # encontrar o gene de variantes que não foram encontradas no ensembl
    crom = chr[i]
    posicao = pos[i]

    # baseado na posiçao da variante, extraimos informação do gene
    server = "http://grch37.rest.ensembl.org"
    ext = f"/overlap/region/human/{crom}:{posicao}-{posicao}?feature=gene"
    r = requests.get(server + ext, headers={"Content-Type": "application/json"})

    if r.ok:
      dados = r.json()
    else:
      gene = "-"

    dados_string = json.dumps(dados)
    pattern_gene = re.compile(r'"external_name":\s*"([^"]*)"')

    # Encontrar a correspondência na string
    match_gene = pattern_gene.search(dados_string)

    # Verificar se foi encontrada uma correspondência
    if match_gene:
      gene_name = match_gene.group(1)
      gene = gene_name
    else:
      gene = "-"

    return gene

def encontrar_id_dbsnp(variant_string, snps, i, nao_encontrados, variant_ann):
    pattern_id = re.compile(r'"id":\s*"(rs\d+)"')
    match_id = pattern_id.findall(variant_string)
    if snps[i] not in match_id:
      # verifica se é um registro duplicado, com o rs do snp anterior
      if snps[i-1] in match_id:
        return False, i, nao_encontrados, variant_ann
      else:
        # então é um registro de um próximo snp da lista, provavelmente
        while snps[i] not in match_id:
          # anota o snp atual como não encontrado e vai para o proximo snp da lista
          gene = encontrar_gene_snps_desconhecidos(chr, pos, i)
          dbsnp_id = "SNP não encontrado"
          fonte = "-"
          a1 = "-"
          freq1 = "-"
          a2 = "-"
          freq2 = "-"
          variant_ann.append([chr[i], pos[i], snps[i], dbsnp_id, gene, fonte, a1, freq1, a2, freq2, dp[i]])
          nao_encontrados.append(snps[i])
          i = i+1

    if elementos_sao_iguais(match_id):
        # o rs do vcf bate com o rs do dbsnp
        dbsnp_id = match_id[0]
    else:
        # o rs do vcf não bate com o rs do dbsnp e precisa ser atualizado
        dbsnp_id = atualizar_rs(match_id, snps[i])

    return dbsnp_id, i, nao_encontrados, variant_ann

def encontrar_gene(variant_string):
    # encontra a informação de gene nos dados vindos do ensembl
    pattern_gene = re.compile(r'"gene_symbol":\s*"([^"]*)"(?:,|})')
    match_gene = pattern_gene.search(variant_string)

    if match_gene:
        gene = match_gene.group(1)
    else:
        gene = "-"

    return gene

def variante_bialelica(alt, i):
    # verifica se a variante é bialelica
    if "," in alt[i]:
        return False
    else:
        return True

def conferir_alelos(ref, alt, i, frequencies_json):
    # verificar se é uma variante bialelica
    bialelica = variante_bialelica(alt, i)

    if bialelica == True:
        # capturar os alelos no vcf
        alelos_snp = [ref[i], alt[i]]
    else:
        # separa os alelos alternativos e adiciona a uma lista
        alelos = alt[i].split(", ")
        alelos_snp = [ref[i]]
        for a in alelos:
            alelos_snp.append(a)

    # encontrar quais alelos estão nos dados de frequencia do Ensembl
    pattern_alelo = re.compile(r'"([TAGC]+)": \{')
    matches_alelo = pattern_alelo.findall(frequencies_json)

    # conferir se o alelo do ensembl está no vcf passado
    # caso tenha só um alelo nos dados ensembl
    if len(matches_alelo) == 1:
        if matches_alelo[0] in alelos_snp:
            alelo = matches_alelo[0]
        else:
            # o alelo do ensembl é diferente do alelo do vcf
            alelo = False
    # caso tenha mais de um alelo, verificar se algum está no vcf
    else:
        alelo = []
        for a in matches_alelo:
            if a in alelos_snp:
                alelo.append(a)
        # se a lista estiver vazia, nenhum alelo do ensembl está no vcf
        if len(alelo) == 0:
            alelo = False
    return alelo

def extrair_freq(frequency_data):
    # primeiro procurar pela freq no 1k gen
    pattern_freq_amr = re.compile(r'"amr":\s*(-?\d+(\.\d+)?([eE][-+]?\d+)?)')
    match_amr = pattern_freq_amr.search(frequency_data)
    if match_amr:
        freq = match_amr.group()
    else:
        # se não encontrar no 1k gen, procurar no gnomead
        pattern_freq_gnomade_amr = re.compile(r'"gnomade_amr":\s*(-?\d+(\.\d+)?([eE][-+]?\d+)?)')
        match_gnomade_amr = pattern_freq_gnomade_amr.search(frequency_data)
        if match_gnomade_amr:
            freq = match_gnomade_amr.group()
        else:
            freq = False
    return freq

def extrair_string_freq(alelo_freq, frequency_data):
    # dentro da string de frequencia, captura apenas as informações para um alelo específico
    pattern_al = re.compile(rf'"{alelo_freq}": \{{[^}}]+\}}')
    match_al = pattern_al.search(frequency_data)
    resultado = match_al.group()
    return resultado

def gerar_termo_freq(a, frequencies_json):
    # procura dados de frequencia para aquele alelo especifico
    frequencia_string = extrair_string_freq(a, frequencies_json)
    # extrai a frequencia e a fonte para aquele alelo
    termo_freq = extrair_freq(frequencia_string)
    return termo_freq

def definir_fonte_freq(termo_freq):
    chave_freq = termo_freq.split(": ")
    # captura apenas a frequencia
    freq = chave_freq[1]
    # captura apenas o banco de dados, amr é a frequencia do 1k Gen, gnomade_amr é a frequencia do gnomade
    banco = chave_freq[0]
    if banco.strip('"') == "amr":
        fonte = "1k Genomes"
    else:
        fonte = "GNOMEAD"

    return freq, fonte

def encontrar_frequencias(variant_string, ref, alt, i):
    # procura a informação de frequencia nos dados json ensembl
    pattern_freq = re.compile(r'"frequencies": \{.*?\}\}', re.DOTALL)
    match_freq = pattern_freq.search(variant_string)

    # caso encontre informações de frequencia
    if match_freq:
        frequencies_json = match_freq.group()

        # confere se o alelo do ensembl é o mesmo do vcf
        # retorna os alelos encontrados em ambos ou False se os alelos forem diferentes
        alelo_freq = conferir_alelos(ref, alt, i, frequencies_json)
        if alelo_freq == False:
            alelo_freq = "Alelo diferente do REF/ALT"
            freq = "-"
            fonte = "-"
            alelo_freq2 = "-"
            freq2 = "-"
        else:
            # se alelo_freq for uma string, só terá uma unica frequencia a ser extraida
            if isinstance(alelo_freq, str):
                termo_freq = extrair_freq(frequencies_json)
                if termo_freq == False:
                  freq = "Frequencia para amr não encontrada"
                  fonte = "-"
                else:
                  freq, fonte = definir_fonte_freq(termo_freq)
                # esses campos ficam em branco pois não temos um segundo alelo
                alelo_freq2 = "-"
                freq2 = "-"
            else:
                # se for uma lista, é porque havia dados de frequencia para mais de um alelo no Ensembl
                # se for uma variante bialelica, teremos só um alelo na lista
                bialelica = variante_bialelica(alt, i)
                if bialelica == True:
                  alelo_freq = alelo_freq[0]
                  termo_freq = gerar_termo_freq(alelo_freq, frequencies_json)
                  freq, fonte = definir_fonte_freq(termo_freq)
                  # não terá frequencia de um segundo alelo pois é variante bialelica
                  alelo_freq2 = "-"
                  freq2 = "-"
                else:
                  # variante multialelica
                  # pode ser que haja dados de freq para apenas um alelo
                  if len(alelo_freq) == 1:
                    termo_freq = gerar_termo_freq(alelo_freq[0], frequencies_json)
                    freq, fonte = definir_fonte_freq(termo_freq)
                    alelo_freq2 = "-"
                    freq2 = "-"
                  else:
                    # se houver dados para mais de um alelo
                    alelo_freq = alelo_freq[0]
                    alelo_freq2 = alelo_freq[1]
                    freqs = []
                    for i in [0,1]:
                      termo_freq = gerar_termo_freq(alelo_freq[i], frequencies_json)
                      freqs.append(termo_freq[0], termo_freq[1])
                    freq = freqs[0]
                    fonte = freqs[1]
                    freq2 = freqs[2]
    # caso não sejam encontrados dados de frequencia
    else:
        fonte = "-"
        alelo_freq = "-"
        freq = "Não há dados"
        alelo_freq2 = "-"
        freq2 = "-"

    return alelo_freq, freq, fonte, alelo_freq2, freq2

# capturar informações de snps sem rs que foram encontrados no ensembl
def info_snps_semrs_1(dados_semrs, chr_sem_rs, pos_sem_rs, id_sem_rs, ref_sem_rs, alt_sem_rs, dp_sem_rs, i):
    dados_semrs = r.json()
    # transformando json em string para usar regex
    dados_semrs_string = json.dumps(dados_semrs)

    # procurando por um rs no registro
    pattern_id = re.compile(r'"id":\s*"(rs\d+)"')
    match_id = pattern_id.findall(dados_semrs_string)

    # se não encontrar, match_id será vazio
    if len(match_id) == 0:
      dbsnp = "."
    elif len(match_id) == 1:
      dbsnp = match_id[0]
    else:
      dbsnp = match_id    # se encontrar mais de um, serão exibidos os dois

    # encontrar o gene a que essa variante pertence
    gene = encontrar_gene(dados_semrs_string)

    # encontrar a frequencia em ao menos um banco populacional
    a1, freq1, fonte, a2, freq2 = encontrar_frequencias(dados_semrs_string, ref_sem_rs, alt_sem_rs, i)

    return [chr_sem_rs[i], pos_sem_rs[i], id_sem_rs[i], dbsnp, gene, fonte, a1, freq1, a2, freq2, dp_sem_rs[i]]

# capturar informação de gene de snps sem rs que não foram encontrados no ensembl
def info_snps_semrs2(chr_sem_rs, pos_sem_rs, id_sem_rs, dp_sem_rs, i):
    gene = encontrar_gene_snps_desconhecidos(chr_sem_rs, pos_sem_rs, i)
    dbsnp_id = "SNP não encontrado"
    fonte = "-"
    a1 = "-"
    freq1 = "-"
    a2 = "-"
    freq2 = "-"
    return [chr_sem_rs[i], pos_sem_rs[i], id_sem_rs[i], dbsnp_id, gene, fonte, a1, freq1, a2, freq2, dp_sem_rs[i]]

# Função principal para recuperar informações dos snps com rs
def main(response_data, chr, pos, snps, ref, alt, dp):

  i = 0
  nao_encontrados = []
  variant_ann = []

  # itera nos registros obtidos por meio da requisição POST
  for variant in response_data:
      variant_string = json.dumps(variant)

      # encontrar o dbsnp_id das variantes e verificar quais não tiveram registros encontrados
      dbsnp_id, i, nao_encontrados, variant_ann = encontrar_id_dbsnp(variant_string, snps, i, nao_encontrados, variant_ann)

      # se o dbsnp_id for False, é um registro duplicado, seguimos para o próximo
      if dbsnp_id == False:
          continue

      # encontrar o gene a que essa variante pertence
      gene = encontrar_gene(variant_string)

      # encontrar a frequencia em ao menos um banco populacional
      a1, freq1, fonte, a2, freq2 = encontrar_frequencias(variant_string, ref, alt, i)

      variant_ann.append([chr[i], pos[i], snps[i], dbsnp_id, gene, fonte, a1, freq1, a2, freq2, dp[i]])

      i = i + 1
      if i >= 100:
        break

  return variant_ann


#arquivo_vcf = './data/NIST.vcf.gz'

# inicia listas vazias para guardar informacoes das variantes, dividindo por variantes com rs e sem rs (.)
chr_rs = []
pos_rs = []
id_rs = []
ref_rs = []
alt_rs = []
dp_rs = []

chr_sem_rs = []
pos_sem_rs = []
id_sem_rs = []
ref_sem_rs = []
alt_sem_rs = []
dp_sem_rs = []

# abre e le o arquivo vcf, adicionando nas listas as informações das variantes
with gzip.open(snakemake.input[0], 'rt') as v:
    for linha in v:
        if linha.startswith('#'):
            continue
        fields = linha.strip().split('\t')
        id = fields[2]

        # pegar a informação de profundidade no campo INFO do vcf
        info = fields[7]
        info_dict = dict(item.split('=') for item in info.split(';') if '=' in item)
        dp = info_dict.get('DP', 'N/A')

        if id != ".":
            chr_rs.append(fields[0])
            pos_rs.append(fields[1])
            id_rs.append(id)
            ref_rs.append(fields[3])
            alt_rs.append(fields[4])
            dp_rs.append(dp)
        else:
            chr_sem_rs.append(fields[0])
            pos_sem_rs.append(fields[1])
            id_sem_rs.append(id)
            ref_sem_rs.append(fields[3])
            alt_sem_rs.append(fields[4])
            dp_sem_rs.append(dp)

## CÓDIGO PARA EXTRAIR INFORMAÇÕES DOS SNPS QUE POSSUEM RS ##

# numero de snps que serão requisitados por vez (máximo 300)
n_snps = 100

# lista para armazenas as listas com informações de cada variante
registros_snps_rs = []

for i in range(0, len(id_rs), n_snps):
    chr = chr_rs[i:i+n_snps]
    pos = pos_rs[i:i+n_snps]
    snps = id_rs[i:i+n_snps]
    ref = ref_rs[i:i+n_snps]
    alt = alt_rs[i:i+n_snps]
    dp = dp_rs[i:i+n_snps]
    snps_dict = {"ids": snps}

    response_data = solicitar_dados_ensemble(snps_dict)

    # Se a solicitação der errado, pula para o proximo bloco de snps
    if response_data == False:
        print("Erro na solicitação")
        continue

    registros = main(response_data, chr, pos, snps, ref, alt, dp)
    registros_snps_rs.extend(registros)


## CÓDIGO PARA EXTRAIR INFORMAÇÕES DOS SNPS QUE NÃO POSSUEM RS ##
registros_snps_sem_rs = []

for i in range(0, len(pos_sem_rs)):
    chr = chr_sem_rs[i]
    pos = pos_sem_rs[i]
    alt = alt_sem_rs[i]
    tam = len(alt)
    server = "http://grch37.rest.ensembl.org"
    ext = f"/vep/human/region/{chr}:{pos}-{pos}:{tam}/{alt}?"
    
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
    
    if r.ok:
        # se a solicitação retornou um registro, extrairemos as informações desse registro
        dados_semrs = r.json()
        infos = info_snps_semrs_1(dados_semrs, chr_sem_rs, pos_sem_rs, id_sem_rs, ref_sem_rs, alt_sem_rs, dp_sem_rs, i)
        registros_snps_sem_rs.append(infos)
    
    else:
        # snp não encontrado, vão ter só informações do gene
        infos = info_snps_semrs2(chr_sem_rs, pos_sem_rs, id_sem_rs, dp_sem_rs, i)
        registros_snps_sem_rs.append(infos)

# juntando os resultados obtidos para marcadores com e sem rs
dados_combinados = registros_snps_rs + registros_snps_sem_rs

# nome das colunas do dataframe
colunas = ['CHR', 'POS', 'ID', 'dbSNP_ID', 'Gene', 'Fonte', 'A1', 'Freq1-AMR', 'A2', 'Freq2-AMR', 'DP']

# criando o dataframe
df = pd.DataFrame(dados_combinados, columns=colunas)

# convertendo 'CHR' e 'POS' para int para garantir a ordenação correta
df['CHR'] = pd.to_numeric(df['CHR'])
df['POS'] = pd.to_numeric(df['POS'])

# convertendo as colunas de freq para float, convertendo valores inválidos em NaN
df['Freq1-AMR'] = pd.to_numeric(df['Freq1-AMR'], errors='coerce')
df['Freq2-AMR'] = pd.to_numeric(df['Freq2-AMR'], errors='coerce')

# ordenando o DataFrame por 'CHR' e 'POS'
df = df.sort_values(by=['CHR', 'POS'])

# resetando o índice do DataFrame
df = df.reset_index(drop=True)

# salvando o DataFrame em um arquivo CSV
df.to_csv(snakemake.output[0], index=False)
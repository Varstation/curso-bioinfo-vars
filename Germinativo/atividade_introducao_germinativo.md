# Curso de introdução à bioinformática com ênfase em diagnóstico - Prática guiada

Este curso tem como objetivo explorar de forma prática as diversas etapas de um pipeline padrão para análises de dados NGS em rotina de testes genéticos. Curso desenvolvido pela equipe de bioinformática do Hospital Israelita Albert Einstein:
* Deyvid Amgarten, Murilo Cervato, Pedro Sebe e Rodrigo Reis

# Dia 1

## Configurações e acesso aos servidores na nuvem
Todo o curso será realizado em servidores Linux (ubuntu) na AWS. Estes servidores possuem todas as ferramentas utilizadas no curso já instaladas, assim como as configurações necessárias. É importante lembrar que na maioria dos casos do dia-a-dia de um bioinformata, estas configurações e instalações podem ser necessárias e de responsabilidade do bioinformata.  

**Usuário e senha**  
Cada aluno receberá uma mensagem no email cadastrado durante a inscrição com os dados de usuário e senha para acesso aos servidores.
Fique atento para as seguintes informações:  
**User  
Password  
IP  
Port**  

**Procedimento de acesso**
* Fazer download e instalar [MOBA](https://mobaxterm.mobatek.net/download-home-edition.html)
* Clique em "Sessão" -> "Nova Sessão" -> "SSH"
* Especificar os dados de acordo com informações recebidas no email: Remote host (IP), username, port 
* Clique em OK para iniciar a sessão

## Manual dos comando de linux
Forma fácil de obter informações sobre comandos específicos do linux:

```
man ls
man cp
man cat
```

## Listar o diretório atual;
```
pwd
```

## Criar diretório de teste
```
mkdir teste
```
## Mudar para o diretório de teste;
```
cd teste
```

## Listar o diretório atual;
```
pwd
```

## Listar o conteúdo diretório atual;
```
ls -l
touch meuArquivo.txt
ls -l
```

## Movimentar para o seu diretório home;
```
cd ~/
```

## Editor de texto básico de linux;
```
vim arq_texto_exemplo.txt

#<Esc>i  		 = inserir texto antes do cursor, até precionar <Esc> </br>
#<Esc>dd  		 = deletar a linha inteira; </br>
#<Esc>:x<Enter>  = sair do vi e SALVAR as alterações;</br>
#<Esc>:wq<Enter> = sair do vi e SALVAR as alterações;</br>
#<Esc>:q<Enter>	 = sair do vi;</br>
#<Esc>:q!<Enter> = forçar sair do vi SEM SALVAR;</br>
```

## Mostrar o conteúdo de um arquivo;
```
cat arq_texto_exemplo.txt 
```

## Contar o número de linhas de um arquivo;
```
wc -l arq_texto_exemplo.txt
```

## Checar as permissões e coverter os BCLs para FASTQ
Os dados que saem do sequenciador são arquivos de imagens em formato BCL. Eles precisam passar pelo processo de _basecalling_, ou chamada de bases. 
**Este procedimento será realizado apenas pelo palestrante**.
```
ls -l /bioinfo/dados/NextSeq_RUN01
chmod -R 775 /bioinfo/dados/NextSeq_RUN01/Files
cd /bioinfo/dados/NextSeq_RUN01/Files
nohup bcl2fastq --no-lane-splitting --barcode-mismatches 1 1>bcl2fastq.log &
```

# DIA 2

## Criar a estrutura de diretórios para trabalhar;
```
mkdir dados
mkdir dados/fastq
mkdir dados/bwa
mkdir dados/picard
mkdir dados/fastqc
mkdir dados/bedtools
mkdir dados/annovar
mkdir dados/freebayes
mkdir dados/gatk

```

## Criar a estrutura de diretório para guardar seu genoma de referencia;
```
mkdir referencia
mkdir referencia/hg19
```

## Listar o diretório atual;
```
pwd
```

## Copiar os FASTQ para sua pasta de análise;
```
cp /bioinfo/dados/NextSeq_RUN01/Files/Data/Intensities/BaseCalls/AMOSTRA01_S1*.fastq.gz dados/fastq/
```

## Listar os arquivos copiados;
```
ls -lh dados/fastq/*
```

## Executar o FASTQC para avaliar a qualidade das sequencias produzidas;
```
time fastqc -o dados/fastqc dados/fastq/AMOSTRA01_S1_R1_001.fastq.gz dados/fastq/AMOSTRA01_S1_R2_001.fastq.gz
```
Manual do [FastQC](https://dnacore.missouri.edu/PDF/FastQC_Manual.pdf).</br>
Exemplo de resultado [BOM](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/good_sequence_short_fastqc.html) e [RUIM](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/bad_sequence_fastqc.html).</br>


## Remover os reads fora do padrão configurado no sequenciamento 75bp e Q20;
```
time cutadapt --minimum-length 75 --maximum-length 75 \
-o dados/fastq/AMOSTRA01_S1_R1_001_cutadapt.fastq \
-p dados/fastq/AMOSTRA01_S1_R2_001_cutadapt.fastq \
dados/fastq/AMOSTRA01_S1_R1_001.fastq.gz \
dados/fastq/AMOSTRA01_S1_R2_001.fastq.gz 
``` 
## Executar o FASTQC para avaliar a qualidade das sequencias produzidas após o cutadapt;
```
time fastqc -o dados/fastqc dados/fastq/AMOSTRA01_S1_R1_001_cutadapt.fastq dados/fastq/AMOSTRA01_S1_R2_001_cutadapt.fastq
```

## Fazer download dos HTMLs gerados com o FastQC e comparar os dois, antes e depois do cutadapt

## Fazer download de um cromossomo para utilizar como referencia
```
# Mover para o diretório do seu genoma de referência
cd referencia/hg19/
pwd
wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr13.fa.gz
gunzip chr13.fa.gz
```

## Criar o índice do BWA
```
# reference.fa​ = chr13.fa

bwa index -a bwtsw reference.fa​ 
```

## Gerar o índice do FASTA (genoma de referência)
```
# reference.fa​  = chr13.fa
samtools faidx reference.fa
```

## Gerar o dicionário das sequências FASTA
```
# reference.fa​  = chr13.fa
# reference.dict = chr13.dict

java -jar /bioinfo/app/picard/picard.jar CreateSequenceDictionary \
REFERENCE=reference.fa \
OUTPUT=reference.dict
```

## Mapear os FASTQ limpos contra o hg19;
```
# Voltar para o HOME
cd ~/
NOME=NOME;
Biblioteca=Biblioteca;
Plataforma=Plataforma;

time bwa mem -M -R "@RG\tID:CAP\tSM:$NOME\tLB:$Biblioteca\tPL:$Plataforma" \
/bioinfo/referencia/hg19/chr1_13_17.fa \
dados/fastq/AMOSTRA01_S1_R1_001_cutadapt.fastq \
dados/fastq/AMOSTRA01_S1_R2_001_cutadapt.fastq >dados/bwa/AMOSTRA01_S1.sam
```

## Utilizar o samtools: fixmate, sort e index
```
time samtools fixmate dados/bwa/AMOSTRA01_S1.sam dados/bwa/AMOSTRA01_S1.bam
time samtools sort -O bam -o dados/bwa/AMOSTRA01_S1_sorted.bam dados/bwa/AMOSTRA01_S1.bam
time samtools index dados/bwa/AMOSTRA01_S1_sorted.bam
```

## Visualizar o BAM com o samtools;
```
time samtools view -H dados/bwa/AMOSTRA01_S1_sorted.bam
time samtools view dados/bwa/AMOSTRA01_S1_sorted.bam
```

## Converter BAM to BED para utilizarmos o BED para analise de cobertura;
```
bamToBed -i dados/bwa/AMOSTRA01_S1_sorted.bam >dados/bedtools/AMOSTRA01_S1_sorted.bed
mergeBed -i dados/bedtools/AMOSTRA01_S1_sorted.bed >dados/bedtools/AMOSTRA01_S1_merged.bed
sortBed -i dados/bedtools/AMOSTRA01_S1_merged.bed >dados/bedtools/AMOSTRA01_S1_merged_sorted.bed
```

## Gerar a cobertura média para cada target do BED
```
coverageBed -a dados/bedtools/AMOSTRA01_S1_merged_sorted.bed \
-b dados/bwa/AMOSTRA01_S1_sorted.bam -mean \
>dados/bedtools/AMOSTRA01_S1_coverageBed.bed
```

## Listar as regiões com cobertura média maior que 30x
```
cat dados/bedtools/AMOSTRA01_S1_coverageBed.bed | \
awk -F "\t" '{if($4>30){print}}' \
> dados/bedtools/AMOSTRA01_S1_coverageBed30x.bed
```

## Chamada de variantes com o Freebayes;

```
# -F --min-alternate-fraction N
#      Require at least this fraction of observations supporting
#      an alternate allele within a single individual in the
#      in order to evaluate the position.  default: 0.05
# -C --min-alternate-count N
#      Require at least this count of observations supporting
#      an alternate allele within a single individual in order
#      to evaluate the position.  default: 2

time /bioinfo/app/freebayes/bin/freebayes -f /bioinfo/referencia/hg19/chr1_13_17.fa \
-F 0.3 -C 15 \
--pooled-continuous dados/bwa/AMOSTRA01_S1_sorted.bam \
>dados/freebayes/AMOSTRA01_S1_sorted.vcf
```

## Chamada de variantes com o GATK;
```
time /bioinfo/app/gatk/gatk-4.1.2.0/gatk HaplotypeCaller -R /bioinfo/referencia/hg19/chr1_13_17.fa \
-I dados/bwa/AMOSTRA01_S1_sorted.bam \
-O dados/gatk/AMOSTRA01_S1_sorted.vcf
```

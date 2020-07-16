# Curso de introdução à bioinformática com ênfase em diagnóstico - Prática guiada

Este curso tem como objetivo explorar de forma prática as diversas etapas de um pipeline padrão para análises de dados NGS em rotina de testes genéticos. Curso desenvolvido pela equipe de bioinformática do Hospital Israelita Albert Einstein:
* Deyvid Amgarten, Murilo Cervato, Pedro Sebe e Rodrigo Reis

# Atividade Pipeline de análises de doenças infecciosas

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


# DIA 2

## Criar a estrutura de diretórios para trabalhar;
```
mkdir viroma
mkdir viroma/fastq
mkdir viroma/bwa
mkdir viroma/fastqc
mkdir viroma/karken2
mkdir reference
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
time fastqc -o viroma/fastqc viroma/fastq/AMOSTRA_JOAQUIM_S17_R1_001.fastq.gz dados/fastq/AMOSTRA_JOAQUIM_S17_R2_001.fastq.gz
```
Manual do [FastQC](https://dnacore.missouri.edu/PDF/FastQC_Manual.pdf).</br>
Exemplo de resultado [BOM](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/good_sequence_short_fastqc.html) e [RUIM](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/bad_sequence_fastqc.html).</br>


## Remover os reads fora do padrão configurado no sequenciamento 75bp e Q20;
```
time cutadapt --minimum-length 75 --maximum-length 200 \

-o viroma/fastq/AMOSTRA_JOAQUIM_S17_R1_001_cutadapt.fastq \
-p viroma/fastq/AMOSTRA_JOAQUIM_S17_R2_001_cutadapt.fastq \
viroma/fastq/AMOSTRA_JOAQUIM_S17_R1_001.fastq.gz \
viroma/fastq/AMOSTRA_JOAQUIM_S17_R1_001.fastq.gz 
``` 
## Executar o FASTQC para avaliar a qualidade das sequencias produzidas após o cutadapt;
```
time fastqc -o dados/fastqc dados/fastq/AMOSTRA01_S1_R1_001_cutadapt.fastq dados/fastq/AMOSTRA01_S1_R2_001_cutadapt.fastq
```

## Fazer download dos HTMLs gerados com o FastQC e comparar os dois, antes e depois do cutadapt

## Fazer download de um cromossomo Humano para utilizar como referencia na remoção de contaminantes
```
# Mover para o diretório do seu genoma de referência
cd referencia/
pwd
cp /dados/... .
```

## Criar o índice do BWA
```
bwa index -a bwtsw chr13.fa
```

## Gerar o índice do FASTA (genoma de referência)
```
# reference.fa​  = chr13.fa
samtools faidx reference.fa
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
viroma/fastq/AMOSTRA_JOAQUIM_S17_R1_001_cutadapt.fastq \
viroma/fastq/AMOSTRA_JOAQUIM_S17_R2_001_cutadapt.fastq >viroma/bwa/AMOSTRA_JOAQUIM_S17.sam
```

## Utilizar o samtools: fixmate, sort e index
```
time samtools fixmate viroma/bwa/AMOSTRA_JOAQUIM_S17.sam viroma/bwa/AMOSTRA_JOAQUIM_S17.bam
time samtools sort -O bam -o viroma/bwa/AMOSTRA_JOAQUIM_S17_sorted.bam viroma/bwa/AMOSTRA_JOAQUIM_S17.bam
time samtools index viroma/bwa/AMOSTRA_JOAQUIM_S17_sorted.bam
```

## Visualizar o BAM com o samtools;
```
time samtools view -H viroma/bwa/AMOSTRA_JOAQUIM_S17_sorted.bam
time samtools view viroma/bwa/AMOSTRA_JOAQUIM_S17_sorted.bam
```

## Gerar FASTQ de reads não mapeados no genoma humano

```

```

## Identificação de patógenos com o Kraken2
```

```

## Geração do relatório de diversidade 
```

```

## Baixar o resultado e ver qual é o patógeno responsável pela patologia do paciente
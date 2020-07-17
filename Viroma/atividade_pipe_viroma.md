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

## Definir qual amostra você vai analisar
```
SAMPLE="AMOSTRA_JOAQUIM_S17"
#SAMPLE="AMOSTRA_MARIA_S19"
#SAMPLE="AMOSTRA_PEDRO_S10"
#SAMPLE="AMOSTRA_MARCELA_S12"
```


## Criar a estrutura de diretórios para trabalhar;
```
mkdir viroma
mkdir viroma/fastq
mkdir viroma/cutadapt
mkdir viroma/bwa
mkdir viroma/fastqc
mkdir viroma/kraken2
mkdir -p viroma/reference/hg19
```

## Listar o diretório atual;
```
pwd
```

## Copiar os FASTQ para sua pasta de análise;
```
rsync --progress /data/FASTQ_VIROMA/${SAMPLE}* viroma/fastq/
```

## Listar os arquivos copiados;
```
ls -lh viroma/fastq/
```

## Executar o FASTQC para avaliar a qualidade das sequencias produzidas;
```
fastqc -o viroma/fastqc viroma/fastq/${SAMPLE}_R1_001.fastq.gz viroma/fastq/${SAMPLE}_R2_001.fastq.gz
```
Manual do [FastQC](https://dnacore.missouri.edu/PDF/FastQC_Manual.pdf).</br>
Exemplo de resultado [BOM](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/good_sequence_short_fastqc.html) e [RUIM](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/bad_sequence_fastqc.html).</br>


## Remover os reads fora do padrão configurado no sequenciamento 75bp e Q20;
```
cutadapt --minimum-length 75 --maximum-length 200 \
-q 20,20 --max-n 0 -u 9 -u -5 -U 9 -U -5 \
-o viroma/cutadapt/${SAMPLE}_R1_001_cutadapt.fastq \
-p viroma/cutadapt/${SAMPLE}_R2_001_cutadapt.fastq \
viroma/fastq/${SAMPLE}_R1_001.fastq.gz \
viroma/fastq/${SAMPLE}_R2_001.fastq.gz 
``` 
## Executar o FASTQC para avaliar a qualidade das sequencias produzidas após o cutadapt;
```
fastqc -o viroma/fastqc viroma/cutadapt/${SAMPLE}_R1_001_cutadapt.fastq viroma/cutadapt/${SAMPLE}_R2_001_cutadapt.fastq
```

## Fazer download dos HTMLs gerados com o FastQC e comparar os dois, antes e depois do cutadapt

## Fazer download de um cromossomo Humano para utilizar como referencia na remoção de contaminantes
```
# Mover para o diretório do seu genoma de referência
cd viroma/reference/hg19
pwd
rsync --progress /data/chr13.fa .
```

## Criar o índice do BWA
```
bwa index -a bwtsw chr13.fa
```

## Gerar o índice do FASTA (genoma de referência)
```
samtools faidx chr13.fa
```

## Mapear os FASTQ limpos contra o hg19;
```
# Voltar para o HOME
cd ~/
NOME=NOME;
Biblioteca=Biblioteca;
Plataforma=Plataforma;

bwa mem -M -R "@RG\tID:CAP\tSM:$NOME\tLB:$Biblioteca\tPL:$Plataforma" \
viroma/reference/hg19/chr13.fa \
viroma/cutadapt/${SAMPLE}_R1_001_cutadapt.fastq \
viroma/cutadapt/${SAMPLE}_R2_001_cutadapt.fastq >viroma/bwa/${SAMPLE}.sam
```

## Utilizar o samtools: fixmate, sort e index
```
samtools fixmate viroma/bwa/${SAMPLE}.sam viroma/bwa/${SAMPLE}.bam
samtools sort -O bam -o viroma/bwa/${SAMPLE}_sorted.bam viroma/bwa/${SAMPLE}.bam
samtools index viroma/bwa/${SAMPLE}_sorted.bam
```

## Visualizar o BAM com o samtools;
```
samtools view -H viroma/bwa/${SAMPLE}_sorted.bam
samtools view viroma/bwa/${SAMPLE}_sorted.bam
```

## Gerar FASTQ de reads não mapeados no genoma humano

```
cd ~/
samtools view -u -f 12 -b viroma/bwa/${SAMPLE}_sorted.bam > viroma/bwa/${SAMPLE}_sorted_unmapped.bam
bamToFastq -i viroma/bwa/${SAMPLE}_sorted_unmapped.bam -fq viroma/bwa/${SAMPLE}_R1_unmapped_human.fastq -fq2 viroma/bwa/${SAMPLE}_R2_unmapped_human.fastq
```

## Identificação de patógenos com o Kraken2
```
kraken2 -db /data/minikraken2_v2_8GB_201904_UPDATE --threads 1 --minimum-base-quality 20 --report viroma/kraken2/${SAMPLE}_kraken_report.txt\
--paired viroma/bwa/${SAMPLE}_R1_unmapped_human.fastq viroma/bwa/${SAMPLE}_R2_unmapped_human.fastq > viroma/kraken2/${SAMPLE}_kraken2_NT.out
```

## Geração do relatório de diversidade 
```
awk '{print $2"\t"$3}' viroma/kraken2/${SAMPLE}_kraken2_NT.out > viroma/kraken2/${SAMPLE}_kraken2_NT_2krona.tab
ktImportTaxonomy viroma/kraken2/${SAMPLE}_kraken2_NT_2krona.tab -o viroma/kraken2/${SAMPLE}_classification.html
```

## Baixar o resultado e ver qual é o patógeno responsável pela patologia do paciente

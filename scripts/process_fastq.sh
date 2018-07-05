mkdir cutadapt
for f in *.fastq; do cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC --match-read-wildcards -m 10 -o cutadapt/$f $f -j 4 ; done

cd cutadapt
mkdir umitools

for f in *.fastq; do umi_tools extract --bc-pattern NNNNNNNN --stdin $f --stdout umitools/$f -L logfile_$f.txt; done

cd umitools
mkdir fastqc

for f in *.fastq; do fastqc $f -o ./fastqc ; done

mkdir star
for f in *.fastq; do sudo STAR --runThreadN 4 --readFilesIn $f --genomeDir ../../../../R64e/genomedir --outSAMtype BAM SortedByCoordinate  --outFilterMismatchNoverReadLmax 0.04 --outFileNamePrefix star/$f --outFilterMultimapNmax 1 --alignIntronMax 1000; done

cd star
for f in *.bam; do samtools index $f; done

mkdir dedup
for f in *.bam; do umi_tools dedup -I $f -S dedup/$f -L $f_log.txt --method="cluster"; done

cd dedup
for f in *; do samtools index $f; done











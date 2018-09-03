# Didymella
Commands used in the assembly and analysis of Didymella (Phoma) genomes

This document details the commands used to assemble and annotate minion sequence data.

Note - all this work was performed in the directory:
/home/groups/harrisonlab/project_files/didymella


# 0. Building of directory structure

### Minion Data


Basecalling from the gridion was used:

```bash
ProjectDir=/home/groups/harrisonlab/project_files/didymella
mkdir -p $ProjectDir
cd $ProjectDir


Organism=D.pinodella
Strain=61B
Date=2018-05-15
RawDatDir=/data/seq_data/minion/2018/61B-20180515
OutDir=$ProjectDir/raw_dna/minion/$Organism/$Strain
mkdir -p $OutDir

cat $RawDatDir/GA30000/*.fastq | gzip -cf > $OutDir/${Strain}_${Date}_albacore_v2.2.7.fastq.gz
```

### MiSeq data

```bash
  RawDatDir=/data/seq_data/miseq/2018/RAW/180716_M04465_0083_000000000-BM8CL/Data/Intensities/BaseCalls
  ProjectDir=/home/groups/harrisonlab/project_files/didymella
  OutDir=$ProjectDir/raw_dna/paired/D.pinodella/61B
  mkdir -p $OutDir/F
  mkdir -p $OutDir/R
  cd $OutDir/F
  cp -s $RawDatDir/61B_S3_L001_R1_001.fastq.gz .
  cd $OutDir/R
  cp -s $RawDatDir/61B_S3_L001_R2_001.fastq.gz .
	cd $ProjectDir
```



## Assembly

### Removal of adapters

Splitting reads and trimming adapters using porechop
```bash
for RawReads in $(ls raw_dna/minion/*/*/*.fastq.gz); do
Organism=$(echo $RawReads| rev | cut -f3 -d '/' | rev)
Strain=$(echo $RawReads | rev | cut -f2 -d '/' | rev)
echo "$Organism - $Strain"
OutDir=qc_dna/minion/$Organism/$Strain
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
qsub $ProgDir/sub_porechop.sh $RawReads $OutDir
done
```

#### QC of MiSeq data

programs:
  fastqc
  fastq-mcf
  kmc

Data quality was visualised using fastqc:
```bash
  for RawData in $(ls raw_dna/paired/*/*/*/*.fastq.gz); do
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc
    echo $RawData;
    qsub $ProgDir/run_fastqc.sh $RawData
  done
```

```bash
  for StrainPath in $(ls -d raw_dna/paired/*/*); do
    ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc
    IlluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/ncbi_adapters.fa
    ReadsF=$(ls $StrainPath/F/*.fastq*)
    ReadsR=$(ls $StrainPath/R/*.fastq*)
    echo $ReadsF
    echo $ReadsR
    qsub $ProgDir/rna_qc_fastq-mcf.sh $ReadsF $ReadsR $IlluminaAdapters DNA
  done
```

# Identify sequencing coverage

For Minion data:
```bash
for RawData in $(ls qc_dna/minion/*/*/*q.gz); do
echo $RawData;
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
GenomeSz=60
OutDir=$(dirname $RawData)
qsub $ProgDir/sub_count_nuc.sh $GenomeSz $RawData $OutDir
done
```

```bash
  for StrainDir in $(ls -d qc_dna/minion/*/* ); do
    Strain=$(basename $StrainDir)
    printf "$Strain\t"
    for File in $(ls $StrainDir/*.txt); do
      echo $(basename $File);
      cat $File | tail -n1 | rev | cut -f2 -d ' ' | rev;
    done | grep -v '.txt' | awk '{ SUM += $1} END { print SUM }'
  done
```
MinION coverage was:
```
61B	140.02
```

For Miseq data:
```bash
	for RawData in $(ls qc_dna/paired/*/*/*/*q.gz | grep 'lactucae'); do
		echo $RawData;
		ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/dna_qc;
		qsub $ProgDir/run_fastqc.sh $RawData;
		GenomeSz=60
		OutDir=$(dirname $RawData)
		qsub $ProgDir/sub_count_nuc.sh $GenomeSz $RawData $OutDir
	done
```

```bash
	for StrainDir in $(ls -d qc_dna/paired/*/* | grep 'lactucae'); do
		Strain=$(basename $StrainDir)
		printf "$Strain\t"
		for File in $(ls $StrainDir/*/*.txt); do
			echo $(basename $File);
			cat $File | tail -n1 | rev | cut -f2 -d ' ' | rev;
		done | grep -v '.txt' | awk '{ SUM += $1} END { print SUM }'
	done
```

Miseq coverage was:
```

```


### Read correction using Canu

```bash
for TrimReads in $(ls qc_dna/minion/*/*/*q.gz); do
Organism=$(echo $TrimReads | rev | cut -f3 -d '/' | rev)
Strain=$(echo $TrimReads | rev | cut -f2 -d '/' | rev)
OutDir=assembly/canu-1.6/$Organism/"$Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/canu
qsub $ProgDir/sub_canu_correction.sh $TrimReads 35m $Strain $OutDir
done
```


### Assembbly using SMARTdenovo

```bash
for CorrectedReads in $(ls assembly/canu-1.6/*/*/*.trimmedReads.fasta.gz); do
Organism=$(echo $CorrectedReads | rev | cut -f3 -d '/' | rev)
Strain=$(echo $CorrectedReads | rev | cut -f2 -d '/' | rev)
Prefix="$Strain"_smartdenovo
OutDir=assembly/SMARTdenovo/$Organism/"$Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/SMARTdenovo
qsub $ProgDir/sub_SMARTdenovo.sh $CorrectedReads $Prefix $OutDir
done
```


Quast and busco were run to assess the effects of racon on assembly quality:

```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/*.dmo.lay.utg); do
  Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
  Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)  
	echo "$Organism - $Strain"
  OutDir=$(dirname $Assembly)
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
  qsub $ProgDir/sub_quast.sh $Assembly $OutDir
	OutDir=gene_pred/busco/$Organism/$Strain/assembly
	BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/ascomycota_odb9)
	ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
	qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```


Error correction using racon:

```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/*.dmo.lay.utg); do
Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
ReadsFq=$(ls qc_dna/minion/*/$Strain/*q.gz)
Iterations=10
OutDir=$(dirname $Assembly)"/racon_$Iterations"
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/racon
qsub $ProgDir/sub_racon.sh $Assembly $ReadsFq $Iterations $OutDir
done
```

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/SMARTdenovo/*/*/racon_10/*.fasta | grep 'round_10'); do
OutDir=$(dirname $Assembly)
echo "" > tmp.txt
ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
$ProgDir/remove_contaminants.py --keep_mitochondria --inp $Assembly --out $OutDir/racon_min_500bp_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
done
```

Quast and busco were run to assess the effects of racon on assembly quality:

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/SMARTdenovo/*/*/racon_10/racon_min_500bp_renamed.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```

```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/racon*/*.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/ascomycota_odb9)
OutDir=gene_pred/busco/$Organism/$Strain/assembly
# OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```
```bash
printf "Filename\tComplete\tDuplicated\tFragmented\tMissing\tTotal\n"
for File in $(ls gene_pred/busco/*/*/assembly/*/short_summary_*.txt); do
FileName=$(basename $File)
Complete=$(cat $File | grep "(C)" | cut -f2)
Duplicated=$(cat $File | grep "(D)" | cut -f2)
Fragmented=$(cat $File | grep "(F)" | cut -f2)
Missing=$(cat $File | grep "(M)" | cut -f2)
Total=$(cat $File | grep "Total" | cut -f2)
printf "$FileName\t$Complete\t$Duplicated\t$Fragmented\t$Missing\t$Total\n"
done
```



# Nanopolish


```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/racon_10/racon_min_500bp_renamed.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
# Step 1 extract reads as a .fq file which contain info on the location of the fast5 files
# Note - the full path from home must be used
ReadDir=raw_dna/nanopolish/$Organism/$Strain
mkdir -p $ReadDir
ReadsFq=$(ls raw_dna/minion/*/$Strain/*.fastq.gz)
CurDir=$PWD
Fast5Dir=$(ls -d /data/seq_data/minion/2018/61B-20180515/GA30000/reads/)
nanopolish index -d $Fast5Dir $ReadsFq
done

for Assembly in $(ls assembly/SMARTdenovo/*/*/racon_10/racon_min_500bp_renamed.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
# Step 1 extract reads as a .fq file which contain info on the location of the fast5 files
# Note - the full path from home must be used
# ReadDir=raw_dna/nanopolish/$Organism/$Strain
# mkdir -p $ReadDir
ReadsFq=$(ls raw_dna/minion/*/$Strain/*.fastq.gz)
OutDir=$(dirname $Assembly)/nanopolish
mkdir -p $OutDir
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/nanopolish
# submit alignments for nanoppolish
qsub $ProgDir/sub_minimap2_nanopolish.sh $Assembly $ReadsFq $OutDir/nanopolish
done
```


 Split the assembly into 50Kb fragments an submit each to the cluster for
 nanopolish correction

```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/racon_10/racon_min_500bp_renamed.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
OutDir=$(dirname $Assembly)/nanopolish
# ReadsFq=$(ls raw_dna/minion/*/$Strain/*.fastq.gz | grep '2017-12-03')
ReadsFq=$(ls raw_dna/minion/*/$Strain/*.fastq.gz)
AlignedReads=$(ls $OutDir/nanopolish/reads.sorted.bam)

NanoPolishDir=/home/armita/prog/nanopolish/nanopolish/scripts
python $NanoPolishDir/nanopolish_makerange.py $Assembly --segment-length 50000 > $OutDir/nanopolish_range.txt

Ploidy=1
echo "nanopolish log:" > $OutDir/nanopolish_log.txt
ls -lh $OutDir/*/*.fa | grep -v ' 0 ' | cut -f8 -d '/' | sed 's/_consensus.fa//g' > $OutDir/files_present.txt
for Region in $(cat $OutDir/nanopolish_range.txt | grep -vwf "$OutDir/files_present.txt"); do
# for Region in $(cat $OutDir/nanopolish_range.txt | grep -e 'contig_6:300000-350200' -e 'contig_6:350000-400200'); do
Jobs=$(qstat | grep 'sub_nanopo' | grep 'qw' | grep -v -e '863720' -e '863721' | wc -l)
while [ $Jobs -gt 1 ]; do
sleep 1m
printf "."
Jobs=$(qstat | grep 'sub_nanopo' | grep 'qw' | grep -v -e '863720' -e '863721' | wc -l)
done		
printf "\n"
echo $Region
echo $Region >> $OutDir/nanopolish_log.txt
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/nanopolish
qsub $ProgDir/sub_nanopolish_variants.sh $Assembly $ReadsFq $AlignedReads $Ploidy $Region $OutDir/$Region
done
done
```

A subset of nanopolish jobs needed to be resubmitted as the ran out of RAM

```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/racon_10/racon_min_500bp_renamed.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
OutDir=$(dirname $Assembly)/nanopolish
# ReadsFq=$(ls raw_dna/minion/*/$Strain/*.fastq.gz | grep '2017-12-03')
ReadsFq=$(ls raw_dna/minion/*/$Strain/*.fastq.gz)
AlignedReads=$(ls $OutDir/nanopolish/reads.sorted.bam)

NanoPolishDir=/home/armita/prog/nanopolish/nanopolish/scripts
# python $NanoPolishDir/nanopolish_makerange.py $Assembly --segment-length 50000 > $OutDir/nanopolish_range.txt

Ploidy=1
echo "nanopolish log:" > $OutDir/nanopolish_high_mem_log.txt
ls -lh $OutDir/*/*.fa | grep -v ' 0 ' | cut -f8 -d '/' | sed 's/_consensus.fa//g' > $OutDir/files_present.txt
for Region in $(cat $OutDir/nanopolish_range.txt | grep -vwf "$OutDir/files_present.txt"); do
echo $Region
echo $Region >> $OutDir/nanopolish_high_mem_log.txt
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/nanopolish
qsub $ProgDir/sub_nanopolish_variants_high_mem.sh $Assembly $ReadsFq $AlignedReads $Ploidy $Region $OutDir/$Region
done
done

```

```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/racon_10/racon_min_500bp_renamed.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
OutDir=assembly/SMARTdenovo/$Organism/$Strain/nanopolish
mkdir -p $OutDir
# cat "" > $OutDir/"$Strain"_nanoplish.fa
NanoPolishDir=/home/armita/prog/nanopolish/nanopolish/scripts
InDir=$(dirname $Assembly)
python $NanoPolishDir/nanopolish_merge.py $InDir/nanopolish/*/*.fa > $OutDir/"$Strain"_nanoplish.fa

echo "" > tmp.txt
ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
$ProgDir/remove_contaminants.py --keep_mitochondria --inp $OutDir/"$Strain"_nanoplish.fa --out $OutDir/"$Strain"_nanoplish_min_500bp_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
done
```

Quast and busco were run to assess the effects of nanopolish on assembly quality:

```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/nanopolish/*_nanoplish_min_500bp_renamed.fasta); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
# Quast
OutDir=$(dirname $Assembly)
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
# Busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/ascomycota_odb9)
OutDir=gene_pred/busco/$Organism/$Strain/assembly
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```

```bash
  for File in $(ls gene_pred/busco/*/*/assembly/*/short_summary_*.txt); do
  Strain=$(echo $File| rev | cut -d '/' -f4 | rev)
  Organism=$(echo $File | rev | cut -d '/' -f5 | rev)
  Complete=$(cat $File | grep "(C)" | cut -f2)
  Single=$(cat $File | grep "(S)" | cut -f2)
  Fragmented=$(cat $File | grep "(F)" | cut -f2)
  Missing=$(cat $File | grep "(M)" | cut -f2)
  Total=$(cat $File | grep "Total" | cut -f2)
  echo -e "$Organism\t$Strain\t$Complete\t$Single\t$Fragmented\t$Missing\t$Total"
  done
```



### Pilon error correction


Assemblies were polished using Pilon

```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/nanopolish/*_nanoplish_min_500bp_renamed.fasta); do
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
IlluminaDir=$(ls -d qc_dna/paired/*/$Strain)
TrimF1_Read=$(ls $IlluminaDir/F/*_trim.fq.gz | head -n1)
TrimR1_Read=$(ls $IlluminaDir/R/*_trim.fq.gz | head -n1)
OutDir=$(dirname $Assembly)/../pilon
Iterations=5
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/pilon
qsub $ProgDir/sub_pilon.sh $Assembly $TrimF1_Read $TrimR1_Read $OutDir $Iterations
done
```

Contigs were renamed
```bash
echo "" > tmp.txt
Assembly=$(ls assembly/SMARTdenovo/*/*/pilon/*.fasta | grep 'pilon_5')
OutDir=$(dirname $Assembly)
ProgDir=~/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/remove_contaminants
$ProgDir/remove_contaminants.py --keep_mitochondria --inp $Assembly --out $OutDir/pilon_min_500bp_renamed.fasta --coord_file tmp.txt > $OutDir/log.txt
```

Quast and busco were run to assess the effects of pilon on assembly quality:

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
for Assembly in $(ls assembly/SMARTdenovo/*/*/pilon/*.fasta | grep 'pilon_min_500bp_renamed.fasta'); do
  Strain=$(echo $Assembly | rev | cut -f2 -d '/' | rev)
  Organism=$(echo $Assembly | rev | cut -f3 -d '/' | rev)  
  OutDir=$(dirname $Assembly)
  qsub $ProgDir/sub_quast.sh $Assembly $OutDir
done
```


```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/pilon/*.fasta ); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
echo "$Organism - $Strain"
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/ascomycota_odb9)
OutDir=gene_pred/busco/$Organism/$Strain/assembly
# OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```


```bash
printf "Filename\tComplete\tDuplicated\tFragmented\tMissing\tTotal\n"
for File in $(ls gene_pred/busco/*/*/assembly/*/short_summary_*.txt); do  
FileName=$(basename $File)
Complete=$(cat $File | grep "(C)" | cut -f2)
Duplicated=$(cat $File | grep "(D)" | cut -f2)
Fragmented=$(cat $File | grep "(F)" | cut -f2)
Missing=$(cat $File | grep "(M)" | cut -f2)
Total=$(cat $File | grep "Total" | cut -f2)
printf "$FileName\t$Complete\t$Duplicated\t$Fragmented\t$Missing\t$Total\n"
done
```

```
short_summary_3018_contigs_unmasked.txt 1292    2       11      12      1315
short_summary_IPT5_contigs_unmasked.txt 1302    3       4       9       1315
short_summary_ITCC_4638_contigs_unmasked.txt    1298    4       1       16      1315
short_summary_JCM_15942_contigs_unmasked.txt    1295    1       7       13      1315
short_summary_CBS_183_55_contigs_unmasked.txt   1290    1       10      15      1315

short_summary_61B_smartdenovo.dmo.lay.txt       1296    0       493     1936    3725

short_summary_61B_smartdenovo_racon_round_1.txt 1114    2       99      102     1315
short_summary_61B_smartdenovo_racon_round_2.txt 1111    3       111     93      1315
short_summary_61B_smartdenovo_racon_round_3.txt 1125    3       100     90      1315
short_summary_61B_smartdenovo_racon_round_4.txt 1134    2       99      82      1315
short_summary_61B_smartdenovo_racon_round_5.txt 1120    3       105     90      1315
short_summary_61B_smartdenovo_racon_round_6.txt 1130    4       95      90      1315
short_summary_61B_smartdenovo_racon_round_7.txt 1117    2       115     83      1315
short_summary_61B_smartdenovo_racon_round_8.txt 1131    3       92      92      1315
short_summary_61B_smartdenovo_racon_round_9.txt 1140    2       90      85      1315
short_summary_61B_smartdenovo_racon_round_10.txt        1125    3       94      96      1315
short_summary_61B_contigs_unmasked.txt  1297    1       7       11      1315
short_summary_61B_nanoplish_min_500bp_renamed.txt       1249    1       29      37      1315

short_summary_pilon_1.txt       1292    1       8       15      1315
short_summary_pilon_2.txt       1296    1       7       12      1315
short_summary_pilon_3.txt       1297    1       7       11      1315
short_summary_pilon_4.txt       1297    1       7       11      1315
short_summary_pilon_5.txt       1297    1       7       11      1315
short_summary_pilon_min_500bp_renamed.txt       1297    1       7       11      1315

```

# Repeat Masking

Repeat masking was performed on the non-hybrid assembly.

```bash
for Assembly in $(ls assembly/SMARTdenovo/*/*/pilon/pilon_min_500bp_renamed.fasta); do
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)  
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
OutDir=repeat_masked/$Organism/"$Strain"/filtered_contigs
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/repeat_masking
qsub $ProgDir/rep_modeling.sh $Assembly $OutDir
qsub $ProgDir/transposonPSI.sh $Assembly $OutDir
done
```

The TransposonPSI masked bases were used to mask additional bases from the
repeatmasker / repeatmodeller softmasked and hardmasked files.

```bash
for File in $(ls repeat_masked/*/*/*/*_contigs_softmasked.fa | grep 'Stocks4'); do
OutDir=$(dirname $File)
TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
OutFile=$(echo $File | sed 's/_contigs_softmasked.fa/_contigs_softmasked_repeatmasker_TPSI_appended.fa/g')
echo "$OutFile"
bedtools maskfasta -soft -fi $File -bed $TPSI -fo $OutFile
echo "Number of masked bases:"
cat $OutFile | grep -v '>' | tr -d '\n' | awk '{print $0, gsub("[a-z]", ".")}' | cut -f2 -d ' '
done
# The number of N's in hardmasked sequence are not counted as some may be present within the assembly and were therefore not repeatmasked.
for File in $(ls repeat_masked/*/*/*/*_contigs_hardmasked.fa | grep 'Stocks4'); do
OutDir=$(dirname $File)
TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
OutFile=$(echo $File | sed 's/_contigs_hardmasked.fa/_contigs_hardmasked_repeatmasker_TPSI_appended.fa/g')
echo "$OutFile"
bedtools maskfasta -fi $File -bed $TPSI -fo $OutFile
done
```


Quast and BUSCO

```bash
for Assembly in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_unmasked.fa); do
Strain=$(echo $Assembly | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Assembly | rev | cut -f4 -d '/' | rev)
OutDir=$(dirname $Assembly)
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/assemblers/assembly_qc/quast
qsub $ProgDir/sub_quast.sh $Assembly $OutDir
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/ascomycota_odb9)
OutDir=$(dirname $Assembly)
qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
done
```


```bash
  for File in $(ls repeat_masked/*/*/filtered_contigs/run_*_contigs_unmasked/short_summary_*.txt); do
  Strain=$(echo $File| rev | cut -d '/' -f4 | rev)
  Organism=$(echo $File | rev | cut -d '/' -f5 | rev)
  Complete=$(cat $File | grep "(C)" | cut -f2)
  Fragmented=$(cat $File | grep "(F)" | cut -f2)
  Missing=$(cat $File | grep "(M)" | cut -f2)
  Total=$(cat $File | grep "Total" | cut -f2)
  echo -e "$Organism\t$Strain\t$Complete\t$Fragmented\t$Missing\t$Total"
  done
```




#### Aligning

Insert sizes of the RNA seq library were unknown until a draft alignment could
be made. To do this tophat and cufflinks were run, aligning the reads against a
single genome. The fragment length and stdev were printed to stdout while
cufflinks was running.

```bash
for Assembly in $(ls repeat_masked/*/*/*/*_contigs_unmasked.fa | grep '61B'); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
for File in $(ls ../../../../../data/scratch/armita/didymella/qc_rna/paired/*/*/*/unpaired/*_trim.fq.gz | grep -v 'PDB'); do
Jobs=$(qstat | grep 'sub_sta' | grep 'qw'| wc -l)
# while [ $Jobs -gt 1 ]; do
# sleep 1m
# printf "."
# Jobs=$(qstat | grep 'sub_sta' | grep 'qw'| wc -l)
# done
printf "\n"
Timepoint=$(echo $File | rev | cut -f4 -d '/' | rev)
echo "$Prefix"
Prefix=$(echo $File | rev | cut -f3 -d '/' | rev)
OutDir=../../../../../data/scratch/armita/didymella/alignment/star/$Organism/$Strain/$Timepoint/$Prefix
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
qsub $ProgDir/sub_star_unpaired.sh $Assembly $File $OutDir
# $ProgDir/sub_star_unpaired.sh $Assembly $File $OutDir
done
done
```



Splitting reads and trimming adapters using porechop

```bash
	for RawReads in $(ls ../../../../../data/scratch/armita/didymella/qc_rna/paired/*/*/*/unpaired/*_trim.fq.gz | grep 'PDB'); do
    Organism=D.rabiei
    Strain=P4
    echo "$Organism - $Strain"
  	OutDir=../../../../../data/scratch/armita/didymella/qc_rna/paired/D.rabiei/P4/PDB/unpaired/split
    mkdir -p $OutDir
    # gunzip -c $RawReads | split -l 10000000 - $OutDir/'PDB_split2_'
    for File in $(ls $OutDir/PDB_split2_* | grep -v 'fq.gz'); do
      cat $File | gzip -cf > ${File}.fq.gz
      rm $File
    done
  done
```


```bash
for Assembly in $(ls repeat_masked/*/*/*/*_contigs_unmasked.fa | grep '61B'); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
for File in $(ls ../../../../../data/scratch/armita/didymella/qc_rna/paired/*/*/*/unpaired/split/PDB_split2_*.fq.gz | grep 'PDB'); do
Jobs=$(qstat | grep 'sub_sta' | grep 'qw'| wc -l)
# while [ $Jobs -gt 1 ]; do
# sleep 1m
# printf "."
# Jobs=$(qstat | grep 'sub_sta' | grep 'qw'| wc -l)
# done
printf "\n"
Timepoint=$(echo $File | rev | cut -f4 -d '/' | rev)
echo "$Prefix"
Prefix=$(echo $File | rev | cut -f3 -d '/' | rev)
OutDir=../../../../../data/scratch/armita/didymella/alignment/star/$Organism/$Strain/$Timepoint/$Prefix
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
# qsub $ProgDir/sub_star_unpaired.sh $Assembly $File $OutDir
$ProgDir/sub_star_unpaired.sh $Assembly $File $OutDir
mv $OutDir/star $OutDir/$(basename ${File%.fq.gz})
done
done
```


Accepted hits .bam file were concatenated and indexed for use for gene model training:


```bash
qlogin -pe smp 4
cd /home/groups/harrisonlab/project_files/didymella
for OutDir in $(ls -d ../../../../../data/scratch/armita/didymella/alignment/star/*/* | grep '61B'); do
  Strain=$(echo $OutDir | rev | cut -d '/' -f1 | rev)
  Organism=$(echo $OutDir | rev | cut -d '/' -f2 | rev)
  echo "$Organism - $Strain"
  # For all alignments
  BamFiles=$(ls $OutDir/*/*/*/*.sortedByCoord.out.bam | tr -d '\n' | sed 's/.bam/.bam /g')
  mkdir -p $OutDir/concatenated
  samtools merge -@ 4 -f $OutDir/concatenated/concatenated.bam $BamFiles
done
logout
```



#### Braker prediction

```bash
	cp /home/armita/.gm_key_64 ~/.gm_key
```

```bash
for Assembly in $(ls repeat_masked/*/*/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep '61B'); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
mkdir -p alignment/$Organism/$Strain/concatenated
OutDir=gene_pred/braker/$Organism/"$Strain"_braker
AcceptedHits=$(ls ../../../../../data/scratch/armita/didymella/alignment/star/$Organism/$Strain/concatenated/concatenated.bam)
GeneModelName="$Organism"_"$Strain"_braker
rm -r /home/armita/prog/augustus-3.1/config/species/"$Organism"_"$Strain"_braker
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/braker1
qsub $ProgDir/sub_braker_fungi.sh $Assembly $OutDir $AcceptedHits $GeneModelName
done
```

** Number of genes predicted:  **


## Supplimenting Braker gene models with CodingQuary genes

Additional genes were added to Braker gene predictions, using CodingQuary in
pathogen mode to predict additional regions.

Firstly, aligned RNAseq data was assembled into transcripts using Cufflinks.

Note - cufflinks doesn't always predict direction of a transcript and
therefore features can not be restricted by strand when they are intersected.

```bash
for Assembly in $(ls repeat_masked/*/*/*/*_contigs_unmasked.fa | grep '61B'); do
Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
OutDir=gene_pred/cufflinks/$Organism/$Strain/concatenated
mkdir -p $OutDir
AcceptedHits=$(ls ../../../../../data/scratch/armita/didymella/alignment/star/$Organism/$Strain/concatenated/concatenated.bam)
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/RNAseq
qsub $ProgDir/sub_cufflinks.sh $AcceptedHits $OutDir
done
```


Secondly, genes were predicted using CodingQuary:

```bash
  for Assembly in $(ls repeat_masked/*/*/*/*_contigs_unmasked.fa | grep '61B'); do
    Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
    Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
    echo "$Organism - $Strain"
    Jobs=$(qstat | grep 'sub_cuffli' | wc -l)
    while [ $Jobs -gt 0 ]; do
    sleep 1m
    printf "."
    Jobs=$(qstat | grep 'sub_cuffli' | wc -l)
    done
    printf "\n"
    OutDir=gene_pred/codingquary/$Organism/$Strain
    CufflinksGTF=$(ls gene_pred/cufflinks/$Organism/$Strain/concatenated/transcripts.gtf)
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
    qsub $ProgDir/sub_CodingQuary.sh $Assembly $CufflinksGTF $OutDir
  done
```

Then, additional transcripts were added to Braker gene models, when CodingQuary
genes were predicted in regions of the genome, not containing Braker gene
models:

Note - only 41 secreted proteins were predicted by CodingQuarry preventing it
being run in pathogen mode



```bash
for BrakerGff in $(ls gene_pred/braker/*/*_braker/*/augustus.gff3 | grep '61B'); do
Strain=$(echo $BrakerGff| rev | cut -d '/' -f3 | rev | sed 's/_braker//g')
Organism=$(echo $BrakerGff | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
Assembly=$(ls repeat_masked/$Organism/$Strain/filtered_contigs/*_softmasked_repeatmasker_TPSI_appended.fa)
CodingQuaryGff=gene_pred/codingquary/$Organism/$Strain/out/PredictedPass.gff3
# PGNGff=gene_pred/codingquary/$Organism/$Strain/out/PGN_predictedPass.gff3
AddDir=gene_pred/codingquary/$Organism/$Strain/additional
FinalDir=gene_pred/final/$Organism/$Strain/final
AddGenesList=$AddDir/additional_genes.txt
AddGenesGff=$AddDir/additional_genes.gff
FinalGff=$AddDir/combined_genes.gff
mkdir -p $AddDir
mkdir -p $FinalDir

bedtools intersect -v -a $CodingQuaryGff -b $BrakerGff | grep 'gene'| cut -f2 -d'=' | cut -f1 -d';' > $AddGenesList
# bedtools intersect -v -a $PGNGff -b $BrakerGff | grep 'gene'| cut -f2 -d'=' | cut -f1 -d';' >> $AddGenesList
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation
$ProgDir/gene_list_to_gff.pl $AddGenesList $CodingQuaryGff CodingQuarry_v2.0 ID CodingQuary > $AddGenesGff
# $ProgDir/gene_list_to_gff.pl $AddGenesList $PGNGff PGNCodingQuarry_v2.0 ID CodingQuary >> $AddGenesGff
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
# -
# This section is edited
$ProgDir/add_CodingQuary_features.pl $AddGenesGff $Assembly > $AddDir/add_genes_CodingQuary_unspliced.gff3
$ProgDir/correct_CodingQuary_splicing.py --inp_gff $AddDir/add_genes_CodingQuary_unspliced.gff3 > $FinalDir/final_genes_CodingQuary.gff3
# -
$ProgDir/gff2fasta.pl $Assembly $FinalDir/final_genes_CodingQuary.gff3 $FinalDir/final_genes_CodingQuary
cp $BrakerGff $FinalDir/final_genes_Braker.gff3
$ProgDir/gff2fasta.pl $Assembly $FinalDir/final_genes_Braker.gff3 $FinalDir/final_genes_Braker
cat $FinalDir/final_genes_Braker.pep.fasta $FinalDir/final_genes_CodingQuary.pep.fasta | sed -r 's/\*/X/g' > $FinalDir/final_genes_combined.pep.fasta
cat $FinalDir/final_genes_Braker.cdna.fasta $FinalDir/final_genes_CodingQuary.cdna.fasta > $FinalDir/final_genes_combined.cdna.fasta
cat $FinalDir/final_genes_Braker.gene.fasta $FinalDir/final_genes_CodingQuary.gene.fasta > $FinalDir/final_genes_combined.gene.fasta
cat $FinalDir/final_genes_Braker.upstream3000.fasta $FinalDir/final_genes_CodingQuary.upstream3000.fasta > $FinalDir/final_genes_combined.upstream3000.fasta


GffBraker=$FinalDir/final_genes_Braker.gff3
GffQuary=$FinalDir/final_genes_CodingQuary.gff3
GffAppended=$FinalDir/final_genes_appended.gff3
cat $GffBraker $GffQuary > $GffAppended
done
```

In preperation for submission to ncbi, gene models were renamed and duplicate gene features were identified and removed.
 * no duplicate genes were identified


 ```bash
 for Gff in $(ls gene_pred/final/*/*/final/final_genes_appended.gff3 | grep '61B'); do
 	Strain=$(echo $Gff | rev | cut -d '/' -f3 | rev)
 	Organism=$(echo $Gff | rev | cut -d '/' -f4 | rev)
 	echo "$Strain - $Organism"
 	cat $Gff | grep -w 'gene' | wc -l
 done
 ```

 ```
 61B - D.pinodella
 10908
 ```


 In preperation for submission to ncbi, gene models were renamed and duplicate gene features were identified and removed.
  * no duplicate genes were identified


 ```bash
 for GffAppended in $(ls gene_pred/final/*/*/final/final_genes_appended.gff3 | grep '61B'); do
 Strain=$(echo $GffAppended | rev | cut -d '/' -f3 | rev)
 Organism=$(echo $GffAppended | rev | cut -d '/' -f4 | rev)
 echo "$Organism - $Strain"
 FinalDir=gene_pred/final_genes/$Organism/$Strain/final
 mkdir -p $FinalDir
 ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
 # $ProgDir/remove_dup_features.py --inp_gff $GffAppended
 # $ProgDir/remove_dup_features.py --inp_gff $GffAppended | grep -A2 'Duplicate gene found' | tail -n1 | cut -f2 -d'=' > $FinalDir/filter_list.tmp
 GffFiltered=$FinalDir/filtered_duplicates.gff
 # cat $GffAppended | grep -v -w -f $FinalDir/filter_list.tmp > $GffFiltered
 # rm $FinalDir/filter_list.tmp
 $ProgDir/remove_dup_features.py --inp_gff $GffAppended --out_gff $GffFiltered
 GffRenamed=$FinalDir/final_genes_appended_renamed.gff3
 LogFile=$FinalDir/final_genes_appended_renamed.log
 ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/codingquary
 $ProgDir/gff_rename_genes.py --inp_gff $GffFiltered --conversion_log $LogFile > $GffRenamed
 rm $GffFiltered

 Assembly=$(ls repeat_masked/$Organism/$Strain/filtered_contigs/*_softmasked_repeatmasker_TPSI_appended.fa)
 $ProgDir/gff2fasta.pl $Assembly $GffRenamed $FinalDir/final_genes_appended_renamed

 # The proteins fasta file contains * instead of Xs for stop codons, these should
 # be changed
 sed -i 's/\*/X/g' $FinalDir/final_genes_appended_renamed.pep.fasta
 done
 ```


## Assessing the Gene space in predicted transcriptomes:

 ```bash
 for Assembly in $(ls gene_pred/final_genes/*/*/final/final_genes_appended_renamed.gene.fasta | grep '61B'); do
 Strain=$(echo $Assembly| rev | cut -d '/' -f3 | rev)
 Organism=$(echo $Assembly | rev | cut -d '/' -f4 | rev)
 echo "$Organism - $Strain"
 ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/busco
 # BuscoDB="Fungal"
 BuscoDB=$(ls -d /home/groups/harrisonlab/dbBusco/ascomycota_odb9)
 OutDir=gene_pred/busco/$Organism/$Strain/genes
 qsub $ProgDir/sub_busco3.sh $Assembly $BuscoDB $OutDir
 done
 ```

 ```bash
 	for File in $(ls gene_pred/busco/*/*/genes/*/short_summary_*.txt); do  
 		echo $File;
 		cat $File | grep -e '(C)' -e 'Total';
 	done
 ```


#Functional annotation

## A) Interproscan

Interproscan was used to give gene models functional annotations.
Annotation was run using the commands below:

Note: This is a long-running script. As such, these commands were run using
'screen' to allow jobs to be submitted and monitored in the background.
This allows the session to be disconnected and reconnected over time.

Screen ouput detailing the progress of submission of interporscan jobs
was redirected to a temporary output file named interproscan_submission.log .

```bash
	screen -a
	ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
	for Genes in $(ls gene_pred/final_genes/*/*/final/final_genes_appended_renamed.pep.fasta | grep '61B'); do
	echo $Genes
	$ProgDir/sub_interproscan.sh $Genes
	done 2>&1 | tee -a interproscan_submisison.log
```

Following interproscan annotation split files were combined using the following
commands:

```bash
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/interproscan
for Proteins in $(ls gene_pred/final_genes/*/*/final/final_genes_appended_renamed.pep.fasta | grep '61B'); do
Strain=$(echo $Proteins | rev | cut -d '/' -f3 | rev)
Organism=$(echo $Proteins | rev | cut -d '/' -f4 | rev)
echo "$Organism - $Strain"
echo $Strain
InterProRaw=gene_pred/interproscan/$Organism/$Strain/raw
$ProgDir/append_interpro.sh $Proteins $InterProRaw
done
```


## B) SwissProt


```bash
for Proteome in $(ls gene_pred/final_genes/*/*/final/final_genes_appended_renamed.pep.fasta | grep '61B'); do
Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
OutDir=gene_pred/swissprot/$Organism/$Strain
SwissDbDir=../../uniprot/swissprot
SwissDbName=uniprot_sprot
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/swissprot
qsub $ProgDir/sub_swissprot.sh $Proteome $OutDir $SwissDbDir $SwissDbName
done
```

## Small secreted proteins

Putative effectors identified within Augustus gene models using a number
of approaches:

* A) From Braker gene models - Signal peptide & small cystein rich protein


### A) From Augustus gene models - Identifying secreted proteins

Required programs:
* SigP
* biopython
* TMHMM


Proteins that were predicted to contain signal peptides were identified using
the following commands:

```bash
 for Proteome in $(ls gene_pred/final_genes/*/*/final/final_genes_appended_renamed.pep.fasta | grep '61B'); do
   SplitfileDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
   ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/signal_peptides
   Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
   Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
   SplitDir=gene_pred/braker_split/$Organism/$Strain
   mkdir -p $SplitDir
   BaseName="$Organism""_$Strain"_braker_preds
   $SplitfileDir/splitfile_500.py --inp_fasta $Proteome --out_dir $SplitDir --out_base $BaseName
   for File in $(ls $SplitDir/*_braker_preds_*); do
   Jobs=$(qstat | grep 'pred_sigP' | wc -l)
   while [ $Jobs -gt '20' ]; do
   sleep 10
   printf "."
   Jobs=$(qstat | grep 'pred_sigP' | wc -l)
   done
   printf "\n"
   echo $File
   qsub $ProgDir/pred_sigP.sh $File signalp-4.1
   done
 done
```

The batch files of predicted secreted proteins needed to be combined into a
single file for each strain. This was done with the following commands:
```bash
 for SplitDir in $(ls -d gene_pred/braker_split/*/* | grep '61B'); do
   Strain=$(echo $SplitDir | cut -d '/' -f4)
   Organism=$(echo $SplitDir | cut -d '/' -f3)
   InStringAA=''
   InStringNeg=''
   InStringTab=''
   InStringTxt=''
   SigpDir=braker_signalp-4.1
   echo "$Organism - $Strain"
   for GRP in $(ls -l $SplitDir/*_braker_preds_*.fa | rev | cut -d '_' -f1 | rev | sort -n); do  
     InStringAA="$InStringAA gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_braker_preds_$GRP""_sp.aa";  
     InStringNeg="$InStringNeg gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_braker_preds_$GRP""_sp_neg.aa";  
     InStringTab="$InStringTab gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_braker_preds_$GRP""_sp.tab";
     InStringTxt="$InStringTxt gene_pred/$SigpDir/$Organism/$Strain/split/"$Organism"_"$Strain"_braker_preds_$GRP""_sp.txt";  
   done
   cat $InStringAA > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_aug_sp.aa
   cat $InStringNeg > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_aug_neg_sp.aa
   tail -n +2 -q $InStringTab > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_aug_sp.tab
   cat $InStringTxt > gene_pred/$SigpDir/$Organism/$Strain/"$Strain"_aug_sp.txt
 done
```

Some proteins that are incorporated into the cell membrane require secretion.
Therefore proteins with a transmembrane domain are not likely to represent
cytoplasmic or apoplastic effectors.

Proteins containing a transmembrane domain were identified:

```bash
 for Proteome in $(ls gene_pred/final_genes/*/*/final/final_genes_appended_renamed.pep.fasta | grep '61B'); do
   Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
   Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
   ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/feature_annotation/transmembrane_helices
   qsub $ProgDir/submit_TMHMM.sh $Proteome
 done
```

Those proteins with transmembrane domains were removed from lists of Signal
peptide containing proteins

```bash
for File in $(ls gene_pred/trans_mem/*/*/*_TM_genes_neg.txt | grep '61B'); do
Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
# echo "$Organism - $Strain"
NonTmHeaders=$(echo "$File" | sed 's/neg.txt/neg_headers.txt/g')
cat $File | cut -f1 > $NonTmHeaders
SigP=$(ls gene_pred/braker_signalp-4.1/$Organism/$Strain/"$Strain"_aug_sp.aa)
OutDir=$(dirname $SigP)
ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
$ProgDir/extract_from_fasta.py --fasta $SigP --headers $NonTmHeaders > $OutDir/"$Strain"_final_sp_no_trans_mem.aa
# echo "Number of SigP proteins:"
TotalProts=$(cat $SigP | grep '>' | wc -l)
# echo "Number without transmembrane domains:"
SecProt=$(cat $OutDir/"$Strain"_final_sp_no_trans_mem.aa | grep '>' | wc -l)
# echo "Number of gene models:"
SecGene=$(cat $OutDir/"$Strain"_final_sp_no_trans_mem.aa | grep '>' | cut -f1 -d't' | sort | uniq |wc -l)
# A text file was also made containing headers of proteins testing +ve
PosFile=$(ls gene_pred/trans_mem/$Organism/$Strain/"$Strain"_TM_genes_pos.txt)
TmHeaders=$(echo $PosFile | sed 's/.txt/_headers.txt/g')
cat $PosFile | cut -f1 > $TmHeaders
printf "$Organism\t$Strain\t$TotalProts\t$SecProt\t$SecGene\n"
done
```

```
  D.pinodella	61B	1062	854	852
```



### C) From Augustus gene models - Effector identification using EffectorP

Required programs:
 * EffectorP.py

```bash
  for Proteome in $(ls gene_pred/final_genes/*/*/final/final_genes_appended_renamed.pep.fasta | grep '61B'); do
    Strain=$(echo $Proteome | rev | cut -f3 -d '/' | rev)
    Organism=$(echo $Proteome | rev | cut -f4 -d '/' | rev)
    BaseName="$Organism"_"$Strain"_EffectorP
    OutDir=analysis/effectorP/$Organism/$Strain
    ProgDir=~/git_repos/emr_repos/tools/seq_tools/feature_annotation/fungal_effectors
    qsub $ProgDir/pred_effectorP.sh $Proteome $BaseName $OutDir
  done
```

Those genes that were predicted as secreted and tested positive by effectorP
were identified:

Note - this doesnt exclude proteins with TM domains or GPI anchors

```bash
  for File in $(ls analysis/effectorP/*/*/*_EffectorP.txt | grep '61B'); do
    Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
    Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
    echo "$Organism - $Strain"
    Headers=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_headers.txt/g')
    cat $File | grep 'Effector' | grep -v 'Effector probability:' | cut -f1 > $Headers
    printf "EffectorP headers:\t"
    cat $Headers | wc -l
    Secretome=$(ls gene_pred/braker_signalp-4.1/$Organism/$Strain/"$Strain"_final_sp_no_trans_mem.aa)
    OutFile=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_secreted.aa/g')
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
    $ProgDir/extract_from_fasta.py --fasta $Secretome --headers $Headers > $OutFile
    OutFileHeaders=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_secreted_headers.txt/g')
    cat $OutFile | grep '>' | tr -d '>' > $OutFileHeaders
    printf "Secreted effectorP headers:\t"
    cat $OutFileHeaders | wc -l
    Gff=$(ls gene_pred/final_genes/$Organism/$Strain/*/final_genes_appended_renamed.gff3)
    EffectorP_Gff=$(echo "$File" | sed 's/_EffectorP.txt/_EffectorP_secreted.gff/g')
    ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
    $ProgDir/extract_gff_for_sigP_hits.pl $OutFileHeaders $Gff effectorP ID > $EffectorP_Gff
  done
```
```
D.pinodella - 61B
EffectorP headers:	1629
Secreted effectorP headers:	116
```

## SSCP

Small secreted cysteine rich proteins were identified within secretomes. These
proteins may be identified by EffectorP, but this approach allows direct control
over what constitutes a SSCP.

```bash
for Secretome in $(ls gene_pred/braker_signalp-4.1/*/*/*_final_sp_no_trans_mem.aa | grep '61B'); do
Strain=$(echo $Secretome| rev | cut -f2 -d '/' | rev)
Organism=$(echo $Secretome | rev | cut -f3 -d '/' | rev)
echo "$Organism - $Strain"
OutDir=analysis/sscp/$Organism/$Strain
mkdir -p $OutDir
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/sscp
$ProgDir/sscp_filter.py --inp_fasta $Secretome --max_length 300 --threshold 3 --out_fasta $OutDir/"$Strain"_sscp_all_results.fa
cat $OutDir/"$Strain"_sscp_all_results.fa | grep 'Yes' > $OutDir/"$Strain"_sscp.fa
printf "number of SSC-rich genes:\t"
cat $OutDir/"$Strain"_sscp.fa | grep '>' | tr -d '>' | cut -f1 -d '.' | sort | uniq | wc -l
printf "Number of effectors predicted by EffectorP:\t"
EffectorP=$(ls analysis/effectorP/$Organism/$Strain/*_EffectorP_secreted_headers.txt)
cat $EffectorP | wc -l
printf "Number of SSCPs predicted by both effectorP and this approach: \t"
cat $OutDir/"$Strain"_sscp.fa | grep '>' | tr -d '>' > $OutDir/"$Strain"_sscp_headers.txt
cat $OutDir/"$Strain"_sscp_headers.txt $EffectorP | cut -f1 | sort | uniq -d | wc -l
echo ""
done
```

```
  D.pinodella - 61B
  % cysteine content threshold set to:	3
  maximum length set to:	300
  No. short-cysteine rich proteins in input fasta:	65
  number of SSC-rich genes:	65
  Number of effectors predicted by EffectorP:	116
  Number of SSCPs predicted by both effectorP and this approach: 	43
```

### C) Identification of MIMP-flanking genes

```bash
for Assembly in $(ls repeat_masked/*/*/filtered_contigs/*_contigs_unmasked.fa | grep '61B'); do
Organism=$(echo "$Assembly" | rev | cut -d '/' -f4 | rev)
Strain=$(echo "$Assembly" | rev | cut -d '/' -f3 | rev)
GeneGff=$(ls gene_pred/final_genes/$Organism/"$Strain"/final/final_genes_appended_renamed.gff3)
OutDir=analysis/mimps/$Organism/$Strain
mkdir -p "$OutDir"
echo "$Organism - $Strain"
ProgDir="/home/armita/git_repos/emr_repos/tools/pathogen/mimp_finder"
$ProgDir/mimp_finder.pl $Assembly $OutDir/"$Strain"_mimps.fa $OutDir/"$Strain"_mimps.gff > $OutDir/"$Strain"_mimps.log
$ProgDir/gffexpander.pl +- 2000 $OutDir/"$Strain"_mimps.gff > $OutDir/"$Strain"_mimps_exp.gff
echo "The number of mimps identified:"
cat $OutDir/"$Strain"_mimps.fa | grep '>' | wc -l
bedtools intersect -u -a $GeneGff -b $OutDir/"$Strain"_mimps_exp.gff > $OutDir/"$Strain"_genes_in_2kb_mimp.gff
echo "The following transcripts intersect mimps:"
MimpProtsTxt=$OutDir/"$Strain"_prots_in_2kb_mimp.txt
MimpGenesTxt=$OutDir/"$Strain"_genes_in_2kb_mimp.txt
cat $OutDir/"$Strain"_genes_in_2kb_mimp.gff | grep -w 'mRNA' | cut -f9 | cut -f1 -d';' | cut -f2 -d'=' | sort | uniq > $MimpProtsTxt
cat $OutDir/"$Strain"_genes_in_2kb_mimp.gff | grep -w 'mRNA' | cut -f9 | cut -f1 -d';' | cut -f2 -d'=' | cut -f1 -d '.'| sort | uniq > $MimpGenesTxt
cat $MimpProtsTxt | wc -l
cat $MimpGenesTxt | wc -l
echo ""
done
```

Extraction of the B-tubulin region of IMV00293 revealed it to be a F. oxysporum
isolate and distinct from the other two isolates.
```
D.pinodella - 61B
The number of mimps identified:
0
The following transcripts intersect mimps:
0
0
```

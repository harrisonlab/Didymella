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
    gunzip -c $RawReads | split -l 100000000 - $OutDir/'PDB_split_'
    for File in $(ls $OutDir/PDB_split_* | grep -v 'fq.gz'); do
      cat $File | gzip -cf > ${File}.fq.gz
      # rm $File
    done
  done
```


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
# qsub $ProgDir/sub_star_unpaired.sh $Assembly $File $OutDir
$ProgDir/sub_star_unpaired.sh $Assembly $File $OutDir
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

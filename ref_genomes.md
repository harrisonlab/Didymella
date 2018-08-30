


from local machine
```bash
scp Did*.zip cluster:/home/groups/harrisonlab/project_files/didymella/assembly/external/.
Phot*.zip cluster:/home/groups/harrisonlab/project_files/didymella/assembly/external/.
```

```bash
cd assembly/external
for File in *.zip; do
unzip $File -d .
done

WorkDir=/home/groups/harrisonlab/project_files/didymella

Organism=D.exigua
Strain=CBS_183_55
OutDir=assembly/external/$Organism/$Strain
mkdir -p $OutDir
mv assembly/external/Didex1 $OutDir/.
cd $OutDir
gunzip */*.gz
cd $WorkDir

Organism=D.zeae_maydis
Strain=3018
OutDir=assembly/external/$Organism/$Strain
mkdir -p $OutDir
mv assembly/external/Didma1 $OutDir/.
gunzip $OutDir/*/*.gz

Organism=P.tracheiphila
Strain=IPT5
OutDir=assembly/external/$Organism/$Strain
mkdir -p $OutDir
mv assembly/external/Photr1 $OutDir/.
gunzip $OutDir/*/*.gz
```

From ncbi:

```bash
WorkDir=/home/groups/harrisonlab/project_files/didymella
Orgnaism=P.herbarum
Strain=JCM_15942
OutDir=assembly/external/$Organism/$Strain/genbank
mkdir -p $OutDir
cd $OutDir
wget ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/BC/GR/BCGR01/BCGR01.1.fsa_nt.gz
cat *.fsa_nt.gz | gunzip -cf > ${Strain}_ncbi.fa
cp -s ${Strain}_ncbi.fa $PWD/genome.ctg.fa
cd $WorkDir
Organism=D.rabiei
Strain=ITCC_4638
OutDir=assembly/external/$Organism/$Strain/genbank
mkdir -p $OutDir
cd $OutDir
wget ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/JY/NV/JYNV01/JYNV01.1.fsa_nt.gz
cat *.fsa_nt.gz | gunzip -cf > ${Strain}_ncbi.fa
cp -s ${Strain}_ncbi.fa $PWD/genome.ctg.fa
cd $WorkDir
```


## RNAseq data


RNAseq data was downloaded from ncbi (one sample per timepoint of Mycelium, 96hpi, 36hpi and 12hpi):
https://www.ncbi.nlm.nih.gov/sra/SRX1080581
https://www.ncbi.nlm.nih.gov/sra/SRX1080577
https://www.ncbi.nlm.nih.gov/sra/SRX1080571
https://www.ncbi.nlm.nih.gov/sra/SRX1080565


```bash
# cd /home/groups/harrisonlab/project_files/N.haematococca
mkdir -p /data/scratch/armita/didymella
cd /data/scratch/armita/didymella
# PDB
OutDir=raw_rna/D.rabiei/P4/PDB
mkdir -p $OutDir
fastq-dump --split-files -A SRR2086542 --gzip --outdir $OutDir
# 96 hpi
OutDir=raw_rna/D.rabiei/P4/chickpea_96hpi
mkdir -p $OutDir
fastq-dump --split-files -A SRR2086538 --gzip --outdir $OutDir
# 36 hpi
OutDir=raw_rna/D.rabiei/P4/chickpea_36hpi
mkdir -p $OutDir
fastq-dump --split-files -A SRR2086532 --gzip --outdir $OutDir
# 12 hpi
OutDir=raw_rna/D.rabiei/P4/chickpea_12hpi
mkdir -p $OutDir
fastq-dump --split-files -A SRR2086526 --gzip --outdir $OutDir
```

```bash
WorkDir=/home/groups/harrisonlab/project_files/didymella
cd $WorkDir
Organism=D.rabiei
Strain=ITCC_4638
OutDir=assembly/external/$Organism/$Strain/genbank
mkdir -p $OutDir
cd $OutDir
wget ftp://ftp.ncbi.nlm.nih.gov/sra/wgs_aux/JY/NV/JYNV01/JYNV01.1.fsa_nt.gz
cat *.fsa_nt.gz | gunzip -cf > ${Strain}_ncbi.fa
cp -s ${Strain}_ncbi.fa $PWD/genome.ctg.fa
cd $WorkDir
```

Trimming of RNAseq data was performed

```bash
cd /home/groups/harrisonlab/project_files/didymella

for ReadsF in $(ls ../../../../../data/scratch/armita/didymella/raw_rna/*/*/*/*.fastq.gz); do    
Strain=$(echo $ReadsF| rev | cut -d '/' -f3 | rev)
Organism=$(echo $ReadsF | rev | cut -d '/' -f4 | rev)
TimePoint=$(echo $ReadsF | rev | cut -d '/' -f2 | rev)
echo "$Organism - $Strain - $TimePoint"
ProgDir=/home/armita/git_repos/emr_repos/tools/seq_tools/rna_qc
IlluminaAdapters=/home/armita/git_repos/emr_repos/tools/seq_tools/ncbi_adapters.fa
# ReadsR=$(echo "$ReadsF" | sed 's/_1.fastq.gz/_2.fastq.gz/g')
OutDir=../../../../../data/scratch/armita/didymella/qc_rna/paired/$Organism/$Strain/$TimePoint
echo $ReadsF
echo $ReadsR
qsub $ProgDir/rna_qc_fastq-mcf_unpaired.sh $ReadsF $IlluminaAdapters RNA $OutDir
done
```



# Repeat Masking

Repeat masking was performed on the downloaded assemblies.

```bash
  for Assembly in $(ls assembly/external/*/*/*/* | grep -e '_AssemblyScaffolds.fasta' -e '_ncbi.fa' | grep 'ITCC_4638'); do
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

for File in $(ls repeat_masked/*/*/*/*_contigs_softmasked.fa); do
OutDir=$(dirname $File)
TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
OutFile=$(echo $File | sed 's/_contigs_softmasked.fa/_contigs_softmasked_repeatmasker_TPSI_appended.fa/g')
echo "$OutFile"
bedtools maskfasta -soft -fi $File -bed $TPSI -fo $OutFile
echo "Number of masked bases:"
cat $OutFile | grep -v '>' | tr -d '\n' | awk '{print $0, gsub("[a-z]", ".")}' | cut -f2 -d ' '
done
# The number of N's in hardmasked sequence are not counted as some may be present within the assembly and were therefore not repeatmasked.
for File in $(ls repeat_masked/*/*/*/*_contigs_hardmasked.fa); do
OutDir=$(dirname $File)
TPSI=$(ls $OutDir/*_contigs_unmasked.fa.TPSI.allHits.chains.gff3)
OutFile=$(echo $File | sed 's/_contigs_hardmasked.fa/_contigs_hardmasked_repeatmasker_TPSI_appended.fa/g')
echo "$OutFile"
bedtools maskfasta -fi $File -bed $TPSI -fo $OutFile
done
```


Quast and busco were run to assess the effects of racon on assembly quality:

```bash
for Assembly in $(ls repeat_masked/*/*/*/*_contigs_unmasked.fa); do
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

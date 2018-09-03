# Commands to build reference databases for Blast searches

species specific effectors need to be identified that show divergence between species. As such, blast databases of reference genomes need to be constructed and searched against.


Key genomes on this list are:

*

blast databases were made in the fusarium project and are found on the cluster
at analysis/diagnostics_blast

## Didymella effectorP analysis

### Prepare queries

```bash
mkdir -p analysis/diagnostics_blast
for File in $(ls analysis/effectorP/D.pinodella/*/*_EffectorP_secreted_headers.txt); do
  Strain=$(echo $File | rev | cut -f2 -d '/' | rev)
  Organism=$(echo $File | rev | cut -f3 -d '/' | rev)
  GeneFa=$(ls gene_pred/final_genes/$Organism/$Strain/final/final_genes_appended_renamed.cdna.fasta)
  cat $File | cut -f1 > tmp.txt
  ProgDir=/home/armita/git_repos/emr_repos/tools/gene_prediction/ORF_finder
  $ProgDir/extract_from_fasta.py --fasta $GeneFa --headers tmp.txt \
    | sed "s/>/>${Strain}|/g"
done > analysis/diagnostics_blast/D.pinodella_effP.fa
```

### Perform BLAST

Build a reference genome database:



Blast vs the reference genome databases:

```bash
for RefGenome in $(ls assembly/external/*/*/*/*_AssemblyScaffolds.fasta); do
Organism=$(echo $RefGenome | rev | cut -f4 -d '/' | rev)
Strain=$(echo $RefGenome | rev | cut -f3 -d '/' | rev)
Prefix="${Organism}_${Strain}"
# Jobs=$(qstat | grep 'run_blast' | grep 'qw'| wc -l)
# while [ $Jobs -gt 1 ]; do
# sleep 20s
# printf "."
# Jobs=$(qstat | grep 'run_blast' | grep 'qw'| wc -l)
# done
# printf "\n"
# Prefix=$(echo $RefGenome | cut -f5,6 -d '/' --output-delimiter '_')
echo $Prefix
OutDir=analysis/diagnostics_blast/vs_ref_genomes/$Prefix
mkdir -p $OutDir
CurDir=$PWD
cd $OutDir
rm ${Prefix}_genome.fa
cp -s $CurDir/$RefGenome ${Prefix}_genome.fa
cd $CurDir
OutDir=analysis/diagnostics_blast/vs_ref_genomes/$Prefix
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
qsub $ProgDir/run_blast2csv.sh analysis/diagnostics_blast/D.pinodella_effP.fa dna $RefGenome $OutDir
done
```

```bash

for RefGenome in $(ls assembly/external/*/*/genbank/genome.ctg.fa); do
Organism=$(echo $RefGenome | rev | cut -f4 -d '/' | rev)
Strain=$(echo $RefGenome | rev | cut -f3 -d '/' | rev)
Prefix="${Organism}_${Strain}"
echo $Prefix
OutDir=analysis/diagnostics_blast/vs_ref_genomes/$Prefix
mkdir -p $OutDir
CurDir=$PWD
cd $OutDir
rm ${Prefix}_genome.fa
cat $CurDir/$RefGenome | sed "s/gi|.*|dbj|//g" | cut -f1 -d'|' | sed 's/\./_/g' > ${Prefix}_genome.fa
cd $CurDir
OutDir=analysis/diagnostics_blast/vs_ref_genomes/$Prefix
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
qsub $ProgDir/run_blast2csv.sh analysis/diagnostics_blast/D.pinodella_effP.fa dna $RefGenome $OutDir
done
```


Blast against themselves

```bash
for RefGenome in $(ls repeat_masked/*/*/*/*_contigs_softmasked_repeatmasker_TPSI_appended.fa | grep '61B'); do
Prefix=$(basename $RefGenome _contigs_softmasked_repeatmasker_TPSI_appended.fa)
OutDir=analysis/diagnostics_blast/vs_own_genomes/$Prefix
# Prefix=$(echo $RefGenome | cut -f5,6 -d '/' --output-delimiter '_')
echo $Prefix
OutDir=analysis/diagnostics_blast/vs_own_genomes/$Prefix
mkdir -p $OutDir
CurDir=$PWD
cd $OutDir
rm ${Prefix}_genome.fa
# cp -s $CurDir/$RefGenome ${Prefix}_genome.fa
cat $CurDir/$RefGenome > ${Prefix}_genome.fa
cd $CurDir
ProgDir=/home/armita/git_repos/emr_repos/tools/pathogen/blast
qsub $ProgDir/run_blast2csv.sh analysis/diagnostics_blast/D.pinodella_effP.fa dna $RefGenome $OutDir
done
```

## Summarise blast hists

```bash
CsvFiles=$(ls analysis/diagnostics_blast/vs_*_genomes/*/*D.pinodella_effP.fa_hits.csv | grep -v 'vs_ref_genomes' | grep -v 'ITCC_4638')
Headers=$(echo $CsvFiles | sed 's&analysis/diagnostics_blast/vs_ref_genomes/&&g' | sed 's&analysis/diagnostics_blast/vs_seq_genomes/&&g' | sed 's&analysis/diagnostics_blast/vs_own_genomes/&&g' | sed -r "s&_D.pinodella_effP.fa_hits.csv&&g")
OutDir=analysis/diagnostics_blast/vs_ref_genomes/extracted
mkdir -p $OutDir
CsvFiles=$(ls analysis/diagnostics_blast/vs_*_genomes/*/*D.pinodella_effP.fa_hits.csv  | grep -v 'vs_ref_genomes' | grep -v 'ITCC_4638')
Genomes=$(ls analysis/diagnostics_blast/vs_*_genomes/*/*_genome.fa | grep -v 'vs_ref_genomes' | grep -v 'ITCC_4638')
ProgDir=/home/armita/git_repos/emr_repos/scripts/fusarium/AHDB_project/blast_searches
$ProgDir/blast_parse_AHDB.py --blast_csv $CsvFiles --headers $Headers --genomes $Genomes --identity 0.50 --evalue 1e-30 --out_prefix $OutDir/D.pinodella_effP
ls $OutDir/D.pinodella_effP.csv
```


A number of interesting genes (as present in JS-169) were identified which showed presence in F. solani, but absent in other Fusarium.

```
JS-169|g297.t1
JS-169|g1364.t1
JS-169|g1548.t1
JS-169|g1634.t1
JS-169|g2349.t1
JS-169|g2358.t1
JS-169|g3150.t1
JS-169|g3156.t1
JS-169|g3907.t1
JS-169|g4181.t1
JS-169|g4671.t1
JS-169|g5415.t1
JS-169|g5689.t1
JS-169|g7048.t1
JS-169|g8051.t1
JS-169|g8329.t1
JS-169|g8399.t1
JS-169|g8720.t1
JS-169|g9452.t1
JS-169|g9606.t1
JS-169|g10403.t1
JS-169|g10978.t1
JS-169|g11008.t1
JS-169|g11028.t1
JS-169|g11116.t1
```



These alignments were downloaded for further study.


```bash
ExtractDir=analysis/diagnostics_blast/vs_ref_genomes/extracted
mkdir -p analysis/diagnostics_blast/selected_loci
GeneList="JS-169|g297.t1 JS-169|g1364.t1 JS-169|g1548.t1 JS-169|g1634.t1 JS-169|g2349.t1 JS-169|g2358.t1 JS-169|g3150.t1 JS-169|g3156.t1 JS-169|g3907.t1 JS-169|g4181.t1 JS-169|g4671.t1 JS-169|g5415.t1 JS-169|g5689.t1 JS-169|g7048.t1 JS-169|g8051.t1 JS-169|g8329.t1 JS-169|g8399.t1 JS-169|g8720.t1 JS-169|g9452.t1 JS-169|g9606.t1 JS-169|g10403.t1 JS-169|g10978.t1 JS-169|g11008.t1 JS-169|g11028.t1 JS-169|g11116.t1"
x=$(echo $GeneList | sed 's/|/_/g')
for Name in $x; do
  ls $ExtractDir/D.pinodella_effP_${Name}_hits.fa
  cp $ExtractDir/D.pinodella_effP_${Name}_hits.fa analysis/diagnostics_blast/selected_loci/.
done

InterPro=$(ls gene_pred/interproscan/D.pinodella/JS-169/JS-169_interproscan.tsv)
x=$(echo $GeneList | sed "s/JS-169|//g")
for Gene in $x; do
  cat $InterPro | grep $Gene
done >> analysis/diagnostics_blast/selected_loci/interpro.tsv

```

```bash
mkdir -p analysis/diagnostics_blast/all_targets
for File in $(ls analysis/diagnostics_blast/vs_*_genomes/*/*_genome.fa); do
cp $File analysis/diagnostics_blast/all_targets/.
done
tar -cz -f analysis/diagnostics_blast/all_targets.tar.gz analysis/diagnostics_blast/all_targets
rm -r analysis/diagnostics_blast/all_targets
```

# Command list for kMetaShot benchmark

## Outline
1.[Benchmark on HMP genomes](#Benchmark on HMP genomes)<br/>
2.[Benchmark on SRR606249 mock](#Benchmark on SRR606249 mock)<br/>
3.[Benchmark on CAMI II data](#Benchmark on CAMI II data)<br/>

### Benchmark on HMP genomes

##### Data Download
The download of Human Microbiome Project (HMP) genome has been executed as follows:
```bash
#!/usr/bin/bash

mkdir all_fna all_rpt

cd all_rpt
wget https://ftp.ncbi.nlm.nih.gov/genomes/HUMAN_MICROBIOM/Bacteria/all.rpt.tar.gz
tar -xf all.rpt.tar.gz
rm all.rpt.tar.gz


cd ../all_fna
wget https://ftp.ncbi.nlm.nih.gov/genomes/HUMAN_MICROBIOM/Bacteria/all.fna.tar.gz
tar -xf all.fna.tar.gz
rm all.fna.tar.gz
for f in $(ls ); do
	cd ${f}
	if [ -e *.scaffold.fna.tgz ]; then
	name=$(basename $(ls | grep .scaffold.fna.tgz) .scaffold.fna.tgz)
	tar -xf *.tgz
	rm *.tgz
	for f1 in $(ls); do
		cat ${f1} >> ${name}.fna
	done
	gzip ${name}.fna
	rm *.fna
	echo ${name} done
	else
	gzip *.fna
	echo $(ls)
	fi
	cd ..
done

rsync --list-only rsync://ftp.ncbi.nlm.nih.gov/genomes/HUMAN_MICROBIOM/Bacteria/ > ftp_hmp_list.txt
```
Then, a composition file was prepared:
```python
import pandas as pd
import os

def df_gen():
    for f in os.listdir('./all_rpt'):
        for f1 in os.path.join('./all_rpt',f):
            if f1.endswith('.rpt'):
                a = pd.read_csv(os.path.join('./all_rpt',f ,f1), 
                                sep=r'[:=] ', 
                                header=None)
                a.set_index(0, inplace=True)
                a = a.T
                a['local_path'] = os.path.join('./all_rpt',f ,f1)
                yield a


ranks = ['superkingdom', 'phylum', 
         'class', 'order',
         'family', 'genus', 
         'species', 'taxid']
assembly_summary = pd.concat(df_gen(), axis=0)
columns = [a.strip().replace(' ','_') for a in assembly_summary.columns]
columns[3] = 'taxid'
assembly_summary.columns = columns
assembly_summary.to_csv('./assembly_summary_HMP.csv')
assum = assembly_summary[['taxid','Taxname','local_path']].drop_duplicates()
kassum = pd.read_hdf('kMetaShot_reference.h5','new_assemblysummary')
kassum = kassum[ranks].drop_duplicates()
final = assum.merge(kassum, on='taxid')
final.to_csv('./HMP_composition_fullpath.csv')
```


##### kMetaShot classification
```bash
kMetaShot_classifier_NV.py -b HMP_simlink_4gtdbtk/ \
                           -r kMetaShot_reference.h5 \
                           -p 10 \
                           -o kMetaShot_classif
```

#### CAMITAX
```bash
nextflow -q \
         -log camitax.log \
         -bg run CAMI-challenge/CAMITAX \
         --db ../../../../camitax_reference/ \
         --i ../HMP_simlink_4gtdbtk/ \
         --x fna
```

#### GTDBtk
```bash
gtdbtk classify_wf --genome_dir HMP_simlink_4gtdbtk/ \
                   --out_dir ./gtdbtk \
                   --cpus 10 \
                   -x gz
```
INFO: Using GTDB-Tk reference data version r89: /home/gdefazio/miniconda2/envs/gtdbtk/share/gtdbtk-1.0.2/db/

Then, the evaluation of classification performances was carried out as follows:
```bash
python3 camitax_transform.py ./camitax_classif/data/camitax.tsv

cd ./gtdbtk_classif
python3 gtdbtk_transform.py
cd ..
              
python3 assignment_evaluation2.py \
  -g ./gtdbtk_classif/gtdbtk_assignment_adjusted.csv \
  -x ./camitax_classif/camitax_analysis_resume_paths.csv \
  -k ./kMetaShot_classif/kMetaShot_classification_resume.csv \
  -a kMetaShot_referece.h5 \
  -n ./assembly_summary_HMP_complete_paths.csv \
  -b 939 \
  -HMP
```

### Benchmark on SRR606249 mock

#### Trimming
```bash
sickle pe -f SRR606249_1.fastq.gz \
          -r SRR606249_2.fastq.gz \
          -t sanger -g \
          -o SRR606249_trm_1.fq.gz \
          -p SRR606249_trm_2.fq.gz \
          -s SRR606249_s.fq.gz -q 32 -l 50
```
#### metaSPAdes assembly benchmark
```commandline
spades.py 
    --meta 
    -1 /home3/gdefazio/kMetaShot_benchmark/mock_SRR606249/SRR606249_trm_1.fq.gz 
    -2 /home3/gdefazio/kMetaShot_benchmark/mock_SRR606249/SRR606249_trm_2.fq.gz 
    -s /home3/gdefazio/kMetaShot_benchmark/mock_SRR606249/SRR606249_s.fq.gz 
    -o /home3/gdefazio/kMetaShot_benchmark/mock_SRR606249/assembly
    -t 40 
    -k 35,57,79,99
```

##### MetaBAT2

```commandline
metabat2 -i assembly/contigs.fasta 
         -o metabat_binning/bin 
         --maxP 99 
         --minS 98 
         -m 1500 
```
##### kMetaShot classification
```commandline
kMetaShot_classifier_NV.py -b metabat_binning/ 
                           -r /home3/gdefazio/kMetaShot_Bacteria_Archaea_ref/clndMashOne4strain/compressed_clndMashOne4strain_bacteria_archaea.h5 
                           -p 10 -o kMetaShot_classif
                           -a 0.1
```

##### CAMITAX
```commandline
nextflow -q 
         -log camitax.log 
         -bg run CAMI-challenge/CAMITAX 
         --db /home3/gdefazio/kMetaShot_benchmark/camitax_reference/ 
         --i ../metabat_binning/ 
         --x fa
```
##### GTDBtk

```commandline
gtdbtk classify_wf 
       --genome_dir ../metabat_binning/ 
       --out_dir ./ 
       --cpus 10 
       -x fa
```

#### MegaHIT assembly benchmark
```bash
megahit -o /home3/gdefazio/kMetaShot_benchmark/mock_SRR606249/megahit_metabat/assembly_35_57_79_99 \
        --k-list 35,57,79,99 \
        -t 20 \
        -1 /home3/gdefazio/kMetaShot_benchmark/mock_SRR606249/sample/SRR606249_trm_1.fq.gz \
        -2 /home3/gdefazio/kMetaShot_benchmark/mock_SRR606249/sample/SRR606249_trm_2.fq.gz \
        -r /home3/gdefazio/kMetaShot_benchmark/mock_SRR606249/sample/SRR606249_s.fq.gz
```

##### MetaBAT2
```bash
metabat2 -i final.contigs.fa \
         -o metabat_binning/bin \
         --maxP 99 \
         --minS 98 \
         -m 1500 
```

##### kMetaShot classification
```bash
kMetaShot_classifier_NV.py -b metabat_binning/ \
                           -r kMetaShot_reference.h5 \
                           -p 10 -o kMetaShot_classif -a 0.1
```

##### CAMITAX
```bash
nextflow -q \
         -log camitax.log \
         -bg run CAMI-challenge/CAMITAX \
         --db ../../../camitax_reference/ \
         --i ../metabat_binning/ \
         --x fa
```
##### GTDBtk

```bash
gtdbtk classify_wf \
       --genome_dir ../metabat_binning/ \
       --out_dir ./ \
       --cpus 10 \
       -x fa 
```

#### MASH distance assessment
This test was carried out to confirm the correctness of kMetaShot classification 
at strain level for bins assigned to species lacking a strain indication in the
composition file of the analyzed mock.

```bash
mash sketch -l to_compare_lst.txt -o to_compare_lst

mash dist -t -l to_compare_lst.msh to_compare_lst.txt > comparison.tsv 
```

### Benchmark on CAMI II data
#### Data download
For plant associated datasets:
```bash
wget https://frl.publisso.de/data/frl:6425521/plant_associated/long_read_pacbio
```

For Gastrointestinal tract datasets:
```bash
java -jar camiClient.jar \
    -l https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/CAMI_Gastrointestinal_tract
```

For Airways datasets:
```bash
java -jar camiClient.jar \
    -l  https://openstack.cebitec.uni-bielefeld.de:8080/swift/v1/CAMI_Airways
```

##### MetaBAT2 binning
For each sample of a specific dataset, contigs pre-assembled by Cami II were used 
to obtain bins by using MetaBAT2 as follows
```bash
metabat2 -i anonymous_gsa.fasta.gz \
         -o metabat_binning/bin \
         --maxP 99 \
         --minS 98 \
         -m 1500 
```

Considering the huge amount of samples to analyze (i.e. 10 Gastrointestinal Short Reads, 
10 Gastrointestinal PacBio, 10 Airways Short Reads, 10 Airways PacBio 
and 20 plant associated ) a script was used to manage this part of benchmark 
as follows. This is an example for Gastrointestinal Short Reads samples:
these commands were slightly adapted for each samples group.
About kMetaShot execution, for Gastrointestinal and Airways samples, MAGs were analyzed by using 
an *ass2ref* at 0.1. Otherwise, a 0.2 *ass2ref* was used for plant associated 
samples.

```bash
#!/bin/bash

kmetashotREF="/home3/gdefazio/kMetaShot_Bacteria_Archaea_ref/clndMashOne4strain/compressed_clndMashOne4strain_bacteria_archaea.h5"
camitaxREF="/home3/gdefazio/kMetaShot_benchmark/camitax_reference/"
shortReadPath="/home3/gdefazio/kMetaShot_benchmark/cami2/GI_tract/short_read"

## CAMI2 automatic GI short reads
for smpl in $( ls $shortReadPath | grep sample_  ); do
  echo "START___${smpl}___"
  samplePath=${shortReadPath}/${smpl}/contigs
  sampleName=$( echo $smpl | sed -e 's/^2018.01.23_11.53.11_//' )
  #   composition adjust and sample name returning
  sampleName=$( conda run -n kmetashot python3 /home3/gdefazio/kMetaShot_HQ/benchmark/cami2_composition_adjust.py \
                             ${samplePath} ${shortReadPath}/metadata.tsv )
  # MetaBAT2 binning
  conda run -n metabat2 /home/gdefazio/miniconda2/envs/metabat2/bin/metabat2 \
                            -i ${samplePath}/anonymous_gsa.fasta.gz \
                            -o ${samplePath}/bins/bin --maxP 99 --minS 98 -m 1500

  # GTDBTK classif
  conda run -n gtdbtk /home/gdefazio/miniconda2/envs/gtdbtk/bin/gtdbtk \
                           classify_wf  --genome_dir ${samplePath}/bins \
                           --out_dir ${samplePath}/gtdbtk_classif --cpus 10 -x fa

  # GTDBTK adjust
  cd ${samplePath}/gtdbtk_classif
  ln -s ${ncbi2gtdbtkB} ./ncbi_vs_gtdb_r89_bacteria.xlsx
  ln -s ${ncbi2gtdbtkA} ./ncbi_vs_gtdb_r89_archaea.xlsx
  conda run -n kmetashot python3 /home3/gdefazio/kMetaShot_HQ/benchmark/gtdbtk_transform.py

  # CAMITAX classif
  mkdir ${samplePath}/camitax_classif
  cd ${samplePath}/camitax_classif
  conda run -n camitax /home/gdefazio/miniconda2/envs/camitax/bin/nextflow \
                            -log ./camitax.log \
                            -bg run CAMI-challenge/CAMITAX \
                            --db ${camitaxREF} \
                            --i ../bins  \
                            --x fa

  while ! [[ -f ./data/camitax.tsv ]]; do
    sleep 60s
  done

  # CAMITAX adjust
  conda run -n kmetashot python3 /home3/gdefazio/kMetaShot_HQ/benchmark/camitax_transform.py \
                                 ./data/camitax.tsv

  # kMetaShot classif
  conda run -n gius \
  /home/gdefazio/miniconda2/envs/gius/bin/kMetaShot_classifier_NV.py \
  -b ${samplePath}/bins/ \
  -r ${kmetashotREF} \
  -p 10 -o ${samplePath}/kMetaShot_classif \
  -a {0.1, 0.2}

  # classif evaluation
  cd ${samplePath}
  conda run -n kmetashot python3 \
  /home3/gdefazio/kMetaShot_HQ/benchmark/assignment_evaluation2.py \
  -g ${samplePath}/gtdbtk_classif/gtdbtk_assignment_adjusted.csv \
  -x ${samplePath}/camitax_classif/camitax_analysis_resume_paths.csv \
  -k ${samplePath}/kMetaShot_classif/kMetaShot_classification_resume.csv \
  -a ${kmetashotREF} \
  -n $( ls ${samplePath} | grep _composition_fullpath.csv ) \
  -b $( ls ${samplePath}/bins | wc -l )

  echo "STOP___${sampleName}___"
done

```
##### CPU, RAM, Time benchmark
For this benchmark we used *cmdbench* Python package embedded in an
*in-house* provided Python script launched as follows:
```bash
./memory_measuring.py
```


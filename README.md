# CS 4775 Final Project

Contributors: Cameron Xu, Jackson Staniec, Alexandra Tarzanin

## Data Proprocessing
<img width="259" alt="image" src="https://github.com/JacksonStaniec/cs4775-project/assets/81320443/98f2dcc0-a7a3-4c50-b0b5-40b3e0d5097c">

Intermediate data from the proprocessing pipeline can be found in BigWig, BED, 
and Fasta formats on [Cornell Box](https://cornell.box.com/s/020dpidze9pz7kkxgdwcj3o960m3ctaf).

The preprocessed data ready for analysis can be found in [`./data`](./data).

## Oligo Analysis
To run oligo analysis on a single fasta file, run the following command:
```
python ./src/oligomer_analysis.py -f [data file] -k [motif length]
```
To repeat the results of oligo analysis used for our report, run the following
command:
```
python ./run_oligo.py
```
The outputs will be located in [`./results/oligo`](./results/oligo).

## Position Analysis
The [`test_motif_analysis`](./src/test_motif_analysis.ipynb) Python notebook gives
examples for how to run position analysis on single and multiple fasta files, as
well as the code used to produce the position analysis results in our report. 

## Motif Analysis
The [`analyze-results`](./results/analyze-results.ipynb) Python notebook walks 
through the process of performing correlation analysis on the discovered motifs
to select motifs for phylogenetic tree reconstruction.

## Phylogenetic Tree Reconstruction
The [`form_newick_tree`](./src/form_newick_tree.py) Python script was used to 
generate phylogenetic trees in Newick format given motif sequences. An 
[online tool](http://etetoolkit.org/treeview/) was used 
to create the Newick tree visualizations used in our report.

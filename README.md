# micit_prediction

Repo to study the metabolic production of species found in the article from [Science](https://science.sciencemag.org/content/368/6494/973), Nejman _et al_., 2020 

Quick explanation of scripts and files

Scripts:

_get_species.py_ is the script used to get the full species identification of hits from S2 and S4 (_get_species.ipynb_ is the same but in jupyter notebook)

_seqs2fasta.ipynb_ is the script used to get the fasta file from the large sequences csv file from the article

_species_frequencies.ipynb_ is the script to get a filtered table with species abundances from the S2 table (removing useless columns)

_exploration.R_ is the R script used to do the data analysis: calculate final species frequencies, calculate micit production by tissue and other stuff

Files:

_final_hits.csv_: file with full specified species and hits where they appear in the tissues

_unique_seq_hits.csv_: full species with reference 16S gene sequences

_species_freqs.xlsx_: file with normalised abundances from script _species_frequencies.ipynb_

_freqs_hits.xlsx_: final file with fixed relative abundances of species both in general and within hits

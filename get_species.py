#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 20 16:23:10 2020

@author: dani
"""

import pandas as pd
import os

path = '/home/dani/Documents/MRC_postdoc/Cancer/micit_prediction'
# change dir to my path
os.chdir(path)


hits = pd.read_excel('./sup_info/aay9189_TableS4.xlsx',
                     sheet_name='Hits',skiprows=2)

cols = ['Breast (N)','Breast (NAT)','Breast (T)','Lung (NAT)','Lung (T)',
        'Melanoma (T)','Pancreas (T)','Ovary (N+NAT)','Ovary (T)','Bone (T)',
        'GBM (T)','Colon (NAT)','Colon (T)']

# sum rows
hits['row_sum'] = hits[cols].sum(axis=1)

is_hit =  hits['row_sum'] >= 1

is_hit
# store only rows as hits
hits_redux = hits[is_hit]


# remove NaN
hits_redux = hits_redux[hits_redux['species'].notnull()]

# which have "species" in their species name
hits_species = hits_redux.species.str.contains('Unknown')

# opposite of previous selection, slopy way I know, but...
hits_species = hits_species == False

hits_final = hits_redux[hits_species]

hits_final.to_csv('final_hits.csv', index = False)

#### load the csv with 16S sequences


seqs = pd.read_csv('./sup_info/curated_taxonomy_and_seq.csv')

seqs.head()

# how many different species
species = list(set(hits_final['species']))

seqs_hit = seqs[seqs.Species.isin(species)]

len(set(seqs_hit['Species']))

seqs_hit.to_csv('seqs_hits_extended.csv',index=False)


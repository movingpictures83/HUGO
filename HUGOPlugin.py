#!/usr/bin/env python
# coding: utf-8

# # PRADA

# ## Prioritization of Regulatory Pathways based on Analysis of RNA Dynamics Alterations

# Dysregulation of RNA stability plays an important role in cancer progression. Key regulators of RNA turnover, such as miRNAs and RNA-binding proteins, have been implicated in a variety of cancers - however, the list of annotated regulatory programs that govern the RNA lifecycle remains incomplete. The development of analytical frameworks for systematic discovery of post-transcriptional regulators is critical for a better understanding of regulatory networks that impact disease progression. For this purpose, we have developed a computational framework, named PRADA, to identify RNA-binding proteins that underlie pathologic gene expression modulations. Using this approach, we uncovered the RNA-binding protein RBMS1 as a novel suppressor of colon cancer progression. Our findings indicate that silencing RBMS1, which is achieved through epigenetic reprogramming, results in increased metastatic capacity in colon cancer cells. Restoring RBMS1 expression, in turn, blunts metastatic capacity. We have shown that RBMS1 functions as a post-transcriptional regulator of RNA stability by binding and stabilizing ~80 target mRNAs. Importantly, our analyses of colon cancer datasets as well as measurements in clinical samples have shown that RBMS1 silencing is associated with disease progression and poor survival. Our findings establish a previously unknown role for RBMS1 in mammalian gene expression regulation and its role in colon cancer metastasis.

# ### 1. Create and clean up input file (log-fold change)

# In this study, we are starting with Illumina arrays from GSE59857, which compares poorly and highly metastatic colon cancer cell lines.

# In[1]:


import sys
import os
import pandas as pd
import re
import numpy as np
import scipy as sp
from collections import defaultdict
from itertools import islice
os.environ['KMP_DUPLICATE_LIB_OK']='True'


import PyPluMA
import PyIO

class HUGOPlugin:
 def input(self, inputfile):
     self.parameters = PyIO.readParameters(inputfile)
 def run(self):
     pass
 def output(self, outputfile):
  exp = pd.read_csv(PyPluMA.prefix()+"/"+self.parameters["RDfile"], sep='\t', header=0, index_col=0)

  motifs = pd.read_csv(PyPluMA.prefix()+"/"+self.parameters["motifmap"], sep='\t', header=0)
  motifs.set_index('RBP', inplace=True, drop=False)


  hgnc_to_ref = pd.read_csv(PyPluMA.prefix()+"/"+self.parameters["hg19"], sep='\t', header=0)
  hgnc_to_ref = hgnc_to_ref.groupby('HGNC')['RefSeq'].apply(lambda x: "%s" % ','.join(x))
  mat = pd.read_csv(PyPluMA.prefix()+"/"+self.parameters["refseq"], sep='\t', header=0, index_col=0)

  for rbp in mat.columns:
    #print(rbp)
    rbp_refs = hgnc_to_ref[rbp].split(',')
    rbp_sum = 0
    rbp_cnt = 0
    rbp_max = 0
    for r in rbp_refs:
        if r in exp.index:
            if (abs(exp.loc[r,'logFC']) > rbp_max):
                motifs.loc[rbp,'diff'] = exp.loc[r,'logFC']
                motifs.loc[rbp,'pval'] = exp.loc[r,'pval']
                rbp_max = abs(exp.loc[r,'logFC'])

  motifs.to_csv(outputfile+'/RBP_motif_diff.txt', sep='\t', index=True, index_label='RefSeq')
  for rbp in mat.columns:
    mat[rbp] = mat[rbp]*motifs.loc[rbp,'diff']
  mat.to_csv(outputfile+'/RBP-v-RefSeq_target_matrix_dExp.txt', sep='\t', index=True, index_label='RefSeq')
  penalties = pd.DataFrame(index=motifs.index)
  penalties['penalties'] = motifs['diff'].apply(lambda x: 1/abs(x))
  penalties.to_csv(outputfile+'/penalties.txt', sep='\t', index=True, index_label='RBP')

  exp_fil = pd.DataFrame(index=mat.index)
  exp_fil['logFC'] = exp.loc[mat.index,'logFC']
  print("HI")
  exp_fil.to_csv(outputfile+'/high-vs-low_metastatic_lines_GSE59857_logFC_refseq_fil.txt', sep='\t', index=True, index_label='RefSeq')





import numpy as np
from os.path import splitext, basename
import os
from os import sys, path
import pandas as pd
import csv
import urllib2
import json
import sys
from IPython.core.debugger import Tracer

if __name__=='__main__':
    gene_info_file = sys.argv[1];
    output_file = sys.argv[2];
    open(output_file, 'w').close()
    #output_file = 'gene_summary.csv';
    iter_csv = pd.read_csv(gene_info_file);
    gene_ids = pd.unique(iter_csv[iter_csv['tax_id']==9606]['gid']);
    chunk_size = 100;
    cn = len(gene_ids)/chunk_size+1
    for i in range(cn):
        chunk_genes = gene_ids[chunk_size*i:np.min([chunk_size*(i+1), len(gene_ids)])];
        gids = ','.join([str(s) for s in chunk_genes])
        print (i+1),'/',cn
        url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id=' + gids + '&retmode=json';
        print url
        data = json.load(urllib2.urlopen(url))
        result = [];
        for g in chunk_genes:
            result.append([g, data['result'][str(g)]['summary'] if str(g) in data['result'] else '']);
        pd.DataFrame(result, columns=['gene_id', 'summary']).to_csv(output_file, index=False, mode='a', header= (i==0))               

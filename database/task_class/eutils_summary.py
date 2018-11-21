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
import util

class EUtilsSummary(object):
    @staticmethod
    def get_gene_info(gene_info_file, output_file, taxid_list):

        open(output_file, 'w').close()
        #print gene_info_file
        #print output_file
        iter_csv = util.read_csv(gene_info_file, skiprows=1, names=["tax_id","gid","Symbol","LocusTag","Synonyms","dbXrefs","chromosome","map_location","description","type_of_gene","Symbol_from_nomenclature_authority","Full_name_from_nomenclature_authority","Nomenclature_status","Other_designations","Modification_date","Feature_type"], sep="\t");
        iter_csv = iter_csv.query('tax_id in [%s]'%','.join(taxid_list))
        gene_ids = pd.unique(iter_csv['gid']);
        print 'Downloading gene summary for %d genes of %d species. Will take a long time.'%(len(gene_ids), len(taxid_list))
        gid2taxid = {iter_csv.iloc[i]['gid']:iter_csv.iloc[i]['tax_id'] for i in range(len(iter_csv))}
        chunk_size = 200;
        cn = len(gene_ids)/chunk_size+1
        for i in range(cn):
            chunk_genes = gene_ids[chunk_size*i:np.min([chunk_size*(i+1), len(gene_ids)])];
            gids = ','.join([str(s) for s in chunk_genes])
            print (i+1),'/',cn
            url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id=' + gids + '&retmode=json';
            print url
            data = json.load(urllib2.urlopen(url))
            #Tracer()()    
            result = [];
            for g in chunk_genes:
                chr_content=[]
                summary=''
                if str(g) in data['result'] and 'genomicinfo' in data['result'][str(g)] and len(data['result'][str(g)]['genomicinfo'])>0:
                    for loc in data['result'][str(g)]['genomicinfo']:
                        try:
                            start = int(str(loc['chrstart'])) +1
                            end = int(str(loc['chrstop'])) + 1
                            if start > end:
                                t =start
                                start = end
                                end = t
                        except:
                            print '%s is not a coordinate of a gene'%str(loc['chrstart']);
                            continue

                        chr_content.append(''.join(['chr',loc['chrloc'], ':', str(start), '-', str(end)]));

                if str(g) in data['result']:
                    summary=data['result'][str(g)]['summary']
                
                #Tracer()()
                if len(chr_content)>0:
                    result.append([g, ';'.join(chr_content), 'genome_location', gid2taxid[g]])

                if len(summary)>0:
                    result.append([g, summary, 'gene_summary', gid2taxid[g]])
                    
            pd.DataFrame(result, columns=['gid', 'content', 'type_name', 'tax_id']).to_csv(output_file, index=False, mode='a', header= (i==0))               

if __name__=='__main__':

    gene_info_file = sys.argv[1];
    output_file = sys.argv[2];
    EUtilsSummary.get_gene_info(gene_info_file,output_file);
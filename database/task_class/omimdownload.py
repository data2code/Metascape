#!/usr/bin/env python
#from .core import *
import numpy as np
import shutil
import urllib
import urlparse
from os.path import splitext, basename
import os
from os import sys, path
from core import *
import util
import pandas as pd
import csv
import re
from pprint import pprint

from IPython.core.debugger import Tracer

#from .. import SyncDB
class OmimDownload(XmlClass):
    def __init__(self, xe=None):
        XmlClass.__init__(self,xe=xe)
        self.fn_dest_omim_term =os.path.join(SyncDB.DOWNLOAD_DIR(),"omim_term.csv")
        self.fn_dest_omim2gene =os.path.join(SyncDB.DOWNLOAD_DIR(),"omim2gene.csv")
        self.fn_dest_omim_annotations =os.path.join(SyncDB.DOWNLOAD_DIR(),"omim_annotation.csv")         
        self.inputs=['ds:omim']
        
    def populate_more(self,root):
        self.outputs = [self.fn_dest_omim_term, self.fn_dest_omim2gene, self.fn_dest_omim_annotations]

    def get_omim_term(self):
        #this file need account, mannually copy now.
        #urllib.urlretrieve("ftp://anonymous:blah%40blah.org@ftp.omim.org/OMIM/genemap2.txt", self.fn_dest_omim_term + '.tmp' )
        df_omim = util.read_csv(os.path.join(SyncDB.DOWNLOAD_DIR(),"geneMap2.txt"), sep=r'\t', comment='#',
                names=['Chromosome','Genomic Position Start','Genomic Position End','Cyto Location','Computed Cyto Location',
                       'Mim Number','Gene Symbols','Gene Name','Approved Symbol','Entrez Gene ID','Ensembl Gene ID',
                       'Comments','Phenotypes','Mouse Gene Symbol/ID'])
        util.rename2(df_omim,{'Mim Number':'mim_num',
                              'Gene Name':'title',
                              'Phenotypes':'disorders',
                              'Comments':'comments'
                            })
        # df_omim = df_omim[df_omim['Mim Number']==155600]
        # pprint(df_omim[:1].to_dict())
        df_omim[['mim_num', 'title', 'disorders','comments']].to_csv(self.fn_dest_omim_term, sep=',', index=False)


    def get_omim2gene(self):
        #urllib.urlretrieve("ftp://anonymous:blah%40blah.org@ftp.omim.org/OMIM/mim2gene.txt", self.fn_dest_omim2gene + '.tmp')
        df_term = util.read_csv(os.path.join(SyncDB.DOWNLOAD_DIR(),"mim2gene.txt"),
                                sep=r'\t',
                                comment='#',
                                names=['MIM Number','Type','Entrez Gene ID (NCBI)',
                                       'Approved Gene Symbol (HGNC)','Ensembl Gene ID (Ensembl)'])
        util.rename2(df_term,{'MIM Number':'term_id',
                              'Entrez Gene ID (NCBI)':'gid',
                            })
        df_term = util.df2sdf(df_term,s_format='%.0f', l_all_str=True)
        df_term = df_term[df_term['Type'].str.contains('gene')]
        df_term = df_term[df_term['gid']!='']
        df_term[['term_id', 'gid']].to_csv(self.fn_dest_omim2gene, sep=',', index=False)

    def disorder_map (self, v):
        v=str(v)
        if v == 'nan':
            return '';
        ret=[];
        for s in str.split(v, ';'):            
            m=re.search('(\s*,\s*\d+\s*)*\(\d+\)$|(\s*,\s*\d+$)', s);
            if (m!=None):
                ret.append(re.sub("\?|\[|\]|\{|\}",'',s[:-len(m.group(0))]));
            else:
                ret.append(re.sub("\?|\[|\]|\{|\}",'',s));
        return '|' + ';'.join(ret)
        
    def value_map(self, v):
        v=str(v)
        return v if v != 'nan' else '';
        
    def get_omim_annotations (self):

        df_omim_terms = util.read_csv(self.fn_dest_omim_term, skiprows = 1, names=['mim_num', 'title', 'disorders', 'comments']);
        df_omim2gene = util.read_csv(self.fn_dest_omim2gene, skiprows = 1, 
                names=['mim_num', 'gid']);
        #Tracer()()
        df_join = pd.merge(df_omim2gene, df_omim_terms, left_on='mim_num', right_on='mim_num', how='left')
        df_join['content'] = "OMIM:" + df_join['mim_num'].map(str) + df_join['disorders'].map(self.disorder_map);
        df_join['annotation_field1'] = df_join['disorders'].map(self.value_map) \
                                       + '|' + df_join['comments'].map(self.value_map)
        df_join['tax_id'] = '9606'
        df_join[['gid','content','annotation_field1', 'tax_id']].to_csv(self.fn_dest_omim_annotations, index=False);
        
        
    def do_update(self):
        self.get_omim2gene()
        self.get_omim_term()
        self.get_omim2gene()
        self.get_omim_annotations()


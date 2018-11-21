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

from IPython.core.debugger import Tracer

#from .. import SyncDB
from gputil import GPUtils


class JaxDownload(XmlClass):
    def __init__(self, xe=None):
        XmlClass.__init__(self,xe=xe)
        self.fn_dest_gene2phenotype =os.path.join(SyncDB.DOWNLOAD_DIR(),"gene2phenotype.csv")
        self.fn_dest_reference =os.path.join(SyncDB.DOWNLOAD_DIR(),"mgi_marker_reference.csv")
        self.fn_dest_jax_annotations =os.path.join(SyncDB.DOWNLOAD_DIR(),"jax_annotation.csv")         
        self.inputs=['ds:jax']
        
    def populate_more(self,root):
        self.outputs = [self.fn_dest_jax_annotations ]
  
    def get_jax_annotations (self):    
        #urllib.urlretrieve("ftp://ftp.informatics.jax.org/pub/reports/HMD_HumanPhenotype.rpt", self.fn_dest_gene2phenotype)
        #urllib.urlretrieve("ftp://ftp.informatics.jax.org/pub/reports/VOC_MammalianPhenotype.rpt", self.fn_dest_reference)
        urllib.urlretrieve('http://www.informatics.jax.org/downloads/reports/HMD_HumanPhenotype.rpt', self.fn_dest_gene2phenotype)
        urllib.urlretrieve('http://www.informatics.jax.org/downloads/reports/VOC_MammalianPhenotype.rpt', self.fn_dest_reference)

        df_gene2phenotype = util.read_csv(self.fn_dest_gene2phenotype, names=['human_symbol', 'gid', 'homolo_gid', 'yes_no','mouse_marker', 'mgi_marker', 'phenotype_ids'], sep=r'\t', index_col=False);
        df_mgi_reference = util.read_csv(self.fn_dest_reference, names=['phenotype_id', 'name', 'description'], sep=r'\t', index_col=False);
        df_gene2phenotype['mgi_marker']= df_gene2phenotype['mgi_marker'].map(str.strip)
        df_mgi_reference['phenotype_id']= df_mgi_reference['phenotype_id'].map(str.strip)
        
        data = [];
        for index,r in df_gene2phenotype.iterrows():
            if r['phenotype_ids']:
                for pid in r['phenotype_ids'].split(' '):
                    data.append({'gid':r['gid'], 'phenotype_id':pid})
                
        df_gene2phenotype = pd.DataFrame(data); 
        df_join = pd.merge(df_gene2phenotype, df_mgi_reference, left_on='phenotype_id', right_on='phenotype_id', how='inner')
        
        data=[]
        for k,g in df_join.groupby('gid', as_index=False):
            data.append({'gid':k, 'content':'; '.join(g['name']), 'annotation_field1':'; '.join([x for x in g['description'] if x is not None]), 'tax_id':'9606'})
        
        pd.DataFrame(data).to_csv(self.fn_dest_jax_annotations, index=False);
               
    def do_update(self):
        self.get_jax_annotations()
        
    def check_inputs(self):
        passed = True
        urls=[]
        print "Checking urls for JAX..."
            
        urls.append("ftp://ftp.informatics.jax.org/pub/reports/HMD_HumanPhenotype.rpt")        
        urls.append("ftp://ftp.informatics.jax.org/pub/reports/VOC_MammalianPhenotype.rpt")

        for url in urls:
            if not GPUtils.check_url(url):
                passed = False;                
                
        return passed                

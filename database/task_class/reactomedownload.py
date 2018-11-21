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
from pprint import  pprint
from IPython.core.debugger import Tracer

#from .. import SyncDB
from gputil import GPUtils


class ReactomeDownload(XmlClass):
    TAX_NAME2ID = {
                 #'Bos taurus': 9913,
                 'Caenorhabditis elegans': 6239,
                 'Danio rerio': 7955,
                 'Drosophila melanogaster': 7227,
                 #'Gallus gallus': 9031,
                 'Homo sapiens': 9606,
                 'Mus musculus': 10090,
                 'Rattus norvegicus': 10116,
                 'Saccharomyces cerevisiae': 4932,
                 #'Xenopus tropicalis': 8364,
                 #'Canis familiaris': -1,
                }
    def __init__(self, xe=None):
        XmlClass.__init__(self,xe=xe)
        self.fn_dest_go_term =os.path.join(SyncDB.DOWNLOAD_DIR(),"reactome_term.csv")
        self.fn_gene_term_pair = os.path.join(SyncDB.DOWNLOAD_DIR(),"reactome_gene_term_pair.csv")
        self.inputs=['ds:reactome']
        

    def populate_more(self,root):
        self.outputs = [self.fn_dest_go_term,
                        self.fn_gene_term_pair]

    def get_go_term_and_go_term2term(self):
        fn_source =os.path.join(SyncDB.DOWNLOAD_DIR(),"NCBI2Reactome_All_Levels.txt")
        urllib.urlretrieve("http://www.reactome.org/download/current/NCBI2Reactome_All_Levels.txt",fn_source)
        df = pd.read_csv(fn_source,sep='\t',names=['gene_id',
                                                   'term_id',
                                                   'url',
                                                   'term_name',
                                                   'evident',
                                                   'species'])
        df.drop(['url', 'evident'], axis=1, inplace=True)
        df = df[ df['species'].isin(ReactomeDownload.TAX_NAME2ID)]
        df[ 'species'] = [ ReactomeDownload.TAX_NAME2ID[x] for x in df['species'] ]
        #term
        df_term = df[['term_id','term_name','species']]
        df_term = df_term.drop_duplicates()
        df_term.to_csv(self.fn_dest_go_term, index=False)
        if False: # to get TAX_NAME2ID
            mydb = db.DB(name=SyncDB.CONNECTION_ID, db=SyncDB.DATABASE)
            df_tmp = mydb.sql_in('''
                            SELECT t.tax_id, t.scientific_name
                            FROM taxid2name t
                            WHERE t.scientific_name IN (''',
                                   ")",
                                   util.unique(df['species']))
            pprint( { r['scientific_name']: r['tax_id'] for i, r in df_tmp.iterrows()})
        #term gene pair
        df.drop(['term_name'], axis=1, inplace=True)
        df = df.drop_duplicates()
        util.rename2(df,{'species':'tax_id',
                      'gene_id':'gid',
                      'term_id':'term_id'})
        df.to_csv(self.fn_gene_term_pair, index=False)

        return


    def do_update(self):
        self.get_go_term_and_go_term2term()    

    def check_inputs(self):
        passed = True
        urls=[]
        print "Checking urls for godownload..."
        urls.append("http://www.reactome.org/download/current/NCBI2Reactome_All_Levels.txt")
        for url in urls:
            if not GPUtils.check_url(url):
                passed = False;                
                
        return passed        
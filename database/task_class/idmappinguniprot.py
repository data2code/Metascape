#!/usr/bin/env python
#from .core import *
import numpy as np
import pandas as pd
import shutil
import urllib
import urlparse
from os.path import splitext, basename
import os
from os import sys, path
import util
from pprint import pprint
import StringIO
import db
from gp import *
from core import *
from idmapping import IdMapping

class IdMappingUniprot(IdMapping):
    def __init__(self, xe):
        #add default col instance to children
        IdMapping.__init__(self,xe=xe)

    def do_update(self):

        #gid,id_type_id,source_id,id_status,ds,tax_id
        IdMapping.do_update(self)
        fn_isoform_idmap = os.path.join(SyncDB.DOWNLOAD_DIR(),self.dependence)
        #annotation_field1,content,gid,tax_id
        df_ensembl = pd.read_csv(fn_isoform_idmap, dtype= str)
        df_ensembl = df_ensembl.query(self.filter_str)
        rows = []
        for i, r in df_ensembl.iterrows():
            lines = r['content'].split('|')
            for line in lines:
                id = {x.split(':')[0]:x.split(':')[1] for x in line.split(',')}['Uniprot_ID']
                if id == '-':
                    continue
                rows.append({
                            'source_id':id,
                             'gid':r['gid'],
                             'tax_id':r['tax_id']
                             })
        df_ensembl = pd.DataFrame(rows)
        df_ensembl[self.type_col ] =  self.type_col_value
        df_ensembl['ds'] = 'ensembl'
        df_ensembl['id_status'] = 'None'
        df_ensembl = df_ensembl.reindex(columns=['gid','id_type_id','source_id','id_status','ds','tax_id'],copy=False)
        df_uniprot = pd.read_csv(self.fn_dest,dtype = str)
        df_all = pd.concat([df_uniprot,df_ensembl])
        df_all = df_all.drop_duplicates(subset=['gid','id_type_id','source_id','id_status','tax_id'])
        df_all.to_csv(self.fn_dest,index=False)





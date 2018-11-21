#!/usr/bin/env python
import numpy as np
import pandas as pd
import shutil
import urllib
import urlparse
from os.path import splitext, basename
import os
from os import sys, path
import util
from core import *
from pprint import pprint
import StringIO
import db
import csv

from IPython.core.debugger import Tracer

class CustomGeneListHelper(UploadCsvConvert):
    def __init__(self, xe,dest,key_col,column_map,base_obj):
        self.base_obj = base_obj
        UploadCsvConvert.__init__(self,xe=xe,dest=dest)
        self.type_col = 'term_category_id'
        self.key_col = key_col
        self.column_map = column_map
        if 'termPrefix' in xe.attrib:
            self.term_prefix = xe.attrib['termPrefix']
        if(key_col == 'gid'):
            self.value_col = 'term_id'
            self.new_value_col = 'term_ids'
        else:
            self.value_col = 'gid'
            self.new_value_col = 'gids'

        self.new_cols = [self.key_col,self.new_value_col,'id_count','ds','tax_id']
        if 'destFile' in xe.attrib:    
            self.use_dest =xe.attrib['destFile'];

    def get_type_col_value_sql(self):
        return 'SELECT t.term_category_id FROM %s.term_category t WHERE t.category_name = ?' % SyncDB.DATABASE

    def get_term_prefix(self,s):
        if s=='term_id':
            if hasattr(self,'term_prefix'):
                return self.term_prefix
            return self.get_type_col_value() + '_'
        return ''

    def get_chunksize(self):
      return None

    def do_one_chunk(self,chunk):
        rows=[]
        data_key_col = self.key_col if self.key_col not in self.column_map else self.column_map[self.key_col];
        data_value_col = self.value_col if self.value_col not in self.column_map else self.column_map[self.value_col];
        tax_id_col = 'tax_id' if 'tax_id' not in self.column_map else self.column_map['tax_id'];
        
        for k,g in chunk.groupby(data_key_col, as_index=False):
            #Tracer()()
            r = {}
            r[self.key_col] = self.get_term_prefix(self.key_col)+str(k)
            allids = self.get_term_prefix(self.value_col) +  g[data_value_col].astype(str)
            if self.value_col == 'gid':
                allids = [x for x in util.unique(allids) if str(x).isdigit()]
            else:
                allids = util.unique(allids)
            r[self.new_value_col] = ','.join(allids)
            r['id_count'] = len(allids);
            r['ds'] = self.ds
            if hasattr(g.iloc[0],tax_id_col):
                all_tax_ids = g[tax_id_col].astype(str)
                r['tax_id']=','.join(util.unique(all_tax_ids))
            else:
                r['tax_id'] = '9606'
            rows.append(r)
                
        return rows
        
    def do_update(self):
        UploadCsvConvert.do_update(self)
        if self.base_obj:
            self.base_obj.create_term_file()


class CustomGeneList(XmlClass):
    def __init__(self, xe):
        XmlClass.__init__(self,xe)
        #Tracer()()
        self.fn_source = os.path.join(SyncDB.DOWNLOAD_DIR(),self.options['source'])
        self.fn_dest = os.path.join(SyncDB.DOWNLOAD_DIR(),self.options['dest'])

        if 'termName' in xe.attrib:
            self.term_name = xe.attrib['termName']
        else: 
            self.term_name = None
            
        self.column_map = {};
        for c in self.children:
            if 'col.Col' in str(c):
                self.column_map[c.options['colName']] = c.options['sourceCol']        
        self.gid2terms = CustomGeneListHelper(xe,dest='gid2terms',key_col='gid', column_map=self.column_map, base_obj=None)
        self.term2gids = CustomGeneListHelper(xe,dest='term2gids',key_col='term_id', column_map=self.column_map, base_obj=self)
        self.term2gids.outputs = [self.fn_dest]                
        self.children = [self.gid2terms,self.term2gids]
        
      
    def create_term_file (self):
        kwargs={}
        if 'oldCols' in self.options:
            kwargs['names']= self.options['oldCols'].split(',')
        if 'read_csv' in self.options:
            for kv_str in  self.options['read_csv'].split(','):
                kv = kv_str.split('=')
                kwargs[kv[0]] = kv[1]
                if kv[1] == 'None':
                    kwargs[kv[0]] = None
                if kv[0].lower() == 'skiprows':
                    kwargs[kv[0]] = int(kv[1])

        iter_csv = util.read_csv(self.fn_source, iterator=True, chunksize=self.get_chunksize(),dtype=str,**kwargs)
        term_id_col = 'term_id' if 'term_id' not in self.column_map else self.column_map['term_id'];
        term_ids = [];
        for chunk in iter_csv:
            term_ids += util.unique(chunk[term_id_col]);            
        term_ids = util.unique(term_ids);

        with open(self.fn_dest, "w") as myfile:
            wr = csv.writer(myfile)
            wr.writerow(['term_id','term_name','term_type']);
            wr.writerows([[term_id, self.term_name if self.term_name else term_id, self.options['typeName']] for term_id in term_ids]);
            
    def get_xmlclass_instances(self):
        r = list(self.children)
        return r
        
    def get_chunksize(self):
      return None

    def do_update(self):
        self.create_term_file();

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
from IPython.core.debugger import Tracer

class TermGenePairHelper(UploadCsvConvert):
    def __init__(self, xe,dest,key_col):
        UploadCsvConvert.__init__(self,xe=xe,dest=dest)
        self.type_col = 'term_category_id'
        self.key_col = key_col
        if 'termPrefix' in xe.attrib:
            self.term_prefix = xe.attrib['termPrefix']
        if(key_col == 'gid'):
            self.value_col = 'term_id'
            self.new_value_col = 'term_ids'
        else:
            self.value_col = 'gid'
            self.new_value_col = 'gids'

        self.new_cols = [self.key_col,self.new_value_col,'id_count','ds','tax_id','term_category_id']

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
        if len(chunk) == 0:
            return rows

        key=[self.key_col]
        if self.key_col == 'term_id' and hasattr(chunk.iloc[0],'tax_id'):
            key.append('tax_id');
            
        for k,g in chunk.groupby(key, as_index=False):

            row = {}
            if self.key_col == 'gid':
                row[self.key_col] = k
            else:
                row[self.key_col] = self.get_term_prefix(self.key_col)+k[0]
            allids = self.get_term_prefix(self.value_col) +  g[self.value_col].astype(str)
            allids = util.unique(allids)
            row[self.new_value_col] = ','.join(allids)
            row['id_count'] = len(allids);
            row['ds'] = self.ds
            if hasattr(g.iloc[0],'tax_id'):
                row['tax_id']=g.iloc[0]['tax_id']
            else:
                row['tax_id'] = None
            row['term_category_id']=self.get_type_col_value()
            rows.append(row)
        return rows


class TermGenePair(XmlClass):
    def __init__(self, xe):
        XmlClass.__init__(self,xe)
        self.gid2terms = TermGenePairHelper(xe,dest='gid2terms',key_col='gid')
        self.term2gids = TermGenePairHelper(xe,dest='term2gids',key_col='term_id')
        self.children = [self.term2gids,self.gid2terms]

    def get_xmlclass_instances(self):
        r = list(self.children)
        return r



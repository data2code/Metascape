#!/usr/bin/env python
#from .core import *
import numpy as np
import pandas as pd
import shutil
import urllib
import urlparse
from os.path import splitext, basename
import os
import util
from pprint import pprint
import StringIO
import db
from core import *
from IPython.core.debugger import Tracer


class Term(UploadCsvConvert):
    def __init__(self, xe):
        #add default col instance to children
        xe.attrib['newCols'] = 'term_id,term_source_id,term_name,term_category_id,description,term_field1,term_field2'
        defined_cols = {x.attrib['colName'] for x in xe}
        UploadCsvConvert.__init__(self,xe=xe,dest='term')
        for c in self.children:
            if c.col_name == 'term_source_id' and ( 'term_source_id' not in defined_cols):
                c.source_col = 'term_id'
        self.type_col = 'term_category_id'
        if 'termPrefix' in xe.attrib:
            self.term_prefix = xe.attrib['termPrefix']
            
        if 'destFile' in xe.attrib:    
            self.use_dest =xe.attrib['destFile'];

    def get_type_col_value_sql(self):
        return 'SELECT t.term_category_id FROM %s.term_category t WHERE t.category_name = ?' % SyncDB.DATABASE

    def get_term_prefix(self):
        if hasattr(self,'term_prefix'):
            return self.term_prefix
        return self.get_type_col_value() + '_'

    def generate_new_row(self,row):
        r = super(Term,self).generate_new_row(row)
        #if the Species is supported
        if r is not None:
            r['term_id'] = self.get_term_prefix() + r['term_id']
        return r


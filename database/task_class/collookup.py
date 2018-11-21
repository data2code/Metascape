#!/usr/bin/env python
#from .core import *
import numpy as np
import pandas as pd
import shutil
import urllib
import urlparse
from os.path import splitext, basename
import os
from core import *
from IPython.core.debugger import Tracer

class ColLookup(Col):
    ALL_LOOKUPS = {}
    def __init__(self, xe):
        Col.__init__(self,xe=xe)
        self.lookup = xe.attrib['lookup']        
        self.lookup_name = xe.attrib['lookupName']
        self.lookup_key = xe.attrib['lookupKey']
        self.lookup_value = xe.attrib['lookupValue']
        self.lookup_table = None
    def get_lookup_table(self):
        if not self.lookup_name in ColLookup.ALL_LOOKUPS:
            fn = os.path.join(SyncDB.DOWNLOAD_DIR(),self.lookup)
            df = pd.read_csv(fn, dtype=str)
            kv = dict(zip(df[self.lookup_key],df[self.lookup_value]))
            ColLookup.ALL_LOOKUPS[self.lookup_name] = kv
        return ColLookup.ALL_LOOKUPS[self.lookup_name]


    def do_work(self,row=None):
        if self.lookup_table is None:
            self.lookup_table = self.get_lookup_table()
        k = str(row[self.source_col])
        return self.lookup_table.get(k, None)

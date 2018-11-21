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
from IPython.core.debugger import Tracer

#from .. import SyncDB
class IdMappingGeneSynonyms(CsvConvert):
    def __init__(self, xe=None):
        CsvConvert.__init__(self,xe=xe)

    def do_post_update(self):
		#The synonyms are concatenated with '|' for each gene. This post_update creates a separate record for each synonym.
        df = util.read_csv(self.fn_dest,dtype='S300')
        h = df.header()
        gid_col_name = h[0]
        synonyms_col_name = h[1]
        tax_id = h[2]
        rows = []
        for i, r in df.iterrows():
            more_rows = [ {gid_col_name: r[gid_col_name], synonyms_col_name:x, tax_id:r[tax_id]} for x in str(r[synonyms_col_name]).split('|')]
            rows.extend(more_rows)
        df = pd.DataFrame(rows)
        df.to_csv(self.fn_dest, index=False)



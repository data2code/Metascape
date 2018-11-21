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
class IdMappingGeneDbXref(CsvConvert):
    def __init__(self, xe=None):
        CsvConvert.__init__(self,xe=xe)

    def do_post_update(self):
		#The synonyms are concatenated with '|' for each gene. This post_update creates a separate record for each synonym.
        df = util.read_csv(self.fn_dest,dtype='S300')
        h = df.header()
        gid_col_name = h[0]
        synonyms_col_name = h[1]
        tax_id_col_name = h[2]
        rows = []
        for i, r in df.iterrows():
            for x in str(r[synonyms_col_name]).split('|'):
                import re
                m = re.search('([^:]*):(.*)',x)
                if m and len(m.groups())==2:
                    db = m.group(1)
                    xref = m.group(2)
                    tax_id  = r[tax_id_col_name]
                    if db in ['MGI', 'WormBase','SGD','FLYBASE','ZFIN','RGD','Araport']:
                        rows.append({
                            gid_col_name: r[gid_col_name],
                            synonyms_col_name:xref,
                            tax_id_col_name:tax_id,
                            'id_status':db
                           })
        df = pd.DataFrame(rows)
        df.to_csv(self.fn_dest, index=False)


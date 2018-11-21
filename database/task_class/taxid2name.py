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


class Taxid2Name(UploadCsvConvert):
    def __init__(self, xe):
        #add default col instance to children
        xe.attrib['newCols'] = 'tax_id,scientific_name,common_name'
        self.xe = xe
        self.dest = 'taxid2name'
        self.id = 'Taxid2Name'
        self.ds = 'NCBI'
        self.taxidList = SyncDB.SUPPORTED_SPECIES;
        self.process()
        UploadCsvConvert.__init__(self, xe=self.xe, dest='taxid2name')

    def populate_more(self,root):
        self.inputs=[]
        self.outputs= [self.dest]

    def process(self):
        if not os.path.exists(os.path.join(SyncDB.DOWNLOAD_DIR(), "names.dmp")):
            urllib.urlretrieve("ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip",os.path.join(SyncDB.DOWNLOAD_DIR(),"taxdmp.zip") )
            util.unix("unzip %s " % os.path.join(SyncDB.DOWNLOAD_DIR(),"taxdmp.zip"))

        df=util.read_csv(os.path.join(SyncDB.DOWNLOAD_DIR(), "names.dmp"), names=['tax_id', 'name', 'uname', 'type_name'], sep="\s*\|\s*", index_col=False).query('type_name in ["scientific name", "genbank common name"]')
        df=df.query('tax_id in [%s]'%','.join(['"'+x+'"' for x in self.taxidList]))
        df=df.sort('tax_id')
        rows = []
        last_taxid=None
        row = None
        for e in df.iterrows():

            if e[1]['tax_id'] != last_taxid:
                if row is not None and 'scientific_name' in row:
                    row['tax_id']=last_taxid
                    if 'common_name' not in row:
                        row['common_name'] = row['scientific_name']
                    rows.append(row)
                row={}
                last_taxid=e[1]['tax_id']

            if 'scientific name' in e[1]['type_name'].strip():
                row['scientific_name'] = e[1]['name'].strip()
            if 'genbank common name' in e[1]['type_name'].strip():
                row['common_name'] = e[1]['name'].strip()

        if row is not None and 'scientific_name' in row:
            rows.append(row)

        if not os.path.exists(os.path.join(SyncDB.UPLOAD_DIR(),self.dest)):
            util.unix('mkdir -p %s'%os.path.join(SyncDB.UPLOAD_DIR(),self.dest))
        pd.DataFrame(rows).to_csv(os.path.join(SyncDB.UPLOAD_DIR(),self.dest) + "/taxid2name.csv", index=False, header=True, columns=['tax_id', 'scientific_name','common_name'])

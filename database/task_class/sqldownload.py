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
import math
#from .. import SyncDB
class SqlDownload(XmlClass):
    def __init__(self, xe=None):
        XmlClass.__init__(self,xe=xe)
        self.query= xe.attrib['query']
        self.fn_dest=os.path.join(SyncDB.DOWNLOAD_DIR(),xe.attrib['dest'])
        if xe.attrib['connection'] == 'ENSEMBL':
            self.con = self.get_ensembl_con()
        else:
            self.con = db.get_con(xe.attrib['connection'])
       
    def get_ensembl_con(self):
        import MySQLdb as mysql
        return mysql.connect('ensembldb.ensembl.org','anonymous','','ensembl_ontology_77')

    def do_work(self):
        df = db.from_sql(self.con,self.query)
        df.to_csv(self.fn_dest, index=False)

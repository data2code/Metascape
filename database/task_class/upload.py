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
from pprint import pprint
import StringIO
import db
import itertools
from core import *
from IPython.core.debugger import Tracer


#print __package__
#from ..core import SyncDB, XmlClass
class Upload(XmlClass):
    def __init__(self, xe):
        XmlClass.__init__(self,xe=xe)
        self.dest = xe.attrib['table']
        self.cols = xe.attrib['cols']

    def get_fn_dest(self):
        return os.path.join(SyncDB.UPLOAD_DIR(),"%s.csv" % self.dest)

    def do_update(self):
        dir = os.path.join(SyncDB.UPLOAD_DIR(),self.dest)
        cmd = "find %s -name '*.csv' -print0  |xargs -0 tail -n +2 -q > %s" % (dir, self.fn_dest)
        print cmd
        s = util.unix(cmd)
        conn_info = db.get_con_info(SyncDB.CONNECTION_ID)
        if conn_info is None:
            print "cannot find connection for " + SyncDB.CONNECTION_ID
            exit()

        #Tracer()()
        print s
        cols = self.cols.split(',')
        if 'tax_id' in cols:
            df = util.read_csv(self.fn_dest, names=self.cols.split(','))
            df['tax_id'] = df['tax_id'].map(str)
            f = lambda x: SyncDB.PARENT_SPECIE_MAP[x] if x in SyncDB.PARENT_SPECIE_MAP else x
            df['tax_id'] = df['tax_id'].map(f)
            df.to_csv(self.fn_dest, index=False, header=False)

        #for interaction table.
        if 'tax_id_A' in cols:
            df = util.read_csv(self.fn_dest, names=self.cols.split(','))
            df['tax_id_A'] = df['tax_id_A'].map(str)
            df['tax_id_B'] = df['tax_id_B'].map(str)
            f = lambda x: SyncDB.PARENT_SPECIE_MAP[x] if x in SyncDB.PARENT_SPECIE_MAP else x
            df['tax_id_A'] = df['tax_id_A'].map(f)
            df['tax_id_B'] = df['tax_id_B'].map(f)
            df.to_csv(self.fn_dest, index=False, header=False)

        mysql_cmd = [
            'mysqlimport',
            '--local',
            '--delete',
            '--host ' + conn_info["HOST"],
            '-u ' + conn_info["USR"],
            '--password=' + conn_info["PWD"],
            '--fields-terminated-by=,',
            '--fields-optionally-enclosed-by="\\""',
            '--ignore-lines=0', 
            '--columns=%s' % self.cols,
            SyncDB.DATABASE,
            self.fn_dest]
        c = " ".join(mysql_cmd)
        print c
        s = util.unix(mysql_cmd)

    def get_id(self):
         return  self.tag+":"+self.dest

    def populate_more(self,root):
        XmlClass.populate_more(self,root);
        dir = os.path.join(SyncDB.UPLOAD_DIR(),self.dest)
        for c in root.all_children:
            if self != c:
                self.inputs.extend([fn for fn in c.outputs if dir in fn])


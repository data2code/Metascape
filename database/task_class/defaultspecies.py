#!/usr/bin/env python
#from .core import *
import numpy as np
import shutil
import urllib
import urlparse
from os.path import splitext, basename
import os
from os import sys, path
import re
from gp import *
from core import *
import math
from IPython.core.debugger import Tracer
from myutils import MyUtils

class DefaultSpecies(object):
    def get_specie_subtree(self, s):
        stack=[s]
        nodes=[]
        while len(stack) != 0:
            n = stack[-1]
            nodes.append(str(n))
            del stack[-1]
            if self.stree[n] is not None:
                stack += self.stree[n]['c']

        return  nodes;

    def buid_sp_tree(self):
        if not os.path.exists(SyncDB.DOWNLOAD_DIR() + "/nodes.dmp"):
            urllib.urlretrieve("ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip", SyncDB.DOWNLOAD_DIR() + "/taxdmp.zip")
            fn =  SyncDB.DOWNLOAD_DIR() + "/taxdmp.zip"
            exdir = SyncDB.DOWNLOAD_DIR()
            util.unix2('unzip -o -d {0} {1}'.format(exdir,fn), l_print=False)


        df = util.read_csv(SyncDB.DOWNLOAD_DIR() + "/nodes.dmp", sep="\s*\|\s*", index_col=False)[["1", "1.1"]]
        tree = [None]*(np.max([np.max(df['1.1']), np.max(df['1'])])+1)
        for row in df.iterrows():
            tree[row[1]['1']] = {'p':row[1]["1.1"], 'c':[]}

        for i in range(len(tree)):
            if tree[i] is not None:
                if tree[tree[i]['p']] is not None:
                    tree[tree[i]['p']]['c'].append(i)

        self.stree = tree

    def __init__(self, xe):
        if 'taxidList' in xe.attrib:
            fn = os.path.join(SyncDB.DOWNLOAD_DIR(),'syncdb_species.pickle')
            if os.path.isfile(fn):
                syncdb_species = MyUtils.load(fn)
                SyncDB.ROOT_SPECIES = syncdb_species['ROOT_SPECIES']
                SyncDB.SUPPORTED_SPECIES = syncdb_species['SUPPORTED_SPECIES']
                SyncDB.PARENT_SPECIE_MAP = syncdb_species['PARENT_SPECIE_MAP']
                SyncDB.SPECIE_SUBTREE = syncdb_species['SPECIE_SUBTREE']
            else:
                self.buid_sp_tree()
                SyncDB.ROOT_SPECIES = filter(None, [x for x in xe.attrib['taxidList'].split(',')])
                SyncDB.SUPPORTED_SPECIES = []
                SyncDB.PARENT_SPECIE_MAP = {}
                SyncDB.SPECIE_SUBTREE = {}
                for s in SyncDB.ROOT_SPECIES:
                    st=self.get_specie_subtree(int(s))
                    print "Taxonomy subtree for %s : [%s]"%(s, ','.join(st))
                    SyncDB.SPECIE_SUBTREE[s] = [];
                    for t in st:
                        SyncDB.PARENT_SPECIE_MAP[t] = s
                        SyncDB.SPECIE_SUBTREE[s].append(t);
                    SyncDB.SUPPORTED_SPECIES += st
                syncdb_species = {'ROOT_SPECIES': SyncDB.ROOT_SPECIES,'SUPPORTED_SPECIES':SyncDB.SUPPORTED_SPECIES,'PARENT_SPECIE_MAP':SyncDB.PARENT_SPECIE_MAP,'SPECIE_SUBTREE':SyncDB.SPECIE_SUBTREE}
                MyUtils.dump_object(syncdb_species,'syncdb_species',SyncDB.DOWNLOAD_DIR())

        else:
            print "taxidList is missing in DefaultSpecies tag"
            exit(-1)
            

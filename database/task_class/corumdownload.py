#!/usr/bin/env python
import numpy as np
import pandas as pd
import util
import re
import db

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
import csv
from pprint import  pprint
from IPython.core.debugger import Tracer

#from .. import SyncDB
class CorumDownload(XmlClass):

    def __init__(self, xe=None):
        XmlClass.__init__(self,xe=xe)
        self.fn_dest_go_term =os.path.join(SyncDB.DOWNLOAD_DIR(),"corum_term.csv")
        self.fn_gene_term_pair = os.path.join(SyncDB.DOWNLOAD_DIR(),"corum_gene_term_pair.csv")
        self.inputs=['ds:corum']
        self.fn_gene_info = os.path.join(SyncDB.DOWNLOAD_DIR(),xe.attrib['source'])
        self.fn_data = os.path.join(SyncDB.DOWNLOAD_DIR(),'allComplexes.txt.zip')
        self.inputs.append(self.fn_gene_info)


    def populate_more(self,root):

        self.outputs = [self.fn_dest_go_term,
                        self.fn_gene_term_pair]

    def do_update(self):
        print '##############################################################'
        download_url = 'http://mips.helmholtz-muenchen.de/corum/download/allComplexes.txt.zip'
        urllib.urlretrieve(download_url, self.fn_data)
        t=pd.read_table(self.fn_data)
        #print t.header()
        #['ComplexID', 'ComplexName', 'Organism', 'Synonyms', 'Cell line', 'subunits(UniProt IDs)', 'subunits(Entrez IDs)', 'Protein complex purification method', 'GO ID', 'GO description', 'FunCat ID', 'FunCat description', 'PubMed ID', 'subunits(Protein name)', 'subunits(Gene name)', 'subunits(Gene name syn)', 'Disease comment', 'Subunits comment', 'Complex comment', 'SWISSPROT organism']
        C_TAX={'Rat':10116, 'Human':9606, 'Mouse':10090}
        c_gene2tax = self.get_gene2tax()

        out_term=[]
        out_gids=[]
        for i in t.index:
            id=t.ix[i,'ComplexID']
            s_go=t.ix[i,'ComplexName']
            s_des=t.ix[i,'Complex comment']
            if s_des=='None':
                s_des=t.ix[i, 'GO description']
                if s_des=='None':
                    s_des=s_go
            gids=t.ix[i,'subunits(Entrez IDs)']
            gids=gids.replace(';',',')
            S=gids.split(',')
            S=[ x for x in S if x in c_gene2tax]
            S_tax=[c_gene2tax[x] for x in S]
            l_new_term=True
            for s_tax in util.unique(S_tax):
                S_gid=[ S[i] for i,x in enumerate(S_tax) if x==s_tax ]
                n=len(S_gid)
                if n<3:
                    #print "Too few proteins: ", id, S
                    continue
                if l_new_term:
                    out_term.append({'term_id':'CORUM:%d' % id, 'term_name':s_go, 'description':s_des})
                    l_new_term=False
                for gid in S_gid:
                    out_gids.append({'gid': gid,
                                     'term_id':'CORUM:%d' % id,
                                     'term_name':s_go,
                                     'type_name':'CORUM',
                                      'tax_id':s_tax})
        #gid, term_id, term_name, type_name, tax_id


        t_term=pd.DataFrame(out_term)
        t_term.to_csv(self.fn_dest_go_term, index=False)
        print "Number of Complexes: %d" % len(t_term)
        t_gids=pd.DataFrame(out_gids)
        t_gids.to_csv(self.fn_gene_term_pair, index=False)
        #print t_gids.header()
        print util.unique_count(t_gids['tax_id'].values)


        #['MINK', 'Rabbit', 'Dog', 'Bovine', 'Rat', 'Mammalia', 'Human', 'Pig', 'Mouse', 'Hamster']

    def get_gene2tax(self):
        # tax_id, gid, Symbol, LocusTag, Synonyms, map_location, description, type_of_gene, Full_name_from_nomenclature_authority, dbXrefs
        # 3702, 814566, rrn26
        dt = pd.read_csv(self.fn_gene_info)
        dt = dt[dt['tax_id'].isin([9606,10090,10116])]
        dt = dt[['gid','tax_id']]
        #dt.to_csv(os.path.join(SyncDB.DOWNLOAD_DIR(),'k.csv'),index=False)
        gid = dt['gid'].astype(dtype=str)
        tax_id = dt['tax_id'].astype(dtype=np.int32)
        gene2tax = dict(zip(gid,tax_id))

        return gene2tax

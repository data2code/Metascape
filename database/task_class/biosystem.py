#!\usr\bin\env python

import numpy as np
import pandas as pd
import shutil
import urllib
import urlparse
from os.path import splitext, basename
import os
from os import sys, path
from core import *
import util
import csv
from gputil import GPUtils
from lxml import etree

from IPython.core.debugger import Tracer

class BioSystem(XmlClass):

    def __init__(self, xe=None):
        XmlClass.__init__(self,xe=xe)
        self.tag = "BioSystem"
        self.fn_terms=os.path.join(SyncDB.DOWNLOAD_DIR(),"biosystem_terms.csv")
        self.fn_term_gene_pair=os.path.join(SyncDB.DOWNLOAD_DIR(),"biosystem_term_gene_pair.csv")
        
    def populate_more(self,root):
        self.outputs = [self.fn_terms,self.fn_term_gene_pair]
        self.inputs = [os.path.join(SyncDB.DOWNLOAD_DIR(),"bsid2info.gz"), os.path.join(SyncDB.DOWNLOAD_DIR(),"biosystems_gene.gz")]
   
    def do_update(self):
        #bsid,source,source_id,name,type,tax_id,description
        taxidList = SyncDB.SUPPORTED_SPECIES;
        for child in self.children:
            if type(child).__name__ == "Species":
                taxidList = child.supported_species
                break;
        taxid_filter = "";
        if len(taxidList) !=0:
            taxid_filter = ['$6==\"\"']+['$6==\"'+ t + '\"' for t in taxidList]
            taxid_filter = "("+"||".join(taxid_filter) + ") &&"
            
		#filter terms of desired species and type.	
        cmd = "time zcat " + SyncDB.DOWNLOAD_DIR() + "/bsid2info.gz | cut -f1,2,3,4,5,7,8 | awk 'BEGIN{FS=\"\\t\"; OFS=\"\\t\"}{ if (" + taxid_filter + "($2==\"KEGG\") && ($5==\"functional_set\" || $5==\"pathway\" || $5==\"structural_complex\")&&($3!~/^ko/)) print $0;}'| sort -t $'\\t' -k1,1 > " + SyncDB.DOWNLOAD_DIR() + "/bsid2info_human.csv";
        print cmd;
        util.unix(cmd);

        #bsid,source_id,type
        cmd = "time cut -f1,3,5 " + SyncDB.DOWNLOAD_DIR() + "/bsid2info_human.csv > " + SyncDB.DOWNLOAD_DIR() + "/bsid2info_ids.csv";
        print cmd;
        util.unix(cmd);
        
		#gene_id to term_id map.
        #bsid,gene_id
        cmd="time zcat " + SyncDB.DOWNLOAD_DIR() + "/biosystems_gene.gz | cut -f1,2 | sort -k1,1 -t $'\\t' > " + SyncDB.DOWNLOAD_DIR() + "/biosystems_gene_selected.csv";        
        print cmd;
        util.unix(cmd);
        
		#join gene_id to term_id data with term info map.
        #bsid,source_id,type,gene_id
        cmd = "time join -1 1 -2 1 -t $'\\t' " + SyncDB.DOWNLOAD_DIR() + "/bsid2info_ids.csv " + SyncDB.DOWNLOAD_DIR() + "/biosystems_gene_selected.csv | sort -k4,4 -t $'\\t' >" + SyncDB.DOWNLOAD_DIR() + "/biosystem_processed.csv";        
        print cmd;
        util.unix(cmd);
        
        if not os.path.exists(os.path.join(SyncDB.DOWNLOAD_DIR(),"gene_info.gz")):
            urllib.urlretrieve("ftp://ftp.ncbi.nih.gov/gene/DATA/gene_info.gz", os.path.join(SyncDB.DOWNLOAD_DIR(),"gene_info.gz"))


        taxid_filter = "";
        if len(taxidList) !=0:
            taxid_filter = ['$1==\"'+ t + '\"' for t in taxidList]
            taxid_filter = "if ("+"||".join(taxid_filter) + ")"
            
		#get gene_id to tax_id map.	
        #gene_id,tax_id
        cmd = "time zcat " + SyncDB.DOWNLOAD_DIR() + "/gene_info.gz | cut -f1,2 | sed 1d | awk 'BEGIN{FS=\"\\t\"; OFS=\"\\t\"}{" + taxid_filter + " print $2,$1;}' | sort -k1,1 -t $'\\t' >" + SyncDB.DOWNLOAD_DIR() + "/geneid2taxid.csv";
        print cmd;
        util.unix(cmd);
        
		# merge gene_id to term data with gene_id to tax_id
        #gene_id,bsid,source_id,type,tax_id
        cmd = "time join -1 4 -2 1 -t $'\\t' " + SyncDB.DOWNLOAD_DIR() + "/biosystem_processed.csv " +  SyncDB.DOWNLOAD_DIR() + "/geneid2taxid.csv > " + self.fn_term_gene_pair;
        print cmd;
        util.unix(cmd);

        cmd = "sed -i '1s/^/gid\\tbsid\\tterm_id\\ttype\\ttax_id\\n/' " + self.fn_term_gene_pair;
        print cmd;
        util.unix(cmd);
        
        #bsid
        cmd = "time cut -f2 " + self.fn_term_gene_pair + " | sort | uniq > " + SyncDB.DOWNLOAD_DIR() + "/bsid_list.csv";
        print cmd;
        util.unix(cmd);

        #bsid,source,source_id,name,type,tax_id,description
        cmd = "time join -1 1 -2 1 -t $'\\t' " + SyncDB.DOWNLOAD_DIR() + "/bsid2info_human.csv " +  SyncDB.DOWNLOAD_DIR() + "/bsid_list.csv > " + self.fn_terms;
        print cmd;
        util.unix(cmd);

        cmd = "sed -i '1s/^/bsid\\tsource\\tsource_id\\tname\\ttype\\ttax_id\\tdescription\\n/' " + self.fn_terms
        print cmd
        util.unix(cmd)

        #fill 0 in tax_id in self.fn_terms, so that mix species term will not be lost.
        df = pd.read_csv(self.fn_terms, sep='\t', dtype=str)
        df['tax_id'] = [x if pd.notnull(x) else '0' for x in df['tax_id']]
        df.to_csv(self.fn_terms, sep='\t', index=False)

    def check_inputs(self):
        return True;

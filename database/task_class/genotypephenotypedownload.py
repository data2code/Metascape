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
from core import *
import util
import csv
import urllib2
import sys
from gputil import GPUtils
import StringIO

from IPython.core.debugger import Tracer
#tew

class GenotypePhenotypeDownload(XmlClass):
    def __init__(self, xe=None):
        XmlClass.__init__(self,xe=xe)
        self.fn_dest =os.path.join(SyncDB.DOWNLOAD_DIR(),xe.attrib['dest']) 
        self.inputs = []
        
    def populate_more(self,root):
        self.outputs = [self.fn_dest]
               
    def get_data(self, type_name):
        if type_name=="dbGap":
            url = 'http://www.ncbi.nlm.nih.gov/projects/gapplus/get_objects.cgi?&src=0&download=1';
        elif type_name=="NHGRI_GWAS":
            url = 'http://www.ncbi.nlm.nih.gov/projects/gapplus/get_objects.cgi?&src=1&download=1';
        
        print url
        tries = 0;
        succeed = False
        while (tries < 10):
            try:
                data = util.read_csv(StringIO.StringIO(urllib2.urlopen(url).read()), sep="\t")
                succeed = True;
                break;
            except:
                print "Loading ", url, " failed. Try again...", tries
                tries += 1;

        if not succeed:
            raise ValueError('Loding url failed')
        result = []
        for index, row in data.iterrows():
            if row["GENE1_ID"] != row["GENE2_ID"]:
                result.append({'gid':row["GENE1_ID"], 'content':row[" Trait"], 'annotation_field1':row['ANALYSIS_NAME']})
                result.append({'gid':row["GENE2_ID"], 'content':row[" Trait"], 'annotation_field1':row['ANALYSIS_NAME']})
            else:
                result.append({'gid':row["GENE1_ID"], 'content':row[" Trait"], 'annotation_field1':row['ANALYSIS_NAME']})           
                    
        out=[]
        for k, g in pd.DataFrame(result).groupby(['gid']):   
            out.append({'gid':k, 'content':";".join(pd.unique(g['content'].tolist())), 'type_name': type_name, 'annotation_field1':";".join(pd.unique(g['annotation_field1'].tolist()))}) 
        
        return out;    
        
    def do_update(self):
        dbGap=self.get_data("dbGap");
        NHGRI=self.get_data("NHGRI_GWAS");
        df=pd.DataFrame(dbGap+NHGRI)
        df['tax_id']='9606'
        df.to_csv(self.fn_dest, index=False, sep=',');                        

    def check_inputs(self):
        passed = True
        urls=[]
        print "Checking urls for ncbi genotype/phenotype..."
            
        urls.append("http://www.ncbi.nlm.nih.gov/projects/gapplus/get_objects.cgi?&src=0&download=1")
        urls.append("http://www.ncbi.nlm.nih.gov/projects/gapplus/get_objects.cgi?&src=1&download=1")
        
        for url in urls:
            if not GPUtils.check_url(url):
                passed = False;                
                
        return passed        

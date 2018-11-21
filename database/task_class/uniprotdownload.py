#!/usr/bin/env python
#from .core import *
from os import sys, path

from gputil import GPUtils

sys.path.insert(0, path.join(path.dirname(path.abspath(__file__)),'../'))
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
import json
import sys
import requests
import eutils
from lxml import etree
import re
import requests

from IPython.core.debugger import Tracer
#tew

class UniProtDownload(XmlClass):
    def __init__(self, xe=None):
        XmlClass.__init__(self,xe=xe)
        self.fn_dest =os.path.join(SyncDB.DOWNLOAD_DIR(),"uniprot.csv")
        self.fn_dest_kinase =os.path.join(SyncDB.DOWNLOAD_DIR(),"uniprot_kinase.csv")
        self.fn_dest_idmapping =os.path.join(SyncDB.DOWNLOAD_DIR(),"uniprot_idmapping.csv")
        self.taxidList = SyncDB.SUPPORTED_SPECIES;
        for child in self.children:
            if type(child).__name__ == "Species":
                self.taxidList = child.supported_species
                break;

    def populate_more(self,root):
        self.outputs = [self.fn_dest,self.fn_dest_idmapping,self.fn_dest_kinase]    
    
    def create_uniprot2gid_map(self):
        df = util.read_csv(os.path.join(SyncDB.DOWNLOAD_DIR(),self.fn_dest_idmapping))
        #df = df[df["tax_id"] == 9606]
        df = df.query('tax_id in [%s]'%','.join(self.taxidList))
        self.uniprot2gid = {};
        for index, el in df.iterrows():
            self.uniprot2gid[el["source_id"]] = el["gid"]
        
    def create_idmapping(self):
        if not os.path.exists(os.path.join(SyncDB.DOWNLOAD_DIR(),"uniprot_idmapping.dat.gz")):
            urllib.urlretrieve("ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz", os.path.join(SyncDB.DOWNLOAD_DIR(),"uniprot_idmapping.dat.gz"))
            print "UniProt idmapping file downloaded...";
        else:
            print "uniprot_idmapping.dat.gz exists. start processing...";

        #keep the lines that contain tax_id and gene_id from the file.
        if not os.path.exists(os.path.join(SyncDB.DOWNLOAD_DIR(),"uniprot_idmapping_cleaned.dat")):
            filter = '\|'.join(['NCBI_TaxID\\t'+x for x in self.taxidList])
            print "extracting rows of desired species with grep command...."
            unix_cmd = "LC_ALL=C zgrep 'NCBI_TaxID\|GeneID' " + os.path.join(SyncDB.DOWNLOAD_DIR(),"uniprot_idmapping.dat.gz") \
                       + " | LC_ALL=C grep --no-group-separator -A 1 $'" + filter + "' > " \
                       + os.path.join(SyncDB.DOWNLOAD_DIR(), "uniprot_idmapping_cleaned.dat");
            print unix_cmd
            util.unix(unix_cmd)
            print "extracting rows of desired species with grep command is done."
            
        df_itr = util.read_csv(os.path.join(SyncDB.DOWNLOAD_DIR(),"uniprot_idmapping_cleaned.dat"), sep="\t", header=None, names=["uniprot_id", "type", "value"],iterator=True, chunksize=50000);
        

        data = []
        chunks = 0;
        cur_uid = None;
        cur_tax_id = None;
        cur_gene_id = None;
        
        for df in df_itr:                       
            for index, el in df.iterrows():
                if cur_uid != el["uniprot_id"]:
                    cur_uid = el["uniprot_id"]
                    #discard previous tax_id and gene_id 
                    cur_tax_id = None
                    cur_gene_id = None
            
                if el["type"] == "GeneID":
                    cur_gene_id = el["value"]
                elif el["type"] == "NCBI_TaxID":
                    cur_tax_id = el["value"]
                else:
                    print "error in idmapping_clean.dat file. expecting 'GeneID' or 'TaxID', seen ", el["type"];
                    
                if cur_gene_id is not None and cur_tax_id is not None:
                    #gene_id and tax_id for the current uid is found
                    if cur_tax_id in [int(x) for x in self.taxidList]:
                        data.append({"source_id":cur_uid, "gid":cur_gene_id, "tax_id":cur_tax_id})
                        #Tracer()()
            chunks += 1;        
            print "processed ", chunks*50000, " found ", len(data);    
            del df;
            #break;
        
        pd.DataFrame(data).to_csv(self.fn_dest_idmapping, index=False)
        print "UniProt idmapping file processed. ", len(data), " ids found";

    def process_uniprot_entry(self, entry):
        row = {"type":[], "gid_xml":[], "acc":[], "content":[], "keyword":[]};

        ORG = entry.find('{http://uniprot.org/uniprot}organism')
        if ORG is None:
            return None;

        for org in ORG.findall('{http://uniprot.org/uniprot}dbReference'):            
            if "type" in org.attrib and org.attrib["type"] == 'NCBI Taxonomy':
                if org.attrib["id"] in self.taxidList:
                    row["tax_id"] = org.attrib["id"]
                    break;
                else:
                    return None;                   
        
        for ac in entry.findall('{http://uniprot.org/uniprot}accession'):
            row ["acc"].append(ac.text);
            
        for dr in entry.findall('{http://uniprot.org/uniprot}dbReference'):
            if "type" in dr.attrib and dr.attrib["type"] == 'GeneID':
                row["gid_xml"].append(dr.attrib["id"]);

        kinase = False
        for feature in entry.findall('{http://uniprot.org/uniprot}feature'):
            if "type" in feature.attrib and "transmembrane" in feature.attrib["type"]:
                row["type"].append("uniprot_transmembrane");

            if "type" in feature.attrib and feature.attrib["type"] == "domain" and "description" in feature.attrib and "kinase" in feature.attrib["description"]:
                row["type"].append("uniprot_kinase_class");
                kinase=True
                    
        if kinase:
            for com in entry.findall('{http://uniprot.org/uniprot}comment'):
                if com.attrib["type"] == "similarity":
                    txt=com.find("{http://uniprot.org/uniprot}text").text;
                    if "protein kinase family" in txt:
                        row["content"].append(txt);
                        
        for keyword in entry.findall('{http://uniprot.org/uniprot}keyword'):
            kt=keyword.text.lower();
            if kt== "secreted":
                row["type"].append("uniprot_secreted");
            if kt== "kinase":
                row["type"].append("uniprot_kinase_class");
            if kt== "transmembrane":
                row["type"].append("uniprot_transmembrane");
            
            if "transmembrane" in kt or "kinase" in kt or "secreted" in kt:
                row["keyword"].append(keyword.text)
        
        if len(row["acc"]):
            row["gid"] = self.uniprot2gid[row["acc"][0]] if row["acc"][0] in self.uniprot2gid else None;

        return row;
        
    def process_uniprot_xml (self):
        if not os.path.exists(os.path.join(SyncDB.DOWNLOAD_DIR(),"uniprot_sprot.xml")):
            urllib.urlretrieve("ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.xml.gz", os.path.join(SyncDB.DOWNLOAD_DIR(),"uniprot_sprot.xml.gz"))
            util.unix("gunzip " + os.path.join(SyncDB.DOWNLOAD_DIR(),"uniprot_sprot.xml.gz"))
            
        #http://stackoverflow.com/questions/7018326/lxml-iterparse-in-python-cant-handle-namespaces
        context = etree.iterparse(os.path.join(SyncDB.DOWNLOAD_DIR(),"uniprot_sprot.xml"), tag='{http://uniprot.org/uniprot}entry')
        out = [];
        i = 0;
        gid_not_found=[];
        multiple_gid=[];
        multiple_type=[];
        missing_type=[];
        missing_keyword=[];
        
        for action, elem in context:
            i = i + 1;
            if i%10000 == 0: print "processed " + str(i) + " uniprot entries";
            if elem.tag=='{http://uniprot.org/uniprot}entry':
                row = self.process_uniprot_entry (elem);
                if row is None: #it's not a desired species.
                    continue
                #Tracer()()    
                
                if len(row['type']) == 0:
                    if len(row['keyword']) != 0:
                        #print "missing_type";
                        #Tracer()();
                        missing_type.append(row);
                    continue;
                    
                if row['gid'] is None:
                    gid_not_found.append(row);
                    #print "gid not found";
                    #Tracer()();                    
                    continue;
                    
                if len(list(set(row['type'])))>1:
                    #print "multiple type";
                    #Tracer()();
                    multiple_type.append(row);

                if len(row['keyword']) == 0:
                    #Tracer()();
                    missing_keyword.append(row);
                    
                for t in list(set(row['type'])):
                    content="yes";
                    if (t == "uniprot_transmembrane"):
                        content = len(row['acc']);
                    if (t == "uniprot_kinase_class"):
                        content = ",".join(row['content']);
                    out.append({"gid":row["gid"], "type_name":t, "content":content, 
                                "annotation_field1":",".join(row["acc"]), "tax_id":row["tax_id"]});
                # processing goes here

            elem.clear()
            # second, delete previous siblings (records)
            while elem.getprevious() is not None:
                del elem.getparent()[0]
            
        
        df=pd.DataFrame(out)

        df.to_csv(self.fn_dest, index=False)
        print "Multiple type:", len(multiple_type);
        print "No GID:", len(gid_not_found);
        print "Missing type:", len(missing_type);
        print "Missing keyword:", len(missing_keyword);

        print "No GID:", gid_not_found;
        return out;      
        
    def parse_uniprot_kinase(self):
        url='http://www.uniprot.org/docs/pkinfam'
        r = requests.post(url)
        if not r.ok:
            util.error_msg('Cannot fetch kinase members from UniProt: %s', r.text)
        S=r.content.split('\n')
        #s_file='pkinfam'
        #f=open(s_file)
        data=[]
        n=len(S)
        i=0
        while i<n:
            line=S[i]
            i+=1
            if re.search(r"^=+", line):
                s_grp=S[i]
                i+=1
                s=S[i]
                i+=1
                if not re.search(r"^=+", line):
                    util.error_msg("Parsing error, expecting: ====")
                s=S[i]
                i+=1
                while re.search('^\W*$', s):
                    s=S[i]
                    i+=1
                    continue
                #if re.search('\w', s):
                #    util.error_msg("Parsing error, expecting a blank line")
                while re.search('\w', s):
                    rslt=re.search('_HUMAN\s+\(<a.+>(\w+)<\/a>\s+\)', s)
                    if rslt is not None:
                        data.append({'annotation_field1': rslt.groups()[0], 'content':s_grp})
                    s=S[i]
                    i+=1
        t_kinase=pd.DataFrame(data)
        t_kinase['gid']=t_kinase.annotation_field1.apply(lambda x: self.uniprot2gid.get(x, 0))
        t_kinase['tax_id'] = '9606'
        t_kinase=t_kinase[t_kinase['gid']>0].copy()
        t_kinase.to_csv(os.path.join(SyncDB.DOWNLOAD_DIR(),self.fn_dest_kinase), index=False)
        print "%d Kinase Proteins Fetched" % len(t_kinase)
               
       
    def do_update(self):
        self.create_idmapping();
        self.create_uniprot2gid_map();
        self.parse_uniprot_kinase();
        self.process_uniprot_xml();

        
    def check_inputs (self):
        passed = True
        urls=[]
        print "Checking urls for godownload..."
            
        urls.append("http://www.uniprot.org/docs/pkinfam")        
        urls.append("ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.xml.gz")
        urls.append("ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz")        
        
        for url in urls:
            if not GPUtils.check_url(url):
                passed = False;                
                
        return passed        

if __name__=="__main__":

    u = UniProtDownload()
    u.do_update()

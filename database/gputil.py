#!/usr/bin/env python
import shlex, subprocess

import numpy as np
import pandas as pd
import math
import os
import shutil
import re
import sys
import util
import urllib
from os import sys, path
sys.path.insert(0, path.join(path.dirname(path.abspath(__file__)),'../'))    
from core import SyncDB

from IPython.core.debugger import Tracer
import urllib2

class GPUtils(object):


    @staticmethod
    def get_sym2gid_map ():
        if hasattr(GPUtils, 'sym2gid_map'):
            return GPUtils.sym2gid_map
            
        if not os.path.exists(path.join(SyncDB.DOWNLOAD_DIR(), "gene_info.gz")):
            urllib.urlretrieve("ftp://ftp.ncbi.nih.gov/gene/DATA/gene_info.gz", os.path.join(SyncDB.DOWNLOAD_DIR(), "gene_info.gz"))           
        gene_info=util.read_csv(os.path.join(SyncDB.DOWNLOAD_DIR(), 'gene_info.gz'), header=None, skiprows=1, sep="\t",
                                names=["tax_id","GeneID","Symbol","LocusTag","Synonyms","dbXrefs","chromosome","map_location","description","type_of_gene","Symbol_from_nomenclature_authority","Full_name_from_nomenclature_authority","Nomenclature_status","Other_designations","Modification_date","Feature_type"]
                                ).query('tax_id in [9606,10090,10116]')[['GeneID', 'Symbol','Synonyms']]
        
        GPUtils.sym2gid_map = {"sym2gid":{}, "gid2sym":{}}
        lookup = {};
        for i in gene_info.index:
            GPUtils.sym2gid_map["sym2gid"][gene_info.at[i, 'Symbol']] = gene_info.at[i, 'GeneID']
            GPUtils.sym2gid_map["gid2sym"][gene_info.at[i, 'GeneID']] = gene_info.at[i, 'Symbol']
            for syn in str(gene_info.at[i, 'Synonyms']).split('|'):
                syn = syn.strip()
                if (len(syn)>1):
                    GPUtils.sym2gid_map["sym2gid"][syn] = gene_info.at[i, 'GeneID']
                    GPUtils.sym2gid_map["gid2sym"][gene_info.at[i, 'GeneID']] = syn
        
        #Tracer()()
        return GPUtils.sym2gid_map

    @staticmethod            
    def get_ensembl2gid_map_ncbi ():
        if hasattr(GPUtils, 'ensembl2gid_map_ncbi'):
            return GPUtils.ensembl2gid_map_ncbi
            
        if not os.path.exists(path.join(SyncDB.DOWNLOAD_DIR(), "gene2ensembl.gz")):
            urllib.urlretrieve("ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz", os.path.join(SyncDB.DOWNLOAD_DIR(), "gene2ensembl.gz"))           
        gene2ens=util.read_csv(os.path.join(SyncDB.DOWNLOAD_DIR(), 'gene2ensembl.gz'), header=None, skiprows=1, sep="\t", names=["tax_id","GeneID","Ensembl_gene_identifier","RNA_nucleotide_accession.version","Ensembl_rna_identifier","protein_accession.version","Ensembl_protein_identifier"]).query('tax_id in [9606]')[['GeneID', 'Ensembl_gene_identifier']]        
        gene2ens = gene2ens.drop_duplicates();
        GPUtils.ensembl2gid_map_ncbi  = {};
        for i in gene2ens.index:
            GPUtils.ensembl2gid_map_ncbi[gene2ens.at[i, 'Ensembl_gene_identifier']] = GPUtils.ensembl2gid_map_ncbi.get(gene2ens.at[i, 'Ensembl_gene_identifier']) or []
            GPUtils.ensembl2gid_map_ncbi[gene2ens.at[i, 'Ensembl_gene_identifier']].append(gene2ens.at[i, 'GeneID'])
        #Tracer()()
        return GPUtils.ensembl2gid_map_ncbi

    @staticmethod
    def get_ensembl_mart_url():
        from ftplib import FTP		
        ftp = FTP('ftp.ensembl.org')
        ftp.login()
        ftp.cwd('pub/current_mysql/')
        ll=ftp.nlst()
        ensembl_mart_file = [l for l in ll if 'ensembl_mart' in l][0]
        ftp.quit()
        return 'ftp://ftp.ensembl.org/pub/current_mysql/' + ensembl_mart_file + '/hsapiens_gene_ensembl__gene__main.txt.gz'
         
    @staticmethod            
    def get_ensembl2gid_map_old():
        if hasattr(GPUtils, 'ensembl2gid_map'):
            return GPUtils.ensembl2gid_map
        hgnc_file = path.join(SyncDB.DOWNLOAD_DIR(), "hgnc_complete_set.txt");     
        ensembl_file = path.join(SyncDB.DOWNLOAD_DIR(), "hsapiens_gene_ensembl__gene__main.txt.gz"); 
        
        if not os.path.exists(hgnc_file):
            urllib.urlretrieve("ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt", hgnc_file)
        if not os.path.exists(ensembl_file):
            urllib.urlretrieve(GPUtils.get_ensembl_mart_url(), ensembl_file)

        hgnc_data=util.read_csv(hgnc_file, sep="\t")[['hgnc_id','symbol', 'ensembl_gene_id','entrez_id']]
        #Tracer()();
        
        ensembl_data=util.read_csv(ensembl_file, header=None, skiprows=1, sep="\t")[[4,5,6]]
        hgnc_lookup = {};
        for i in hgnc_data.index:
            hgnc_lookup[hgnc_data.at[i, 'hgnc_id']] = hgnc_data.at[i, 'entrez_id']

        out=[];
        lookup = {};
        for i in ensembl_data.index:        
            hgnc_id = ensembl_data.at[i, 4];
            l = hgnc_id.find("HGNC:")
            if  l < 0:
                #print "HGNC: was not found in ", hgnc_id;
                continue
            hgnc_id = hgnc_id[l:-1]    
            
            if hgnc_id in hgnc_lookup:
                lookup[ensembl_data.at[i, 6]] = hgnc_lookup[hgnc_id]
                out.append({'ensembl_gene_id':ensembl_data.at[i, 6], 'gene_id':hgnc_lookup[hgnc_id]})
        
        pd.DataFrame(out).sort(['ensembl_gene_id']).to_csv('ebi_ensembl_map.csv', index=False);
        #Tracer()()
        GPUtils.ensembl2gid_map = lookup
        return lookup

    @staticmethod        
    def get_ensembl2gid_map():
        if hasattr(GPUtils, 'ensembl2gid_map'):
            return GPUtils.ensembl2gid_map                   
        GPUtils.ensembl2gid_map = GPUtils.get_ensembl2gid_map_ensembl()
        
        #ebi_map = GPUtils.get_ensembl2gid_map_ebi()        
        #Tracer()()
        #for g in ebi_map:
        #    if g not in GPUtils.ensembl2gid_map:
        #        GPUtils.ensembl2gid_map[g] = ebi_map[g];
        
        ncbi_map = GPUtils.get_ensembl2gid_map_ncbi()        

        for g in ncbi_map:
            if g not in GPUtils.ensembl2gid_map:
                GPUtils.ensembl2gid_map[g] = np.unique(ncbi_map[g]);
            else:
                GPUtils.ensembl2gid_map[g] = np.unique(ncbi_map[g] + GPUtils.ensembl2gid_map[g]);
        
        return GPUtils.ensembl2gid_map
        
    @staticmethod            
    def get_ensembl2gid_map_ebi():
        if hasattr(GPUtils, 'ensembl2gid_map_epi'):
            return GPUtils.ensembl2gid_map_ebi

        hgnc_file = path.join(SyncDB.DOWNLOAD_DIR(), "hgnc_complete_set.txt");     
        
        if not os.path.exists(hgnc_file):
            urllib.urlretrieve("ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt", hgnc_file)

        hgnc_data=util.read_csv(hgnc_file, sep="\t")[['hgnc_id','symbol', 'ensembl_gene_id','entrez_id']]
        #Tracer()();
        
        GPUtils.ensembl2gid_map_ebi = {};
        for i in hgnc_data.index:
            GPUtils.ensembl2gid_map_ebi[hgnc_data.at[i, 'ensembl_gene_id']] = GPUtils.ensembl2gid_map_ebi.get(hgnc_data.at[i, 'ensembl_gene_id']) or [];
            GPUtils.ensembl2gid_map_ebi[hgnc_data.at[i, 'ensembl_gene_id']].append(hgnc_data.at[i, 'entrez_id']);
        #Tracer()()
        
        return GPUtils.ensembl2gid_map_ebi;
        
    @staticmethod            
    def get_ensembl2gid_map_ensembl ():
        if hasattr(GPUtils, 'ensembl2gid_map_ensembl'):
            return GPUtils.ensembl2gid_map_ensembl
        
        ensembl_file = path.join(SyncDB.DOWNLOAD_DIR(), "ensembl_genes_info.csv"); 
        
        if not os.path.exists(ensembl_file):
            cmd='wget -O ' + ensembl_file + ' - \'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0"\
            encoding="UTF-8"?><!DOCTYPE Query><Query virtualSchemaName = "default"\
            formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "hsapiens_gene_ensembl" interface = "default"\
            ><Attribute name = "ensembl_gene_id"/><Attribute name = "entrezgene"/></Dataset></Query>\'';
            
            util.unix(cmd);

        ensembl_data=util.read_csv(ensembl_file, sep="\t", header=None, names=['ensembl_gene_id','gene_id']);
        ensembl_data=ensembl_data[ensembl_data['gene_id'].notnull()];
        ensembl_data[['gene_id']]=ensembl_data[['gene_id']].astype(int)
        GPUtils.ensembl2gid_map_ensembl = {};
        for i in ensembl_data.index:        
            GPUtils.ensembl2gid_map_ensembl[ensembl_data.at[i, 'ensembl_gene_id']] = GPUtils.ensembl2gid_map_ensembl.get(ensembl_data.at[i, 'ensembl_gene_id']) or [];
            GPUtils.ensembl2gid_map_ensembl[ensembl_data.at[i, 'ensembl_gene_id']].append(ensembl_data.at[i, 'gene_id']);
        #Tracer()()
        return GPUtils.ensembl2gid_map_ensembl
        
    @staticmethod 
    def check_inputs():
        passed = True
        urls=[]
        print "Checking urls in gputil..."
            
        urls.append("ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt")        
        urls.append(GPUtils.get_ensembl_mart_url())
        urls.append("ftp://ftp.ncbi.nih.gov/gene/DATA/gene_info.gz")
        for url in urls:
            if not GPUtils.check_url(url):
                passed = False;                                
        return passed

    @staticmethod
    def check_url3(url):
        """ Check a URL """
        import urllib2
        try:
            urllib.urlopen(url)
        except:
            #print(e.code)
            print("error in fetching url:", url)
            return False
        return True

    @staticmethod
    def check_url2(url):
        """ Check a URL """
        from IPython.core.debugger import Tracer

        import httplib2
        try:
            h = httplib2.Http()
            resp = h.request(urllib.quote(url, safe="%/:=&?~#+!$,;'@()*[]"), 'HEAD')
        except:
            return False

        if int(resp[0]['status']) < 400:
            return True
        return False


    @staticmethod
    def check_url1(url):
        """ Check a URL """
        from IPython.core.debugger import Tracer
        try:
            resp = urllib2.urlopen(HeadRequest(urllib.quote(url, safe="%/:=&?~#+!$,;'@()*[]")))
            return True
        except:
            return False

        return False

    @staticmethod
    def check_url(url):
        """ Check a URL """
        from IPython.core.debugger import Tracer
        tries = 1;
        print "checking url:", url
        while tries <= 10:
            if GPUtils.check_url1(url):
                return True
            elif GPUtils.check_url2(url):
                return True
            else:
                print "***error in fetching url:", url, " Try again...", tries
                tries += 1

        return False

class HeadRequest(urllib2.Request):
    def get_method(self):
        return "HEAD"

        
if __name__=='__main__':
    a = GPUtils.get_sym2gid_map()
    b = a["sym2gid"]
    if False:
        sys.path.insert(0, path.join(path.dirname(path.abspath(__file__)),'../'))    
        map=GPUtils.get_ensembl2gid_map();
        Tracer()()
        exit()
        ebi=GPUtils.get_ensembl2gid_map_ebi();
        ens=GPUtils.get_ensembl2gid_map_ensembl();
        #exit();
        missing_ens = [];
        diff_ens = [];
        missing_ebi=[];
        
        for g in ebi:
            if g not in ens:
                missing_ens.append({'ens_id':g, 'gene_id':ebi[g]});
            elif ens[g] != ebi[g]:
                diff_ens.append({'ens_id':g, 'ebi_gene_id':ebi[g], 'ens_gid':ens[g]});

        for g in ens:
            if g not in ebi:
                missing_ebi.append({'ens_id':g, 'gene_id':ens[g]});
                
        Tracer()();
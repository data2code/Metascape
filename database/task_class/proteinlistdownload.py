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
import json
import sys
import requests
import eutils
import Entity.ClinVar as cv
from lxml import etree

from IPython.core.debugger import Tracer
#tew
from gputil import GPUtils


class ProteinListDownload(XmlClass):
    def __init__(self, xe=None):
        XmlClass.__init__(self,xe=xe)
        self.fn_dest_targets =os.path.join(SyncDB.DOWNLOAD_DIR(),"target_drugbank.csv")   
        self.fn_dest_tissue_specific =os.path.join(SyncDB.DOWNLOAD_DIR(),"tissue_specificity_tiger.csv")   
        self.fn_dest_clinvar =os.path.join(SyncDB.DOWNLOAD_DIR(),"pathogenic_lof_clinvar.csv")   
        self.inputs=['ds:ncbi']
        
    def populate_more(self,root):
        self.outputs = [self.fn_dest_targets, self.fn_dest_tissue_specific, self.fn_dest_clinvar]

    def tissue_specific(self):
        #http://bioinfo.wilmer.jhu.edu/tiger/download/ref2tissue-Table.txt
        url='http://bioinfo.wilmer.jhu.edu/tiger/download/ref2tissue-Table.txt'
        # RefSeq Tissue(s)
        # http://www.biomedcentral.com/1471-2105/9/271
        r = requests.post(url)
        if not r.ok:
            util.error_msg('Cannot fetch tissue specific data from JHU: %s', r.text)
        S=r.content.split('\n')
        for i,s in enumerate(S):
            S[i]=s.replace('\t', ',', 1).replace('\t', ' ')
        import cStringIO
        return pd.read_csv(cStringIO.StringIO("\n".join(S)))
        
    def get_tissue_specific(self):
        if not path.isfile(os.path.join(SyncDB.DOWNLOAD_DIR(), "gene2refseq.gz")):
            urllib.urlretrieve("ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2refseq.gz", os.path.join(SyncDB.DOWNLOAD_DIR(), "gene2refseq.gz"))

        ref2gene = util.read_csv(os.path.join(SyncDB.DOWNLOAD_DIR(), "gene2refseq.gz"), skiprows=1, header=None, sep='\t', names=["tax_id","GeneID","status","RNA_nucleotide_accession.version","RNA_nucleotide_gi","protein_accession.version","protein_gi","genomic_nucleotide_accession.version","genomic_nucleotide_gi","start_position_on_the_genomic_accession","end_position_on_the_genomic_accession","orientation","assembly","mature_peptide_accession.version","mature_peptide_gi","Symbol"]).query('tax_id in [9606]');
        
        #Tracer()()       
        self.ref2gene_map = {}
        for i in ref2gene.index:                
            if ref2gene.at[i, 'RNA_nucleotide_accession.version'] != '-':
                self.ref2gene_map[ref2gene.at[i, 'RNA_nucleotide_accession.version'].split('.')[0]] = ref2gene.at[i, 'GeneID']

            if ref2gene.at[i, 'protein_accession.version'] != '-':
                self.ref2gene_map[ref2gene.at[i, 'protein_accession.version'].split('.')[0]] = ref2gene.at[i, 'GeneID']

            if ref2gene.at[i, 'genomic_nucleotide_accession.version'] != '-':
                self.ref2gene_map[ref2gene.at[i, 'genomic_nucleotide_accession.version'].split('.')[0]] = ref2gene.at[i, 'GeneID']
                

        t_tissue=self.tissue_specific()
        t_tissue.rename2({'Tissue(s)':'Tissue'})
        t_tissue['gene_id'] = t_tissue.RefSeq.apply(lambda x: self.ref2gene_map.get(x,0))
        t_tissue = t_tissue.query('gene_id > 0');
        data=[]
        for k,t_v in t_tissue.groupby('gene_id'):
            if k==0: continue
            S=[x for x in t_v['Tissue'] if not pd.isnull(x)]
            s=" ".join(S)
            S=util.unique(s.split(" "))
            data.append({'gene_id':k, 'Tissues':";".join(S)})
        t_tissue=pd.DataFrame(data)
        t_tissue['tax_id'] = '9606'
        #Tracer()()        
        t_tissue.to_csv(self.fn_dest_tissue_specific, index=False)
        print "%d Tissue-specific Genes Fetched" % len(t_tissue)
        

    def create_uniprot2gid_map(self):
        df = util.read_csv(os.path.join(SyncDB.DOWNLOAD_DIR(),"uniprot_idmapping.csv"))
        df = df[df["tax_id"] == 9606]
        self.uniprot2gid = {};
        for index, el in df.iterrows():
            self.uniprot2gid[el["source_id"]] = el["gid"]

    def process_drugbank_entry(self, entry):
        row = {"name":"", "targets":[]};

        name = entry.find('{http://www.drugbank.ca}name')
        if name is None:
            return None;
        
        row["name"] = name.text;
        
        targets = entry.find('{http://www.drugbank.ca}targets')
        if targets is None:
            return None;

        for target in targets.findall('{http://www.drugbank.ca}target'):   
            organism = target.find('{http://www.drugbank.ca}organism')
            if organism is None or organism.text != "Human":
                continue;
                
            polyp = target.find('{http://www.drugbank.ca}polypeptide')
            if polyp is None:
                continue

            gid = None
            acc_id = None
            general_function = polyp.find('{http://www.drugbank.ca}general-function')
            general_function = general_function.text if general_function is not None else ""
            
            extids = polyp.find('{http://www.drugbank.ca}external-identifiers')
            if extids is None:
                continue;
            for extid in extids.findall('{http://www.drugbank.ca}external-identifier'):
                resource = extid.find('{http://www.drugbank.ca}resource')
                id = extid.find('{http://www.drugbank.ca}identifier')                
                if resource is not None and resource.text == "UniProtKB" and id is not None:
                    gid = self.uniprot2gid.get(id.text, None)
                    acc_id = id.text
                    if gid is not None:
                        break;
            
            if gid is not None:
                row["targets"].append({"gid":gid, "function":general_function, "uniprot_id":acc_id})      
        return row;

    def get_drug_target(self):
    
        self.create_uniprot2gid_map();
        if not os.path.exists(os.path.join(SyncDB.DOWNLOAD_DIR(),"drugbank_all_full_database.xml.zip")):
            err_msg = '!!!!!!!!!!!Cannot find drugbank_all_full_database.xml.zip, we mannually download this file'
            print err_msg
            SyncDB.ERROR_LOG.write(err_msg)

        #http://stackoverflow.com/questions/7018326/lxml-iterparse-in-python-cant-handle-namespaces
        from zipfile import ZipFile
        zz = ZipFile(open(os.path.join(SyncDB.DOWNLOAD_DIR(),"drugbank_all_full_database.xml.zip")))
        context = etree.iterparse(zz.open(zz.namelist()[0]), tag='{http://www.drugbank.ca}drug')
        
        data = [];
        i = 0;
       
        for action, elem in context:
            i = i + 1;
            if i%10000 == 0: print "processed " + str(i) + " uniprot entries";
            if elem.tag=='{http://www.drugbank.ca}drug' and elem.getparent().tag == '{http://www.drugbank.ca}drugbank':
                row = self.process_drugbank_entry(elem);
                if row is not None:
                    for target in row["targets"]:
                        if target["function"] is not None:
                            data.append({"gid":target["gid"], "content": row["name"] + " (" + target["function"] + ")", "uniprot_id": target["uniprot_id"]});
                        else:
                            data.append({"gid":target["gid"], "content": row["name"], "uniprot_id": target["uniprot_id"]});

        df = pd.DataFrame(data);
        out = [];
        for gid,gr in df.groupby('gid'):
            out.append({"gid":gid, "content":";".join(gr["content"]), "annotation_field1":",".join(list(set(gr["uniprot_id"])))})

        df = pd.DataFrame(out);
        df['tax_id'] = '9606'
        df.to_csv(self.fn_dest_targets, index=False, encoding='utf-8')
        
        print "%d Drug-Targets entries processed" % len(out)
        
    def get_clinvar(self):
        eu=eutils.EUtils()

        # Filters: (Pathogenic or likely pathogenic) and (frameshif, missense or nonsense)
        id,args=eu.esearch({'db':'clinvar', 'term':'(clinsig+pathogenic[Filter] OR clinsig+likely+path[Filter]) AND (mol+cons+frameshift[Filter] OR mol+cons+missense[Filter] OR mol+cons+nonsense[Filter])'})
        # This gives ~28000 records (variants)

        # id,args=eu.esearch({'db':'clinvar', 'term':'(clinsig+pathogenic[Filter] OR clinsig+likely+path[Filter]) AND (mol+cons+frameshift[Filter] OR mol+cons+missense[Filter] OR mol+cons+nonsense[Filter]) AND (var+deletion[Filter] OR var+indel[Filter]'})
        # The above query gave ~5000 records, but it missed missense mutations because of deletion and indel filters

        args['db']='clinvar'
        out=eu.esummary(args,count = args['count'])
        x = cv.SummaryList(out)
        records = x.to_list()

        # Combine results into DataFrame

        rows = []
        for entry in records:
          data = [entry['variant_id']]
          data.extend(entry['variant'])
          for gene in entry['genes']:
            record = data[:]
            record.extend(gene) 
            for key in entry['trait'].keys():
                row = record[:] 
                row.extend([key,"; ".join(entry['trait'][key])])  
                rows.append(row)        

        df = pd.DataFrame(rows, columns = ['id','VariantType','VarinatName','Pathogenic','Symbol','GeneID','Trait','Source'])      
        df.Trait = [trait.encode('ascii','ignore') for trait in df.Trait.values]
        data=[]
        for k,t_v in df.groupby('GeneID'):
            S1=[x for x in util.unique(t_v['VariantType']) if not pd.isnull(x)]
            S2=[x[0] + "(" + x[1] + ")" for x in zip(t_v['Source'], t_v['Trait']) if not pd.isnull(x[0]) and not pd.isnull(x[1]) ]
            #print S1[:], S2[:], S3[:]
            data.append({'gene_id':k, 'variant_type':";".join(S1), 'source_trait':";".join(S2)})
            
        df=pd.DataFrame(data).query('source_trait != ""')
        df['tax_id'] = '9606'
        df.to_csv(self.fn_dest_clinvar, index=False)
        print "%d clinvar entries processed."%len(df);
        
    def do_update(self):
        #return;
        self.get_drug_target()
        self.get_tissue_specific()
        self.get_clinvar()
        
    def check_inputs (self):
        passed = True
        urls=[]
        print "Checking urls for protein list download..."
            
        urls.append("ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2refseq.gz")        
        urls.append("http://bioinfo.wilmer.jhu.edu/tiger/download/ref2tissue-Table.txt")

        for url in urls:
            if not GPUtils.check_url(url):
                passed = False;                
                
        return passed        
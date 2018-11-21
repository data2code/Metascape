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

from IPython.core.debugger import Tracer
#tew

class ProteinListDownload(XmlClass):
    def __init__(self, xe=None):
        XmlClass.__init__(self,xe=xe)
        self.fn_dest_secreted =os.path.join(SyncDB.DOWNLOAD_DIR(),"uniprot_secreted.csv")
        self.fn_dest_transmembrane =os.path.join(SyncDB.DOWNLOAD_DIR(),"uniprot_transmembrane.csv")
        self.fn_dest_targets =os.path.join(SyncDB.DOWNLOAD_DIR(),"target_drugbank.csv")   
        self.fn_dest_tissue_specific =os.path.join(SyncDB.DOWNLOAD_DIR(),"tissue_specificity_tiger.csv")   
        self.fn_dest_clinvar =os.path.join(SyncDB.DOWNLOAD_DIR(),"pathogenic_lof_clinvar.csv")   
        self.inputs=['ds:hithub']
           
    def populate_more(self,root):
        self.outputs = [self.fn_dest_secreted, self.fn_dest_transmembrane, self.fn_dest_targets, self.fn_dest_tissue_specific, self.fn_dest_clinvar]

    def prepare (self):
        self.dbcon=db.get_con('HITHUB')
        s_uniprot2gene="SELECT u.uniprot_accn AS accn,e.gene_id FROM ensembl.uniprot u JOIN ensembl.entrez e USING(ensembl_gene_id) JOIN ensembl.gene_info g USING(ensembl_gene_id) WHERE g.tax_id=9606"
        t_map=self.fetch(s_uniprot2gene)
        self.c_uniprot2gene={ t_map.at[i, 'accn']:t_map.at[i, 'gene_id'] for i in t_map.index }
        print "%d UniProt-Entrez Gene Mapping Fetched" % len(self.c_uniprot2gene)

        
    def fetch(self, s_sql):
        t=db.from_sql(self.dbcon, s_sql)
        if 'gene_id' not in t.header():
            t['gene_id']=t.accn.apply(lambda x: self.c_uniprot2gene.get(x, 0))
            t=t[t['gene_id']>0].copy()
        return t

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

        ref2gene = util.read_csv(os.path.join(SyncDB.DOWNLOAD_DIR(), "gene2refseq.gz"), skiprows=1, header=None, sep='\t', names=["tax_id","GeneID","status","RNA_nucleotide_accession.version","RNA_nucleotide_gi","protein_accession.version","protein_gi","genomic_nucleotide_accession.version","genomic_nucleotide_gi","start_position_on_the_genomic_accession","end_position_on_the_genomic_accession","orientation","assembly","mature_peptide_accession.version","mature_peptide_gi","Symbol"]).query('tax_id in [9606,10090,10116]');
        
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
        
        #Tracer()()        
        t_tissue.to_csv(self.fn_dest_tissue_specific, index=False)
        print "%d Tissue-specific Genes Fetched" % len(t_tissue)
        
    def get_secreted(self):
        s_secreted="SELECT DISTINCT accn,'secreted' As type FROM uniprot.uniprotkb kb JOIN uniprot.entry_keyword  USING(accn) WHERE kb.tax_id = 9606 AND origin='uniprot/swiss_prot' AND keyword='Secreted'"
        # AND amino_acid_count<=300
        t_secreted=self.fetch(s_secreted)
        t_secreted['type'] = 'uniprot_secreted'
        t_secreted.rename2({'type':'term_id', 'gene_id':'gid'})
        t_secreted.to_csv(self.fn_dest_secreted, index=False)
        print "%d Secreted Proteins Fetched" % len(t_secreted)
        #c_se={ t_secreted.at[i, 'accn']:True for i in t_secreted.index }
            
    def get_transmembrane(self):
        s_transmem="SELECT accn, COUNT(1) AS tmcount FROM uniprot.feature f JOIN uniprot.uniprotkb USING(accn) WHERE origin = 'uniprot/swiss_prot' AND f.name = 'transmem' AND tax_id = 9606 GROUP BY accn"
        t_transmem=self.fetch(s_transmem)
        t_transmem['term_id'] = 'uniprot_transmembrane';
        t_transmem.rename2({'gene_id':'gid'})
        t_transmem.to_csv(self.fn_dest_transmembrane, index=False)
        print "%d Transmembrane Proteins Fetched" % len(t_transmem)

    def get_target(self):
        s_target="SELECT d.name As drug,px.identifier As accn,p.name,p.general_function from drugbank.drug_target t JOIN drugbank.drugs d USING(drugbank_id) JOIN drugbank.partner_external_identifiers px USING(partner_id) JOIN drugbank.partners p USING(partner_id) where p.species_uniprot_taxon_id=9606 AND px.resource='UniProtKB'"
        t_target=self.fetch(s_target)
        data=[]
        for k,t_v in t_target.groupby('gene_id'):
            #Tracer()()
            S1=[x[0] + "(" + x[1] + ")" for x in zip(t_v['drug'], t_v['general_function']) if not pd.isnull(x[0]) and not pd.isnull(x[1]) ]
            data.append({'gene_id':k, 'known_drug':";".join(S1)})            
        t_target=pd.DataFrame(data).query('known_drug != ""')        
        t_target.to_csv(self.fn_dest_targets, index=False)
        print "%d Drug-Targets Pairs Fetched" % len(t_target)
        
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
        df.to_csv(self.fn_dest_clinvar, index=False)
        
    def do_update(self):
        #return;
        self.prepare()
        self.get_secreted()
        self.get_transmembrane()
        self.get_target()
        self.get_tissue_specific()
        self.get_clinvar()
        
    def check_inputs (self):
        dbcon=db.get_con('HITHUB')
        print "Check database connection for HITHUB..."
        try:            
            db.from_sql(dbcon, "SELECT * FROM ensembl.uniprot limit 10")            
        except:
            print "Error in connecting to the HITHUB database"
            dbcon.close()
            return False            
        dbcon.close()            
        return True

#!/usr/bin/env python
#from .core import *
from os import sys, path

from task_class.ensembldownload import EnsemblDownload

sys.path.insert(0, path.join(path.dirname(path.abspath(__file__)),'../mylib'))
sys.path.insert(0, path.join(path.dirname(path.abspath(__file__)),'../'))
import numpy as np
import pandas as pd
import shutil
import urllib
import urlparse
from os.path import splitext, basename
import os
from core import *
import util
import urllib2
import json
import sys
from gputil import GPUtils

from IPython.core.debugger import Tracer

class IsoformsIdMap(XmlClass):
    def __init__(self, xe=None):
        #return #DEBUG
        XmlClass.__init__(self,xe=xe)
        self.fn_dest =os.path.join(SyncDB.DOWNLOAD_DIR(),"isoform_idmap.csv")
        self.taxidList = SyncDB.SUPPORTED_SPECIES;
        for child in self.children:
            if type(child).__name__ == "Species":
                self.taxidList = child.supported_species
                break;

        util.unix('mkdir -p ' + SyncDB.DOWNLOAD_DIR() + '/ensembl_files');
        
    def populate_more(self,root):
        self.outputs = [self.fn_dest]    
           
    def uniprot_idmapping_data(self):
        if not os.path.exists(os.path.join(SyncDB.DOWNLOAD_DIR(),"idmapping_selected.tab.gz")):
            urllib.urlretrieve("ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz", os.path.join(SyncDB.DOWNLOAD_DIR(),"idmapping_selected.tab.gz"))
            print "idmapping_selected.tab.gz file downloaded...";
        else:
            print "idmapping_selected.tab.gz exists. start processing...";
           
        #keep the lines that contain tax_id and gene_id from the file.
        if not os.path.exists(os.path.join(SyncDB.DOWNLOAD_DIR(),"uniprot_isoform_idmapping.dat")):
            print "extracting rows of desired species with awk command...."
            unix_cmd = "zcat " + os.path.join(SyncDB.DOWNLOAD_DIR(),"idmapping_selected.tab.gz") + " | awk -F'\t' '{if ($13 == \"9606\" || $13==\"10090\" || $13==\"10116\") print $1,$3,$4,$13,$19,$20,$21}' > " + os.path.join(SyncDB.DOWNLOAD_DIR(),"uniprot_isoform_idmapping.dat");
            print unix_cmd
            util.unix(unix_cmd)
            print "extracting rows of desired species with grep command is done."
            
        df = util.read_csv(os.path.join(SyncDB.DOWNLOAD_DIR(),"uniprot_isoform_idmapping.dat"), sep="\t", header=None, names=["Uniprot_ID","GeneID","Refseq_Prot","Tax_ID","Ensembl_Gene","Ensembl_Trans","Ensembl_Prot"]);
        df["Source"] = "Uniprot";
        df=df.iloc[np.where(df["GeneID"].notnull())[0]]
        print "Uniprot ID mapping data processed..."

        return df;
        
    def ncbi_idmapping_data(self):
        if not os.path.exists(path.join(SyncDB.DOWNLOAD_DIR(), "gene2refseq.gz")):
            urllib.urlretrieve("ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2refseq.gz", os.path.join(SyncDB.DOWNLOAD_DIR(), "gene2refseq.gz"))

        if not os.path.exists(path.join(SyncDB.DOWNLOAD_DIR(), "gene2ensembl.gz")):
            urllib.urlretrieve("ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2ensembl.gz", os.path.join(SyncDB.DOWNLOAD_DIR(), "gene2ensembl.gz"))

        if not os.path.exists(path.join(SyncDB.DOWNLOAD_DIR(), "gene_refseq_uniprotkb_collab.gz")):
            urllib.urlretrieve("ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_refseq_uniprotkb_collab.gz", os.path.join(SyncDB.DOWNLOAD_DIR(), "gene_refseq_uniprotkb_collab.gz"))
        
        gene2ref=util.read_csv(SyncDB.DOWNLOAD_DIR()+'/gene2refseq.gz', sep='\t', header=False, skiprows=1, names=['Tax_ID', 'GeneID', 'status','Refseq_RNA','RNA_nucleotide_gi','Refseq_Prot','protein_gi','genomic_nucleotide_accession.version','genomic_nucleotide_gi','start_position_on_the_genomic_accession','end_position_on_the_genomic_accession','orientation','assembly','mature_peptide_accession.version','mature_peptide_gi','Symbol'])[['Tax_ID', 'GeneID', 'Refseq_RNA','Refseq_Prot','Symbol', 'status']].query('status != "SUPPRESSED" and Refseq_RNA != "-" and Tax_ID == [9606,10090,10116]').drop_duplicates(); 

        gene2ref=gene2ref.iloc[np.where(gene2ref["GeneID"].notnull())[0]];
        gene2ref=gene2ref.iloc[np.where(gene2ref["Refseq_RNA"].notnull())[0]];
        gene2ref = gene2ref[gene2ref['Refseq_RNA'].str.contains("X")==False] #Remove computed Refseq
        ref2uniprot=util.read_csv(SyncDB.DOWNLOAD_DIR()+'/gene_refseq_uniprotkb_collab.gz', sep='\t', header=False, skiprows=1, names=['Refseq_Prot','Uniprot_ID']); 
        ref2uniprot=ref2uniprot.iloc[np.where(ref2uniprot["Refseq_Prot"].notnull())[0]].drop_duplicates(subset='Refseq_Prot', take_last=True).fillna('-')
        
        gene2ens=util.read_csv(SyncDB.DOWNLOAD_DIR()+'/gene2ensembl.gz', sep='\t', header=False, skiprows=1, names=['Tax_ID','GeneID','Ensembl_Gene','Refseq_RNA','Ensembl_Trans','Refseq_Prot','Ensembl_Prot']).query('Tax_ID == [9606,10090,10116]')[['Ensembl_Gene','Refseq_RNA','Ensembl_Trans','Ensembl_Prot']];
        gene2ens=gene2ens.iloc[np.where(gene2ens["Refseq_RNA"].notnull())[0]].drop_duplicates(subset='Refseq_RNA', take_last=True).fillna('-')
        
        
        f=lambda x: x.split('.')[0];
        gene2ref['Refseq_Prot']=gene2ref['Refseq_Prot'].map(f);
        gene2ref['Refseq_RNA']=gene2ref['Refseq_RNA'].map(f);
        gene2ens['Refseq_RNA']=gene2ens['Refseq_RNA'].map(f);        
        df = pd.merge(gene2ref, ref2uniprot, on='Refseq_Prot', how='left');
        df = pd.merge(df, gene2ens, on='Refseq_RNA', how='left');                     
        df["Source"] = "NCBI";
        df=df.fillna('-')
        print "NCBI ID mapping data processed..."
        return df;
        
    def ucsc_idmapping_data_by_taxid(self, tax_id):
        if tax_id == 9606:
            db_name="hg38";
            fname="human_ucsc_isoform_map.csv";
        elif tax_id == 10090:
            db_name="mm10";
            fname="mouse_ucsc_isoform_map.csv";
        else:
            print "species not supported";
            return None;
    
        if not os.path.exists(path.join(SyncDB.DOWNLOAD_DIR(), fname)):
            cmd = 'mysql --host=genome-mysql.cse.ucsc.edu --user=genomep --password=password --database=' + db_name + ' -e "select a.*, refLink.protAcc as Refseq_Prot, knownToEnsembl.value as Ensembl_Trans from (SELECT kgXref.refseq as Refseq_RNA, knownToLocusLink.value as GeneID, knownGene.name as UCSC_ID, knownGene.proteinID as Uniprot_ID from kgXref,knownGene,knownToLocusLink WHERE knownToLocusLink.name = knownGene.name and kgXref.kgID=knownGene.name) a left join refLink on refLink.mrnaAcc = a.Refseq_RNA left join knownToEnsembl on a.UCSC_ID=knownToEnsembl.name" > ' + path.join(SyncDB.DOWNLOAD_DIR(), fname);
            util.unix(cmd);
            print "UCSC Query is done..."

        df=util.read_csv(path.join(SyncDB.DOWNLOAD_DIR(), fname), sep='\t');
        df["Tax_ID"]=tax_id;
        df=df.iloc[np.where(df["GeneID"].notnull())[0]].drop_duplicates();
        print "UCSC ID mapping data processed..."
        return df;
        
    def ucsc_idmapping_data(self):
        ense_trans2gid_map = self.ensembl_trans2gene_map();
        df= pd.concat([self.ucsc_idmapping_data_by_taxid(9606),self.ucsc_idmapping_data_by_taxid(10090)]);
        df.fillna('-', inplace=True);
        f=lambda x: x.split('.')[0];
        df["Ensembl_Trans"] = df["Ensembl_Trans"].map(f);
        df["Ensembl_Gene"]=df["Ensembl_Trans"].map(lambda x: ense_trans2gid_map[x]["Ensembl_Gene"] if x in ense_trans2gid_map is not None else '-');
        df["Ensembl_Prot"]=df["Ensembl_Trans"].map(lambda x: ense_trans2gid_map[x]["Ensembl_Prot"] if x in ense_trans2gid_map is not None else '-');
        df["Ensembl_Trans"]=df["Ensembl_Trans"].map(lambda x: x if "ENS" in x else '-');
        
        return df;
        
    def create_ucsc_idmapping(self):
        if not os.path.exists(path.join(SyncDB.DOWNLOAD_DIR(), "ucsc_refseq_to_id_map.csv")):
            cmd = 'mysql --host=genome-mysql.cse.ucsc.edu --user=genomep --password=password --database=hg38 -e "SELECT refLink.mrnaAcc as Refseq_RNA, kgXref.kgID as UCSC_ID from refLink, kgXref WHERE kgXref.refseq=refLink.mrnaAcc" > ' + path.join(SyncDB.DOWNLOAD_DIR(), "ucsc_refseq_to_id_map.csv");
            util.unix(cmd);
        
        if not os.path.exists(path.join(SyncDB.DOWNLOAD_DIR(), "ucsc_id_to_ensemblid_map.csv")):
            cmd = 'mysql --host=genome-mysql.cse.ucsc.edu --user=genomep --password=password --database=hg38 -e "SELECT knownToEnsembl.name as UCSC_ID, knownToEnsembl.value as  Ensembl_Trans from knownToEnsembl" > ' + path.join(SyncDB.DOWNLOAD_DIR(), "ucsc_id_to_ensemblid_map.csv");
            util.unix(cmd);
                   
        
        df=util.read_csv(SyncDB.DOWNLOAD_DIR()+'/ucsc_id_to_ensemblid_map.csv', sep='\t');
        global ensemblTrans2ucscID;
        ensemblTrans2ucscID = {}
        for index, el in df.iterrows():
            ensemblTrans2ucscID[el["Ensembl_Trans"]] = el["UCSC_ID"]

        df=util.read_csv(SyncDB.DOWNLOAD_DIR()+'/ucsc_refseq_to_id_map.csv', sep='\t');
        global refseq2ucscID;
        refseq2ucscID = {}
        for index, el in df.iterrows():
            refseq2ucscID[el["Refseq_RNA"]] = el["UCSC_ID"]

    def ensembl_trans2gene_map_by_taxid(self, tax_id):     
        if tax_id == 9606:
            db_name="hsapiens_gene_ensembl";
            fname="ensembl_genes_human_trans2gid_map.csv";
        elif tax_id == 10090:
            db_name="mmusculus_gene_ensembl";
            fname="ensembl_genes_mouse_trans2gid_map.csv";
        elif tax_id == 10116:
            db_name="rnorvegicus_gene_ensembl";
            fname="ensembl_genes_rat_trans2gid_map.csv";

        util.unix('mkdir -p ' + SyncDB.DOWNLOAD_DIR() + '/ensembl_files');       
        ensembl_file = SyncDB.DOWNLOAD_DIR() + "/ensembl_files/%s"%fname; 
        
        if not os.path.exists(ensembl_file ):
            cmd = 'wget -O ' + ensembl_file + ' \'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0"\
            encoding="UTF-8"?><!DOCTYPE Query><Query virtualSchemaName = "default"\
            formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "' + db_name + '" interface = "default">\
            <Attribute name = "ensembl_gene_id"/><Attribute name = "ensembl_transcript_id"/><Attribute name = "ensembl_peptide_id"/>\
            </Dataset></Query>\'';
            util.unix(cmd);
            
        
        df=util.read_csv(ensembl_file, sep='\t', names=['Ensembl_Gene','Ensembl_Trans','Ensembl_Prot']);
        df=df.iloc[np.where(df["Ensembl_Trans"].notnull())[0]]
        df=df.drop_duplicates();
        
        df["Tax_ID"]=tax_id;
        
        return df;
    
    def ensembl_trans2gene_map(self):             
        df = pd.concat([self.ensembl_trans2gene_map_by_taxid(9606),self.ensembl_trans2gene_map_by_taxid(10090),self.ensembl_trans2gene_map_by_taxid(10116)]);
        map = {};
        for index, el in df.iterrows():
            map[el["Ensembl_Trans"]] = {"Ensembl_Gene":el["Ensembl_Gene"], "Ensembl_Prot":el["Ensembl_Prot"]}
        
        print "Ensembl ID mapping data created..."
        return map;
    
    def ensembl_idmapping_by_taxid(self, tax_id):     
        if tax_id == 9606:
            db_name="hsapiens_gene_ensembl";
            fname="ensembl_genes_human_idmap";
            ucsc_attr='<Attribute name = "ucsc"/>'
        elif tax_id == 10090:
            db_name="mmusculus_gene_ensembl";
            fname="ensembl_genes_mouse_idmap";
            ucsc_attr='<Attribute name = "ucsc"/>'
        elif tax_id == 10116:
            db_name="rnorvegicus_gene_ensembl";
            fname="ensembl_genes_rat_idmap";
            ucsc_attr='';


        ensembl_file = SyncDB.DOWNLOAD_DIR() + "/ensembl_files/%s"%fname; 
        
        if not os.path.exists(ensembl_file + '_p1.csv'):
            cmd = 'wget -O ' + ensembl_file + '_p1.csv \'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0"\
            encoding="UTF-8"?><!DOCTYPE Query><Query virtualSchemaName = "default"\
            formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "' + db_name + '" interface = "default">\
            <Attribute name = "ensembl_gene_id"/><Attribute name = "ensembl_transcript_id"/><Attribute name = "ensembl_peptide_id"/><Attribute name = "refseq_mrna"/><Attribute name = "refseq_peptide"/>' +  ucsc_attr + '</Dataset></Query>\'';
            util.unix(cmd);
            
        if not os.path.exists(ensembl_file + '_p2.csv'):        
            cmd = 'wget -O ' + ensembl_file + '_p2.csv \'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0"\
            encoding="UTF-8"?><!DOCTYPE Query><Query virtualSchemaName = "default"\
            formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "' + db_name + '" interface = "default">\
            <Attribute name = "ensembl_transcript_id"/><Attribute name = "uniprot_swissprot"/><Attribute name = "entrezgene"/>\
            </Dataset></Query>\'';
            util.unix(cmd);        
       
        if ucsc_attr == '':
            df1=util.read_csv(ensembl_file + '_p1.csv', sep='\t', names=['Ensembl_Gene','Ensembl_Trans','Ensembl_Prot','Refseq_RNA','Refseq_Prot']);
        else:
            df1=util.read_csv(ensembl_file + '_p1.csv', sep='\t', names=['Ensembl_Gene','Ensembl_Trans','Ensembl_Prot','Refseq_RNA','Refseq_Prot','UCSC_ID']);
            df1["UCSC_ID"] = '-';
        df2=util.read_csv(ensembl_file + '_p2.csv', sep='\t', names=['Ensembl_Trans','Uniprot_ID','GeneID']).drop_duplicates(subset='Ensembl_Trans', take_last=True);
        df = pd.merge(df1, df2, on='Ensembl_Trans', how='left');
        df=df.iloc[np.where(df["GeneID"].notnull())[0]]
        df["Tax_ID"]=tax_id;
        df["Source"] = "Ensembl";

        return df.drop_duplicates();
        
    def ensembl_idmapping_data(self):             
        df = pd.concat([self.ensembl_idmapping_by_taxid(9606),self.ensembl_idmapping_by_taxid(10090),self.ensembl_idmapping_by_taxid(10116)]);
        print "Ensembl ID mapping data processed..."
        return df;


    def get_ensembl_connection(self, db):
        import MySQLdb as mysql
        return mysql.connect('ensembldb.ensembl.org', 'anonymous', '', db)

    def ensembl_mapping(self, tax_id):
        ens_db = EnsemblDownload.get_ens_dbname_by_taxid(tax_id)
        if ens_db is None:
            print "Unsupported organism: %s" % tax_id
            return None;

        ensembl_file = SyncDB.DOWNLOAD_DIR() + "/ensembl_files/ensembl_isoform_id_map_%s.csv"%tax_id;

        print "Creating ensembl map for %s"%tax_id
        if not os.path.exists(ensembl_file):
            con = self.get_ensembl_connection(ens_db)
            try:
                #t_gt = db.from_sql(con, "select t1.ENSG, t1.ENST, t1.appris APPRIS, t2.stable_id ENSP  from (SELECT g.stable_id ENSG,t.stable_id ENST, t.transcript_id tid, ta.value appris from gene g,transcript t, transcript_attrib ta,attrib_type at where g.gene_id=t.gene_id and t.transcript_id=ta.transcript_id and ta.attrib_type_id=at.attrib_type_id and at.code='appris' ) t1 left join  translation t2 on t1.tid = t2.transcript_id")

                t_gt = db.from_sql(con, "select t1.ENSG, t1.ENST, t2.stable_id ENSP  from (SELECT g.stable_id ENSG,t.stable_id ENST, t.transcript_id tid from gene g,transcript t where g.gene_id=t.gene_id ) t1 left join  translation t2 on t1.tid = t2.transcript_id")

                t_appris = db.from_sql(con, "SELECT distinct t.stable_id ENST, ta.value APPRIS from transcript t, transcript_attrib ta,attrib_type at where t.transcript_id=ta.transcript_id and ta.attrib_type_id=at.attrib_type_id and at.code='appris'")

                # t_gt=pd.read_csv('gene_trans.csv', header=False, names=['ENSG','ENST'])
                t_g = db.from_sql(con, """
                          SELECT g.stable_id ENSG,x.dbprimary_acc GeneID FROM object_xref oxr, xref x,external_db xdb, gene g WHERE oxr.xref_id=x.xref_id AND x.external_db_id=xdb.external_db_id AND oxr.ensembl_id=g.gene_id AND oxr.ensembl_object_type='Gene' and x.dbprimary_acc is not null AND (xdb.db_name='EntrezGene')""")
                t_g['GeneID'] = t_g['GeneID'].astype(str)
                t_t = db.from_sql(con,
                    "SELECT transcript.stable_id ENST, xref.display_label ACCESSION,external_db.db_name DB FROM transcript, object_xref, xref,external_db WHERE transcript.transcript_id = object_xref.ensembl_id AND object_xref.ensembl_object_type = 'Transcript' AND object_xref.xref_id = xref.xref_id AND xref.external_db_id = external_db.external_db_id")
                t_ucsc = t_t[t_t.DB == 'UCSC'].copy()
                t_ucsc.drop(['DB'], axis=1, inplace=True)
                t_ucsc.rename2({'ACCESSION': 'UCSC_ID'})

                t_refseq_rna = t_t[ (t_t.DB == 'RefSeq_mRNA') | (t_t.DB == 'RefSeq_ncRNA') | (t_t.DB == 'RefSeq_rna') | (t_t.DB == 'RefSeq_mRNA_predicted') | (t_t.DB == 'RefSeq_mRNA_predicted')| (t_t.DB == 'RefSeq_ncRNA_predicted')| (t_t.DB == 'RefSeq_rna_predicted')].copy()
                t_refseq_rna.drop(['DB'], axis=1, inplace=True)
                t_refseq_rna.rename2({'ACCESSION': 'Refseq_RNA'})
                t_refseq_rna['Refseq_RNA'] = t_refseq_rna['Refseq_RNA'].map(lambda x: x.split('.')[0])


                t_p = db.from_sql(con, "SELECT transcript.stable_id ENST, xref.display_label ACCESSION, external_db.db_name DB FROM transcript, translation, object_xref, xref,external_db WHERE translation.translation_id=object_xref.ensembl_id and transcript.transcript_id = translation.transcript_id AND object_xref.xref_id = xref.xref_id AND xref.external_db_id = external_db.external_db_id and external_db.db_name in ('Uniprot/SWISSPROT','RefSeq_peptide')")
                '''
                t_refseq_prot = t_p[t_p.DB == 'RefSeq_peptide'].copy()
                t_refseq_prot.drop(['DB'], axis=1, inplace=True)
                t_refseq_prot.rename2({'ACCESSION': 'Refseq_Prot'})
                '''
                t_uniprot = t_p[(t_p.DB == 'Uniprot/SWISSPROT')].copy()
                t_uniprot.drop(['DB'], axis=1, inplace=True)
                t_uniprot.rename2({'ACCESSION': 'Uniprot_ID'})

                t_gt = t_gt.merge(t_appris, left_on='ENST', right_on='ENST', how='left')
                t_g = t_g.merge(t_gt, left_on='ENSG', right_on='ENSG', how='left')
                t_g = t_g.merge(t_ucsc, left_on='ENST', right_on='ENST', how='left')
                t_g = t_g.merge(t_refseq_rna, left_on='ENST', right_on='ENST', how='left')
                t_g = t_g.merge(t_uniprot, left_on='ENST', right_on='ENST', how='left')
                cur_refseq = self.refseq_idmaping.query('tax_id in [%s]'%','.join(SyncDB.SPECIE_SUBTREE[str(tax_id)]))
                t_g = cur_refseq.merge(t_g, left_on=['GeneID','Refseq_RNA'], right_on=['GeneID','Refseq_RNA'], how='outer')
                t_g.fillna('-', inplace = True)
                t_g = t_g.query("Refseq_RNA!='-' or ENSG !='-' or ENST !='-' or UCSC_ID !='-' or Uniprot_ID !='-' ")
                #t_g = t_g.merge(t_refseq_prot, left_on='ENST', right_on='ENST', how='left')
                t_g['tax_id'] = tax_id;
                t_g.to_csv(ensembl_file, index=False)
                print 'completed ', ensembl_file;
                return t_g

            except:
                print "Error in creating ensembl map for %s"%tax_id
                return None
        else:
            df = pd.read_csv(ensembl_file)
            return df

    def do_update(self):
        def get_rank(a):
            rank = '0' if a['Refseq_RNA'] == '-' else '1'
            rank += '0' if a['status'] == 'MODEL' else '1' if a['status'] == 'WGS' else '2' if a['status'] == 'INFERRED' else '3' if a['status'] == 'PREDICTED' else '3' if a['status'] == 'PROVISIONAL' else '4' if a['status'] == 'VALIDATED' else '5' if a['status'] == 'REVIEWED' else '0';
            rank += '0' if a['APPRIS'] == 'alternative2' else '1' if a['APPRIS'] == 'alternative1' else '2' if a['APPRIS'] == 'principal5' else '3' if a['APPRIS'] == 'principal4' else '4' if a['APPRIS'] == 'principal3' else '5' if a['APPRIS'] == 'principal2' else '6' if a['APPRIS'] == 'principal1' else '0';
            rank += '0' if a['ENSG'] == '-' else '1'
            rank += '0' if a['ENST'] == '-' else '1'
            return  rank;

        def compare(a, b):
            key_a = get_rank(a)
            key_b = get_rank(b)
            if key_a < key_b:
                return 1
            elif key_a > key_b:
                return -1
            else:
                return 0

        self.refseq_idmaping = util.read_csv(SyncDB.DOWNLOAD_DIR() + "/CsvInMem_gene2refseq")[['GeneID','status','RNA_nucleotide_accession_version','protein_accession_version', 'tax_id']].drop_duplicates()
        self.refseq_idmaping.rename2({'RNA_nucleotide_accession_version':'Refseq_RNA', 'protein_accession_version':'Refseq_Prot'})
        self.refseq_idmaping = self.refseq_idmaping.query("status!='SUPPRESSED' and tax_id in [%s]"%','.join(self.taxidList))
        self.refseq_idmaping = self.refseq_idmaping.fillna('-')
        self.refseq_idmaping['GeneID'] = self.refseq_idmaping['GeneID'].astype(str)
        self.refseq_idmaping['Refseq_RNA'] = self.refseq_idmaping['Refseq_RNA'].map(lambda x: x.split('.')[0])
        self.refseq_idmaping['Refseq_Prot'] = self.refseq_idmaping['Refseq_Prot'].map(lambda x: x.split('.')[0])
        #Tracer()()

        dfa=[]
        stat={}
        #self.taxidList = ['559292'] #DEBUG
        for tax_id in self.taxidList:
            map = self.ensembl_mapping(tax_id)
            if map is None:
                continue
            stat[tax_id]={'gid':0,'iso':0, 'trans':0, 'refseq':0, 'uniprot':0}
            dfa.append(map)

        df = pd.concat(dfa)
        df=df.fillna('-')
        data=[]

        for g, row in df.groupby('GeneID'):
            iso = sorted(row.T.to_dict().values(),cmp=compare)
            content=[]
            tax_id=str(row.get('tax_id').values[0])
            stat[tax_id]['gid'] += 1
            stat[tax_id]['iso'] += len(iso)
            for el in iso:
                content.append(','.join([
                                'Refseq_RNA:'+el['Refseq_RNA'],
                                'Refseq_Prot:'+el['Refseq_Prot'],
                                'Ensembl_Trans:'+el['ENST'],
                                'Ensembl_Gene:'+el['ENSG'],
                                'Ensembl_Prot:'+el['ENSP'],
                                'UCSC_ID:' + el['UCSC_ID'],
                                'Uniprot_ID:'+el['Uniprot_ID'],
                                'APPRIS:' + el['APPRIS'],
                                'Status:' + el['status']
                                ]));
            data.append({'gid':g, 'content':'|'.join(content), 'annotation_field1':str(len(content)), 'tax_id':tax_id})

        for t in stat:
            print "%s-Genes:%d, ISO:%d"%(t, stat[t]['gid'], stat[t]['iso'])

        pd.DataFrame(data).to_csv(self.fn_dest, index=False, sep=',');

    def do_update_old(self):
        ucsc_data = self.ucsc_idmapping_data().fillna('-');
        print 'UCSC data done...'
    
        ncbi_data = self.ncbi_idmapping_data().fillna('-');
        print 'NCBI data done...'    

        ensembl_data=self.ensembl_idmapping_data().fillna('-');
        print 'ENSEMBL data done...'

        merged_data = {};
        ens_trans = set([])
        refsq_rna = set([])
        isoforms = 0;
        for index, el in ucsc_data.iterrows():
            if el["GeneID"] not in merged_data:
                merged_data[el["GeneID"]] = [];

            merged_data[el["GeneID"]].append({
                                'Refseq_RNA':el['Refseq_RNA'],
                                'Refseq_Prot':el['Refseq_Prot'],
                                'UCSC_ID':el['UCSC_ID'],
                                'Ensembl_Trans':el['Ensembl_Trans'],
                                'Ensembl_Gene':el['Ensembl_Gene'],
                                'Ensembl_Prot':el['Ensembl_Prot'],
                                'Uniprot_ID':el['Uniprot_ID'],
                                'Tax_ID':el['Tax_ID'],
                                'Source':'UCSC',
                                });

            ens_trans.add(el['Ensembl_Trans']);
            refsq_rna.add(el['Refseq_RNA']);

        print 'Merging ENSEMBL Data';
        for index, el in ensembl_data.iterrows():

            if el["GeneID"] not in merged_data:               
                refsq_rna.add(el['Refseq_RNA']);
                merged_data[el["GeneID"]]= [{
                                    'Refseq_RNA':el['Refseq_RNA'],
                                    'Refseq_Prot':el['Refseq_Prot'],
                                    'UCSC_ID':'-',
                                    'Ensembl_Trans':el['Ensembl_Trans'],
                                    'Ensembl_Gene':el['Ensembl_Gene'],
                                    'Ensembl_Prot':el['Ensembl_Prot'],
                                    'Uniprot_ID':el['Uniprot_ID'],
                                    'Tax_ID':el['Tax_ID'],
                                    'Source':'ENSEMBL',                                    
                                    }];      
            else:
                if el["Ensembl_Trans"] not in ens_trans:
                    refsq_rna.add(el['Refseq_RNA']);
                    merged_data[el["GeneID"]].append({
                                'Refseq_RNA':el['Refseq_RNA'],
                                'Refseq_Prot':el['Refseq_Prot'],
                                'UCSC_ID':'-',
                                'Ensembl_Trans':el['Ensembl_Trans'],
                                'Ensembl_Gene':el['Ensembl_Gene'],
                                'Ensembl_Prot':el['Ensembl_Prot'],
                                'Uniprot_ID':el['Uniprot_ID'],
                                'Tax_ID':el['Tax_ID'],
                                'Source':'ENSEMBL',                                                                    
                                });
                                
        print 'Merging NCBI Data';                        
        for index, el in ncbi_data.iterrows():
            if el["GeneID"] not in merged_data:               
                merged_data[el["GeneID"]]= [{
                                    'Refseq_RNA':el['Refseq_RNA'],
                                    'Refseq_Prot':el['Refseq_Prot'],
                                    'UCSC_ID':'-',
                                    'Ensembl_Trans':el['Ensembl_Trans'],
                                    'Ensembl_Gene':el['Ensembl_Gene'],
                                    'Ensembl_Prot':el['Ensembl_Prot'],
                                    'Uniprot_ID':el['Uniprot_ID'],
                                    'Tax_ID':el['Tax_ID'],
                                    'Source':'NCBI',                                    
                                    }];

            else:
                if el["Refseq_RNA"] not in refsq_rna:
                    merged_data[el["GeneID"]].append({
                                'Refseq_RNA':el['Refseq_RNA'],
                                'Refseq_Prot':el['Refseq_Prot'],
                                'UCSC_ID':'-',
                                'Ensembl_Trans':'-',
                                'Ensembl_Gene':'-',
                                'Ensembl_Prot':'-',
                                'Uniprot_ID':el['Uniprot_ID'],
                                'Tax_ID':el['Tax_ID'],
                                'Source':'NCBI',                                                                    
                                });

        data=[]
        isoforms = {"9606":0, "10090":0, "10116":0}
        def compare(a, b):
            key_a = '0' if a['Refseq_RNA'] =='-' else '1'
            key_a += '0' if a['Ensembl_Trans'] =='-' else '1'
            key_a += '0' if a['UCSC_ID'] =='-' else '1'
            key_a += '0' if a['Uniprot_ID'] =='-' else '1'

            key_b = '0' if b['Refseq_RNA'] =='-' else '1'
            key_b += '0' if b['Ensembl_Trans'] =='-' else '1'
            key_b += '0' if b['UCSC_ID'] =='-' else '1'
            key_b += '0' if b['Uniprot_ID'] =='-' else '1'
            if key_a < key_b:
                return 1
            elif key_a > key_b:
                return -1
            else:
                return 0
                
        for gid in merged_data:
            content=[];
            merged_data[gid] = sorted(merged_data[gid],cmp=compare)
            for el in merged_data[gid]:
                tax_id = el['Tax_ID']
                content.append(','.join([
                                'UCSC_ID:'+el['UCSC_ID'],                
                                'Refseq_RNA:'+el['Refseq_RNA'],
                                'Refseq_Prot:'+el['Refseq_Prot'],
                                'Ensembl_Trans:'+el['Ensembl_Trans'],
                                'Ensembl_Gene:'+el['Ensembl_Gene'],
                                'Ensembl_Prot:'+el['Ensembl_Prot'],
                                'Uniprot_ID:'+el['Uniprot_ID'],
                                'Source:'+el['Source']]));
                isoforms[str(tax_id)] += 1;                
            if (len(content) != 0):
                data.append({'gid':gid, 'content':'|'.join(content), 'annotation_field1':'', 'tax_id':tax_id});
            
        pd.DataFrame(data).to_csv(self.fn_dest, index=False, sep=',');
        print 'Human isoforms: ', isoforms["9606"];
        print 'Mouse isoforms: ', isoforms["10090"];
        print 'Rat isoforms: ', isoforms["10116"];
        
    def check_inputs (self):
        return True;
                
if __name__=='__main__':

    iso = IsoformsIdMap()
    iso.fn_dest = SyncDB.DOWNLOAD_DIR() + "/isoform.csv"
    iso.do_update()
    exit();
    tmpdir = '../id_files'
    #UCSC
    if not os.path.exists(tmpdir + '/ucsc_ids_test.csv'):
    
        cmd = 'mysql --host=genome-mysql.cse.ucsc.edu --user=genomep --password=password --database=hg38 -e "SELECT  knownGene.*, kgXref.refseq as Refseq_RNA, knownToLocusLink.value as GeneID, knownGene.name as UCSC_ID, knownGene.proteinID as Uniprot_ID, knownGene.alignID as Ensembl_Trans from kgXref,knownGene,knownToLocusLink WHERE knownGene.name=knownToLocusLink.name and kgXref.kgID=knownGene.name;" > ' + path.join(tmpdir, "ucsc_ids_test.csv");
        util.unix(cmd);        
    ucsc_data=util.read_csv(tmpdir+'/ucsc_ids_test.csv', sep='\t');
    f=lambda x: x.split('.')[0];
    ucsc_data['Ensembl_Trans']=ucsc_data['Ensembl_Trans'].map(f);
    print "UCSC IDs done";
    
    #NCBI
    if not os.path.exists(tmpdir + '/ncbi_id_data.csv'):
        #Tracer()()
        gene_info = pd.read_csv(SyncDB.DOWNLOAD_DIR()+"/gene_info.csv")[["gid","type_of_gene"]]           
        gene2ref=util.read_csv(SyncDB.DOWNLOAD_DIR()+'/gene2refseq.gz', sep='\t', header=False, skiprows=1, names=['Tax_ID', 'GeneID', 'status','Refseq_RNA','RNA_nucleotide_gi','Refseq_Prot','protein_gi','genomic_nucleotide_accession.version','genomic_nucleotide_gi','start_position_on_the_genomic_accession','end_position_on_the_genomic_accession','orientation','assembly','mature_peptide_accession.version','mature_peptide_gi','Symbol']).query('Tax_ID == 9606 and Refseq_RNA != "-"')[['GeneID', 'Refseq_RNA','Refseq_Prot']].drop_duplicates();        
        gene2ref = gene2ref[gene2ref["Refseq_RNA"].str.contains("X")==False] #Include only the validated entries. "NM_xx and NR_xxx"
        gene2ref=pd.merge(gene2ref, gene_info, how='left', left_on='GeneID', right_on='gid').query('type_of_gene!=["pseudo","other"]')
        
        ref2uniprot=util.read_csv(SyncDB.DOWNLOAD_DIR()+'/gene_refseq_uniprotkb_collab.gz', sep='\t', header=False, skiprows=1, names=['Refseq_Prot','Uniprot_ID']); 
        gene2ens=util.read_csv(SyncDB.DOWNLOAD_DIR()+'/gene2ensembl.gz', sep='\t', header=False, skiprows=1, names=['Tax_ID','GeneID','Ensembl_Gene','Refseq_RNA','Ensembl_Trans','Refseq_Prot','Ensembl_Prot']).query('Tax_ID == 9606 and Refseq_RNA != "-"')[['Ensembl_Gene','Refseq_RNA','Ensembl_Trans']].drop_duplicates();

        f=lambda x: x.split('.')[0];
        gene2ref['Refseq_Prot']=gene2ref['Refseq_Prot'].map(f);
        gene2ref['Refseq_RNA']=gene2ref['Refseq_RNA'].map(f);
        gene2ens['Refseq_RNA']=gene2ens['Refseq_RNA'].map(f);

        ncbi_data=pd.merge(gene2ref, gene2ens, how='left', on='Refseq_RNA');    
        ncbi_data=pd.merge(ncbi_data, ref2uniprot, how='left', on='Refseq_Prot');    
        ncbi_data=ncbi_data.iloc[np.where(ncbi_data["GeneID"].notnull())[0]]
        ncbi_data.to_csv(tmpdir+'/ncbi_id_data.csv', index=False, sep=',');
    else:
        ncbi_data = util.read_csv(tmpdir+'/ncbi_id_data.csv');    
    print "NCBI IDs done";
    
    #ENSEMBL
    if not os.path.exists(tmpdir + '/ensembl_id_data.csv'):
        db_name="hsapiens_gene_ensembl";
        fname="ensembl_genes_human_idmap";
        ensembl_file = SyncDB.DOWNLOAD_DIR() + "/ensembl_files/%s"%fname;     
        if not os.path.exists(ensembl_file + '_p1.csv'):
            cmd = 'wget -O ' + ensembl_file + '_p1.csv \'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0"\
            encoding="UTF-8"?><!DOCTYPE Query><Query virtualSchemaName = "default"\
            formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "' + db_name + '" interface = "default">\
            <Attribute name = "ensembl_gene_id"/><Attribute name = "ensembl_transcript_id"/><Attribute name = "ensembl_peptide_id"/><Attribute name = "refseq_mrna"/><Attribute name = "refseq_ncrna"/><Attribute name = "ucsc"/>\
            </Dataset></Query>\'';
            util.unix(cmd);
        if not os.path.exists(ensembl_file + '_p2.csv'):        
            cmd = 'wget -O ' + ensembl_file + '_p2.csv \'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0"\
            encoding="UTF-8"?><!DOCTYPE Query><Query virtualSchemaName = "default"\
            formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "' + db_name + '" interface = "default">\
            <Attribute name = "ensembl_gene_id"/><Attribute name = "ensembl_transcript_id"/><Attribute name = "uniprot_swissprot"/><Attribute name = "entrezgene"/><Attribute name = "hgnc_symbol"/>\
            </Dataset></Query>\'';
            util.unix(cmd);        
            
        df1=util.read_csv(ensembl_file + '_p1.csv', sep='\t', names=['Ensembl_Gene','Ensembl_Trans','Ensembl_Prot','Refseq_RNA','Refseq_Prot','UCSC_ID']);
        df2=util.read_csv(ensembl_file + '_p2.csv', sep='\t', names=['Ensembl_Trans','Uniprot_ID','GeneID','Gene_Symbol']);
        ens_data = pd.merge(df1, df2, on='Ensembl_Trans', how='left');
        ens_data=ens_data.iloc[np.where(ens_data["GeneID"].notnull())[0]]
        ens_data.to_csv(tmpdir+'/ensembl_id_data.csv', index=False, sep=',');
    else:
        ens_data = util.read_csv(tmpdir+'/ensembl_id_data.csv');
    print "ENSEBL IDs done";    
    
    '''
    ucsc_uniq_values = {}
    for c in ucsc_data.columns: 
        ucsc_uniq_values[c] = pd.unique(ucsc_data[c])
        
    ens_uniq_values = {}
    for c in ens_data.columns: 
        ens_uniq_values[c] = pd.unique(ens_data[c])

    ncbi_uniq_values = {}
    for c in ncbi_data.columns: 
        ncbi_uniq_values[c] = pd.unique(ncbi_data[c])
        
        
    for c in ["GeneID", "Refseq_RNA", "Ensembl_Trans", "UCSC_ID"]:
        print c
        print "\tUCSC Counts:%d"%(ucsc_uniq_values[c].size)
        print "\tENSE Counts:%d"%(ens_uniq_values[c].size)
        if c != 'UCSC_ID':
            print "\tNCBI Counts:%d"%(ncbi_uniq_values[c].size)
        print "\tUCSC & ENSE:%d"%(len(set(ucsc_uniq_values[c]) & set(ens_uniq_values[c])))
        if c != 'UCSC_ID':
            print "\tUCSC & NCBI:%d"%(len(set(ucsc_uniq_values[c]) & set(ncbi_uniq_values[c])))
            print "\tNCBI& ENSE:%d"%(len(set(ncbi_uniq_values[c]) & set(ens_uniq_values[c])))
            print "\tUCSC & ENSE & NCBI:%d"%(len(set(ucsc_uniq_values[c]) & set(ens_uniq_values[c])& set(ncbi_uniq_values[c])))
    '''
    gene_id_trans_per_source = {};
    
    for g, r in ucsc_data.groupby('GeneID'):
        gene_id_trans_per_source[g] = {"UCSC": len(r), "NCBI":0, "ENSEMBL":0};
        
    for g, r in ncbi_data.groupby('GeneID'):
        if g in gene_id_trans_per_source:
            gene_id_trans_per_source[g]["NCBI"] = len(r);
        else:
            gene_id_trans_per_source[g] = {"UCSC":0 , "NCBI":len(r), "ENSEMBL":0};
            
    for g, r in ens_data.groupby('GeneID'):
        if g in gene_id_trans_per_source:
            gene_id_trans_per_source[g]["ENSEMBL"] = len(r);
        else:
            gene_id_trans_per_source[g] = {"UCSC":0 , "NCBI":0, "ENSEMBL":len(r)};
        
    gene_rank = {}
    
    def calculate_rank(t):
        t1=sorted([{"s":x, "c":t[x]} for x in ["UCSC","NCBI","ENSEMBL"]], key=lambda k:-k["c"]);
        #Tracer()()
        rank={};
        if t1[0]["c"] != t1[1]["c"]:
            rank[t1[0]["s"]] = 1;
            if t1[1]["c"] != t1[2]["c"]:
                rank[t1[1]["s"]] = 2;
                rank[t1[2]["s"]] = 3;
            else:
                rank[t1[1]["s"]] = 2.5;
                rank[t1[2]["s"]] = 2.5;            
        else:
            if t1[1]["c"] != t1[2]["c"]:
                rank[t1[0]["s"]] = 1.5;
                rank[t1[1]["s"]] = 1.5;
                rank[t1[2]["s"]] = 3;
            else:
                rank[t1[0]["s"]] = 2;            
                rank[t1[1]["s"]] = 2;
                rank[t1[2]["s"]] = 2;            
        
        return rank;
        
    for g in gene_id_trans_per_source:
        gene_id_trans_per_source[g]["total"] = gene_id_trans_per_source[g]["UCSC"] +gene_id_trans_per_source[g]["NCBI"] +gene_id_trans_per_source[g]["ENSEMBL"] ;
        gene_rank[g] = calculate_rank(gene_id_trans_per_source[g]);

    sdata =sorted([gene_id_trans_per_source[x] for x in gene_id_trans_per_source], key=lambda k:-k["total"]);    
    
    rank_array = pd.DataFrame([gene_rank[x] for x in gene_rank]);
    trans_count = pd.DataFrame(sdata);
    print "Gene covered (UCSC):",len(trans_count.query('UCSC!=0'));
    print "Gene covered (NCBI):",len(trans_count.query('NCBI!=0'));
    print "Gene covered (ENSEMBL):",len(trans_count.query('ENSEMBL!=0'));
    
    print "Average Transcript Count (UCSC):",np.mean(trans_count["UCSC"]);
    print "Average Transcript Count (NCBI):",np.mean(trans_count["NCBI"]);
    print "Average Transcript Count (ENSEMBL):",np.mean(trans_count["ENSEMBL"]);
    
    print "Average Transcript Rank (UCSC):",np.mean(rank_array["UCSC"]);
    print "Average Transcript Rank (NCBI):",np.mean(rank_array["NCBI"]);
    print "Average Transcript Rank (ENSEMBL):",np.mean(rank_array["ENSEMBL"]);
    
    print "Number of genes with most transcript (UCSC):",len(rank_array.query('UCSC==1'));
    print "Number of genes with most transcript (NCBI):",len(rank_array.query('NCBI==1'));
    print "Number of genes with most transcript (ENSEMBL):",len(rank_array.query('ENSEMBL==1'));    
    
    Tracer()()
    '''
    ind = np.arange(len(sdata))
    import matplotlib.pyplot as plt

    plt.clf()
    p1 = plt.bar(ind, (x["UCSC"] for x in sdata), color='r')
    p2 = plt.bar(ind, (x["NCBI"] for x in sdata), color='g')
    p3 = plt.bar(ind, (x["ENSEMBL"] for x in sdata), color='b')
    plt.savefig(tmpdir + "/isoform_count_plot.png")

    Tracer()()
    '''
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
import db

from IPython.core.debugger import Tracer

class EnsemblDownload(XmlClass):
    ENSEMBL_VERSION = '93'
    def __init__(self, xe=None):
        XmlClass.__init__(self,xe=xe)
        EnsemblDownload.set_ENSEMBL_VERSION()
        self.tag = "EnsemblDownload"
        self.fn_annotations=os.path.join(SyncDB.DOWNLOAD_DIR(),"ensembl_annotations.csv")
        #self.fn_phenotype=os.path.join(SyncDB.DOWNLOAD_DIR(),"ensembl_phenotype_annotations.csv")
        self.fn_ensembl_id_map=os.path.join(SyncDB.DOWNLOAD_DIR(),"ensembl_id_map.csv")
        self.taxidList = SyncDB.SUPPORTED_SPECIES;
        for child in self.children:
            if type(child).__name__ == "Species":
                self.taxidList = child.supported_species
                break;

        self.chrList = [str(x+1) for x in range(22)];
        self.chrList += ['X','Y','MT']
        util.unix('mkdir -p ' + SyncDB.DOWNLOAD_DIR() + '/ensembl_files');
        self.mart_ftp = "ftp://ftp.ensembl.org/pub/release-84/mysql/ensembl_mart_84/"
        #self.taxidList=["9606"]

    @staticmethod
    def set_ENSEMBL_VERSION():
        import MySQLdb as mysql
        con = mysql.connect('ensembldb.ensembl.org', 'anonymous', '')
        cursor = con.cursor()
        cursor.execute('''
                        SELECT 
                        SCHEMA_NAME FROM INFORMATION_SCHEMA.SCHEMATA 
                        WHERE 
                        SCHEMA_NAME LIKE 'homo_sapiens_core_%';
                            ''')
        tables = cursor.fetchall()
        tables_names = [x[0] for x in tables]
        all_versions = []
        for x in tables_names:
            import re
            m = re.search('homo_sapiens_core_(\d*)_*', x)
            if m:
                all_versions.append(int(m.group(1)))
        EnsemblDownload.ENSEMBL_VERSION = str(max(all_versions))


    def populate_more(self,root):
        self.outputs = [self.fn_annotations,self.fn_ensembl_id_map]

    def mapping(self, tax_id):
        ens_db = EnsemblDownload.get_ens_dbname_by_taxid(tax_id)
        if ens_db is None:
            print "Unsupported organism: %s" % tax_id
            return None;

        print "Creating ensembl map for %s"%tax_id
        con = self.get_ensembl_connection(ens_db)
        try:
            #Get all the isoforms available in Ensembl (ENSG,ENST,ENSP map.)
            t_gt = db.from_sql(con, "select t1.ENSG, t1.ENST, t2.stable_id ENSP  from (SELECT g.stable_id ENSG,t.stable_id ENST, t.transcript_id tid from gene g,transcript t where g.gene_id=t.gene_id) t1 left join  translation t2 on t1.tid = t2.transcript_id")
			
            #ENSG to EntrezGeneID map.
            t_g = db.from_sql(con, """
        SELECT g.stable_id ENSG,syn.synonym SYMBOL,xdb.db_name DB,'SYN' TYPE FROM object_xref oxr, xref x, external_db xdb,external_synonym syn,gene g WHERE oxr.xref_id=x.xref_id AND x.xref_id=syn.xref_id AND x.external_db_id=xdb.external_db_id AND oxr.ensembl_id=g.gene_id AND oxr.ensembl_object_type='Gene' and syn.synonym is not null
        UNION
        SELECT g.stable_id ENSG,x.display_label SYMBOL,xdb.db_name DB,'DIS' TYPE FROM object_xref oxr, xref x,external_db xdb, gene g WHERE oxr.xref_id=x.xref_id AND x.external_db_id=xdb.external_db_id AND oxr.ensembl_id=g.gene_id AND oxr.ensembl_object_type='Gene' and x.display_label is not null
        UNION
        SELECT g.stable_id ENSG,x.dbprimary_acc SYMBOL,xdb.db_name DB,'PRI' Type FROM object_xref oxr, xref x,external_db xdb, gene g WHERE oxr.xref_id=x.xref_id AND x.external_db_id=xdb.external_db_id AND oxr.ensembl_id=g.gene_id AND oxr.ensembl_object_type='Gene' and x.dbprimary_acc is not null AND (xdb.db_name='EntrezGene')""")
            t_g = t_g[(t_g.DB == 'EntrezGene') & (t_g.TYPE == 'PRI')].copy()
            t_g.drop(['DB', 'TYPE'], axis=1, inplace=True)
            t_g.rename2({'SYMBOL': 'GeneID'})
            t_g['GeneID'] = t_g['GeneID'].astype(str)

            #UCSC to ENST map.
            t_t = db.from_sql(con,
                "SELECT transcript.stable_id ENST, xref.display_label ACCESSION,external_db.db_name DB FROM transcript, object_xref, xref,external_db WHERE transcript.transcript_id = object_xref.ensembl_id AND object_xref.ensembl_object_type = 'Transcript' AND object_xref.xref_id = xref.xref_id AND xref.external_db_id = external_db.external_db_id")
            t_ucsc = t_t[t_t.DB == 'UCSC'].copy()
            t_ucsc.drop(['DB'], axis=1, inplace=True)
            t_ucsc.rename2({'ACCESSION': 'UCSC'})
            t_g = t_g.merge(t_gt, left_on='ENSG', right_on='ENSG', how='left')
            t_g = t_g.merge(t_ucsc, left_on='ENST', right_on='ENST', how='left')
			
            #t_g is a complete map between ENSG,ENST,ENSP,UCSC and EntrezGeneID
            return t_g

        except:
            print "Error in creating ensembl map for %s"%tax_id
            return None

    def get_dbname_by_taxid(self, tax_id):
        map={
        "9606":
            "hsapiens_gene_ensembl",
        "10090":
            "mmusculus_gene_ensembl",
        "10116":
            "rnorvegicus_gene_ensembl",
        "7955":
            "drerio_gene_ensembl",
        "9031":
            "ggallus_gene_ensembl",
        "9544":
            "mmulatta_gene_ensembl",
        "7227":
            "dmelanogaster_gene_ensembl",
        "6239":
            "celegans_gene_ensembl",
        "4932":
            "scerevisiae_gene_ensembl"
        }

        if tax_id in map:
            return map[tax_id]
        return None

    @staticmethod
    def get_ens_dbname_by_taxid(tax_id):
        map={
        "9606":
            "homo_sapiens_core_{0}_".format(EnsemblDownload.ENSEMBL_VERSION), #38
        "10090":
            "mus_musculus_core_{0}_".format(EnsemblDownload.ENSEMBL_VERSION),
        "10116":
            "rattus_norvegicus_core_{0}_".format(EnsemblDownload.ENSEMBL_VERSION),
        "7955":
            "danio_rerio_core_{0}_".format(EnsemblDownload.ENSEMBL_VERSION),
        "9031":
            "gallus_gallus_core_{0}_".format(EnsemblDownload.ENSEMBL_VERSION),
        "9544":
            "macaca_mulatta_core_{0}_".format(EnsemblDownload.ENSEMBL_VERSION),
        "7227":
            "drosophila_melanogaster_core_{0}_".format(EnsemblDownload.ENSEMBL_VERSION),
        "6239":
            "caenorhabditis_elegans_core_{0}_".format(EnsemblDownload.ENSEMBL_VERSION),
        "4932":
            "saccharomyces_cerevisiae_core_{0}_".format(EnsemblDownload.ENSEMBL_VERSION)
        }

        if tax_id not in map:
            return None
        prefix =  map[tax_id]
        return  EnsemblDownload.get_ensembl_latest_version(prefix)

    @staticmethod
    def get_ensembl_latest_version(prefix):
        import MySQLdb as mysql
        con = mysql.connect('ensembldb.ensembl.org', 'anonymous', '')
        cursor = con.cursor()
        cursor.execute('show databases')
        tables = cursor.fetchall()
        return [x[0] for x in tables if prefix in x[0]][0]

    #deleted by Ali Begin
    def get_ensembl2gid_df_web (self, tax_id, type):
        #mmusculus_gene_ensembl   10090
        #rnorvegicus_gene_ensembl 10116
        print "Get %s to gene id for %s"%(type, tax_id)
        attr  = '<Attribute name = "'+type+'" />'
        db_name = self.get_dbname_by_taxid(tax_id)
        if db_name is None:
            return None

        fname='ensembl2gid_%s_%s'%(type, tax_id)

        util.unix('mkdir -p ' + SyncDB.DOWNLOAD_DIR() + '/ensembl_files')
        valid_files= []
        for chr in self.chrList:
            ensembl_file = SyncDB.DOWNLOAD_DIR() + "/ensembl_files/%s_chr%s"%(fname, chr)
            print "downloading %s from %s for chr %s..."%(type,db_name,chr)
            if not os.path.exists(ensembl_file):
                cmd = 'wget -O ' + ensembl_file + ' \'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0"\
                encoding="UTF-8"?><!DOCTYPE Query><Query virtualSchemaName = "default"\
                formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "' + db_name + '" interface = "default"><Filter name="chromosome_name" value="' + chr + '" filter_list=""/>' + attr + '<Attribute name = "entrezgene"\
                /></Dataset></Query>\'';
                util.unix(cmd);
            try:
                tdf = util.read_csv(ensembl_file, sep="\t", header=None, names=['source_id', 'gid'], nrows=1)
            except Exception as exp:
                tdf = pd.DataFrame()

            if len(tdf) != 0:
                valid_files.append(ensembl_file)

        if len(valid_files) == 0:
            return None


        cmd = 'cat %s >> %s'%(' '.join(valid_files), SyncDB.DOWNLOAD_DIR() + "/ensembl_files/" + fname)
        print cmd
        util.unix (cmd)
        print "downloading %s from %s done."%(type,db_name)
        ensembl_data=util.read_csv(SyncDB.DOWNLOAD_DIR() + "/ensembl_files/" + fname, sep="\t", header=None, names=['source_id','gid'])
        ensembl_data=ensembl_data[ensembl_data['gid'].notnull()]
        ensembl_data=ensembl_data[ensembl_data['source_id'].notnull()]
        ensembl_data[['gid']]=ensembl_data[['gid']].astype(int)
        ensembl_data['tax_id']=tax_id
        ensembl_data['type_name']=type
        return ensembl_data

    def get_idmap_source_file(self, type, tax_id):
        org_name = self.get_dbname_by_taxid(tax_id)
        if org_name is None:
            return None

    def get_ensembl2gid_df_not_used (self, tax_id, type):
        print "Get %s to gene id for %s"%(type, tax_id)
        source_file_name = self.get_idmap_source_file(type, tax_id)
        if source_file_name is None:
            return None

        source_file_name = SyncDB.DOWNLOAD_DIR() + "/ensembl_files/" + source_file_name

        out_file=SyncDB.DOWNLOAD_DIR() + "/ensembl_files/" + 'ensembl2gid_%s_%s'%(type, tax_id)

        if not os.path.exists(source_file_name):
            urllib.urlretrieve(self.mart_ftp + "/" + self.get_idmap_source_file(type, tax_id), source_file_name)

        ensembl_data=util.read_csv(source_file_name, sep="\t", header=None, names=['source_id','gid']);
        ensembl_data=ensembl_data[ensembl_data['gid'].notnull()];
        ensembl_data=ensembl_data[ensembl_data['source_id'].notnull()];
        ensembl_data[['gid']]=ensembl_data[['gid']].astype(int)
        ensembl_data['tax_id']=tax_id;
        ensembl_data['type_name']=type;
        return ensembl_data

    def get_idmap_data_old(self):
        
        h_map = GPUtils.get_ensembl2gid_map()
        data=[]
        for g in h_map:
            for gid in h_map[g]:
                data.append({'source_id':g, 'gid':gid, 'tax_id':9606, 'type_name':'ensembl_gene_id'})
        dfa = [pd.DataFrame(data)];
        for tax_id in self.taxidList:
            if tax_id != '9606':
                df=self.get_ensembl2gid_df(tax_id, 'ensembl_gene_id')
                if df is not None:
                    dfa.append(df)

        
        for tax_id in self.taxidList:
            for type in ['ensembl_peptide_id','ensembl_transcript_id', 'ucsc']:
                df=self.get_ensembl2gid_df(tax_id, type)
                if df is not None:
                    df['id_status'] = None
                    dfa.append(df)

        df = pd.concat(dfa).drop_duplicates(subset=['source_id','gid','tax_id'])
        df.to_csv(self.fn_ensembl_id_map, index=False)

    #deleted by Ali End
    def get_idmap_data(self):
        dfa=[]
        for tax_id in self.taxidList:
            #get ENSG,ENST,ENSP,UCSC to GeneID map for the given tax_id.
            map = self.mapping(tax_id)
            if map is None:
                continue
            ensg = map[["GeneID", "ENSG"]].drop_duplicates()
            ensg.rename2({"GeneID":"gid", "ENSG":"source_id"})
            ensg['type_name'] = 'ensembl_gene_id'
            ensg['tax_id'] = tax_id

            enst = map[["GeneID", "ENST"]].drop_duplicates()
            enst.rename2({"GeneID":"gid", "ENST":"source_id"})
            enst['type_name'] = 'ensembl_transcript_id'
            enst['tax_id'] = tax_id

            ensp = map[["GeneID", "ENSP"]].drop_duplicates()
            ensp.rename2({"GeneID":"gid", "ENSP":"source_id"})
            ensp=ensp[ensp['source_id'].notnull()]
            ensp['type_name'] = 'ensembl_peptide_id'
            ensp['tax_id'] = tax_id

            ucsc = map[["GeneID", "UCSC"]].drop_duplicates()
            ucsc.rename2({"GeneID":"gid", "UCSC":"source_id"})
            ucsc['type_name'] = 'ucsc'
            ucsc['tax_id'] = tax_id
            dfa.append(pd.concat([ensg, enst, ensp, ucsc]))

        df = pd.concat(dfa)
        df = df[df['gid'].notnull()]
        df.to_csv(self.fn_ensembl_id_map, index=False)

    #deleted by Ali begin
    def get_phenotype_data_old(self):
        dfa=[]
        for tax_id in self.taxidList:
            df = self.get_phenotype_data_by_taxid_old(tax_id)
            if df is not None:
                dfa.append(df)

        pd.concat(dfa).to_csv(self.fn_phenotype, index=False)

    def get_phenotype_data_by_taxid_old(self, tax_id):
        print "Get phenotype for %s"%tax_id
        db_name = self.get_dbname_by_taxid(tax_id)
        if db_name is None:
            return None

        out_file = path.join(SyncDB.DOWNLOAD_DIR(), "ensembl_phenotype_%s.csv"%tax_id);

        valid_files= []
        for chr in self.chrList:
            ensembl_file = SyncDB.DOWNLOAD_DIR() + "/ensembl_files/ensembl_phenotype_%s_chr%s"%(tax_id, chr);
            print "downloading ensembl_phenotype_%s_chr%s..."%(tax_id, chr);
            if not os.path.exists(ensembl_file):
                util.unix('wget -O ' + ensembl_file + ' - \'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0"\
                encoding="UTF-8"?><!DOCTYPE Query><Query virtualSchemaName = "default"\
                formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "'+db_name + '" interface = "default"\
                ><Filter name="chromosome_name" value="' + chr + '" filter_list=""/><Attribute name = "ensembl_gene_id" /><Attribute name = "entrezgene" /><Attribute name = "phenotype_description"\
                /><Attribute name = "source_name" /><Attribute name = "study_external_id" /></Dataset></Query>\'');

            try:
                tdf = util.read_csv(ensembl_file, sep="\t", header=None, names=['ensembl_gene_id', 'gid', 'phenotype_description', 'source_name','study_external_id'], nrows=1);
            except Exception as exp:
                tdf = pd.DataFrame()

            if len(tdf) != 0:
                valid_files.append(ensembl_file)

        if len(valid_files) == 0:
            return None

        cmd = 'cat %s >> %s'%(' '.join(valid_files), out_file)
        print  cmd
        util.unix (cmd)

        try:
            ensembl_data=util.read_csv(out_file, sep=r"\t", header=None, names=['ensembl_gene_id', 'gid', 'phenotype_description', 'source_name','study_external_id']);
        except Exception as exp:
            print "Phenotype is not available for %s"%tax_id
            return  None

        ensembl_data=ensembl_data[ensembl_data['gid'].notnull()];
        data=[]
        h_map = GPUtils.get_ensembl2gid_map()
        for k, row in ensembl_data.groupby(['source_name','gid']):
            cnt=[];
            for i in row.index:
                v1 = self.tointstr(row.at[i,'study_external_id']);
                if v1:
                    cnt.append('[' + v1 +'] '+row.at[i,'phenotype_description']);
                else:
                    cnt.append(row.at[i,'phenotype_description']);
            content = ';'.join(cnt)
            data.append({'gid':k[1], 'content':content,'annotation_field1':'', 'type_name':k[0].replace(' ','_'), 'tax_id':tax_id})
            
        return self.collapse_by_gid(pd.DataFrame(data))
        #Tracer()()        

	#deleted by Ali end
    def tointstr (self, i):
        if i is None: return '';
        try:
            return str(int(i))
        except ValueError:
            s = str(i)
            if s!='nan':
                return s;
            return '';
    
    def collapse_by_gid(self, df):
        if len(df) == 0:
            return pd.DataFrame()
        data=[]
        for k, grow in df.groupby(['gid','type_name']):
            if k[1]=='TRANSMEMBRANE_ENSEMBL':
                data.append({'gid':k[0], 'content': 'Yes','annotation_field1': ';'.join(pd.unique(grow['annotation_field1'].values)), 'type_name':k[1], 'tax_id':grow['tax_id'].values[0]})
            else:        
                data.append({'gid':k[0], 'content': ';'.join(pd.unique(grow['content'].values)),'annotation_field1': ';'.join(pd.unique(grow['annotation_field1'].values)), 'type_name':k[1], 'tax_id':grow['tax_id'].values[0]})
        return pd.DataFrame(data);
    #deleted by Ali begin
    def get_annotation_helper(self, attributes, type, is_boolean, tax_id):
        db_name = self.get_dbname_by_taxid(tax_id)
        if db_name is None:
            return []

        util.unix('mkdir -p ' + SyncDB.DOWNLOAD_DIR() + '/ensembl_files');
        out_file = SyncDB.DOWNLOAD_DIR() + "/ensembl_files/ensembl_annotations_%s_%s.csv"%(type,tax_id);
        
        attr_str = '<Attribute name = "ensembl_gene_id" /><Attribute name = "entrezgene" />'
        for attr in attributes:
            attr_str += '<Attribute name = "' + attr + '" />';
        print "Processing annotation for %s specie %s"%(type,tax_id);

        valid_files= []
        for chr in self.chrList:
            ensembl_file = SyncDB.DOWNLOAD_DIR() + "/ensembl_files/ensembl_annotations_%s_%s_chr%s.csv"%(type, tax_id, chr);
            print "downloading ensembl_phenotype_%s_chr%s..."%(tax_id, chr);
            if not os.path.exists(ensembl_file):
                util.unix('wget -O ' + ensembl_file + ' - \'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0"\
                encoding="UTF-8"?><!DOCTYPE Query><Query virtualSchemaName = "default"\
                formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "%s" interface = "default"\
                >' % db_name + '<Filter name="chromosome_name" value="' + chr + '" filter_list=""/>' + attr_str + '</Dataset></Query>\'');

            try:
                tdf = util.read_csv(ensembl_file, sep="\t", header=None, names=['ensembl_gene_id', 'gid']+ attributes, nrows=1);
            except Exception as exp:
                tdf = pd.DataFrame()

            if len(tdf) != 0:
                valid_files.append(ensembl_file)

        if len(valid_files) == 0:
            return []

        cmd = 'cat %s >> %s'%(' '.join(valid_files), out_file)
        print  cmd
        util.unix (cmd)

        try:
            ensembl_data=util.read_csv(out_file, sep="\t", header=None, names=['ensembl_gene_id', 'gid'] + attributes).drop_duplicates();
        except Exception as exp:
            print "No data available for %s specie %s" % (type, tax_id);
            return []

        ensembl_data=ensembl_data[ensembl_data['gid'].notnull()];
        data=[]
        #Tracer()()
        for k, grow in ensembl_data.groupby(['gid']):
                #Tracer()()
                cnt=[]
                for i in grow.index:
                    
                    v1 = self.tointstr(grow.at[i,attributes[0]]);
                    #print i, v1;
                    if len(attributes) > 1:
                        v2 = self.tointstr(grow.at[i,attributes[1]])
                        if v1 != '':
                            v1 = '[' + v1 + '] ' + v2;
                        else:
                            v1 = v2;
                    if v1 != '':         
                        cnt.append(v1)

                cnt = pd.unique(cnt);        
                content = "";
                if is_boolean:
                    if len(cnt)>0:
                        content="Yes"                        
                else:
                    content = ';'.join(cnt)
                    
                if content != '':
                    data.append({'gid':k, 'content':content,'annotation_field1':grow.at[i,'ensembl_gene_id'], 'type_name':type, 'tax_id':tax_id})
        
        return data;
    #deleted by Ali end

    def get_annotation_mart_query(self, type, db):
        if type == "GOA_ENSEMBL":
            return "select goa.description_1074 description, goa.dbprimary_acc_1074 term, ent.dbprimary_acc_1074 as gid, main.stable_id_1023 as gene from %s__translation__main main, %s__ox_goslim_goa__dm goa, %s__ox_EntrezGene__dm ent where goa.translation_id_1068_key = main.translation_id_1068_key and ent.gene_id_1020_key = main.gene_id_1020_key and goa.dbprimary_acc_1074 is not NULL and ent.dbprimary_acc_1074 is not Null"%(db,db,db)
        elif type == "INTERPRO_ENSEMBL":
            return "select intp.description_1074 description, intp.interpro_ac_1026 term, ent.dbprimary_acc_1074 as gid, main.stable_id_1023 as gene from %s__translation__main main, %s__interpro__dm intp, %s__ox_EntrezGene__dm ent where intp.translation_id_1068_key = main.translation_id_1068_key and ent.gene_id_1020_key = main.gene_id_1020_key and intp.interpro_ac_1026 is not NULL and ent.dbprimary_acc_1074 is not Null"%(db,db,db)
        elif type == "APPRIS_ENSEMBL":
            return "select NULL description, appr.value_1065 term, ent.dbprimary_acc_1074 as gid, main.stable_id_1023 as gene from %s__transcript__main main, %s__tra_appris__dm appr, %s__ox_EntrezGene__dm ent where appr.transcript_id_1064_key = main.transcript_id_1064_key and ent.gene_id_1020_key = main.gene_id_1020_key and appr.value_1065 is not NULL and ent.dbprimary_acc_1074 is not Null" % (
            db, db, db)
        elif type in ["Orphanet", "GOA", "DDG2P"]:
            return "select phen.description_20125 description, phen.external_id_20125 term, ent.dbprimary_acc_1074 as gid, phen.stable_id_1023 as gene  from %s__phenotype__dm phen, %s__ox_EntrezGene__dm ent where phen.gene_id_1020_key = ent.gene_id_1020_key and phen.external_id_20125 is not NULL and phen.source_20125='%s' and ent.dbprimary_acc_1074 is not Null" % (
                db, db, type)
        elif type == "TRANSMEMBRANE_ENSEMBL":
            return "select Null description, hmm.hit_name_1048 term, ent.dbprimary_acc_1074 as gid, main.stable_id_1023 as gene from %s__translation__main main, %s__protein_feature_tmhmm__dm hmm, %s__ox_EntrezGene__dm ent where hmm.translation_id_1068_key = main.translation_id_1068_key and ent.gene_id_1020_key = main.gene_id_1020_key and hmm.hit_name_1048 is not NULL and ent.dbprimary_acc_1074 is not Null"%(db,db,db)

        return None;

    def get_annotation_martdb(self, tax_id, a_type, is_boolean=False):
        import math
        db_name=self.get_dbname_by_taxid(tax_id)
        if db_name is None:
            return None;
        file = SyncDB.DOWNLOAD_DIR() + "/ensembl_files/biomart_%s_%s.csv"%(a_type, tax_id);
        print "Running query to get %s for %s from martdb"%(a_type, tax_id)


        query = self.get_annotation_mart_query(a_type, db_name)

        con = self.get_biomart_connection()
        try:
            if os.path.exists(file):
                df = util.read_csv(file);
            else:
                df = db.from_sql(con, query).drop_duplicates()
                df.to_csv(file, index=False);
        except:
            print "error in getting %s data for %s"%(a_type, tax_id)
            return None;

        data = []
        #Tracer()()
        for k, grow in df.groupby(['gid']):
            # Tracer()()
            cnt = []
            for i in grow.index:
                v1 = grow.at[i, "term"];
                v2 = grow.at[i, "description"]
                try:
                    if type(v1) is str or not math.isnan(v1):
                        try:
                            if type(v2) is str or not math.isnan(v2):
                                cnt.append('[%s] %s'%(str(v1),str(v2)))
                            else:
                                cnt.append(str(v1))
                        except:
                            cnt.append(str(v1))
                except:
                    pass

            cnt = pd.unique(cnt);
            content = ''
            if is_boolean:
                if len(cnt) > 0:
                    content = "Yes"
            else:
                content = ';'.join(cnt)

            if content != '':
                data.append({'gid': k, 'content': content, 'annotation_field1': grow.at[i, 'gene'],
                             'type_name': a_type, 'tax_id': tax_id})

        return data;
    #deleted by Ali begin
    def get_annotation_data_old(self):
        data = [];
        for tax_id in self.taxidList:
            data += self.get_annotation_helper(['transmembrane_domain'], "TRANSMEMBRANE_ENSEMBL", True, tax_id);
            data += self.get_annotation_helper(['goslim_goa_accession','goslim_goa_description'], "GOA_ENSEMBL", False, tax_id);
            data += self.get_annotation_helper(['interpro','interpro_description'], "INTERPRO_ENSEMBL", False, tax_id);
            data += self.get_annotation_helper(['family','family_description'], "FAMILY_ENSEMBL", False, tax_id);
            data += self.get_annotation_helper(['transcript_appris'], "APPRIS_ENSEMBL", False, tax_id);

        #data += self.get_annotation_helper(['ncoils'], "NCOILS_ENSEMBL");
        #data += self.get_annotation_helper(['signal_domain'], "SIGNAL_DOMAIN_ENSEMBL");
        #data += self.get_annotation_helper(['low_complexity'], "LOW_COMPLEXITY_ENSEMBL");

        data += self.get_variations();
        self.collapse_by_gid(pd.DataFrame(data)).to_csv(self.fn_annotations, index=False)
    #deleted by Ali end
    def get_annotation_data(self):
        data = [];
        for tax_id in self.taxidList:
            for type in ["APPRIS_ENSEMBL", "GOA_ENSEMBL","INTERPRO_ENSEMBL", "Orphanet", "GOA", "DDG2P"]:
                d = self.get_annotation_martdb(tax_id, type);
                if d is not None:
                    data += d
            d = self.get_annotation_martdb(tax_id, "TRANSMEMBRANE_ENSEMBL", is_boolean=True);
            if d is not None:
                data += d;

        data += self.get_variations();
        self.collapse_by_gid(pd.DataFrame(data)).to_csv(self.fn_annotations, index=False)

    def get_ensembl_connection(self, db):
        import MySQLdb as mysql
        return mysql.connect('ensembldb.ensembl.org', 'anonymous', '', db)

    def get_biomart_connection(self):
        import MySQLdb as mysql
        return mysql.connect(host='martdb.ensembl.org', port=5316, user='anonymous', passwd='', db='ensembl_mart_84')

    def get_variations(self):
        ensembl_file = SyncDB.DOWNLOAD_DIR() + "/ensembl_files/ensembl_variations.csv";
        print "Processing variations"
        
        if os.path.exists(ensembl_file):
            t=util.read_csv(ensembl_file);
        else:    
            con = self.get_ensembl_connection(EnsemblDownload.get_ensembl_latest_version('homo_sapiens_variation_{0}_'.format(EnsemblDownload.ENSEMBL_VERSION)))
            query = "select distinct pf.object_id as variation_name,p.description,v.clinical_significance,vg.gene_name, s.name as source_name from source s, phenotype_feature pf, phenotype p, variation v, variation_genename vg where pf.type ='Variation' and pf.phenotype_id = p.phenotype_id and v.name=pf.object_id and v.variation_id=vg.variation_id and v.source_id=s.source_id and v.clinical_significance in ('likely pathogenic','pathogenic','risk factor','association','drug response')"
            t = db.from_sql(con,query, params=[])
            t.to_csv(ensembl_file, index=False);
               
        map = GPUtils.get_sym2gid_map()["sym2gid"];
        data=[]
        for gene, row in t.groupby(['gene_name']):
            if gene in map:
                #Tracer()()
                content = ['['+r[1]['variation_name']+'] '+r[1]['description']+'{'+r[1]['clinical_significance']+'}('+r[1]['source_name']+')' for r in row.iterrows()];
                data.append({'gid':map[gene],'content':';'.join(content),'annotation_field1':gene, 'type_name':'VARIATIONS_ENSEMBL', 'tax_id':'9606'});

        return data;
                
                    
    def do_update(self):
        self.get_idmap_data()
        self.get_annotation_data()
        
    def check_inputs(self):
        return True;

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

class UCSCDownload(XmlClass):

    def __init__(self, xe=None):
        XmlClass.__init__(self,xe=xe)
        self.tag = "UCSCDownload"
        self.fn_id_map=os.path.join(SyncDB.DOWNLOAD_DIR(),"ucsc_id_map.csv")
        self.fn_annotations=os.path.join(SyncDB.DOWNLOAD_DIR(),"ucsc_annotations.csv")
        self.taxidList = SyncDB.SUPPORTED_SPECIES;
        for child in self.children:
            if type(child).__name__ == "Species":
                self.taxidList = child.supported_species
                break;

    def populate_more(self,root):
        #self.outputs = [self.fn_id_map,self.fn_annotations]
        self.outputs = [self.fn_annotations]

    def get_dsname_by_taxid(self, tax_id):
        map={
            "9606":
                "hg38",
            "10090":
                "mm10",
            "10116":
                "rn6",
            "7955":
                "danRer10",
            "9031":
                "galGal4",
            "9544":
                "rheMac8",
            "7227":
                "dm6",
            "6239":
                "ce11",
            "4932":
                "sacCer3"
        }

        if tax_id in map:
            return map[tax_id]
        print 'warning: specie %s is not support in our UCSC build.'%tax_id
        return None

    def get_ucsc2gid_df (self, tax_id):
        #mmusculus_gene_ensembl   10090
        #rnorvegicus_gene_ensembl 10116
        db_name = self.get_dsname_by_taxid(tax_id)
        if db_name is None:
                return pd.DataFrame()

        util.unix('mkdir -p ' + SyncDB.DOWNLOAD_DIR() + '/ucsc_files');       
        file = SyncDB.DOWNLOAD_DIR() + "/ucsc_files/ucscid2gid_%s"%tax_id;
        
        if not os.path.exists(file):
            con=db.get_con('UCSC')
            try:
                df=db.from_sql(con,"select name as source_id, value as gid from %s.knownToLocusLink"%db_name)
            except Exception as exp:
                return pd.DataFrame()
            df.to_csv(file, index=False);

        data=util.read_csv(file);
        data=data[data['gid'].notnull()];
        data=data[data['source_id'].notnull()];
        data[['gid']]=data[['gid']].astype(int)
        data['tax_id']=tax_id;
        return data

            
    def get_idmap_data(self):
        df=pd.DataFrame();
        for tax_id in self.taxidList:
            df=df.append(self.get_ucsc2gid_df(tax_id))
                
        df.drop_duplicates().to_csv(self.fn_id_map, index=False)

    def get_chromosome_data(self):
        data = [];
        for tax_id in self.taxidList:
            data += self.get_chromosome_data_by_taxid(tax_id)
        return data

    def get_chromosome_data_by_taxid(self, tax_id):
        db_name = self.get_dsname_by_taxid(tax_id)
        if db_name is None:
                return []
        print 'Retrieving chromosome info for %s'%tax_id
        file = path.join(SyncDB.DOWNLOAD_DIR(), "ucsc_genome_location_info_%s.csv"%tax_id);
        print '****************************'
        import sys
        from pprint import pprint
        pprint(sys.path)
        if not os.path.exists(file):
            con=db.get_con('UCSC')
            try:
                if db_name == 'sacCer3':
                    query  = """SELECT rl.mrnaAcc AS name, rl.locusLinkId AS gid, rg.chrom, rg.txStart AS start, rg.txEnd AS end
                                        FROM sacCer3.sgdGene rg, hgFixed.refLink rl
                                        WHERE rl.name=rg.name
                                        """
                else:
                    query  = """select rg.name, rl.locusLinkId as gid, rg.chrom, rg.txStart as start, rg.txEnd as end
                                            from %s.refGene rg, %s.refLink rl
                                            where rl.mrnaAcc=rg.name
                                        """ % (db_name, 'hgFixed')
                df=db.from_sql(con,query)
            except Exception as exp:
                print 'Error in getting gene locations for %s'%tax_id
                import traceback
                print traceback.format_exc()
                return []

            df.to_csv(file, index=False);

        loc_data=util.read_csv(file);
        data=[]
        for gid, row in loc_data.groupby('gid'):
            content = pd.unique([r[1]['chrom']+':'+str(r[1]['start']+1)+'-'+str(r[1]['end']) for r in row.iterrows()]);            
            data.append({'gid':gid, 'content':';'.join(content),'annotation_field1':';'.join(row['name'].values), 'type_name':'genome_location', 'tax_id':tax_id})
                
        return data

    def get_annotation_data(self):
        data = [];       
        data += self.get_chromosome_data();                   
        pd.DataFrame(data).drop_duplicates().to_csv(self.fn_annotations, index=False)                         
                    
    def do_update(self):
        # we will get UCSC ID->Gene ID from Ensembl
        #self.get_idmap_data()
        self.get_annotation_data()
        
    def check_inputs(self):
        return True;
#!/usr/bin/env python
#from .core import *
import numpy as np
import shutil
import urllib
import urlparse
from os.path import splitext, basename
import os
from core import *
import util
import xml.etree.ElementTree as ET
from pprint import pprint
#from .. import SyncDB
class MsigdbDownload(XmlClass):
    MSIGDB_DICT = {
        'H:': 'Hallmark Gene Sets',        
        'C1:': 'Positional Gene Sets',
        'C2:CGP':'Chemical And Genetic Perturbations',
        'C2:CP':'Canonical Pathways',
        'C2:CP:BIOCARTA':'Biocarta Gene Sets',
        'C2:CP:REACTOME':'Reactome Gene Sets',
        'C3:MIR':'Microrna Targets',
        'C3:TFT':'Transcription Factor Targets',
        'C4:CGN':'Cancer Gene Neighborhoods',
        'C4:CM':'Cancer Modules',
        'C6:':'Oncogenic Signatures',
        'C7:':'Immunologic Signatures'}
        
    ORGNAME2TAXID = {
        'Danio rerio':'7955',
        'Homo sapiens':'9606',
        'Macaca mulatta':'9544',
        'Mus musculus':'10090',
        'Rattus norvegicus':'10116',
    }
    
    def __init__(self, xe=None):
        XmlClass.__init__(self,xe=xe)
        self.fn_term =os.path.join(SyncDB.DOWNLOAD_DIR(),"msigdb_term.csv") 
        self.fn_term_gene_pair =os.path.join(SyncDB.DOWNLOAD_DIR(),"msigdb_term_gene_pair.csv") 
        self.inputs=['ds:msigdb']

    def populate_more(self,root):
        self.outputs = [self.fn_term,self.fn_term_gene_pair]
    def make_gid2taxid(self):
        taxidList = SyncDB.SUPPORTED_SPECIES;
        for child in self.children:
            if type(child).__name__ == "Species":
                taxidList = child.supported_species
                break;
        if not os.path.exists(os.path.join(SyncDB.DOWNLOAD_DIR(), "geneid2taxid.csv")):
            taxid_filter = "";
            if len(taxidList) != 0:
                taxid_filter = ['$1==\"' + t + '\"' for t in taxidList]
                taxid_filter = "if (" + "||".join(taxid_filter) + ")"
            # gene_id,tax_id
            cmd = "time zcat " + SyncDB.DOWNLOAD_DIR() + "/gene_info.gz | cut -f1,2 | sed 1d | awk 'BEGIN{FS=\"\\t\"; OFS=\"\\t\"}{" + taxid_filter + " print $2,$1;}' | sort -k1,1 -t $'\\t' >" + SyncDB.DOWNLOAD_DIR() + "/geneid2taxid.csv";
            print cmd;
            util.unix(cmd);

        df = util.read_csv(SyncDB.DOWNLOAD_DIR() + "/geneid2taxid.csv", names=["gid", "tax_id"], sep=r'\t')
        self.gid2taxid = {str(df.ix[i, 'gid']):str(df.ix[i, 'tax_id']) for i in df.index};

    def do_update(self):
        self.make_gid2taxid()
        tree=ET.parse(os.path.join(SyncDB.DOWNLOAD_DIR(),"msigdb.xml"))
        root=tree.getroot()
        cols = ['STANDARD_NAME','SYSTEMATIC_NAME','ORGANISM','EXTERNAL_DETAILS_URL','CHIP',
                'CATEGORY_CODE','SUB_CATEGORY_CODE','CONTRIBUTOR','CONTRIBUTOR_ORG',
                'DESCRIPTION_BRIEF','MEMBERS_EZID','DESCRIPTION_FULL']
        #cols = ['STANDARD_NAME','CATEGORY_CODE']
        convert = [{attrib:xe.attrib[attrib].encode('utf-8') for attrib in cols} for xe in root]
        df = pd.DataFrame(columns=cols)
        df = df.append(convert,ignore_index=True)
        #df['tax_id'] = [MsigdbDownload.ORGNAME2TAXID[x] for x in df['ORGANISM']]
        util.rename2(df,{'STANDARD_NAME':'term_name','DESCRIPTION_BRIEF':'description','SYSTEMATIC_NAME':'term_id'})
        df = df[(df['CATEGORY_CODE']!='C5') & (df['SUB_CATEGORY_CODE']!='CP:KEGG')]
        col_type = [MsigdbDownload.MSIGDB_DICT.get(r['CATEGORY_CODE'] + ":" + r['SUB_CATEGORY_CODE'],'Error') for index, r in df.iterrows()]
        df['type_name'] = col_type

        def term_name_map(l):
            #return ' '.join([('' if len(x) == 0 else x[0].upper())+('' if len(x) <= 1 else x[1:].lower()) for x in l.split('_')])
            return l.replace('_', ' ');
        df['term_name']=df['term_name'].map(term_name_map);
        
        df_term = df.loc[:,('term_id','description','type_name','term_name')]
        df_term['tax_id']='0'
        df_term.to_csv(self.fn_term,index=False)
        rows=[]
        for index,r in df.iterrows():
            for gid in r['MEMBERS_EZID'].split(','):
                tax_id = self.gid2taxid[gid] if gid in self.gid2taxid else None;
                if tax_id is not None:
                    rows.append({'term_id':r['term_id'],
                                 'gid': gid,
                                 'type_name': r['type_name'],
                                 'term_name':r['term_name'],
                                 'tax_id':tax_id})
        df_term_gene_pair = pd.DataFrame(rows)
        df_term_gene_pair.to_csv(self.fn_term_gene_pair,index=False)

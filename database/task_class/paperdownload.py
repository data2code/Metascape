#!/usr/bin/env python
import numpy as np
import shutil
import urllib
import urlparse
import os
from core import *
import util
from pprint import pprint
import pandas as pd

class PaperDownload(XmlClass):

    def __init__(self, xe=None):
        XmlClass.__init__(self,xe=xe)
        self.dest = xe.attrib['dest']
        self.s_file_obo = os.path.join(SyncDB.DOWNLOAD_DIR(),'hp.obo')
        self.s_file_gene2hpo = os.path.join(SyncDB.DOWNLOAD_DIR(),'genes_to_phenotype.txt')
        self.fn_hpo_ann = os.path.join(SyncDB.DOWNLOAD_DIR(),'hpo_ann.csv')
        self.fn_trrust_rawdata_human = os.path.join(SyncDB.DOWNLOAD_DIR(),'trrust_rawdata.human.tsv')
        self.fn_trrust_rawdata_mouse = os.path.join(SyncDB.DOWNLOAD_DIR(),'trrust_rawdata.mouse.tsv')
        self.fn_DisGeNET_source = os.path.join(SyncDB.DOWNLOAD_DIR(), 'curated_gene_disease_associations.tsv')
        self.fn_DisGeNET_ann = os.path.join(SyncDB.DOWNLOAD_DIR(), 'disgenet_ann.csv')
        self.fn_trrust_human_term = os.path.join(SyncDB.DOWNLOAD_DIR(),'trrust_human.csv')
        self.fn_trrust_mouse_term = os.path.join(SyncDB.DOWNLOAD_DIR(),'trrust_mouse.csv')
        self.fn_symbol = os.path.join(SyncDB.UPLOAD_DIR(),'gid2source_id','symbol.csv')
        self.fn_synonym = os.path.join(SyncDB.UPLOAD_DIR(),'gid2source_id','gene_synonym.csv')
        self.fn_description = os.path.join(SyncDB.UPLOAD_DIR(),'annotation','gene_description.csv')
        self.inputs=['ds:paper',self.fn_symbol,self.fn_synonym, self.fn_description]
        self.fn_trrust_term =  os.path.join(SyncDB.DOWNLOAD_DIR(), 'trrust_term.csv')
        self.fn_trrust_term_pair =  os.path.join(SyncDB.DOWNLOAD_DIR(), 'trrust_term_pair.csv')
        self.fn_disgenet_term =  os.path.join(SyncDB.DOWNLOAD_DIR(), 'disgenet_term.csv')
        self.fn_disgenet_term_pair =  os.path.join(SyncDB.DOWNLOAD_DIR(), 'disgenet_term_pair.csv')
        self.fn_PaGenBase_term =  os.path.join(SyncDB.DOWNLOAD_DIR(), 'PaGenBase_term.csv')
        self.fn_PaGenBase_term_pair =  os.path.join(SyncDB.DOWNLOAD_DIR(), 'PaGenBase_term_pair.csv')



    def get_fn_dest(self):
        return os.path.join(SyncDB.DOWNLOAD_DIR(),self.dest)

    def populate_more(self,root):
        XmlClass.populate_more(self,root)
        self.outputs.extend([self.fn_hpo_ann,
                             self.fn_trrust_term,
                             self.fn_trrust_term_pair,
                             self.fn_disgenet_term,
                             self.fn_disgenet_term_pair,
                             self.fn_DisGeNET_ann,
                             self.fn_PaGenBase_term,
                             self.fn_PaGenBase_term_pair,
                             ])

    def do_update(self):
        self.get_parse_PaGenBase()
        self.parse_disgenet()
        t_term_human, t_term_pair_human = self.parse_trrust(self.fn_trrust_rawdata_human, 9606, start_id=1)
        t_term_mouse, t_term_pair_mouse = self.parse_trrust(self.fn_trrust_rawdata_mouse, 10090,start_id=len(t_term_human)+1)
        t_term = pd.concat([t_term_human,t_term_mouse])
        t_term_pair = pd.concat([t_term_pair_human,t_term_pair_mouse])
        t_term.to_csv(self.fn_trrust_term, index=False)
        t_term_pair.to_csv(self.fn_trrust_term_pair, index=False)
        parent_child = self.parse_hp(self.s_file_obo)
        # print(parent_child['HP:0000001'])
        pheno_level = self.get_level(parent_child)
        # print(pheno_level['HP:0012823'], pheno_level['HP:0000001'])
        self.parse_gp(self.s_file_gene2hpo, pheno_level)

    def get_parse_PaGenBase(self):
        count_start = 0
        S_term = []
        S_pair = []

        S_file = [
            ('hotisp.txt', 'Tissue-specific', 9606),
            ('hocesp.txt', 'Cell-specific', 9606),
            ('mutisp.txt', 'Tissue-specific', 10090),
            ('mucesp.txt', 'Cell-specific', 10090),
            ('ratisp.txt', 'Tissue-specific', 10116),
            ('drtisp.txt', 'Tissue-specific', 7227)
        ]

        for fn in S_file:
            fn = fn[0]
            urllib.urlretrieve('http://bioinf.xmu.edu.cn/PaGenBase/browse/{0}'.format(fn),
                               os.path.join(SyncDB.DOWNLOAD_DIR(), fn))

        for (s_file, s_ann, tax_id) in S_file:
            s_file = os.path.join(SyncDB.DOWNLOAD_DIR(), s_file)
            t_term, t_pair, count_start = self.parse_PaGenBase(s_file, s_ann, tax_id, count_start)
            S_term.append(t_term)
            S_pair.append(t_pair)

        t_term = pd.concat(S_term, ignore_index=True)
        t_pair = pd.concat(S_pair, ignore_index=True)
        t_term.to_csv(self.fn_PaGenBase_term, index=False)
        t_pair.to_csv(self.fn_PaGenBase_term_pair, index=False)

        pass

    def parse_PaGenBase(self, s_file, s_ann, tax_id, count_start=0):
        t = pd.read_table(s_file, skiprows=7)
        t.rename2({'Gene Symbol': 'Symbol'})
        t = t[['Symbol', 'Sample']].copy()
        S_term = util.unique(t.Sample)
        data = []
        c_id = {}



        for x in S_term:
            count_start += 1
            term_id = 'PGB:%05d' % count_start
            term_name = s_ann + ': ' + x
            data.append({'term_id': term_id, 'term_name': term_name, 'description': term_name})
            c_id[x] = term_id
        t_term = pd.DataFrame(data)
        t_pair = t[['Symbol', 'Sample']].copy()
        t_pair.rename2({'Sample': 'term_name'})
        t_pair['term_id'] = t_pair.term_name.apply(lambda x: c_id[x])
        t_pair['term_name'] = t_pair.term_name.apply(lambda x: s_ann + ': ' + x)
        t_pair['tax_id'] = tax_id
        t_pair['type_name'] = 'PaGenBase'
        t_pair.drop_duplicates(['term_id', 'Symbol'], inplace=True)

        #convert symbol to gid
        dt = pd.read_csv(self.fn_symbol)
        dt = dt[dt['tax_id']==tax_id]
        symbol2gene_id = dict(zip(dt.source_id, dt.gid.astype(str)))
        dt = pd.read_csv(self.fn_synonym)
        dt = dt[dt['tax_id']==tax_id]
        symbol2gene_id.update(dict(zip(dt.source_id, dt.gid.astype(str))))
        t_pair['gid'] = t['Symbol'].apply(lambda x: symbol2gene_id.get(x, ''))
        t_pair = t_pair[t_pair.gid != ''].copy()

        return (t_term, t_pair, count_start)


    def parse_trrust(self, s_file, tax_id, start_id):
        dt = pd.read_csv(self.fn_symbol)
        dt = dt[dt['tax_id']==tax_id]
        symbol2gene_id = dict(zip(dt.source_id, dt.gid.astype(str)))
        dt = pd.read_csv(self.fn_synonym)
        dt = dt[dt['tax_id']==tax_id]
        symbol2gene_id.update(dict(zip(dt.source_id, dt.gid.astype(str))))
        dt = pd.read_csv(self.fn_description)
        dt = dt[dt['tax_id']==tax_id]
        gene_id2description = dict(zip(dt.gid.astype(str), dt.content))

        t = pd.read_table(s_file, header=None, names=['TF', 'Target', 'Drection', 'PMID'])
        t['gid_TF'] = t['TF'].apply(lambda x: symbol2gene_id.get(x, ''))
        t['gid_Target'] = t['Target'].apply(lambda x: symbol2gene_id.get(x, ''))
        t = t[(t.gid_TF != '') & (t.gid_Target != '')].copy()
        t['term_name'] = t['TF'].apply(lambda x: 'Regulated by: ' + x)
        util.rename2(t,{'gid_Target':'gid'})
        all_tf = {}
        term_ids = []
        for i, row in t.iterrows():
            tf = row['TF']
            if tf not in all_tf:
                all_tf[tf] = 'TRR{0:05d}'.format(start_id+len(all_tf))
            term_ids.append(all_tf[tf])
        t['term_id'] = term_ids
        t['tax_id'] = str(tax_id)
        t['type_name'] = 'TRRUST'
        t_term_pair = t[['gid','tax_id','term_id','term_name','type_name']]

        t = t[['gid_TF','term_id','term_name']]
        t = t.drop_duplicates()
        t['description'] = [x['term_name']
                             + '; '
                             + gene_id2description.get(x['gid_TF'],'')[:60]
                             for i, x in t[['gid_TF','term_name']].iterrows()]
        t_term = t[['term_id','term_name','description']]
        return (t_term,t_term_pair)

    def parse_disgenet(self):
        t = pd.read_table(self.fn_DisGeNET_source)
        #term_id,term_name,description
        #geneId  geneSymbol      diseaseId       diseaseName            score                   NofPmids        NofSnps source
        #10      NAT2            C0005695        Bladder Neoplasm        0.245871429880008       5               0       CTD_human
        util.rename2(t,{'diseaseId':'term_id',
                        'diseaseName':'description',
                        'geneId':'gid'})
        t['term_name'] = t['description']
        t_term = t[['term_id','description','term_name']]
        t_term = t_term.drop_duplicates()
        t_term.to_csv(self.fn_disgenet_term,index=False)

        t_term_pair = t[['gid','term_id','term_name']]
        t_term_pair['tax_id'] = '9606'
        t_term_pair['type_name'] = 'DisGeNET'
        t_term_pair.to_csv(self.fn_disgenet_term_pair,index=False)

        #generate annotation file for DisGenNet
        t = pd.read_csv(self.fn_DisGeNET_source, delimiter="\t")
        c_cnt = util.unique_count(t.diseaseId)
        t['diseaseCnt'] = t['diseaseId'].apply(lambda x: c_cnt[x])
        t.sort_values(['geneId', 'diseaseCnt', 'diseaseName'], inplace=True)
        data = []
        for k, t_v in t.groupby('geneId'):
            s_ann = "; ".join(t_v.diseaseName)
            data.append({'gid': k, 'tax_id': 9606, 'content': s_ann})
        t = pd.DataFrame(data)
        t.to_csv(self.fn_DisGeNET_ann, index=False)

    def parse_hp(self, s_file):
        #parses the human phenotype file (assumes that there's a .txt file called "hp.txt")
        pc_dict = dict() #parent-child dictionary (key: parent, value: children)
        S=util.read_list(s_file)
        i=0
        n=len(S)
        while i<n:
            while i<n and S[i].rstrip()!='[Term]':
                i+=1
            if i>=n: break
            i+=1
            while S[i].rstrip()!='':
                contents = S[i].rstrip().split(': ')
                if contents[0] == 'id':
                    id = contents[1]
                elif contents[0] == 'is_a':
                    parent = contents[1].split(' ! ')[0]
                    if parent not in pc_dict:
                        pc_dict[parent] = []
                    pc_dict[parent].append(id)
                i+=1
        return pc_dict
    
    def get_level(self, pc_dict): #, pl_dict, phenotype, level_dict):
        #finds level
        #pc_dict: parent-child dictionary, pl_dict: phenotype-level dictionary (key: phenotype, value: level)
        #phenotype: phenotype that we want to know the level of, level_dict: dictionary with level as key and set of phenotypes as values
        root='HP:0000001'
        pl_dict = {root:0}
        Q=[root] #processing queue
        Q2=[] # child queue
        level=0
        while len(Q):
            for x in Q:
                for y in pc_dict.get(x, []):
                    if y in pl_dict: continue
                    pl_dict[y]=level+1
                    Q2.append(y)
            Q=Q2
            Q2=[]
            level+=1
        return pl_dict
    
    def parse_gp(self, s_file_gene2hpo, pheno_level):
        #parses genes_to_phenotype file (assumes that there's a .txt file called "genes_to_phenotype.txt")
        t=pd.read_table(s_file_gene2hpo, skiprows=1, names=['GeneID','Symbol','TermName','HP_ID'])
        c_cnt=util.unique_count(t.HP_ID)
        c_id2name=dict(zip(t.HP_ID, t.TermName))
        data=[]
        for k,t_v in t.groupby('GeneID'):
            S=t_v.HP_ID.tolist()
            S = sorted(S, key=lambda x: (c_cnt[x], - pheno_level[x]))
            data.append({'Gene':k,
                         'Phenotype': "; ".join([c_id2name[x] for x in S]),
                         'tax_id':9606})
        t_ann=pd.DataFrame(data)
        t_ann.to_csv(self.fn_hpo_ann, index=False)
        return t_ann


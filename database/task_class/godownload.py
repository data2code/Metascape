#!/usr/bin/env python
#from .core import *
import numpy as np
import shutil
import urllib
import urlparse
from os.path import splitext, basename
import os
from os import sys, path
from core import *
import util
import pandas as pd
import csv

from IPython.core.debugger import Tracer

#from .. import SyncDB
from gputil import GPUtils


class GoDownload(XmlClass):
    def __init__(self, xe=None):
        XmlClass.__init__(self,xe=xe)
        self.fn_dest_go_term =os.path.join(SyncDB.DOWNLOAD_DIR(),"go_term.csv")
        self.fn_source_by_tax =os.path.join(SyncDB.DOWNLOAD_DIR(),"gene2go_9606_extended.csv")
        self.fn_gene2go_org =os.path.join(SyncDB.DOWNLOAD_DIR(),"gene2go_9606_org.csv")                
        self.fn_dest_term2term =os.path.join(SyncDB.DOWNLOAD_DIR(),"go_term2term.csv") 
        self.fn_dest_go_annotations =os.path.join(SyncDB.DOWNLOAD_DIR(),"go_annotation.csv")         
        self.inputs=['ds:go']
        

    def populate_more(self,root):
        self.outputs = [self.fn_dest_go_term,self.fn_source_by_tax,self.fn_dest_term2term,self.fn_dest_go_annotations]

    def get_go_term_and_go_term2term(self):
        fn_source =os.path.join(SyncDB.DOWNLOAD_DIR(),"go_daily-termdb-tables.tar.gz") 

        dir = os.path.join(SyncDB.DOWNLOAD_DIR(),"go_daily-termdb-tables")
        urllib.urlretrieve("http://archive.geneontology.org/latest-termdb/go_daily-termdb-tables.tar.gz",fn_source)
        cmd = "tar -zxvf %s -C %s go_daily-termdb-tables/term.txt go_daily-termdb-tables/graph_path.txt go_daily-termdb-tables/term_definition.txt" % (fn_source,SyncDB.DOWNLOAD_DIR())
        util.unix(cmd)
        fn_term = os.path.join(dir,"term.txt")
        df_term = util.read_csv(fn_term, sep=r'\t', skiprows = 1, 
                names=['id','term_name','term_type','term_id','is_obsolete','is_root','is_relation'])
        df_term = df_term[(df_term.is_root == 0) & (df_term.is_relation == 0) ]
        df_term['term_type'] = df_term['term_type'].map({'biological_process':'BP',
            'molecular_function':'MF',
            'cellular_component':'CC'})
        fn_term_description = os.path.join(dir,"term_definition.txt")
        #remove bad lines.
        util.unix("awk 'BEGIN{FS=\"\\t\"} {if(NF==5) print $0}' " + fn_term_description + " > " + fn_term_description + ".tmp")
        util.unix("mv " + fn_term_description + ".tmp " + fn_term_description)
        df_term_description = util.read_csv(fn_term_description, sep=r'\t', skiprows = 1, 
                names=['id','description','dbxref_id','term_comment','reference'],
                usecols=['id','description'])
        df_term_final = pd.merge(df_term, df_term_description, on='id', how='left')
        df_term_final.drop(['id','is_obsolete','is_root','is_relation'],axis=1,inplace=True)
        df_term_final = df_term_final[df_term_final['term_id'].str.contains(r'^GO:')] 
        df_term_final = df_term_final.reindex(columns=['term_id','term_name','term_type','description'])

        fn_gp = os.path.join(dir,"graph_path.txt")
        df_gp = util.read_csv(fn_gp, sep=r'\t', skiprows = 1, 
                names=['tt_id','term1_id','term2_id','relationship_type_id','distance','relation_distance'])
        df_gp = pd.merge(df_gp,df_term,left_on='term1_id', right_on='id',how='left')
        util.rename2(df_gp, {'term_id':'parent_term_id'})
        df_gp = df_gp.loc[:,['term2_id','distance','parent_term_id']]

        df_gp = pd.merge(df_gp,df_term,left_on='term2_id', right_on='id',how='left')
        df_gp = df_gp.loc[:,['parent_term_id','distance','term_id']]

        df_gp = df_gp[df_gp['parent_term_id'].str.contains(r'^GO:') & df_gp['term_id'].str.contains(r'^GO:')] 
        df_gp.drop_duplicates(inplace=True)

        df_gp.to_csv(self.fn_dest_term2term,index=False)

        df_term_final.to_csv(self.fn_dest_go_term, index = False)
    
    def get_parents(self, term):
        return self.parent_lookup[term] if term in self.parent_lookup and len(self.parent_lookup[term])!=0 else None;
        
	#returns all the ancestors of t. 
    def extend_go_term(self, t):
        extended = [t];
        p = self.get_parents(t)
        if p:
            extended += p;
        return np.unique(extended).tolist();
       
    def build_parent_lookup(self):
        self.parent_lookup ={}
        term2term = util.read_csv(self.fn_dest_term2term, sep=r',', skiprows = 1, 
                names=['parent_term_id','distance','term_id'])

        for k,g in term2term.groupby('term_id', as_index=False):
            self.parent_lookup[k] = np.unique(g['parent_term_id']).tolist();
            if k in self.parent_lookup[k]:
                self.parent_lookup[k].remove(k);                
        
	#Creates gene_id to go_term map. Map each gene_id to all the ancestors of the immediate term_id.
	
    def get_gene2go(self):
        terms = util.read_csv(self.fn_dest_go_term, sep=r',', skiprows = 1).values        
        #Tracer()()
        term_type_index = {terms[i][0]:terms[i][2] for i in range(len(terms))}
        term_description_index = {terms[i][0]:terms[i][1] for i in range(len(terms))}
                
        fn_source =os.path.join(SyncDB.DOWNLOAD_DIR(),"gene2go.gz") 
        urllib.urlretrieve("ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz",fn_source)

        taxidList = SyncDB.SUPPORTED_SPECIES;
        for child in self.children:
            if type(child).__name__ == "Species":
                taxidList = child.supported_species
                break;
        taxid_filter = "";
        if len(taxidList) !=0:
            taxid_filter = ['$1=="'+ t + '"' for t in taxidList]
            taxid_filter = "if ("+"||".join(taxid_filter) + ") "
        
        #cmd = "zegrep '^9606' %s |cut -f2,3,1,8,6 > %s" % ( fn_source, self.fn_gene2go_org)
        cmd = ('zcat %s | awk \'BEGIN {FS="\t"}{' + taxid_filter + ' print $0;}\' |cut -f2,3,1,8,6 > %s') % ( fn_source, self.fn_gene2go_org)
        print 'aaaaaaaaaaaaa'
        print cmd
        print 'aaaaaaaaaaaaa'
        util.unix(cmd)
        gene2go = util.read_csv(self.fn_gene2go_org, sep=r'\t', 
                names=['tax_id','gene_id','term_id','type', 'description'])
        extended_terms = [];
        
        for k,g in gene2go.groupby('gene_id', as_index=False):
            extended_term_ids = [];
			#map gene_id to all the ancestors of t. self.extend_go_term(t) returns list of ancestors of t which includes t itself.
            for t in g['term_id']:
                extended_term_ids += self.extend_go_term(t);
                    
            extended_term_ids = np.unique(extended_term_ids);

            taxId = g['tax_id'].values[0];
            extended_terms +=[[taxId,k, t, term_type_index[t] if t in term_type_index else '', term_description_index[t] if t in term_description_index else ''] for t in extended_term_ids];

        extended_terms_df = pd.DataFrame(extended_terms);
        extended_terms_df.columns = gene2go.columns;
        extended_terms_df.to_csv(self.fn_source_by_tax, index=False, sep='\t', header=False)
    
    def get_term_score(self, term):
        if term in self.term_score_map:
            return self.term_score_map[term];
        parents = self.get_parents(term);
        if parents is None:
            #self.term_score_map[term] = -np.log((self.org_goterm_count[term] if term in self.org_goterm_count else 1)/ float(self.extended_goterm_count[term]));
            if term not in self.extended_goterm_count:
                self.term_score_map[term] = 0;
            else:    
                self.term_score_map[term] = -np.log(float(self.extended_goterm_count[term])/50000);            
            return self.term_score_map[term];

        scores = [];
        for p in parents:
            parent_score = self.get_term_score(p);
            if term not in self.extended_goterm_count:
                scores.append(parent_score)
            elif p not in self.extended_goterm_count:
                scores.append(1.1*parent_score - np.log(self.extended_goterm_count[term]/float(50000)));
            else:
                scores.append(1.1*parent_score - np.log(self.extended_goterm_count[term]/float(self.extended_goterm_count[p])));
        
        self.term_score_map[term] = max(scores);
        return self.term_score_map[term];
        
    def build_term_score_map(self):
        self.org_goterm_count = self.build_go_term_count(self.fn_gene2go_org)
        self.extended_goterm_count = self.build_go_term_count(self.fn_source_by_tax)
        self.term_score_map = {};
        #self.get_term_score('GO:0003756')
        path_len = [];
        for t in self.org_goterm_count:
            self.term_score_map[t] = self.get_term_score(t);
            path_len.append(self.path_to_root_lengh(t));
            
        #Tracer()();
        
    def build_go_term_count(self, file):
        goterm_count_map = {}
        gene2go = util.read_csv(file, sep=r'\t', 
                names=['tax_id','gene_id','term_id','type','description'])
          
        for k,g in gene2go.groupby('term_id', as_index=False):
            goterm_count_map[k] = len(util.unique(g['gene_id'].values))
        
        return goterm_count_map;
    
    def path_to_root_lengh (self, term):
        count = 0;
        p = self.get_parents(term);
        while p is not None:            
            count +=1;
            p = self.get_parents(p[0]);
            
        return count;
    
    def print_path_to_root (self, term, level, visited):
        visited.update({term:True});
        print '--'*level, term, "total:", (self.extended_goterm_count[term] if term in self.extended_goterm_count else -1.0), ' self:', (self.org_goterm_count[term] if term in self.org_goterm_count else -1.0), ' score:', self.term_score_map[term]                    
        parents = self.get_parents(term);
        if parents is not None:
            for p in parents:
                if p not in visited:
                    self.print_path_to_root(p, level+1, visited);


    def get_go_annotations(self):
        self.build_term_score_map();       
        gene2go = util.read_csv(self.fn_source_by_tax, sep=r'\t', 
                names=['tax_id','gene_id','term_id','type','description'])
        output = [];
        for k,g in gene2go.groupby('gene_id', as_index=False):
            terms = {};
            for go_term in g.values:
                if go_term[2] is None or go_term[3] is None:
                    continue;
                if go_term[3] not in terms:
                    terms[go_term[3]] = [];                    
                terms[go_term[3]].append({'score':self.term_score_map[go_term[2]], 'count': self.extended_goterm_count[go_term[2]], 'content':go_term[2] + ' ' +  go_term[4] if go_term[4] is not None else '-', 'term':go_term[2]});
            
            for t in terms:
                sorted_terms = sorted(terms[t], key=lambda k:-k['score']);

                output.append([';'.join([sorted_terms[i]['content'] for i in range(min(3,len(sorted_terms)))]), k, str(t), str(g.values[0][0]), 'go' ])

        with open(self.fn_dest_go_annotations, "w") as myfile:
            wr = csv.writer(myfile)
            wr.writerow(['content','gid','type', 'tax_id', 'ds']);
            wr.writerows(output);
        
    def do_update(self):
        self.get_go_term_and_go_term2term()
        self.build_parent_lookup()
        self.get_gene2go()    
        self.get_go_annotations()

    def check_inputs(self):
        passed = True
        urls=[]
        print "Checking urls for godownload..."
            
        urls.append("http://archive.geneontology.org/latest-termdb/go_daily-termdb-tables.tar.gz")        
        urls.append("ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz")
        for url in urls:
            if not GPUtils.check_url(url):
                passed = False;                
                
        return passed        
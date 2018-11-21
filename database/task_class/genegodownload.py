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
import itertools as it
import csv
import urllib2
import json
import sys
import requests
import eutils

from IPython.core.debugger import Tracer
#tew

class GeneGoDownload(XmlClass):
    def __init__(self, xe=None):
        XmlClass.__init__(self,xe=xe)
        self.fn_dest_genego_annotations = os.path.join(SyncDB.DOWNLOAD_DIR(),"genego_annotations.csv")             
        self.fn_dest_genego_terms = os.path.join(SyncDB.DOWNLOAD_DIR(),"genego_terms.csv")     
        self.fn_dest_genego_term_gene_pair = os.path.join(SyncDB.DOWNLOAD_DIR(),"genego_term_gene_pair.csv")             
        self.inputs=['ds:genego']

        self.taxidList = SyncDB.SUPPORTED_SPECIES;
        for child in self.children:
            if type(child).__name__ == "Species":
                self.taxidList = child.supported_species
                break;
        self.dir = SyncDB.DOWNLOAD_DIR() + '/genego_files';
        util.unix('mkdir -p ' + self.dir);


    def populate_more(self,root):
        self.outputs = [self.fn_dest_genego_terms, self.fn_dest_genego_term_gene_pair,self.fn_dest_genego_annotations]

    def prepare (self):
        self.dbcon=db.get_con('GENEGO')
        
    def fetch(self, s_sql):
        t=db.from_sql(self.dbcon, s_sql)
        return t
    
    def get_gene_disease_association(self):
        print 'Getting GeneGo disease association data';
        #df = self.fetch("select distinct a.ref as gene_id, a.disid, a.disname, a.note from (select ga.note_id, g17.ref, d.disid, d.disname, d.note from disease_associations_all_v ga, GeneDBS_17 g17, diseases d, geneorgs go1 where d.disid=ga.dis_id and ga.GENE_ID=g17.gene and go1.gene=ga.gene_id and go1.org=1) a");
        #Tracer()()
        fn = self.dir + "/gene_disease_association.csv"
        if not os.path.exists(fn):
            df = self.fetch("select distinct d.disid TERM_ID, d.disname as TERM_NAME, gdb.ref as GID, d.note as description, orgs.taxonomyid as tax_id from gene_netw gn, genediss gd, diseases d, genedbs gdb, geneorgs o, orgs where gn.gene = gd.gene and gd.dis = d.disid and gn.gene=gdb.gene and gn.gene=o.gene and o.org=orgs.orgid and orgs.taxonomyid in (" +','.join(self.taxidList)+ ") and gdb.db=17 and d.rtyp > 0 and d.disname not like 'By %' and gd.dis in (select distinct dismbr from disrelflat where disgrp = -1173899567 and dismbr <> -1173899567)");
            df['TERM_ID'] = 'gDIS' + df['TERM_ID'].map(str);
            df.rename2({"TERM_ID":"term_id", "TERM_NAME":"term_name", "GID":"gid", "DESCRIPTION":"description", "TAX_ID":"tax_id"})
            df['type_name'] = 'GeneGo Disease Association';
            #remove disease which has more than 500 genes
            df = df.drop(df.index[list(it.chain.from_iterable([ g for k, g in df.groupby('term_id').groups.items() if len(g)>=500]))])
            df.to_csv(fn, index=False)
        else:
            df = util.read_csv(fn)
        df1=pd.DataFrame(df.copy())[['gid','term_id','term_name','type_name', 'tax_id']]
        df2=pd.DataFrame(df.copy())[['term_id','term_name','type_name', 'description']]
        df2=df2.drop_duplicates();
        self.disease_gid2term = df1;        
        self.disease_terms = df2;
        self.disease_done = True;
        
        data=[]
        for k,t_v in df.groupby('gid'):
            S=util.unique([x for x in t_v['term_name'] if not pd.isnull(x)])
            data.append({'gid':k, 'content':"; ".join(S), 'type_name':t_v['type_name'].values[0], 'annotation_field1':len(S), 'tax_id':str(int(t_v['tax_id'].values[0]))})
        
        self.disease_annotations=pd.DataFrame(data)

        print 'GeneGo disease association data captured';
        
    def get_gene_pathway_map(self):
        print 'Getting GeneGo pathway data';
        #df = self.fetch("select distinct g17.ref as gid, m.MAPNAME as term_name, m.imid as term_id from genemaps gm, imagemap_table m, genedbs_17 g17, geneorgs go, genes g where gm.im=m.imid and gm.gene=g17.gene and go.gene = gm.gene and go.org=1 and g.geneid=go.gene");
        fn = self.dir + "/pathway.csv"
        if not os.path.exists(fn):
            df = self.fetch("select distinct i.imid as term_id, i.imagename url, i.mapname term_name, d.ref as gid, orgs.taxonomyid as tax_id from pw_imagemap_shapes s, pw_imagemap_class c, imagemap_table i, gene_netw n, genedbs d, geneorgs o, orgs where s.id=c.shape_id and i.imid=s.im and n.id=c.object_id and d.gene=n.gene and o.gene=n.gene and i.publish=1 and o.org=orgs.orgid and orgs.taxonomyid in (%s) and d.db=17"%(','.join(self.taxidList)));
            #Tracer()()
            df['TERM_ID'] = 'gMAP' + df['TERM_ID'].map(str);
            df.rename2({"TERM_ID":"term_id", "TERM_NAME":"term_name", "GID":"gid", "URL":"term_field1", "TAX_ID":"tax_id"})
            df['type_name'] = 'GeneGo Pathway';
            df.to_csv(fn, index=False)
        else:
            df = util.read_csv(fn)


        self.pathway_gid2term = df;

        df2=pd.DataFrame(df.copy())[['term_id','term_name','type_name']]
        df2['description']=df2['term_name']        
        df2=df2.drop_duplicates();
        self.pathway_terms=df2        
        self.pathway_done = True;              
        print 'GeneGo pathway data captured';

    def get_go_processes(self):
        print 'Getting GeneGo GO Processes';
        #df = self.fetch("select distinct g17.ref as gid, m.MAPNAME as term_name, m.imid as term_id from genemaps gm, imagemap_table m, genedbs_17 g17, geneorgs go, genes g where gm.im=m.imid and gm.gene=g17.gene and go.gene = gm.gene and go.org=1 and g.geneid=go.gene");
        fn = self.dir + "/processes.csv"
        if not os.path.exists(fn):
            df = self.fetch("select distinct b.id as term_id, b.name as term_name, d.ref as gid, orgs.taxonomyid as tax_id from posobjs p, gene_netw g, listblocks l, blocks b, genedbs d, geneorgs o, orgs where p.obj=g.id and l.block= p.blck  and b.id = l.block and g.gene=d.gene and o.gene=g.gene and l.lst=-900 and  o.org=orgs.orgid and orgs.taxonomyid in (%s) and d.db=17"%(','.join(self.taxidList)));
            df['TERM_ID'] = 'gNWT' + df['TERM_ID'].map(str);
            #Tracer()()
            df.rename2({"TERM_ID":"term_id", "TERM_NAME":"term_name", "GID":"gid", "TAX_ID":"tax_id"})
            df['type_name'] = 'GeneGo GO Processes';
            df.to_csv(fn, index=False)
        else:
            df = util.read_csv(fn)

        self.go_processes_gid2term = df;

        term2gidcount = {};
        for k,t_v in df.groupby('term_id'):
            term2gidcount[k] = len(t_v['gid'])
            
        df2=pd.DataFrame(df.copy())[['term_id','term_name','type_name', 'tax_id']]
        df2['description']=df2['term_name']        
        df2=df2.drop_duplicates();
        self.go_processes_terms=df2

        data=[]
        for k,t_v in self.go_processes_gid2term.groupby('gid'):
            #Tracer()()                                
            S = [{'term_name':t_v['term_name'].values[i], 'count':term2gidcount[t_v['term_id'].values[i]]} for i in range(len(t_v['term_name']))]
            len_s = len(S);
            d = {}            
            #to do unique
            for obj in S:
                d[obj['term_name']] = obj
            S = d.values()            
            S = sorted(S, key=lambda k:k['count'])[0:3];
            T = [x['term_name'] for x in S];
            data.append({'gid':k, 'content':"; ".join(T), 'type_name':t_v['type_name'].values[0], 'annotation_field1':len_s, 'tax_id':str(int(t_v['tax_id'].values[0]))})
        self.go_processes_annotations=pd.DataFrame(data)
        
        self.go_processes_done = True;              
        print 'GeneGo GO Processes data captured';
        
    def get_gene_drug_target(self):
        print 'Getting GeneGo drug target data';
        fn = self.dir + "/drug_targets.csv"
        if not os.path.exists(fn):
            df = self.fetch("select g17.ref as gene_id, dt.drug_id as term_id, dt.drug_name as term_name, dt.disease_name, dt.status, orgs.taxonomyid as tax_id from orgs, geneorgs go, genedbs_17 g17, gene_netw gt, (  \
                              select distinct u.id as drug_id, u.name as drug_name, druglinks.target as target_id, ro.name as target_name, d.disid as disease_id, d.disname as disease_name, drug.status\
                              from units u, diseases d, (select unit, dis, ds.dscr as status\
                              from noticediss d, notices_table n, noticequant q, quantity_comp_mv qc, comps c, notice_drugstatus nd, knd_drugstatus ds\
                              where d.note=n.id and n.trust>=0 and n.type_id=4 and n.id=q.note and nd.note = n.id and ds.id = nd.drugtype and q.quant=qc.quantity and qc.op in(43,80) and qc.comp=c.id ) drug \
                              left join (select distinct c.unit as drug, r.id2 as target, dl.link_id, dl.disease_id as dis\
                            from dis_links dl, regulation_rels r, regulation_comp rc, comps c\
                              where dl.link_id = r.link_id and r.id1 = rc.id and rc.object_id = c.id) druglinks on druglinks.dis = drug.dis and druglinks.drug = drug.unit\
                              left join regulation_objects ro on ro.id = druglinks.target\
                              where drug.unit=u.id and drug.dis = d.disid and target is not NULL) dt where gt.id=dt.target_id and g17.gene = gt.gene and go.gene=g17.gene and go.org=orgs.orgid and orgs.taxonomyid in (" + ','.join(self.taxidList) + ")");

            df['TERM_ID'] = 'gDRG' + df['TERM_ID'].map(str);
            df.to_csv(fn, index=False)
        else:
            df = util.read_csv(fn)

        #Tracer()()                
        #self.gene_drug_target_gid2terms = pd.DataFrame(df.copy())[['DRUG_ID','DRUG_NAME', 'GENE_ID']]
        self.gene_drug_target_gid2terms = df[['TERM_ID','TERM_NAME', 'GENE_ID', 'TAX_ID']]
        self.gene_drug_target_gid2terms.rename2({"TERM_ID":"term_id", "TERM_NAME":"term_name", "GENE_ID":"gid", "TAX_ID":"tax_id"})
        self.gene_drug_target_gid2terms['type_name'] = 'GeneGo Drug Target';
        self.gene_drug_target_gid2terms['description'] = "Disease:" + df['DISEASE_NAME'] +  '; Status:' + df['STATUS'];               

        
        data=[]
        for k,t_v in df.groupby('TERM_ID'):
            S=util.unique([x for x in t_v['DISEASE_NAME'] if not pd.isnull(x)])
            data.append({'term_id':k, 'term_name':t_v['TERM_NAME'].values[0], 'type_name':'GeneGo Drug Target', 'description':";".join(S)})        
        self.gene_drug_target_terms=pd.DataFrame(data)

        

        data=[]
        for k,t_v in self.gene_drug_target_gid2terms.groupby('gid'):
            S=util.unique([x for x in t_v['term_name'] if not pd.isnull(x)])
            data.append({'gid':k, 'content':"; ".join(S), 'type_name':t_v['type_name'].values[0], 'annotation_field1':len(S), 'tax_id':str(int(t_v['tax_id'].values[0]))})
        self.gene_drug_target_annotations=pd.DataFrame(data)

        self.drug_target_done = True;
                
        print 'GeneGo drug target captured';

    def get_annotations(self):
    
        print 'Getting GeneGo annotations';
        fn = self.dir + "/annotations.csv"
        if not os.path.exists(fn):
            df = self.fetch("select distinct d.ref as gid, ic.id as class_id, ic.name as functional_class, lc.locid location_id, lc.locname as location, orgs.taxonomyid as tax_id \
                            from orgs, gene_netw g, object_icons ri, regulation_comp rc, comps c, genedbs d,geneorgs o, icons ic, locs lc \
                            where g.id = ri.id (+) and ri.type='N' and g.id=rc.id (+) and rc.object_id = c.id (+) and ri.icon_id=ic.id and c.loc=lc.locid and g.gene=d.gene and o.gene=g.gene and o.org=orgs.orgid and orgs.taxonomyid in (%s) and d.db=17 order by d.ref"%(','.join(self.taxidList)));
            df.to_csv(fn, index=False)
        else:
            df = util.read_csv(fn)

        fc = df[['GID','FUNCTIONAL_CLASS', 'CLASS_ID', 'TAX_ID']].drop_duplicates();
        data=[]
        brief_class_data = []
        for k,t_v in fc.groupby('GID'):
            S=util.unique([x for x in t_v['FUNCTIONAL_CLASS'] if not pd.isnull(x)])
            data.append({'gid':k, 'content':";".join(S), 'type_name':'GeneGo Functional Class', 'annotation_field1':len(S),'tax_id':str(int(t_v['TAX_ID'].values[0]))})
            bclass = [];            
            if len(set([7,64]).intersection(t_v['CLASS_ID'])) != 0:
                bclass.append('Transcription factors')
            if len(set([10,11,15,64]).intersection(t_v['CLASS_ID'])) != 0:
                bclass.append('Receptors')
            if len(set([63]).intersection(t_v['CLASS_ID'])) != 0:
                bclass.append('Ligands')
            if len(set([15,36,61,62]).intersection(t_v['CLASS_ID'])) != 0:
                bclass.append('Kinases')
            if len(set([5,6]).intersection(t_v['CLASS_ID'])) != 0:
                bclass.append('Proteases')
            if len(set([32, 59, 60]).intersection(t_v['CLASS_ID'])) != 0:
                bclass.append('Phosphatases')
            if len(set([2,4]).intersection(t_v['CLASS_ID'])) != 0:
                bclass.append('Enzymes')
            if len(bclass)!=0:
                brief_class_data.append({'gid':k, 'content':";".join(bclass), 'type_name':'GeneGo Brief Class', 'annotation_field1':len(bclass), 'tax_id':str(int(t_v['TAX_ID'].values[0]))})

        self.functional_class_annotation=pd.DataFrame(data)
        self.brief_class_annotation=pd.DataFrame(brief_class_data)
        
        loc = df[['GID','LOCATION', 'LOCATION_ID', 'TAX_ID']].drop_duplicates();
        data=[]
        for k,t_v in loc.groupby('GID'):
            S=util.unique([x for x in t_v['LOCATION'] if not pd.isnull(x)])
            data.append({'gid':k, 'content':";".join(S), 'type_name':'GeneGo Location', 'annotation_field1':len(S),'tax_id':str(int(t_v['TAX_ID'].values[0]))})
        self.location_annotation=pd.DataFrame(data)
        
        self.annotations_done = True;
                
        print 'GeneGo annotations captured';

    def merge_data(self):
        open(self.fn_dest_genego_annotations, 'w').close()    
        open(self.fn_dest_genego_terms, 'w').close()
        open(self.fn_dest_genego_term_gene_pair, 'w').close()
        header_annotation = True;
        header_terms = True
        dfs_terms=[]
        dfs_gid2terms=[]
        dfs_annotations = []
        if hasattr(self, 'pathway_done') and self.pathway_done:
            dfs_terms.append(self.pathway_terms)
            dfs_gid2terms.append(self.pathway_gid2term)

        if hasattr(self, 'disease_done') and self.disease_done:
            dfs_terms.append(self.disease_terms)
            dfs_gid2terms.append(self.disease_gid2term)
            dfs_annotations.append(self.disease_annotations)

        if hasattr(self, 'drug_target_done') and self.drug_target_done:
            dfs_terms.append(self.gene_drug_target_terms)
            dfs_gid2terms.append(self.gene_drug_target_gid2terms)
            dfs_annotations.append(self.gene_drug_target_annotations)

        if hasattr(self, 'go_processes_done') and self.go_processes_done:
            dfs_terms.append(self.go_processes_terms)
            dfs_gid2terms.append(self.go_processes_gid2term)
            dfs_annotations.append(self.go_processes_annotations)

        if hasattr(self, 'annotations_done') and self.annotations_done:
            dfs_annotations.append(self.functional_class_annotation)
            dfs_annotations.append(self.brief_class_annotation)
            dfs_annotations.append(self.location_annotation)

        pd.concat(dfs_terms).reindex(columns=['term_id','term_name','type_name','description']).to_csv(self.fn_dest_genego_terms, index=False, header= True)
        gid2terms = pd.concat(dfs_gid2terms).reindex(columns=['gid','term_id','term_name','type_name', 'tax_id'])
        gid2terms['tax_id'] = gid2terms['tax_id'].map(int)
        gid2terms['tax_id'] = gid2terms['tax_id'].map(str)
        gid2terms.to_csv(self.fn_dest_genego_term_gene_pair, index=False, header= True)
        pd.concat(dfs_annotations).reindex(columns=['gid','content','type_name', 'annotation_field1', 'tax_id']).to_csv(self.fn_dest_genego_annotations, index=False, header= True)
            
    def do_update(self):
        self.prepare()
        self.get_gene_disease_association()
        self.get_gene_pathway_map()
        self.get_go_processes()
        self.get_gene_drug_target()
        self.get_annotations()
        self.merge_data()
        
    def check_inputs (self):
        self.prepare()
        print "Check database connection for GeneGO..."
        try:
            df = self.fetch("select * from geneorgs o where ROWNUM < 10");
        except:
            print "Error in connecting to the GeneGo database"
            return False            
        return True
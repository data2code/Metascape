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
from pprint import pprint
import StringIO
import db
import itertools
from core import *
import util
from IPython.core.debugger import Tracer
import re

#print __package__
#from ..core import SyncDB, XmlClass
class PostProcess():
    tmp_directory="postprocess_tmp_files"
    @staticmethod       
    def load_table(table_name, data):
        fn = PostProcess.tmp_directory + "/" + table_name + ".csv";
        data.to_csv(fn,index=False, header=False)
        conn_info = db.get_con_info(SyncDB.CONNECTION_ID)        
        mysql_cmd = [
            'mysqlimport',
            '--local',
            '--delete',
            '--host ' + conn_info["HOST"],
            '-u ' + conn_info["USR"],
            '--password=' + conn_info["PWD"],
            '--fields-terminated-by=,',
            '--fields-optionally-enclosed-by="\\""',
            '--ignore-lines=0', 
            SyncDB.DATABASE,
            fn]
        c = " ".join(mysql_cmd)
        print c
        s = util.unix(mysql_cmd)
        
    @staticmethod        
    def gene_id_cleanup():

        if not path.exists(path.join(SyncDB.DOWNLOAD_DIR(), "gene_info.csv")):
            print 'PostProcess: required file gene_info.csv does not exist'
            sys.exit()
    
        con = db.get_con(SyncDB.CONNECTION_ID, db=SyncDB.DATABASE)
        #delete gid which doesn't have symbol and history
        to_deleted_gids = db.from_sql(con, '''
                                    SELECT t.gid FROM gid2source_id t
                                    WHERE NOT EXISTS(
                                    SELECT 1 FROM gid2source_id a
                                    WHERE a.gid = t.gid
                                    AND a.id_type_id IN (1,5))
                                    ''')['gid'].tolist()
        to_deleted_gids = ['26148','100507739','101930400','728127','728013','554223','80761','142937','390760','100289290','105369230','112268369','112268338','112268368','112268369','112268368','105369230']
        db.sql_in(s_sql_left = 'DELETE FROM gid2source_id  WHERE gid IN (',
                  s_sql_right= ')',
                  S_id=to_deleted_gids,
                  con = con)
        #print 'cleaning up gene_history';
        query="Select id_type_id from {0}.id_type where id_type_name='Gene_History'";
        query = query.format(SyncDB.DATABASE)
        t = db.from_sql(con,query)
        gene_history_id_type_id = t.ix[0,0].astype(str)
        #gene_history_gids = pd.unique(db.from_sql(con,"Select distinct gid from {0}.gid2source_id where id_type_id={1}".format(SyncDB.DATABASE, gene_history_id_type_id))['gid'])
        #invalid_genes =set(gene_history_gids).difference(PostProcess.valid_genes)
        #if (len(invalid_genes) != 0):
        #    db.sql_in("Delete from {0}.gid2source_id where id_type_id={1} and gid in (".format(SyncDB.DATABASE, gene_history_id_type_id), ")", invalid_genes, con=con);
            
        gene_history = db.from_sql(con,"Select * from {0}.gid2source_id where id_type_id={1}".format(SyncDB.DATABASE, gene_history_id_type_id))
        gene_history_map={}
        for row in gene_history.iterrows():
            gene_history_map[row[1]['source_id']] = row[1]['gid'];

        print 'cleaning up interaction';        
        interations = db.from_sql(con,"Select * from {0}.interaction".format(SyncDB.DATABASE))
        new_rows = [];
        for row in interations.iterrows():
            gidA = str(row[1]['gid_A'])
            if gidA in gene_history_map:
                row[1]['gid_A'] = gene_history_map[gidA]            
                gidB = str(row[1]['gid_B'])
                if gidB in gene_history_map:
                    row[1]['gid_B'] = gene_history_map[gidB]                
                    new_rows.append(row[1])
            
        clean_interactions = pd.DataFrame(new_rows).query("gid_A != gid_B").drop_duplicates()
        PostProcess.load_table("interaction", clean_interactions)

        print 'cleaning up annotaitons';        
        annotaitons = db.from_sql(con,"Select * from {0}.annotation".format(SyncDB.DATABASE))
        new_rows = [];
        for row in annotaitons.iterrows():
            gid = str(row[1]['gid'])
            if gid in gene_history_map:
                row[1]['gid'] = gene_history_map[gid]
                new_rows.append(row[1])
            
        clean_annotations = pd.DataFrame(new_rows).drop_duplicates();
        PostProcess.load_table("annotation", clean_annotations);
        
        print 'cleaning up gid2terms';        
        gid2terms = db.from_sql(con,"Select * from {0}.gid2terms".format(SyncDB.DATABASE))
        new_rows = [];
        for row in gid2terms.iterrows():
            gid = str(row[1]['gid'])
            if gid in gene_history_map:
                row[1]['gid'] = gene_history_map[gid]
                new_rows.append(row[1])
            
        clean_gid2terms = pd.DataFrame(new_rows).drop_duplicates();
        PostProcess.load_table("gid2terms", clean_gid2terms);
              
        print 'cleaning up term2gids';        
        term2gids = db.from_sql(con,"Select * from {0}.term2gids".format(SyncDB.DATABASE))
        new_rows = [];
        for row in term2gids.iterrows():
            new_gids = pd.unique(map(str,[gene_history_map[x] for x in row[1]['gids'].split(',') if x in gene_history_map]));
            row[1]['gids'] = ','.join(new_gids)
            row[1]['id_count'] = len(new_gids);
            if len(new_gids)!=0:
                new_rows.append(row[1])
            
        clean_term2gids = pd.DataFrame(new_rows).drop_duplicates();
        PostProcess.load_table("term2gids", clean_term2gids);
        
        print 'cleaning up homologene table';        
        homologene = db.from_sql(con,"Select * from {0}.homologene".format(SyncDB.DATABASE))
        new_rows = [];
        for row in homologene.iterrows():
            gid = str(row[1]['gid'])
            tax_id = row[1]['tax_id']
            #if tax_id not in [9606,10090,10116]:
            #    new_rows.append(row[1])
            if gid in gene_history_map:
                row[1]['gid'] = gene_history_map[gid]
                new_rows.append(row[1])
            
        clean_homologene = pd.DataFrame(new_rows).drop_duplicates();
        PostProcess.load_table("homologene", clean_homologene);

    @staticmethod
    def add_pf_alias():
        x=db.DB(name=SyncDB.CONNECTION_ID, db=SyncDB.DATABASE)
        #delete old
        x.from_sql("DELETE FROM gid2source_id  WHERE id_type_id = 6 AND id_status='post_pf_alias'")

        #get fn
        url = 'http://plasmodb.org/common/downloads/Current_Release/Pfalciparum3D7/txt/'

        import re
        import urllib2

        pattern = '<a href=".*?Pfalciparum3D7_GeneAliases.txt">(.*?)</a>'
        response = urllib2.urlopen(url).read()
        for filename in re.findall(pattern, response):
            fn_GeneAliases = filename
            break
        print filename

        url = url + fn_GeneAliases
        fn_source =os.path.join(SyncDB.DOWNLOAD_DIR(),fn_GeneAliases)
        urllib.urlretrieve(url,fn_source)
        # one row per gene
        S=util.read_list(fn_source)
        # ticket sent to PlasmoDB, reply says they will add old-style names to the alias file (Build30)
        # i.e., ^PF[A-Z]_\d\d\d\d[cw] is considered an old id
        data=[]
        for i,s in enumerate(S):
            S_id=s.split("\t")
            for s_id in S_id:
                data.append({'index':i, 'name':s_id})
        t_plasmodb=pd.DataFrame(data)

        # obtain name - gene id mapping from NCBI
        t_ncbi=x.from_sql("select gid,source_id name from gid2source_id where id_type_id in (1,6,15) and tax_id=5833")
        t_ncbi.drop_duplicates(inplace=True)
        corrected=[]
        c_seen={}
        for i in t_ncbi.index:
            s_id=t_ncbi.ix[i, 'name']
            c_seen[s_id]=True
            # some genes have bad names in NCBI, e.g., PFA0820w is showed as PFA_0820w
            # ticket sent to NCBI Ticket #28045-111721
            # title: Gene name for P. falciparum
            # waiting for reply
            if re.search(r'^PF[A-Z]_\d\d\d\d[cw]', s_id):
                s_id2=re.sub('_', '', s_id)
                if s_id2 in c_seen: continue
                # add the entry for the correct name
                corrected.append({ 'gid': t_ncbi.ix[i,'gid'], 'name':s_id2})
                c_seen[s_id2]=True
        if len(corrected):
            t=pd.DataFrame(corrected)
            t_ncbi=pd.concat([t_ncbi, t], ignore_index=True)

        t=t_plasmodb.merge(t_ncbi, left_on='name', right_on='name', how='left')
        data=[]
        for k,t_v in t.groupby('index'):
            t_v=t_v.copy()
            S=util.unique(t_v.gid.values, is_rm_na=True)
            if len(S)==0:
                # no name in this group can be mapped to an Entrez Gene, so the whole group is dropped
                #print "==============="
                #print t_v
                continue
            elif len(S)>1:
                if len(S)==len(t_v):
                    # if each symbol has their own Gene ID, let them be, even PlasmoDB maps them into one gene
                    # as we would not know which gene ID to use (we could add both, but let it be for now)
                    data.append(t_v)
                    continue
                t0=t_v[~pd.isnull(t_v.gid)] # names mapped to Entrez Gene, let them be the reference entries
                data.append(t0)
                t1=t_v[pd.isnull(t_v.gid)]  # names cannot be mapped
                # we map the names to each of the found Entrez Gene ID
                # e.g., PF3D7_0209400 is linked to both PFB0423c,8444970  and PFB0425c,812686
                # therefore, we add PF3D7_0209400 as alias to both gene id: 8444970 and 812686
                for i in range(len(S)):
                    t1=t1.copy()
                    t1['gid']=S[i]
                    data.append(t1)
            elif len(S)==1:
                # A unique Entrez Gene ID is found, all names will use that ID
                t_v['gid']=S[0]
                data.append(t_v)
        t=pd.concat(data, ignore_index=True)
        t['gid']=t.gid.astype(int)
        t.drop('index', axis=1, inplace=True)
        # just in case
        t.drop_duplicates(inplace=True)
        # we can exclude those aliases that are already in Metascape
        t_exist=x.from_sql("select gid,source_id name from gid2source_id where tax_id=5833")
        c={ t_exist.ix[i,'name']:t_exist.ix[i,'gid'] for i in t_exist.index }
        I=[ t.ix[i, 'gid']!=c.get(t.ix[i, 'name'],0) for i in t.index ]
        #print len(t)
        t=t[I].copy()
        t['ds']='PlasmoDB'
        #print len(t)
        insert_sql = '''
                INSERT INTO gid2source_id (  gid,  id_type_id,  source_id,  id_status,  ds,  tax_id)
                VALUES(?,?,?,?,?,?)
        '''
        for i, r in t.iterrows(): #name,gid,ds
            x.from_sql(insert_sql,
                       [
                    r['gid'],
                    6, #id_type_id,
                    r['name'], #''source_id',
                    'post_pf_alias',#''id_status',
                    'PlasmoDB',#''ds',
                    '5833',#'tax_id'
                       ])

    @staticmethod
    def orthoMCL():
        #delete old one
        x=db.DB(name=SyncDB.CONNECTION_ID, db=SyncDB.DATABASE)
        #delete old
        x.from_sql("DELETE FROM homologene WHERE comments = 'orthoMCL'")
        fn_source_gz = os.path.join(SyncDB.DOWNLOAD_DIR(), 'aa_deflines_OrthoMCL-5.txt.gz')
        fn_source = os.path.join(SyncDB.DOWNLOAD_DIR(), 'aa_deflines_OrthoMCL-5.txt')
        fn_dest = os.path.join(SyncDB.DOWNLOAD_DIR(), 'OrthoMCL.csv')
        if True:
            url = 'http://orthomcl.org/common/downloads/release-5/aa_deflines_OrthoMCL-5.txt.gz'
            urllib.urlretrieve(url,fn_source_gz)
            util.ungzip_it(fn_source_gz)
        #get c_pf_name2id
        df = x.from_sql('''
                        select gid,source_id name from gid2source_id where tax_id=5833
                        ''')
        c_pf_name2id = { x['name']:int(x['gid']) for i, x in df.iterrows()}
        x=db.DB(name=SyncDB.CONNECTION_ID, db=SyncDB.DATABASE)
        t_ensm = x.from_sql("select gid,source_id name from gid2source_id where id_type_id in (10) and tax_id=9606")
        c_ensm = {t_ensm.ix[i, 'name']: t_ensm.ix[i, 'gid'] for i in t_ensm.index}

        S = util.read_list(fn_source)
        out = []
        c = {'>hsap': '9606', '>pfal': '5833'}
        for s in S:
            if re.search(r'^>(hsap|pfal)', s):
                X = s.split("|")
                out.append([X[2].strip(), X[1].strip(), c.get(X[0])])
        t = pd.DataFrame(out, columns=['grp_id', 'name', 'tax_id'])
        out = []
        for k, t_v in t.groupby('grp_id'):
            if k == 'no_group': continue
            if len(t_v) < 2: continue
            t_v = t_v.copy()
            t_v['gid'] = t_v['name'].apply(lambda x: c_ensm.get(x, 0) or c_pf_name2id.get(x, 0))
            t_v = t_v[t_v.gid > 0].copy()
            S_tax_id = set(t_v.tax_id)
            if '9606' in S_tax_id and '5833' in S_tax_id:
                out.append(t_v)
        t = pd.concat(out, ignore_index=True)
        t.drop_duplicates(inplace=True)
        t = t[['tax_id', 'grp_id', 'gid']]
        insert_sql = '''
                INSERT INTO homologene (  tax_id,  homologene_id,  gid,  comments)
                VALUES (?,?,?,'orthoMCL')
        '''
        for i, r in t.iterrows():  # name,gid,ds
            x.from_sql(insert_sql,
                       [
                           r['tax_id'],
                           r['grp_id'],  # ''homologene_id',
                           r['gid']
                       ])


    @staticmethod
    def add_yeast_alias():
        # download alias file from:
        # http://downloads.yeastgenome.org/curation/chromosomal_feature/SGD_features.tab
        # one row per gene
        # 1.   Primary SGDID (mandatory)
        # 2.   Feature type (mandatory)
        # 3.   Feature qualifier (optional)
        # 4.   Feature name (optional)
        # 5.   Standard gene name (optional)
        # 6.   Alias (optional, multiples separated by |)
        # 7.   Parent feature name (optional)
        # 8.   Secondary SGDID (optional, multiples separated by |)
        # 9.   Chromosome (optional)
        # 10.  Start_coordinate (optional)
        # 11.  Stop_coordinate (optional)
        # 12.  Strand (optional)
        # 13.  Genetic position (optional)
        # 14.  Coordinate version (optional)
        # 15.  Sequence version (optional)
        # 16.  Description (optional)

        x=db.DB(name=SyncDB.CONNECTION_ID, db=SyncDB.DATABASE)
        #delete old
        x.from_sql("DELETE FROM gid2source_id  WHERE id_type_id = 6 AND id_status='post_yeast_alias'")
        x.from_sql("DELETE FROM annotation  WHERE annotation_type_id = 11 AND annotation_field1='post_yeast_summary'")

        url = 'http://downloads.yeastgenome.org/curation/chromosomal_feature/SGD_features.tab'
        fn_source =os.path.join(SyncDB.DOWNLOAD_DIR(),'SGD_features.tab')
        urllib.urlretrieve(url,fn_source)

        t=pd.read_table(fn_source, names=['SGDID','TYPE','QUALIFIER','NAME','GENE_NAME','ALIAS','PARENT','SGDID2','CHROMOSOME','START','STOP','STRAND','POSITION','COORDINATE','SEQUENCE','DESCRIPTION'])
        #print t[2:]
        t=t[t.TYPE=='ORF'].copy()

        t_ncbi=x.from_sql("select gid,source_id NAME from gid2source_id where id_type_id in (1,6,15) and tax_id=4932")
        t_ncbi.drop_duplicates(inplace=True)

        t_sgd=t[['SGDID','NAME','GENE_NAME','DESCRIPTION']]
        #t_sgd.to_csv('t.csv')

        # missing entries are discontinued entries according to NCBI
        #t_sgd=t_sgd.merge(t_ncbi, left_on='NAME', right_on='NAME', how='left')
        #t_sgd=t_sgd[ pd.isnull(t_sgd.gid) ]
        #t_sgd.to_csv('t.csv')
        t_sgd=t_sgd.merge(t_ncbi, left_on='NAME', right_on='NAME')
        t_sgd['ds']='SGD'
        t_exist=x.from_sql("select gid,source_id name from gid2source_id where tax_id=4932")
        c={ t_exist.ix[i,'name']:t_exist.ix[i,'gid'] for i in t_exist.index }
        data=[]
        for i in t_sgd.index:
            s_name1=t_sgd.ix[i, 'NAME']
            s_name2=t_sgd.ix[i, 'GENE_NAME']
            sgdid=t_sgd.ix[i, 'SGDID']
            gid=t_sgd.ix[i, 'gid']
            if gid!=c.get(sgdid, 0):
                data.append({'gid':gid, 'name':sgdid})
            if not pd.isnull(s_name1) and gid!=c.get(s_name1, 0):
                data.append({'gid':gid, 'name':s_name1})
            if not pd.isnull(s_name2) and gid!=c.get(s_name2, 0):
                data.append({'gid':gid, 'name':s_name2})
            if s_name1=='' or s_name2=='':
                print t_sgd[i]
        t=pd.DataFrame(data)
        t.drop_duplicates(inplace=True)
        t['ds']='SGD'
        insert_sql = '''
                INSERT INTO gid2source_id (  gid,  id_type_id,  source_id,  id_status,  ds,  tax_id)
                VALUES(?,?,?,?,?,?)
        '''
        for i, r in t.iterrows(): #name,gid,ds
            x.from_sql(insert_sql,
                       [
                    r['gid'],
                    6, #id_type_id,
                    r['name'], #''source_id',
                    'post_yeast_alias',#''id_status',
                    'SGD',#''ds',
                    4932,#'tax_id'
                       ])

        t=t_sgd[['gid','DESCRIPTION','ds']]
        t.drop_duplicates(inplace=True)
        t['DESCRIPTION']=t['DESCRIPTION'].apply(lambda x: re.sub('[\r\n\t]', ' ', x))

        insert_sql = '''
                INSERT INTO annotation (  gid,  annotation_type_id,  content,  annotation_field1,  ds,  tax_id)
                VALUES  (    ?,    11,    ?,    'post_yeast_summary',    'SGD',  4932  );
        '''
        for k, g in t.groupby('gid'):
            gid = k
            description = ';'.join(g['DESCRIPTION'])
            x.from_sql(insert_sql,
                       [
                    gid,
                    description
                       ])


    @staticmethod
    def add_yeast_microarray_id():

        x=db.DB(name=SyncDB.CONNECTION_ID, db=SyncDB.DATABASE)
        #delete old
        x.from_sql("DELETE FROM gid2source_id  WHERE id_type_id = 14 AND id_status='post_yeast_micorarray' ")

        fn_gz =os.path.join(SyncDB.DOWNLOAD_DIR(),'yeast2.db_3.2.3.tar.gz')
        fn_yeast = os.path.join(SyncDB.DOWNLOAD_DIR(),'yeast2.db/inst/extdata','yeast2.sqlite')
        cmd = 'tar xvf {0} -C {1}'.format(fn_gz,SyncDB.DOWNLOAD_DIR())
        util.unix(cmd)
        con_sqlite = db.open_sqlite(fn_yeast)
        df_probe_id = db.from_sql(con=con_sqlite,
                         sql='''select probe_id, sgd_id from probes
                                where sgd_id is not null
                         ''' )
        #print df_probe_id
        df_sgid = x.sql_in('''SELECT t.source_id, t.gid
                                FROM gid2source_id t
                                where t.source_id in (
                            ''',
                           ')',
                           df_probe_id['sgd_id'])
        df = df_probe_id.merge(df_sgid,
                               how='inner',
                               left_on='sgd_id',
                               right_on='source_id')
        #print df_sgid
        #print df

        insert_sql = '''
                INSERT INTO gid2source_id (  gid,  id_type_id,  source_id,  id_status,  ds,  tax_id)
                VALUES(?,?,?,?,?,?)
        '''
        for i, r in df.iterrows(): #name,gid,ds
            x.from_sql(insert_sql,
                       [
                    r['gid'],
                    14, #id_type_id,
                    r['probe_id'], #''source_id',
                    'post_yeast_micorarray',#''id_status',
                    'bioconduct',#''ds',
                    4896,#'tax_id'
                       ])



    @staticmethod
    def do_postprocess():
        if not os.path.exists(PostProcess.tmp_directory):
            util.unix("mkdir "+PostProcess.tmp_directory)
        if True:
            PostProcess.gene_id_cleanup()
            # generate additional synonyms (we need to add a new data source called PlasmoDB)
            PostProcess.add_pf_alias()
            #YeastAlias.csv and YeastSummary.csv
            #first file add additional synonyms (need to add a new data source called SGD, second file can be uploaded as annotation.summary for yeast).
            PostProcess.add_yeast_alias()

        PostProcess.add_yeast_microarray_id()
        #add more orthoMCL from to homologene table
        PostProcess.orthoMCL()





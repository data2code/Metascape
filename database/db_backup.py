#!/usr/bin/env python
from os import sys, path
import socket

import shlex
from IPython.core.debugger import Tracer

sys.path.insert(0, path.join(path.dirname(path.abspath(__file__)),'mylib'))
    
import pandas as pd
import util
import time
import glob
import sys
import db
import argparse as arg
from pprint import pprint
import os
import xml.etree.ElementTree as ET
import itertools
from  core import SyncDB
from pprint import pprint


class DatabaseDump():
    def __init__(self, args):
        self.CONNECTION_ID = args.connection
        self.DATABASE = args.database
        self.ACTION = args.action
        self.OUTPUT = args.output
        self.INPUT = args.input        
        self.RESTORE_DB = args.restoredb
        self.SYNC_CATEGORY = args.sync_category_data
        
    def run(self):
        if self.SYNC_CATEGORY !='':
            self.sync_category_tables(self.SYNC_CATEGORY);
        elif self.ACTION == 'schema':
            self.dump_schema()
        elif self.ACTION == 'backup':
            self.dump_all()
        elif self.ACTION == 'restore':
            self.restore()
        elif self.ACTION == 'getdata':
            self.dump_data_tables()
        elif self.ACTION == 'restore_genego':
            self.insert_genego_to_term_category()

    def restore(self):
        print 'restoring...'
        if len(self.RESTORE_DB) == 0:
            print "You need to specify a database name to be restored using -r option."
            exit()

        if len(self.INPUT) == 0:
            print "You need to specify a .sql file name -i option."
            exit()
            
        conn_info = db.get_con_info(self.CONNECTION_ID)
        if conn_info is None:
            print "cannot find connection for " + self.CONNECTION_ID
            exit()
        cmd = [
            'mysql',
            '--host ' + conn_info["HOST"],
            '-u ' + conn_info["USR"],
            '-p' + conn_info["PWD"],
            '-e',
            '"create database {0}"'.format(self.RESTORE_DB),
            ]
        print ' '.join(cmd)
        util.unix(cmd)
 
        cmd = [
            'mysql',
            '--host ' + conn_info["HOST"],
            '-u ' + conn_info["USR"],
            '-p' + conn_info["PWD"],
            self.RESTORE_DB,
            '<' + self.INPUT]
        print ' '.join(cmd)
        util.unix(cmd)
        
    def dump_all(self):
        print 'backup...'
        conn_info = db.get_con_info(self.CONNECTION_ID)
        if conn_info is None:
            print "cannot find connection for " + self.CONNECTION_ID
            exit()

        tables = ['annotation',
                  'annotation_type',
                  'gid2source_id',
                  'gid2terms',
                  'homologene',
                  'id_type',
                  'interaction',
                  'interaction_type',
                  'taxid2name',
                  'term',
                  'term_category',
                  'term2gids',
                  'term2term',
                  'statistics',
                  ]
        cmd = [
            'mysqldump',
            '--host ' + conn_info["HOST"],
            '-u ' + conn_info["USR"],
            '-p' + conn_info["PWD"],
            self.DATABASE,
	        ' '.join(tables),
            '>' + self.OUTPUT]
        print cmd
        util.unix(cmd)


    def insert_genego_to_term_category(self):
        with db.DB(self.CONNECTION_ID) as mydb:
            q = '''
            INSERT  INTO term_category(`term_category_id`,`category_name`,`category_group`,`description`,`ds`,`ds_name`,`url`,`default_selected`,`display_order`,`used_in_enrichment`,`category_group_name_membership`,`used_in_membership`,`display_order_membership`) 
                VALUES 
                (27,'GeneGo Pathway',2,'Various curated terms from GeneGo, GO and KEGG.  Here we first enter keywords and find the terms of interest, then Membership Analysis will flag a gene as Y, as long as it falls into any of the selected terms.  For memberships such as secreted proteins, transmembrane proteins, we probably would prefer to use precompiled Custom Gene Sets instead of reconstructing them via keyword searches.','GeneGo','GeneGo','local','N',1,'Y','4_Pathway','N',2),
                (28,'GeneGo Disease Association',1,'GeneGo disease-gene associations (disease markers)','GeneGo','GeneGo','local','N',2,'Y','3_Functional Set','Y',8),
                (29,'GeneGo Drug Target',4,'Human-curated lists of drug targets.','GeneGo','GeneGo',NULL,'N',3,'Y','5_Signature Module','Y',5),
                (31,'GeneGo GO Processes',2,'Human-curated lists of pathways, shown as pathway maps.','GeneGo','GeneGo','local','N',4,'Y','3_Functional Set','N',3)
            '''
            mydb.from_sql(q)


    def dump_schema(self):
        tables_with_data = ['annotation','annotation_type','gid2source_id','gid2terms','homologene',
                            'id_type','interaction','interaction_type','server_status','session','statistics',
                            'taxid2name','term','term2gids','term2term','term_category']
        conn_info = db.get_con_info(self.CONNECTION_ID)
        if conn_info is None:
            print "cannot find connection for " + self.CONNECTION_ID
            exit()
        metascape_db_in_xxx = 'gp_metascape'
        #todo uncomment next line
        if len(tables_with_data):
            cmd = [
                'mysqldump',
                '--host ' + conn_info["HOST"],
                '-u ' + conn_info["USR"],
                '-p' + conn_info["PWD"],
                metascape_db_in_xxx,
                ' '.join(tables_with_data),
                '>' + self.OUTPUT + "_P2"]

            print " ".join(cmd)
            util.unix(cmd)
            util.unix("cat " + self.OUTPUT + "_P2 >> " + self.OUTPUT)
            util.unix("rm " + self.OUTPUT + "_P2")

    def dump_data_tables(self):
        conn_info = db.get_con_info(self.CONNECTION_ID)
        if conn_info is None:
            print "cannot find connection for " + self.CONNECTION_ID
            exit()
        
        tables = ['annotation', 'gid2source_id', 'gid2terms', 'homologene', 'term', 'term2gids', 'term2term']

        cmd = [
            'mysqldump',
            '--host ' + conn_info["HOST"],
            '-u ' + conn_info["USR"],
            '-p' + conn_info["PWD"],
            self.DATABASE,
            ' '.join(tables),
            '>' + self.OUTPUT]

        util.unix(cmd)        
        print " ".join(cmd)
        
    def sync_category_tables(self, dest_db):
        print 'synchronizing category tables of the gp_devel with the production gp.'
        con_prod = db.get_con_info("MYSQL01_RW")
        cmd = [
            'mysqldump',
            '--host ' + con_prod["HOST"],
            '-u ' + con_prod["USR"],
            '-p' + con_prod["PWD"],
            'gp',
            ' '.join(['annotation_type', 'id_type', 'term_category']),
            '>category_tables.sql']
            
        print ' '.join(cmd)            
        util.unix(cmd)            

    
if __name__=='__main__':
    util.no_warning()
    print '--------------------------------------------------------------------------------------------------------------------'
    print '--------------------------------------------------------------------------------------------------------------------'
    opt=arg.ArgumentParser(description='Database dump')
    opt.add_argument('-c','--connection', default='GP_DEVEL', help='')    
    opt.add_argument('-o','--output', default='dump.sql', help='output file')
    opt.add_argument('-i','--input', default='', help='input sql file for restoring a database')    
    opt.add_argument('-d','--database', default='gp', help='database name')
    opt.add_argument('-r','--restoredb', default='', help='database name to be restored. (this is required if action is restore)')    
    opt.add_argument('-a','--action', default='schema', help='action')
    opt.add_argument('-s','--sync_category_data', default='', help='copy the contents of the category tables (i.e. term_category, annotation_type, id_type) from go to gp_devel')

    argString = None
    args = opt.parse_args(shlex.split(argString) if argString else None)

    dbd = DatabaseDump(args)
    dbd.run()
    print '---------------'
    

    exit()
    

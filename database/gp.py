#!/usr/bin/env python
import sys
from os import path

sys.path.insert(0, path.join(path.dirname(path.abspath(__file__)),'mylib'))
from pprint import pprint
#pprint(sys.path)


import shlex
import urllib
import socket
from IPython.core.debugger import Tracer

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
from meta_statistics import Statistics
from  postprocess import PostProcess
from myutils import  MyUtils

if __name__=='__main__':
    util.no_warning()
    print '--------------------------------------------------------------------------------------------------------------------'
    print '--------------------------------------------------------------------------------------------------------------------'
    opt=arg.ArgumentParser(description='Gene Database sync')
    opt.add_argument('-n','--connection', default='', help='')    
    opt.add_argument('-c','--config', default='config.xml', help='')
    opt.add_argument('-r','--rebuild', default='Y', help='rebuild or sync, defualt rebuild=Y')
    opt.add_argument('-s','--source', default='', help='source, seperated by,')
    opt.add_argument('-d','--database', default='gp', help='database name')
    opt.add_argument('-a','--debug_node_ids', default='', help='node ids, rerun these nodes only')
    opt.add_argument('-a1','--debug_node_ids_only', default=False, action='store_true', help='only debug some node')
    opt.add_argument('-a2','--debug', default=False, action='store_true', help='node ids, rerun these nodes only')

    opt.add_argument('-H','--HOME', default='', help='Analysis HOME directory')
    opt.add_argument('-t','--history', action='store_true',help='statistic history')
    opt.add_argument('-v','--validate', action='store_true',help='validate urls and local files')    
    opt.add_argument('-sp','--skippostprocess', action='store_true',help='skip post process (cleaning retired gene ids)')
    opt.add_argument('-sv','--skipvalidation', action='store_true',help='skip datasource validation process (urls and database connections)')    
    opt.add_argument('-p','--report',  default='' ,help='produce HTML build reports in provided path')
    opt.add_argument('-el','--errorlog',  default='error_log.txt',help='error logs.')
    opt.add_argument('-rf','--recoveryfolder',  default='',help='A floder that contains last correct build to recover data from in case of error in the current build.')
    opt.add_argument('-sr','--skiprecovery', action='store_true',help='do not recover in case of error')
    opt.add_argument('-d1','--makestatisticshtml', default=False,action='store_true',help='make statistics html')

    import  datetime
    timenow_str = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    argString = '{0}'.format(debug_str)
    #--debug_node_ids_only this is just a flag. we must have node ids,  which split,by,comma,space,prefix,\,
    argString = '{0} --debug_node_ids CorumDownload,TermGenePair:CORUM,TermGene:CORUM'.format(debug_str,timenow_str)



    #todo (debug-done) uncomment one line below in production
    argString = None

    args = opt.parse_args(shlex.split(argString) if argString else None)

    if not args.HOME:
        import os.path
        path = os.path.split(os.path.abspath(sys.argv[0]))[0]
        args.HOME=path
    SyncDB.set_HOME_DIR(args.HOME)
    SyncDB.DONT_RECOVER = args.skiprecovery
    SyncDB.RECOVERY_FOLDER=args.recoveryfolder
    SyncDB.ERROR_LOG=open(args.errorlog,"w");
    SyncDB.CONNECTION_ID = args.connection
    SyncDB.SUPPORTED_SPECIES = {'9606','10090','10116'}
    SyncDB.DEBUG_NODE_IDS_ONLY = args.debug_node_ids_only
    if args.database:
        SyncDB.DATABASE = args.database

##########################################
    if args.debug:
        s=Statistics()
        s.generate_statistic_html()
        #PostProcess.add_pf_alias()
        #PostProcess.do_postprocess()
        exit()

    #Begin some util commands. Each must have exit
    if args.makestatisticshtml:
        #SyncDB.history()
        s=Statistics()
        s.generate_statistic_html()
        exit()
    #End some util commands




    #useless old update logging
    if args.report:
        SyncDB.stat_js(args.report)
        exit()
        SyncDB.report(args.report)
        exit()
        
    if args.history:
        SyncDB.history()
        exit()

    if args.rebuild:
        SyncDB.IS_REBUILD = args.rebuild == 'Y'
        

    if args.debug_node_ids:
        SyncDB.DEBUG_NODE_IDS = args.debug_node_ids.split(',')
        SyncDB.IS_REBUILD = False
        
    if args.source:
        SyncDB.TO_SYNC_DS = args.source.split(',')
                
    cr = SyncDB.load_xml_file(args.config)
    if not args.skipvalidation and not cr.check_inputs():
        print "Some input files are missing or urls are broken!"
        SyncDB.ERROR_LOG.close()
        exit()
        
    if args.validate and not args.skipvalidation:
        print "Validation complete."
        SyncDB.ERROR_LOG.close()
        exit()
        
    cr.do_work()

    if not args.skippostprocess:
        PostProcess.do_postprocess()



    print "gp database was successfully built";
    SyncDB.ERROR_LOG.close()
    exit()
    

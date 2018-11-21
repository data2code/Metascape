#!/usr/bin/env python3
from __future__ import absolute_import
from __future__ import print_function
import os
import sys
from six.moves import range
from six.moves import zip
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../mylib'))
import ms.biolists as bl
import ms.msobject as msobj
import pandas as pd
import numpy as np
import util
import os
import entrez as ez
import go as go
import ppi as ppi
import ms.circos as circos
import ms.report as rpt
import re
import glob
import stats
import ms.evidence_regression as er
import cluster
import xgmml
import ms.report as report
import time

class MetaEngine:

    def __init__(self, s_output, S_session_files=None, l_WEB_MODE=False):
        if s_output is None:
            util.error_msg('Output folder must be specified!')
        self.mkdir(s_output)

        self.session={
            's_out':s_output,
            #'s_out2':s_output+"/PPI",
            's_cache':s_output+"/Cache",
            # global settings
            'tax_id':9606,
            'myez':None,
            'mygo':None,
            'myppi':None,
            'CPU':4,
            # gene lists
            'n_list':None,
            'lists':None,
            'list_benchmark':None,
            'S_background':None,
            'S_color':None,
            'R_weight':None, #Weights for gene lists is a list/array of +/-1, 1 means hits are desirable, -1 means hits are undesirable (counter assay), 0 means do not care
            # analysis switches
            'l_PPI':True,
            'l_GO':True,
            'l_SIMPLE':False, #if l_SIMPLE, both l_PLOT and l_REPORT will be considered as False
            'l_PLOT':True,
            'l_EXPORT':True,
            'l_REPORT':True,
            # for metascape web site
            'l_WEB_MODE':False,
            'l_EVIDENCE_ONLY':False,
            # GO analysis settings
            'max_clusters':20,
            'max_member_per_go_cluster':10,
            'min_list_size_for_go':3,
            'max_list_size_for_go':3000,
            'max_nof_enriched_go':2000, # change from 1000
            'S_go_category':[11,19,23,6,33,24],
            'cluster_similarity':0.3,
            'max_nodes_in_go_network':200,
            'min_overlap':3,
            'min_enrichment':0,
            'p_cutoff':0.01,
            'l_background_by_ontology':False,
            'l_go_selective':False,
            'l_go_piechart':False,
            'l_L1k_only':False,
            # circos
            'circos_link_logp':-3,
            'circos_link_go_size':100,
            'circos_link_enrichment':2.0,
            'circos_show_symbol':False,
            'circos_max_go_terms':250, # sample no more than this amount for edge display, prevent too many edges
            # PPI analysis settings
            'l_indirect_PPI':True, #?
            'max_size_for_indirect_ppi':10,
            'min_ppi_size':3,
            'max_ppi_size':1000,
            'max_list_size_for_ppi':3000,
            'ppi_datasource':None,
            'l_merge_for_ppi':True,
            'l_connect_in_merge':False, # MCODE, connect interactions cross MCODE clusters
            'l_exclude_mcode_evidence':True,
            'max_ppi_size_for_plot':500, # don't plot network if it's too complicated
            # GPEC analysis settings
            'l_GPEC':True,
            'min_evi_term_size':3,
            'max_evi_term_size':100,
            'logp_evi_term':-2.0,
            'min_benchmark_for_ML':10, # at least 10 known genes to start machine learning the evidence weights
            'l_BYPASS_ML':False, # force to skip ML, even if there are enough benchmark genes
            'max_lists_ML_individual':5, # if <=5, we regress weights for individual evidence lines (15 total), if more lists
                                        # we first combine evidence lines into E_hit, E_go, E_ppi, then regress these 3
            'max_GPEC_iteration':10, # max iteration for RSA-GO
            'gpec_ppi_target_size':250,
            'max_go_annotation_per_mcode':3,
            'golists':None,
            'full_nets':None,
            't_gene_mem':None,
            'gpec_gonet_fn':None,
            'gpec_gonet_nnode':None,
            'gpec_gonet_ncluster':None,
            't_evi_ppi':None,
            't_evi':None,
            'gpec_net': None,
            'search_words':'',
            'CIRCOS_BIN':None,
            'CYTOSCAPE_ON_WINDOWS': False,
            'CYTOSCAPE_HOST':'localhost',
            'CYTOSCAPE_PORT':1234,
            'ADD_LOGO': True
        }

        S=s_output.strip('/').split('/')
        s_session_id=S[-2] if (S[-1]=='output' and len(S)>1) else S[-1]
        util.MSG_PREFIX="%s:%s|" % (s_session_id, os.getpid())

        if S_session_files is not None:
            if type(S_session_files) is dict:
                self.puts(S_session_files)
            elif type(S_session_files) is str or type(S_session_files) is list:
                if type(S_session_files) is str:
                    self.load_session(S_session_files)
                else:
                    for x in S_session_files:
                        self.load_session(x)
            else:
                util.error_msg('Unrecognize S_session_files argument type: %s!' % str(type(S_session_files)))
        else:
            if l_WEB_MODE:
                self.restore_session()
        self.put('l_WEB_MODE', l_WEB_MODE)
        self.is_okay() # check if settings are consistent
        self.sw=util.StopWatch("MetaEngine")

    def reset_gene_lists(self, genelists, S_color=None, R_weight=None):
        self.puts({'lists': genelists, 'S_color':S_color, 'R_weight':R_weight, 'n_list':len(genelists)})

    def save_session(self, s_out, S_names=None):
        s_cache=self.get('s_cache')
        self.mkdir(s_cache)
        if S_names is None:
            S_names=[x for x in self.session.keys() if x not in ('myez','mygo','myppi')]
        else:
            S_names=[x for x in S_names if x in self.session]
        my_session={k:self.session[k] for k in S_names}
        msobj.MSObject.dump_object(my_session, s_out, s_cache_dir=s_cache)

    def read_session(self, s_session_file=None):
        if s_session_file is None: return {}
        s_file=os.path.join(self.get('s_cache'), s_session_file)
        if not os.path.exists(s_file):
            #util.warn_msg('Session file: %s does not exist, skip!' % s_file)
            return {}
        else:
            c=msobj.MSObject.load(s_file)
            if c is None:
                util.warn_msg('Session file: %s failed to load, skip!' % s_file)
                return {}
            return c

    def load_session(self, s_session_file=None):
        s_file=os.path.join(self.get('s_cache'), s_session_file)
        if not os.path.exists(s_file): return
        l_base=s_session_file=='genelists_overlap.pickle'
        n_try = 3 if l_base else 1 # base file might be affected by another process, so may need a retry
        i=1
        while i<=n_try:
            c=self.read_session(s_session_file)
            if (l_base and 'lists' in c) or i>=n_try: break
            time.sleep(3)
            i+=1
        if len(c):
            self.session.update(c)
        #print("++++++++++++++++++++++++", s_session_file, "+++++++++++++")
        #S=list(c.keys())
        #print(S)
        #if 'lists' in S: print(c['lists'])
        #if 'golists' in S: print(c['golists'])
        #print("++++++++++++++++++++++++", s_session_file, "+++++++++++++")


    def update_session(self, s_session_file):
        """rewrite session file, as values may have changed"""
        c=self.read_session(s_session_file)
        self.save_session(s_session_file, list(c.keys()))

    def restore_session(self):
        """Expect .pickle files under s_out/s_cache/*pickle"""
        for x in ['genelists_overlap.pickle', 'golists.pickle', 'ppi.pickle', 'evidence.pickle','gpec_go.pickle', 'gpec_ppi.pickle']:
            self.load_session(s_session_file=x)

    def put(self, name, value):
        if name=='S_background':
            if value is not None:
                if type(value) is not set:
                    value=set(value)
                if type(list(value)[0])!=str: value={str(x) for x in value}
                if len(value)==0:
                    value=None
            if value is not None:
                if self.session['lists'] is not None:
                    for x in self.get('lists').values:
                        value|=x.data
                value=list(value)
        elif name=='S_go_category':
            if value is None or len(value)==0:
                util.warn_msg('Providing zero length S_go_category, ignored!')
                return
        self.session[name]=value

    def puts(self, data):
        for k,v in data.items():
            self.put(k, v)
        #self.session.update(data)

    def get(self, name):
        value=self.session.get(name, None)
        if value is not None: return value
        if name=='myez':
            value=ez.EntrezGene(tax_id=self.get('tax_id'), l_use_GPDB=True)
        elif name=='mygo':
            value=go.GO(l_use_GPDB=True, tax_id=self.get('tax_id'), entrez=self.get('myez'))
        elif name=='mygo_L1k':
            value=go.GO(l_use_GPDB=True, tax_id=9606, entrez=self.get('myez'), l_L1k=True)
        elif name=='myppi':
            value=ppi.PPI(tax_id=self.get('tax_id'), S_DB=self.get('ppi_datasource'), l_use_GPDB=True)
        elif name=='S_name':
            value=self.get('lists').names
        elif name=='n_list':
            value=len(self.get('lists'))
        elif name=='S_color':
            if self.get('lists') is None: return None
            value=bl.CytoPlot.get_qualitative_colors(len(self.get('lists')))
        elif name=='R_weight':
            if self.get('lists') is None: return None
            n_list=len(self.get('lists'))
            value=np.ones(n_list)
        elif name=='s_out':
            value='Output'
        elif name=='s_out2':
            value=self.get('s_out')+"/PPI"
        elif name=="s_cache":
            value=self.get('s_out')+"/Cache"
        if value is not None:
            self.put(name, value)
        return value

    def is_okay(self):
        if self.get('lists') is not None:
            n_list=len(self.get('lists'))
            if self.session.get('S_color', None) is not None:
                if len(self.get('S_color'))!=n_list:
                    util.error_msg('Length of color must match length of gene lists: %d vs %d!' % (len(self.get('S_color')), len(self.get('lists'))))
            if self.session.get('R_weight', None) is not None:
                if len(self.session.get('R_weight'))!=n_list:
                    util.error_msg('Length of weight must match the length of gene lists: %d vs %d!' % (len(self.get('R_weight')), len(self.get('lists'))))


    def mkdir(self, s_out):
        if s_out is None:
            util.error_msg('Cannot create None folder!')
        if s_out in ('.', ''):
            return s_out
        n_try=3
        i_try=1
        while True:
            os.makedirs(s_out, exist_ok=True)
            if os.path.exists(s_out):
                return s_out
            if i_try>=n_try:
                util.error_msg('Fail to create folder %s after %d attempts: %s' % (s_out, i_try))
                break
            else:
                i_try+=1
                time.sleep(0.5)
        return s_out

    def get_subfolder(self, s_out2):
        s_out=self.get('s_out')
        s_out2=s_out+"/"+s_out2
        return self.mkdir(s_out2)

    def genelists_overlap_analysis(self):
        #if self.session['lists'] is None and l_WEB_MODE:
        #    self.load_session('genelists.pickle')
        l_WEB_MODE=self.get('l_WEB_MODE')
        s_out=self.get('s_out')
        s_out2=self.get_subfolder('Overlap_circos')
        lists=self.get('lists')
        n_list=self.get('n_list')
        S_color=self.get('S_color')
        S_name=self.get('S_name')
        R_weight=self.get('R_weight')
        CIRCOS_BIN=self.get('CIRCOS_BIN')
        l_PLOT=self.get('l_PLOT')
        t_gene_mem=lists.membership_table(l_gene_as_index=True)
        self.sw.check("Membership")
        s_symbol=None
        if self.get('circos_show_symbol'):
            myez=self.get('myez')
            t_gene_mem['Symbol']=t_gene_mem.Gene.apply(lambda x: myez.C_GENENEW.get(x, x))
            s_symbol="Symbol"
        #print(">>>>>>>>>>>>>", l_PLOT, n_list)
        if l_PLOT and n_list>1:
            cir=circos.Circos(BIN=CIRCOS_BIN)
            cir.plot_membership(t_gene_mem, outputdir=s_out2, outputfile="CircosOverlapByGene", S_color=S_color, s_symbol="Symbol")
        t_gene_mem.to_csv(s_out2+'/MembershipByGene.csv', index=False)

        # Machine Learning Component
        gl_benchmark=self.get('list_benchmark')
        pos=neg=0
        if gl_benchmark is not None:
            S_gene=set(t_gene_mem.Gene.tolist())
            pos=len([x for x in S_gene if x in gl_benchmark.data])
            neg=len(S_gene)-pos

        R_weight=self.get('R_weight')
        c_sign={}
        if R_weight is not None:
            if len(R_weight)!=n_list:
                util.warn_msg('Lenght of R_w is not the same as num of gene lists!!!')
            c_sign={x:R_weight[i] for i,x in enumerate(lists.names) }
        self.sw.check("genelists_overlap_analysis")
        self.puts({'t_gene_mem': t_gene_mem, 'benchmark_pos':pos, 'benchmark_neg':neg, 'c_sign':c_sign})
        if l_WEB_MODE:
            self.save_session('genelists_overlap.pickle') #, ['lists','n_list','list_benchmark','S_name','R_weight','S_color','t_gene_mem','S_background','benchmark_pos', 'benchmark_neg', 'c_sign'])

    def go_analysis(self):
        #if self.session['lists'] is None and l_WEB_MODE:
        #    self.load_session('genelists_overlap.pickle')

        def go_analysis_L1k(lists, n_CPU, max_list_size, S_go_category, SRC_GENE, min_overlap, p_cutoff, min_enrichment):
            S_go_category_L1k=[x for x in S_go_category if x>90]
            S_go_category=[x for x in S_go_category if x<90]
            #print(S_go_category, S_go_category_L1k)
            golists=golists_L1k=None
            if len(S_go_category)>0:
                mygo=self.get('mygo')
                golists=lists.go_analysis(mygo, n_CPU=CPU, max_list_size=max_list_size, S_go_category=S_go_category, SRC_GENE=S_background, min_overlap=min_overlap, p_cutoff=p_cutoff, min_enrichment=min_enrichment)
            if len(S_go_category_L1k):
                mygo=self.get('mygo_L1k')
                golists_L1k=lists.go_analysis(mygo, n_CPU=CPU, max_list_size=max_list_size, S_go_category=S_go_category_L1k, SRC_GENE=S_background, min_overlap=min_overlap, p_cutoff=p_cutoff, min_enrichment=min_enrichment)
            if golists_L1k is None:
                return golists
            if golists is None:
                return golists_L1k
            golists.combine(golists_L1k)
            return golists

        l_WEB_MODE=self.get('l_WEB_MODE')
        s_out=self.get('s_out')
        lists=self.get('lists')
        n_list=self.get('n_list')
        #mygo=self.get('mygo')
        max_list_size_for_go=self.get('max_list_size_for_go')
        CPU=self.get('CPU')
        S_go_category=self.get('S_go_category')
        S_background=self.get('S_background')
        min_evi_term_size=self.get('min_evi_term_size')
        max_evi_term_size=self.get('max_evi_term_size')
        min_overlap=self.get('min_overlap')
        min_enrichment=self.get('min_enrichment')
        p_cutoff=self.get('p_cutoff')
        logp_evi_term=self.get('logp_evi_term')

        s_out2=self.get_subfolder("Enrichment_GO")
        golists=go_analysis_L1k(lists, n_CPU=CPU, max_list_size=max_list_size_for_go, S_go_category=S_go_category, SRC_GENE=S_background, min_overlap=min_overlap, p_cutoff=p_cutoff, min_enrichment=min_enrichment)
        t_evi_go=golists.gene_evidence_table(min_size=min_evi_term_size, max_size=max_evi_term_size, logp=logp_evi_term)
        t_go_mem=golists.membership_table(S_matrix="LogP")
        # as the order of gene list in t_go_mem is different, we need to get a new color and label list
        # _MEMBER_* are sorted, also some genelists may not have enriched go_terms in t_go_mem
        l_go_selective=self.get('l_go_selective')
        if l_go_selective and n_list>1:
            S_logp=[x for x in t_go_mem.header() if x.startswith('_LogP_')]
            X=np.clip(np.abs(t_go_mem[S_logp].values), 1.0, None) # set min to 1.0
            t_go_mem['GiniIndex']=np.apply_along_axis(stats.gini, 1, X)
        else:
            t_go_mem['GiniIndex']=0
        self.sw.check("GO membership: t_go_mem %d x %d" % (len(t_go_mem), len(t_go_mem.header())))
        self.sw.check("GO Enrichment")
        self.puts({'t_evi_go': t_evi_go, 'golists':golists, 't_go_mem':t_go_mem})
        if l_WEB_MODE:
            # no longer needed
            #self.update_session('genelists_overlap.pickle')
            self.save_session('golists.pickle', ['golists','t_evi_go','t_go_mem','max_list_size_for_go','S_go_category','S_background','min_evi_term_size','max_evi_term_size','log[_evi_Term' ])
        golists.save(s_out2+'/GO_AllLists.csv')
        util.df2sdf(t_go_mem).to_csv(s_out2+'/GO_membership.csv')

    def ppi_analysis(self):
        #if self.session['lists'] is None and l_WEB_MODE:
        #    self.load_session('genelists_overlap.pickle')
        l_WEB_MODE=self.get('l_WEB_MODE')
        s_out=self.get('s_out')
        #s_out2=self.get('s_out2')
        lists=self.get('lists')
        n_list=self.get('n_list')
        s_out2=self.get_subfolder("Enrichment_PPI")
        s_out3=self.get_subfolder("Enrichment_PPI/xgmml")
        s_merge_name="MERGE"
        myppi=self.get('myppi')
        l_MCODE=True #not self.get('l_SIMPLE')
        if self.get('l_EVIDENCE_ONLY'): l_MCODE=False
        self.sw.check('Start PPI Analysis')
        min_ppi_size=self.get('min_ppi_size')
        max_ppi_size=self.get('max_ppi_size')
        l_merge_for_ppi=self.get('l_merge_for_ppi')
        l_indirect_ppi=self.get('l_indirect_ppi')
        max_size_for_indirect_ppi=self.get('max_size_for_indirect_ppi')
        max_list_size_for_ppi=self.get('max_list_size_for_ppi')
        CPU=self.get('CPU')
        (all_nets, t_mcode, t_over, full_nets, t_evi_ppi)=lists.ppi_analysis(myppi, s_name="MERGE", l_MCODE=l_MCODE, min_size=min_ppi_size, max_size=max_ppi_size, l_overconnect=False, l_propagation=False, opt_over=None, n_CPU=CPU, l_MCODE_optimize=False, l_individual=True, l_merge=False, l_indirect_PPI=False, indirect_size=max_size_for_indirect_ppi, max_list_size=max_list_size_for_ppi)
        mcode_nets=[x for x in all_nets if x.name.endswith('_MCODE_ALL')]
        self.puts({'t_evi_ppi': t_evi_ppi, 't_mcode':t_mcode, 'full_nets': full_nets, 'mcode_nets':mcode_nets})
        if l_WEB_MODE:
            #self.update_session('genelists_overlap.pickle')
            self.save_session('ppi.pickle', ['t_evi_ppi','t_mcode','full_nets','mcode_nets','min_ppi_size','max_ppi_size'])
        if mcode_nets is not None:
            for x in mcode_nets:
                s_name=x.name
                x.to_xgmml(s_out3+"/"+s_name)
        if full_nets is not None:
            for k,v in full_nets.items():
                if v is not None:
                    v.to_xgmml(s_out3+"/"+k)
        if t_mcode is not None and len(t_mcode):
            t_mcode.to_csv(s_out2+"/MCODE.csv", index=False)

    def GPEC_evidence(self):
        l_WEB_MODE=self.get('l_WEB_MODE')
        l_GO=self.get('l_GO')
        l_PPI=self.get('l_PPI')
        #if l_WEB_MODE:
        #    if self.session['lists'] is None: self.load_session('genelists_overlap.pickle')
        #    if l_GO and self.session['golists'] is None: self.load_session('golists.pickle')
        #    if l_PPI and self.session['full_nets'] is None: self.load_session('ppi.pickle')
        s_out=self.get('s_out')
        lists=self.get('lists')
        n_list=self.get('n_list')
        t_gene_mem=self.get('t_gene_mem')
        t_evi_hit=t_gene_mem.copy()
        t_evi_go=self.get('t_evi_go')
        t_evi_ppi=self.get('t_evi_ppi')
        l_GPEC=self.get('l_GPEC')
        c_sign=self.get('c_sign')
        #print t_evi_hit.header()
        #print t_evi_go.header()
        #print t_evi_ppi.header()
        #exit()
        c_rename={}
        S=['Gene']
        for x in t_evi_hit.header():
            if x.startswith('_MEMBER_'):
                c_rename[x]='EvidenceHit:%s' % re.sub(r'^_MEMBER_','',x)
                S.append(c_rename.get(x, x))
        t_evi_hit.rename2(c_rename)
        t_evi=t_evi_hit[S]
        #print t_evi.col_types(l_dict=True)
        if t_evi_go is not None and len(t_evi_go)>0:
            #print t_evi_go.col_types(l_dict=True)
            t_evi=t_evi.merge(t_evi_go, left_on='Gene', right_on='Gene', how='outer')
        if t_evi_ppi is not None and len(t_evi_ppi)>0:
            #print t_evi_ppi.col_types(l_dict=True)
            l_exclude_mcode_evidence=self.get('l_exclude_mcode_evidence')
            if l_exclude_mcode_evidence:
                S_mcode=[x for x in t_evi_ppi.header() if x.startswith('EvidenceMCODE')]
                if len(S_mcode):
                    t_evi_ppi.drop(S_mcode, axis=1, inplace=True)
            t_evi=t_evi.merge(t_evi_ppi, left_on='Gene', right_on='Gene', how='outer')
        t_evi=t_evi.fillna(0)
        S_evi=[x for x in t_evi.header() if x.startswith('Evidence')]
        n_evi=len(S_evi)
        R_w=np.ones(n_evi)
        for i,s in enumerate(S_evi):
            s=re.sub(r'.+:', '', s)
            if s not in c_sign:
                util.warn_msg('Gene list name in Evidence table not recognized: %s!' % s)
            R_w[i]=c_sign.get(s, 1)
        X=t_evi[S_evi].values
        t_evi=t_evi.astype(int)
        t_evi['Gene']=t_evi.Gene.astype(str)
        t_evi['TotalEvidenceCounts']=(X*R_w).sum(axis=1)
        t_evi=t_evi.move_column('TotalEvidenceCounts', 1)

        # Machine Learning Component
        gl_benchmark=self.get('list_benchmark')
        #print "======================="
        #print gl_benchmark
        min_benchmark_for_ML=self.get('min_benchmark_for_ML')
        #print l_GPEC, min_benchmark_for_ML
        accuracy=None
        l_BYPASS_ML=self.get('l_BYPASS_ML')
        s_out2=self.get_subfolder("Evidence")
        if l_GPEC and gl_benchmark is not None:
            t_evi['TruePositives']=[ 1 if x in gl_benchmark.data else 0 for x in t_evi.Gene ]
            y=t_evi.TruePositives.values
            #print ">>>>>>>>>>>>>>>>>>>>", sum(y), sum(1-y)
            if not l_BYPASS_ML and (sum(y)>=min_benchmark_for_ML and sum(1-y)>=min_benchmark_for_ML): # machine learning
                n_max_lists=self.get('max_lists_ML_individual')
                l_evi_merge=False
                if n_list > n_max_lists and np.all(R_w>0.01): # too many evidence lines, aggregate them into three
                    X=[]
                    S_col=[]
                    for s in ['Hit','GO','PPI']:
                        S_col2=[x for x in t_evi.header() if x.startswith("Evidence"+s)]
                        if len(S_col2):
                            X.append(t_evi[S_col2].sum(axis=1).reshape(-1,1))
                            S_col.append(s)
                            print(X[-1].shape)
                    R_w=np.ones(len(X))
                    X=np.hstack(X)/n_list
                    l_evi_merge=True
                else:
                    S_col=[re.sub(r'^Evidence', '', x) for x in S_evi]
                print(S_col)
                self.sw.check('Optimizing evidence weights ...')
                # find optimal weights
                import traceback
                try:
                    R_w2=np.array([ np.sign(x) for x in R_w ])
                    (R, accuracy, y_pred, prob)=er.FindWeight.evidence_weight_wrapper(X, y, signs=R_w2, lb_penalty=-5, ub_penalty=5, num_penalty=11, l_auto_weight=True, X_test=X)
                    R_w=R
                    print("Successful machine learning")
                    #y_pred=(X*R_w[1:]).sum(axis=1)
                    t_evi['TotalEvidenceCounts']=y_pred
                    #y_pred+=R_w[0]
                    t_evi['Probability']=prob #1.0/(1.0+np.exp(-y_pred))
                    #S_col=[re.sub(r'^Evidence', '', x) for x in S_evi]
                    S_file=[s_out2+'/EvidenceWeight.png', s_out2+'/EvidenceWeight.pdf']
                    Weights=er.FindWeight.auto_weight(y)
                    er.FindWeight.plot(S_col, R_w, y, y_pred, t_evi.Probability.values, S_file, Weights)
                    t_w=pd.DataFrame({'Evidence':['_OFFSET_']+S_col, 'Weight':R_w})
                    t_w['RelativeWeight']=np.abs(R_w)/np.max(R_w[1:])
                    t_w['Sign']=np.sign(R_w)
                    t_w.ix[0, 'RelativeWeight']=0
                    t_w.to_csv(s_out2+'/EvidenceWeight.csv', index=False)
                except:
                    accuracy=None
                    util.warn_msg('Failed to machine learn the weights!')
                    print(traceback.format_exc())
            else:
                if l_BYPASS_ML:
                    util.warn_msg('By pass machine learning!')
                else:
                    util.warn_msg('Too few benchmark genes for machine learning: %d' % sum(y))
        print("Final Weights: ", R_w)
        t_evi.sort_values('TotalEvidenceCounts', ascending=False, inplace=True)
        t_evi.index=list(range(len(t_evi)))

        min_list_size_for_go=self.get('min_list_size_for_go')
        max_list_size_for_go=self.get('max_list_size_for_go')
        max_GPEC_iteration=self.get('max_GPEC_iteration')

        # old pandas version does not support keep=last in drop_duplicates, so we reverse it
        tmp=t_evi.TotalEvidenceCounts[::-1]
        R_idx=np.array(tmp.drop_duplicates().index[::-1])
        R_idx=R_idx[R_idx >= min_list_size_for_go] # at least three genes
        if max_list_size_for_go>0:
            R_idx=R_idx[R_idx< max_list_size_for_go]
        # pick at most max_iter cutoffs to run RSA
        if not l_GPEC: # only keep one option
            R_idx=R_idx[-1:]
        if len(R_idx)>max_GPEC_iteration:
            R_pick=cluster.cluster_array_to_k_groups(R_idx.reshape(len(R_idx),1), max_GPEC_iteration)
            R_idx=np.array(sorted([ R_idx[int(x[0])] for x in R_pick ]))
        R_evi=t_evi.ix[R_idx, 'TotalEvidenceCounts'].values

        self.puts({'t_evi':t_evi, 'R_idx':R_idx, 'R_evi':R_evi, 'S_evi':S_evi, 'R_GPEC_Weight':R_w, 'ML_accuracy':accuracy })
        if l_WEB_MODE:
            #self.update_session('genelists_overlap.pickle')
            self.save_session('evidence.pickle', ['t_evi', 'R_idx', 'R_evi', 'l_GPEC', 'l_GO', 'l_PPI', 'S_evi', 'R_GPEC_Weight','ML_accuracy'])
        if t_evi is not None and len(t_evi)>0:
            t_evi.to_csv(s_out2+'/Evidence.csv', index=False)

    def GPEC_go_analysis(self):
        #if l_WEB_MODE:
        #    if self.session['lists'] is None: self.load_session('genelists_overlap.pickle')
        #    if self.session['golists'] is None: self.load_session('golists.pickle')
        #    if self.session['t_evi'] is None: self.load_session('evidence.pickle')
        l_WEB_MODE=self.get('l_WEB_MODE')
        s_out=self.get('s_out')
        R_idx=self.get('R_idx')
        R_evi=self.get('R_evi')
        golists=self.get('golists')
        t_evi=self.get('t_evi')
        l_GPEC=self.get('l_GPEC')
        l_PLOT=self.get('l_PLOT')
        n_list=self.get('n_list')
        l_skip_GPEC=(not l_GPEC) and n_list<=1
        S_background=self.get('S_background')
        CPU=self.get('CPU')
        min_overlap=self.get('min_overlap')
        min_enrichment=self.get('min_enrichment')
        p_cutoff=self.get('p_cutoff')
        l_background_by_ontology=self.get('l_background_by_ontology')
        max_list_size_for_go=self.get('max_list_size_for_go')
        S_go_category=self.get('S_go_category')
        t_go=None
        s_out2=self.get_subfolder("Enrichment_GO")
        s_out3=self.get_subfolder("Enrichment_heatmap")
        sw2=util.StopWatch("MetaEngine::GPEC_go_analysis")

        if l_skip_GPEC:
            golists=self.get('golists')
            if len(golists)==1:
                t_go=golists.values[0].data
        else:
            S_go=[]
            for x in golists.values:
                S_go.extend(list(x.data.CategoryID.astype(str)+"_"+x.data.GO))
            S_go=util.unique(S_go)
            c_hitlist={}
            c_go={}
            #mygo=self.get('mygo')
            # RSA-style go analysis
            for i,idx in enumerate(R_idx):
                if idx>max_list_size_for_go: continue
                s_name=str(R_evi[i])
                c_hitlist[s_name] = t_evi.Gene.tolist()[:idx+1]
                c_go[s_name]=S_go
                #print s_name, len(c_hitlist[s_name]), len(S_go)
            if len(c_go)>0:
                S_go_category_L1k=[x for x in S_go_category if x>90]
                S_go_category=[x for x in S_go_category if x<90]
                if len(S_go_category)>0:
                    mygo=self.get('mygo')
                    t_go=mygo.analysis(c_hitlist, S_go=c_go, SRC_GENE=S_background, S_go_category=S_go_category, n_CPU=CPU, min_overlap=min_overlap, min_enrichment=min_enrichment, p_cutoff=p_cutoff, l_background_by_ontology=l_background_by_ontology)
                if len(S_go_category_L1k):
                    mygo=self.get('mygo_L1k')
                    t_go_L1k=mygo.analysis(c_hitlist, S_go=c_go, SRC_GENE=S_background, S_go_category=S_go_category_L1k, n_CPU=CPU, min_overlap=min_overlap, min_enrichment=min_enrichment, p_cutoff=p_cutoff, l_background_by_ontology=l_background_by_ontology)
                    if t_go is None:
                        t_go=t_go_L1k
                    elif t_go_L1k is not None:
                        t_go=pd.concat([t_go, t_go_L1k], ignore_index=True)
                        t_go=t_go.sort_values(['LogP','Enrichment','#GeneInGOAndHitList'], ascending=[True,False,False])
            #t_go.to_csv(s_out+"/_ALL_GPEC_GO.csv", index=False)
        if not l_GPEC and n_list>1:
            data=[]
            if t_go is not None: data.append(t_go)
            for x in golists.values:
                if len(x): data.append(x.data)
            if len(data):
                t_go=pd.concat(data, ignore_index=True)
                if 'Name' in t_go.header():
                    t_go.drop('Name', axis=1, inplace=True)
            #print t_go.header()

        gpec_golist=t_gpec_go=None
        l_go_selective=self.get('l_go_selective')
        s_xgmml=n_node=n_cluster=None

        if t_go is not None:
            t_go_mem=self.get('t_go_mem')
            t_go_mem2=t_go_mem.drop('Description', axis=1)
            t_go=t_go_mem2.merge(t_go, left_on='GO', right_on='GO')
            t_go.sort_values(['GO','CategoryID','GiniIndex','LogP'], inplace=True)
            #util.df2sdf(t_go).to_csv(s_out+'/t_go.csv', index=False)
            t_go.drop_duplicates(['GO','CategoryID'], inplace=True)
            #t_go.sort(['GiniIndex','LogP'], inplace=True)
            max_nof_enriched_go=self.get('max_nof_enriched_go')
            #if len(t_go)>max_nof_enriched_go:
            #    t_go=t_go[:max_nof_enriched_go]
            if 'Name' in t_go.header(): # Name comes from GPEC analysis
                t_go['EvidenceCutoff']=t_go.Name.astype(float)
                t_go.drop('Name', axis=1, inplace=True)
            else:
                t_go['EvidenceCutoff']=0
            S_background=self.get('S_background')
            golist=bl.GOList(t_go, 'Merged', genelist=None, go_opts={'SRC_GENE':S_background, 'min_overlap':min_overlap, 'min_enrichment':min_enrichment, 'p_cutoff':p_cutoff, 'n_CPU':CPU, 'S_go_category':S_go_category, 'l_background_by_ontology':l_background_by_ontology})
            #golist.data=t_go
            cluster_similarity=self.get('cluster_similarity')
            t_gpec_go=golist.cluster(n_CPU=CPU, similarity=cluster_similarity, max_terms=max_nof_enriched_go, l_go_selective=l_go_selective)
            max_clusters=self.get('max_clusters')
            n_CPU=self.get('CPU')
            t_gene_go_mem=golist.membership(max_cluster=max_clusters, n_CPU=n_CPU, l_cluster=True)
            header = util.header(t_gene_go_mem)
            for c in header:
                if c.startswith('GRP'):
                    v = float(c.split('|')[2])
                    r = 0
                    if v<-100:
                        r=0.7
                    elif v<-20:
                        r=0.6
                    elif v<-10:
                        r=0.5
                    elif v<-6:
                        r=0.4
                    elif v<-4:
                        r=0.3
                    elif v<-3:
                        r=0.2
                    elif v<-2:
                        r=0.1
                    t_gene_go_mem[c] = [ r if x>0.5 else 0 for x in t_gene_go_mem[c].tolist()]
            t_gene_go_mem.to_csv(s_out2+'/GeneGo_membership.csv', index=False)
            gpec_golist=golist

            golist.summary_plot(s_out3+'/HeatmapSelectedGO', t_go_mem, max_clusters=max_clusters)
            max_member_per_go_cluster=self.get('max_member_per_go_cluster')
            max_nodes_in_go_network=self.get('max_nodes_in_go_network')
            l_go_selective=self.get('l_go_selective')
            net=golist.network(max_clusters=max_clusters, max_members=max_member_per_go_cluster, max_nodes=max_nodes_in_go_network, l_go_selective=l_go_selective)
            print("Add node count ...")
            t_gene_mem=self.get('t_gene_mem')
            bl.net_annotation(net, t_gene_mem)
            s_xgmml=s_out2+'/GONetwork'
            net.to_xgmml(s_xgmml)
            n_node=net.nof_nodes()
            n_cluster=int(net.T_node.GROUP_ID.max())

        if gpec_golist is not None and l_PLOT:
            self.put('gpec_golist', gpec_golist)
            self.plot_circos_go()
            self.put('gpec_golist', None) # wipe it out, as it contains distance matrix from cluster, increase pickle size considerably
        self.puts({'gpec_gonet_fn':s_xgmml, 'gpec_gonet_nnode':n_node, 'gpec_gonet_ncluster':n_cluster}) #'gpec_golist':gpec_golist})
        if l_WEB_MODE:
            #self.update_session('genelists_overlap.pickle')
            self.save_session('gpec_go.pickle', ['R_idx','R_evi','gpec_gonet_fn','gpec_gonet_nnode','gpec_gonet_ncluster']) #'gpec_golist'])
        if t_gpec_go is not None and len(t_gpec_go):
            t_gpec_go.to_csv(s_out2+'/_FINAL_GO.csv', index=False)

    def GPEC_ppi_analysis(self):
        #if l_WEB_MODE:
        #    if self.session['lists'] is None: self.load_session('genelists_overlap.pickle')
        #    if self.session['t_evi_ppi'] is None: self.load_session('ppi.pickle')
        #    if self.session['t_evi'] is None: self.load_session('evidence.pickle')

        def find_net(S_hits, c_evi_score, target_size=200, max_list_size=1000):
            """Search within R_idx, evidence_score cutoff, find the network that is closest to target_size"""
            if len(S_hits)==0: return (None, 0)
            myppi=self.get('myppi')
            full_net=myppi.subnetwork(S_hits)
            n_full_nodes=full_net.nof_nodes()
            #print(">>>>>>>>>>", n_full_nodes)
            if n_full_nodes==0: return (None, n_full_nodes)
            full_net.add_a_node_attr('EvidenceScore', c_evi_score, s_NULL=0.0)
            t=full_net.T_edge
            t['EvidenceScore_A']=t.Gene_A.apply(lambda x: c_evi_score.get(x, 0.0))
            t['EvidenceScore_B']=t.Gene_B.apply(lambda x: c_evi_score.get(x, 0.0))
            t['EvidenceScore']=t[['EvidenceScore_A','EvidenceScore_B']].min(axis=1)
            #return (full_net, n_full_nodes)
            t.sort_values('EvidenceScore', ascending=False, inplace=True)
            #t.to_csv('t.csv', index=False)
            #full_net.to_xgmml('full_net')
            n=len(t)
            t.index=range(n)
            iB=iE=0
            I_node=np.zeros(n)
            genes=set()
            for i in range(1, n+1):
                if i>=n or t.EvidenceScore[i]!=t.EvidenceScore[i-1]:
                    iE=i-1
                    genes|=set(t.ix[iB:iE, 'Gene_A'])|set(t.ix[iB:iE, 'Gene_B'])
                    I_node[iB:(iE+1)]=len(genes)
                    iB=i
            print(util.unique(I_node))
            I_node=I_node[::-1]
            idx=n-np.argmin(np.abs(I_node-target_size))
            #print(idx)
            #print(len(S_hits))
            S_hits=list(set(set(t.Gene_A[:idx])|set(t.Gene_B[:idx])))
            #print(len(S_hits), n_full_nodes)
            if len(S_hits)<n_full_nodes:
                opt_net=full_net.subnetwork(S_hits)
            else:
                opt_net=full_net
            if opt_net.nof_nodes()>max_list_size: return (None, n_full_nodes)
            return (opt_net, n_full_nodes)

        l_WEB_MODE=self.get('l_WEB_MODE')
        s_out=self.get('s_out')
        t_evi=self.get('t_evi')
        myppi=self.get('myppi')
        #R_idx=self.get('R_idx')
        R_evi=self.get('R_evi')
        l_GPEC=self.get('l_GPEC')
        n_list=self.get('n_list')
        t_gene_mem=self.get('t_gene_mem').copy()

        s_out2=self.get_subfolder("Enrichment_PPI")
        s_out3=self.get_subfolder("Enrichment_PPI/xgmml")
        s_merge_name="_FINAL"
        gpec_net=gpec_mcode_net=t_gpec_mcode=None
        gpec_ppi_target_size=self.get('gpec_ppi_target_size')
        min_ppi_size=self.get('min_ppi_size')
        max_ppi_size=self.get('max_ppi_size')
        max_list_size_for_ppi=self.get('max_list_size_for_ppi')

        sw2=util.StopWatch("MetaEngine::GPEC_ppi_analysis")
        # counts initialized
        if not l_GPEC and n_list<=1:
            #idx=len(R_idx)-1
            pass
        else:
            c_evi_score=dict(zip(t_evi.Gene, t_evi.TotalEvidenceCounts))
            (gpec_net, n_full_nodes)=find_net(t_evi.Gene.tolist(), c_evi_score, target_size=gpec_ppi_target_size, max_list_size=max_list_size_for_ppi)
            if gpec_net is not None:
                n_nodes=gpec_net.nof_nodes()
                if n_nodes<min_ppi_size or n_nodes>max_ppi_size:
                    util.warn_msg('Merged GPEC network is too small or too large: nodes=%d' % n_nodes)
                    gpec_net=None

        #print idx, gpec_net.nof_nodes(), n_full_nodes
        # single list, no need to generate network if GPEC net is the same as the full network
        if n_list<=1 and gpec_net is not None and gpec_net.nof_nodes()==n_full_nodes:
            util.warn_msg('GPEC PPI is skipped, as it is the same as the original network!')
            gpec_net=None
        else:
            if gpec_net is not None:
                n_nodes=gpec_net.nof_nodes()
                print("Use GPEC Network> GPEC size=%d, MERGE size=%d" % (n_nodes, n_full_nodes))
                gpec_net.name=s_merge_name
                glist=bl.GeneList('_FINAL', gpec_net.nodes())
                l_MCODE_optimize=self.get('l_MCODE_optimize')
                CPU=self.get('CPU')
                l_connect_in_merge=self.get('l_connect_in_merge')
                (gpec_mcode_net, t_gpec_mcode)=bl.GeneList.mcode_analysis(gpec_net, n_CPU=CPU, l_MCODE_optimize=False, min_size=min_ppi_size, max_size=max_ppi_size, l_connect_in_merge=l_connect_in_merge)

                S=[x for x in t_gene_mem.header() if x.startswith('_MEMBER_') ]
                nets=[gpec_net]
                if gpec_mcode_net is not None:
                    nets.append(gpec_mcode_net)
                for i,x in enumerate(S):
                    c_dict={}
                    # _MEMBER_* columns might have been modified by Circos plot into 0.3333
                    for j in t_gene_mem.index:
                        c_dict[str(t_gene_mem.ix[j,'Gene'])]=1 if t_gene_mem.ix[j, x]>=0.99 else 0
                    for y in nets:
                        y.add_a_node_attr('#COUNT_%03d' % (i+1), c_dict)

        full_nets=self.get('full_nets') or {}
        full_nets=list(full_nets.values())
        if gpec_net is not None and not gpec_net.is_empty():
            full_nets.append(gpec_net)
        mcode_nets=self.get('mcode_nets') or []
        if gpec_mcode_net is not None and not gpec_mcode_net.is_empty():
            mcode_nets.append(gpec_mcode_net)

        c_hit={}
        c_hit={x.name:x.nodes() for x in full_nets}
        for x in mcode_nets:
            c_hit[x.name]=x.nodes()

        t_mcode=pd.DataFrame() if self.get('t_mcode') is None else self.get('t_mcode')
        if t_gpec_mcode is not None and len(t_gpec_mcode)>0:
            t_mcode=pd.concat([t_mcode, t_gpec_mcode], ignore_index=True)
        if len(t_mcode):
            #Plot PPI
            for k,t_v in t_mcode.groupby(['Network','Cluster']):
                c_hit["%s_MCODE_%d" % (k[0], k[1])]=list(t_v.Gene)

        def sort_by_name(t, s_col='Name'):
            S_X=t[s_col].tolist()
            S_Y=[0]*len(t)
            for i,s in enumerate(S_X):
                m=re.search(r'(?P<net>.*)_(?P<id>\d+)$', s)
                if m:
                    #print(t.ix[i, s_col], m.groups())
                    net_name=m.group('net')
                    if re.search(r'_SUB\d+_MCODE$', net_name):
                        net_name=re.sub(r'\d+_MCODE$', '', net_name)
                    S_X[i]=net_name
                    S_Y[i]=int(m.group('id'))
            t['_X']=S_X
            t['_Y']=S_Y
            t['_Z']=range(len(t)) # keep old order
            t.sort_values(by=['_X', '_Y', '_Z'], inplace=True)
            t.drop(['_X','_Y','_Z'], axis=1, inplace=True)

        sw2.check("***Start GO Analysis***")
        t_mcode_go=t_mcode_ann=None
        if len(c_hit):
            mygo=self.get('mygo')
            S_background=self.get('S_background')
            CPU=self.get('CPU')
            S_go_category=self.get('S_go_category')
            max_go_annotation_per_mcode=self.get('max_go_annotation_per_mcode')
            S_go_category=[x for x in S_go_category if x < 90] # remove L1k categories
            t_mcode_go=mygo.analysis(c_hit, SRC_GENE=S_background, n_CPU=CPU, S_go_category=S_go_category)
            sw2.check("Finish GO analysis")
            if t_mcode_go is not None and len(t_mcode_go):
                sort_by_name(t_mcode_go, 'Name')
                util.df2sdf(t_mcode_go).to_csv(s_out2+"/GO_MCODE.csv", index=False)
                data=[]
                S_net=[]
                S_enrh=[]
                for k,t_v in t_mcode_go.groupby('Name'):
                    t_v=t_v[:max_go_annotation_per_mcode].copy()
                    t_v['Description']=t_v['Description'].apply(lambda x: re.sub(';', ' ', x))
                    S_net.append(k)
                    S=t_v.GO+"|"+t_v.Description+"|"+t_v.LogP.apply(lambda x: '%.1f' % x)
                    S_enrh.append(";".join(list(S)))
                t_mcode_ann=pd.DataFrame(data={'Network':S_net, 'Annotation':S_enrh}, columns=["Network", "Annotation"])
                sort_by_name(t_mcode_ann, 'Network')
                t_mcode_ann.to_csv(s_out2+"/GO_MCODE_Top3.csv", index=False)
        sw2.check("*****GO Analysis within PPI Analysis")

        #print R_idx, len(R_idx), R_evi, len(R_evi), idx
        gpec_ppi_evidence_cutoff=gpec_ppi_size=None
        if gpec_net is not None:
            gpec_ppi_evidence_cutoff=gpec_net.T_node.EvidenceScore.min()
            gpec_ppi_size=sum(t_evi.TotalEvidenceCounts>= gpec_ppi_evidence_cutoff)
        self.puts({'gpec_net':gpec_net, 'gpec_mcode_net':gpec_mcode_net, 't_gpec_mcode':t_gpec_mcode, 'GPEC_PPI_Evidence_Cutoff':gpec_ppi_evidence_cutoff, 'GPEC_PPI_Size': gpec_ppi_size, 't_mcode_go':t_mcode_go, 't_mcode_ann':t_mcode_ann})
        if l_WEB_MODE:
            #self.update_session('genelists_overlap.pickle')
            self.save_session('gpec_ppi.pickle', ['gpec_net','gpec_mcode_net','t_gpec_mcode','GPEC_PPI_Evidence_Cutoff','GPEC_PPI_Size','t_mcode_go','t_mcode_ann'])
        if gpec_net is not None:
            gpec_net.to_xgmml(s_out2+'/xgmml/'+gpec_net.name)
        if gpec_mcode_net is not None:
            gpec_mcode_net.to_xgmml(s_out2+'/xgmml/'+gpec_mcode_net.name)
        if t_gpec_mcode is not None and len(t_gpec_mcode):
            t_gpec_mcode.to_csv(s_out2+'/_FINAL_MCODE.csv', index=False)


    def plot_circos_go(self):
        # this method must be executed within GPEC_go_analysis(), when self.get('gpec_golist') is temporarily set for this method
        #if l_WEB_MODE:
        #    if self.session['t_gene_mem'] is None: self.load_session('genelists_overlap.pickle')
        #    if self.session['gpec_golist'] is None: self.load_session('gpec_go.pickle')
        l_WEB_MODE=self.get('l_WEB_MODE')
        CIRCOS_BIN=self.get('CIRCOS_BIN')
        circos_max_go_terms=self.get('circos_max_go_terms')
        l_go_selective=self.get('l_go_selective')
        max_clusters=self.get('max_clusters')
        max_member_per_go_cluster=self.get('max_member_per_go_cluster')

        s_out2=self.get_subfolder("Overlap_circos")
        n_list=self.get('n_list')
        if n_list<2: return
        golist=self.get('gpec_golist')
        if golist is None or len(golist.data)==0: return
        t_gene_mem=self.get('t_gene_mem')
        S_color=self.get('S_color')
        # expanded circos plot
        tmp_go=golist.data.copy()
        circos_link_logp=self.get('circos_link_logp')
        circos_link_go_size=self.get('circos_link_go_size')
        circos_link_enrichment=self.get('circos_link_enrichment')
        tmp_go=tmp_go[(tmp_go.LogP<=circos_link_logp) & (tmp_go['#GeneInGO']<=circos_link_go_size) & (tmp_go.Enrichment>=circos_link_enrichment)]
        if 'GROUP_ID' in tmp_go.header() and circos_max_go_terms>0 and len(tmp_go)>circos_max_go_terms:
            # trim so go terms
            S_go=go.GO_Cluster.sample_rows(tmp_go, max_clusters=max_clusters, max_members=max_member_per_go_cluster, max_nodes=circos_max_go_terms, l_go_selective=l_go_selective)
            tmp_go=tmp_go[tmp_go.GO.apply(lambda x: x in S_go)]
        tmp_go=tmp_go[['GeneID']]
        cir=circos.Circos(BIN=CIRCOS_BIN)
        #t_go=golist.data
        #t_go=t_go[(t_go.LogP<=circos_link_logp) & (t_go['#GeneInGO']<=circos_link_go_size) & (t_go.Enrichment>=circos_link_enrichment)]
        #t_go=t_go[['GeneID']]
        s_symbol=None
        if 'Symbol' in t_gene_mem.header():
            s_symbol='Symbol'
        cir.plot_membership(t_gene_mem, t_go=tmp_go, outputdir=s_out2, outputfile="CircosOverlapByGO", S_color=S_color, s_symbol=s_symbol)
        self.sw.check("Circos GO Overlap")

    def plot_GPEC_go(self):
        #if l_WEB_MODE:
        #    if self.session['t_gene_mem'] is None: self.load_session('genelists_overlap.pickle')
        #    if self.session['golists'] is None: self.load_session('golists.pickle')
        #    if self.session['gpec_golist'] is None: self.load_session('gpec_go.pickle')
        l_WEB_MODE=self.get('l_WEB_MODE')
        s_out=self.get('s_out')
        s_out2=self.get_subfolder("Enrichment_GO")
        #golist=self.get('gpec_golist')
        #if golist is None: return
        s_xgmml=self.get('gpec_gonet_fn')
        if s_xgmml is None: return
        t_go_mem=self.get('t_go_mem')
        #t_gene_mem=self.get('t_gene_mem')
        n_list=self.get('n_list')
        S_name=self.get('S_name')
        S_color=self.get('S_color')
        #max_clusters=self.get('max_clusters')
        #max_member_per_go_cluster=self.get('max_member_per_go_cluster')
        #max_nodes_in_go_network=self.get('max_nodes_in_go_network')
        #l_go_selective=self.get('l_go_selective')
        #golist.summary_plot(s_out+'/HeatmapSelectedGO', t_go_mem, max_clusters=max_clusters)
        self.sw.check('Start Plot GO')
        cp=bl.CytoPlot(host=self.get('CYTOSCAPE_HOST'), port=self.get('CYTOSCAPE_PORT'), on_windows=self.get('CYTOSCAPE_ON_WINDOWS'), add_logo=self.get('ADD_LOGO'))

        n_node=self.get('gpec_gonet_nnode')
        n_cluster=self.get('gpec_gonet_ncluster')
        s_session=os.path.abspath(s_xgmml+".cys")
        if self.get('CYTOSCAPE_ON_WINDOWS'):
            s_session=util.format_path(s_session, OS="win")

        S_idx=bl.get_color_label(t_go_mem, S_name)
        S_color2=[ S_color[i] for i in S_idx if i>=0]
        S_label2=[ S_name[i] for i in S_idx if i>=0]
        #print S_color2, S_label2
        cp.plot_go_network(s_xgmml+".xgmml", n_list, n_cluster, s_out_dir=s_out2, S_fmt=["png","pdf"], s_session=s_session, S_label=S_label2, S_color=S_color2, n_node=n_node)
        print("Save session file: "+s_session)
        self.sw.check("CytoPlot")

    def plot_GPEC_ppi(self):
        #if l_WEB_MODE:
        #    if self.session['t_gene_mem'] is None: self.load_session('genelists_overlap.pickle')
        #    if self.session['full_nets'] is None: self.load_session('ppi.pickle')
        #    if self.session['gpec_net'] is None: self.load_session('gpec_ppi.pickle')
        l_WEB_MODE=self.get('l_WEB_MODE')
        max_ppi_size_for_plot=self.get('max_ppi_size_for_plot')
        s_out=self.get('s_out')
        s_out2=self.get_subfolder('Enrichment_PPI')
        s_out3=self.get_subfolder("Enrichment_PPI/xgmml")
        full_nets=self.get('full_nets') or {}
        full_nets=list(full_nets.values())
        mcode_nets=self.get('mcode_nets') or []
        gpec_net=self.get('gpec_net')
        if gpec_net is not None and not gpec_net.is_empty():
            full_nets.append(gpec_net)
        gpec_mcode_net=self.get('gpec_mcode_net')
        if gpec_mcode_net is not None and not gpec_mcode_net.is_empty():
            mcode_nets.append(gpec_mcode_net)

        if len(full_nets)==0: return

        lists=self.get('lists')
        n_list=self.get('n_list')
        S_name=self.get('S_name')
        S_color=self.get('S_color')

        s_merge_name="_FINAL"

        S_xgmml=[]
        S_node=[]
        c_legend_prefix={}

        cp=bl.CytoPlot(host=self.get('CYTOSCAPE_HOST'), port=self.get('CYTOSCAPE_PORT'), on_windows=self.get('CYTOSCAPE_ON_WINDOWS'), add_logo=self.get('ADD_LOGO'))

        for x in mcode_nets:
            s_name=re.sub('_MCODE_ALL', '', x.name)
            S_xgmml.append(s_out3+"/"+x.name+".xgmml")
            S_node.append(x.nof_nodes())
            n_cluster=int(max(x.T_node.MCODE_CLUSTER_ID))
            c_legend_prefix[x.name]=["MCODE%d" % (i+1) for i in range(n_cluster)]
            c_legend_prefix[s_name]=c_legend_prefix[x.name]
        n_mcode_cluster=1
        max_mcode_cluster=20
        t_mcode=pd.DataFrame() if self.get('t_mcode') is None else self.get('t_mcode')
        t_gpec_mcode=self.get('t_gpec_mcode')
        if t_gpec_mcode is not None and len(t_gpec_mcode)>0:
            t_mcode=pd.concat([t_mcode, t_gpec_mcode], ignore_index=True)
        if len(t_mcode):
            n_mcode_cluster=max(t_mcode.Cluster)

        S_full_network=[]
        S_node_full=[]
        for x in full_nets:
            if x.nof_nodes()<=max_ppi_size_for_plot:
                S_full_network.append(s_out3+"/"+x.name+".xgmml")
                S_node_full.append(x.nof_nodes())
        S_by_count=[]
        # just in case there were bugs
        S_xgmml=util.unique(S_xgmml)
        #S_full_network=util.unique(S_full_network)
        if len(S_xgmml+S_full_network)>0:
            #S_node=[ x.nof_nodes() for x in best_out if x.name.endswith('_MCODE_ALL') ]
            if n_list>1:
                S_by_count=[x for x in S_xgmml if x.startswith(s_out3+"/"+s_merge_name+'_')]
            if len(S_full_network):
                S_xgmml+=S_full_network
                S_node+=S_node_full
                if n_list>1:
                    S_by_count+=[x for x in S_full_network if x==s_out3+"/"+s_merge_name+'.xgmml']
            s_session2=os.path.abspath(s_out2+"/MCODE_PPI.cys")
            if self.get('CYTOSCAPE_ON_WINDOWS'):
                s_session2=util.format_path(s_session2, OS="win")
            self.sw.check('Start Plot PPI')
            #print S_by_count, S_xgmml, n_list, n_mcode_cluster, S_name, S_color, S_node, c_legend_prefix
            cp.plot_ppi_network(S_xgmml, n_list, n_mcode_cluster, s_out_dir=s_out2, S_fmt=['png','pdf'], s_session=s_session2, max_cluster=max_mcode_cluster, S_by_count=S_by_count, S_label=S_name, S_color=S_color, S_node=S_node, s_legend_prefix=c_legend_prefix)
            self.sw.check("Plot PPI")

    def pptx_report(self):
        #if l_WEB_MODE:
        #    if self.session['t_gene_mem'] is None: self.load_session('genelists_overlap.pickle')
        # PPTX report
        l_WEB_MODE=self.get('l_WEB_MODE')
        s_out=self.get('s_out')
        s_out_ppi=s_out+"/Enrichment_PPI" # self.get_subfolder('Enrichment_PPI'), don't create the folder, if no PPI analysis
        n_list=self.get('n_list')
        lists=self.get('lists')
        S_name=self.get('S_name')
        search_words=self.get('search_words')
        S_cnt=[len(x) for k,x in lists]
        t_lists=pd.DataFrame(data={'Name':S_name, 'Total':S_cnt, 'Unique':S_cnt})

        c_table={}
        S_ppi=glob.glob(s_out_ppi+'/*.png')
        S_colorby=[x for x in S_ppi if x.endswith('PPIColorByCluster.png')]
        S_countby=[x for x in S_ppi if x.endswith('PPIColorByCounts.png')]
        if os.path.exists(s_out_ppi+'/GO_MCODE_Top3.csv'):
            c_table['GO_MCODE']=s_out_ppi+'/GO_MCODE_Top3.csv'
        acc=self.get('ML_accuracy')
        l_ML=acc is not None and acc>0
        if n_list>1:
            if len(t_lists)>100:
                t_lists=t_lists[:100]
                util.warn_msg('Too many gene lists, only keep the first 100!')
            c_img={'membership': s_out+'/Membership_PieChart/membership.png', 'circos': s_out+'/Overlap_circos/CircosOverlapByGene.png', 'circos_go': s_out+'/Overlap_circos/CircosOverlapByGO.png', 'GOBargraph': s_out+'/Enrichment_hetamap/HeatmapSelectedGO.png', 'GOColorByCluster':s_out+'/Enrichment_GO/ColorByCluster.png', 'GOColorByPValue':s_out+'/Enrichment_GO/ColorByPValue.png', 'GOColorByCounts': s_out+'/Enrichment_GO/ColorByCounts.png'}
            if len(S_colorby):
                c_img['PPIColorByCluster']=S_colorby
            if len(S_countby):
                c_img['PPIColorByCounts']=S_countby
            if l_ML:
                c_img['EvidenceWeight']=s_out+'/Evidence/EvidenceWeight.png'
            rpt.pptx_multiple(
                s_out+'/AnalysisReport.pptx',
                t_lists,
                c_txt={'SearchTerm': search_words},
                c_img=c_img,
                c_table=c_table
            )
        elif n_list==1:
            c_img={'membership': s_out+'/membership.png', 'GOBargraph': s_out+'/Enrichment_GO/HeatmapSelectedGO.png', 'GOColorByCluster':s_out+'/Enrichment_GO/ColorByCluster.png', 'GOColorByPValue': s_out+'/Enrichment_GO/ColorByPValue.png'}
            if len(S_colorby):
                c_img['PPIColorByCluster']=S_colorby
            if l_ML:
                c_img['EvidenceWeight']=s_out+'/Evidence/EvidenceWeight.png'
            rpt.pptx_single(
                s_out+'/AnalysisReport.pptx',
                c_txt={'SearchTerm': search_words, '#Total':t_lists.Total[0], '#Unique':t_lists.Unique[0]},
                c_img=c_img,
                c_table=c_table
            )

    def plot_go_pie(self):
        s_out=self.get('s_out')
        l_WEB_MODE=self.get('l_WEB_MODE')
        t_go=None
        if os.path.exists(s_out+'/Enrichment_GO/_FINAL_GO.csv'):
            t_go=pd.read_csv(s_out+'/Enrichment_GO/_FINAL_GO.csv')
            t_go=t_go[(t_go.GROUP_ID<=self.get('max_clusters')) & (t_go.FirstInGroupByLogP==1)]
        if t_go is None or len(t_go)==0: return
        t_go.sort_values('GROUP_ID', inplace=True)
        S_go=t_go.GO.tolist()

        lists=self.get('lists')
        s_out2=self.get_subfolder('Enrichment_PieChart')

        S_background=self.get('S_background')

        MAX_LISTS=20 # plot at most 20 lists
        S_color=self.get('S_color')[:MAX_LISTS]
        S_name=self.get('S_name')[:MAX_LISTS]
        c_idx={ x: i for i,x in enumerate(S_name) }

        # this does not work for L1000 categories
        mygo=self.get('mygo')
        l_L1k=False
        if sum(t_go['CategoryID']>90)>0:
            l_L1k=True
            mygo_L1k=self.get('mygo_L1k')
            if S_background is not None:
                S_has_GO_L1k=set(mygo_L1k.ALL_GENE & S_background)
            else:
                S_has_GO_L1k=mygo_L1k.ALL_GENE
                N_L1k=len(S_has_GO_L1k)
        S_has_GO=set(mygo.GENE_GO.keys())
        N=len(S_has_GO)
        c_hit={}
        for s_name in S_name:
            S_hit=lists[s_name].data & S_has_GO
            if S_background is not None:
                S_hit=S_hit & S_background
            c_hit[s_name]=S_hit
        S_n2=[len(c_hit[x]) for x in S_name]
        for i in t_go.index:
            s_go=str(t_go.at[i, 'CategoryID'])+'_'+t_go.at[i, 'GO']
            if t_go.at[i, 'CategoryID']<90:
                S_inGO=set(mygo.GO_GENE[s_go])
                n1=len(S_has_GO & S_inGO)
                N_total=N
            else:
                S_inGO=set(mygo_L1k.GO_GENE[s_go])
                n1=len(S_has_GO_L1k & S_inGO)
                N_total=N_L1k
            S_n=[len(S_inGO & c_hit[x]) for x in S_name]
            bl.PiePlot.plot(N_total, n1, S_n2, S_n, s_out2+'/EnrichmentPie_Cluster%02d' % t_go.at[i, 'GROUP_ID'])
        return



    def reset_analysis(self):
        """Remove all results after overlap analysis, users might have changed S_background, GO categories, etc. So all analyses will be redone"""
        s_out=self.get('s_out')
        l_WEB_MODE=self.get('l_WEB_MODE')

        def remove(s_file):
            if os.path.exists(s_file): os.remove(s_file)

        def remove_pattern(s_pat):
            for x in glob.glob(s_out+'/'+s_pat):
                os.remove(x)

        if l_WEB_MODE:
            s_cache=self.get('s_cache')
            for x in ['golists.pickle', 'ppi.pickle', 'evidence.pickle', 'gpec_go.pickle', 'gpec_ppi.pickle']:
                remove(s_cache+"/"+x)

        import shutil
        for x in ['Overlap_circos','Enrichment_GO', 'Enrichment_PPI', 'Evidence', 'Enrichment_heatmap', 'Enrichment_PieChart', 'icon' ]:
            s_out2=s_out+"/"+x
            if os.path.exists(s_out2): shutil.rmtree(s_out2)

        for x in ['all.zip', 'AnalysisReport.pptx', 'AnalysisReport.html']:
            remove(s_out+"/"+x)

        #for x in ['ColorBy*', 'GONetwork.*', 'HeatmapSelectedGO.*']:
        #    remove_pattern(x)
        self.puts({'t_evi_go':None, 'golists':None, 't_go_mem':None, 't_evi_ppi':None, 't_mode':None, 'full_nets':None, 'mcode_nets':None, 't_evi':None, 'R_idx':None, 'R_evi':None, 'S_evi':None, 'R_GPEC_Weight':None, 'R_evi':None, 'gpec_gonet_fn':None,'gpec_gonet_nnode':None,'gpec_gonet_ncluster':None,'gpec_net':None, 'gpec_mcode_net':None, 't_gpec_mcode':None, 't_mcode_go':None, 't_mcode_ann':None, 'GPEC_PPI_Evidence_Cutoff':None, 'GPEC_PPI_Size':None, 'ML_accuracy':None})

    def analysis_L1000_only(self):
        mygo_L1k=self.get('mygo_L1k')
        s_out=self.get('s_out')
        lists=self.get('lists')
        n_list=self.get('n_list')
        max_list_size_for_go=self.get('max_list_size_for_go')
        CPU=self.get('CPU')
        S_go_category=[91,92,93]
        min_overlap=1
        min_enrichment=1.5
        p_cutoff=0.01
        S_background=None
        min_evi_term_size=self.get('min_evi_term_size')
        max_evi_term_size=self.get('max_evi_term_size')
        print(lists)
        golists=lists.go_analysis(mygo_L1k, n_CPU=CPU, max_list_size=max_list_size_for_go, S_go_category=S_go_category, SRC_GENE=S_background, min_overlap=min_overlap, p_cutoff=p_cutoff, min_enrichment=min_enrichment)
        self.puts({'t_evi_go': None, 'golists':golists, 't_go_mem':None})
        s_out2=self.get_subfolder("Enrichment_GO")
        golists.save(s_out2+'/GO_AllLists.csv')

    def analysis(self):
        #self.get('mygo')
        #self.get('myppi')
        sw=util.StopWatch("MetaEngine::analysis")
        s_out=self.get('s_out')
        l_SIMPLE=self.get('l_SIMPLE')
        l_PLOT=self.get('l_PLOT')
        if l_SIMPLE: l_PLOT=False
        l_REPORT=self.get('l_REPORT')
        l_GO=self.get('l_GO')
        l_PPI=self.get('l_PPI')
        l_GPEC=self.get('l_GPEC')
        lists=self.get('lists')
        n_list=self.get('n_list')
        l_WEB_MODE=self.get('l_WEB_MODE')
        l_EVIDENCE_ONLY=self.get('l_EVIDENCE_ONLY')
        data=[]
        if lists.is_empty():
            util.warn_msg('No gene in input')
            return
        self.mkdir(s_out)
        self.genelists_overlap_analysis()
        dur=sw.check("******* OVERLAP")
        data.append(['Overlap', dur])

        #print(util.memory_usage())
        #mygo=self.get('mygo')
        #print(util.memory_usage())
        #myppi=self.get('myppi')
        #print(util.memory_usage())

        if l_GO:
            self.get('mygo')
            self.go_analysis()
            dur=sw.check("***** GO ANALYSIS")
            data.append(['GO Analysis', dur])

        if l_PPI:
            self.get('myppi')
            self.ppi_analysis()
            dur=sw.check("***** PPI ANALYSIS")
            data.append(['PPI Ananlysis', dur])

        self.GPEC_evidence()
        dur=sw.check("*****  EVIDENCE")
        data.append(['Evidence Matrix', dur])

        if l_EVIDENCE_ONLY:
            return pd.DataFrame(data, columns=['Step', 'Duration'])

        if l_GO:
            self.GPEC_go_analysis()
            dur=sw.check("****** GPEC GO")
            data.append(['GPEC GO', dur])

        if l_PPI:
            self.GPEC_ppi_analysis()
            dur=sw.check("****** GPEC PPI")
            data.append(['GPEC PPI', dur])

        if l_PLOT:
            if l_GO:
                dur=sw.check("****** PLOT CIRCOS")
                data.append(['Circos Plot', dur])

                self.plot_GPEC_go()
                dur=sw.check("****  PLOT GO")
                data.append(['GO Plot', dur])

                if self.get('l_go_piechart'):
                    self.plot_go_pie()
                    dur=sw.check("****  PLOT GO PieCharts")
                    data.append(['GO PiePlot', dur])

            if l_PPI:
                self.plot_GPEC_ppi()
                dur=sw.check("***** PLOT PPI")
                data.append(['PPI Plot', dur])

        if l_REPORT:
            self.pptx_report()

        if not l_WEB_MODE:
            import ms.analysis_report as anarpt
            myez=self.get('myez')
            dbcon=myez.db.con()
            t_go=None
            if os.path.exists(s_out+'/Enrichment_GO/_FINAL_GO.csv'):
                t_go=pd.read_csv(s_out+'/Enrichment_GO/_FINAL_GO.csv')
            t_lists=pd.DataFrame({'Name':lists.names, 'Unique':[len(x) for x in lists.values]})
            img_circos_gene=img_circos_go=None
            if os.path.exists(s_out+'/Overlap_circos/CircosOverlapByGene.png'):
                img_circos_gene='Overlap_circos/CircosOverlapByGene.png'
            if os.path.exists(s_out+'/Overlap_circos/CircosOverlapByGO.png'):
                img_circos_go='Overlap_circos/CircosOverlapByGO.png'
            bm=0
            if self.get('list_benchmark'):
                bm=len(self.get('list_benchmark'))
            words=self.get('search_words')
            img_membership=None
            if words is not None and words!='':
                if os.path.exists(s_out+'/Membership_PieChart/membership.png'):
                    img_membership='Membership_PieChart/membership.png'
            img_go_heatmap=img_go_nw_clstr=img_go_nw_pval=img_go_nw_cnt=""
            if os.path.exists(s_out+'/Enrichment_heatmap/HeatmapSelectedGO.png'):
                img_go_heatmap='Enrichment_heatmap/HeatmapSelectedGO.png'
            if os.path.exists(s_out+'/Enrichment_GO/ColorByCluster.png'):
                img_go_nw_clstr='Enrichment_GO/ColorByCluster.png'
            if os.path.exists(s_out+'/Enrichment_GO/ColorByPValue.png'):
                img_go_nw_pval='Enrichment_GO/ColorByPValue.png'
            if os.path.exists(s_out+'/Enrichment_GO/ColorByCounts.png'):
                img_go_nw_cnt='Enrichment_GO/ColorByCounts.png'
            bg=0
            if self.get('S_background'):
                bg=len(self.get('S_background'))
            ppi_img_fn=[os.path.basename(x) for x in glob.glob(s_out+'/Enrichment_PPI/*.png')]
            ppi_ds=self.get('ppi_datasource')
            if ppi_ds is None:
                if self.get('tax_id')==9606:
                    ppi_ds=["BioGrid","InWeb_IM","OmniPath"]
                else:
                    ppi_ds=['BioGrid']
            myopt={
                'URL_BASE':'.',
                'LISTS':n_list,
                'ENRICHMENT':t_go,
                'DATAFRAME_LISTS':t_lists,
                'IMG_CIRCO_BY_GENE':img_circos_gene,
                'IMG_CIRCO_BY_GO':img_circos_go,
                'MEMBERSHIP':words,
                'IMG_MEMBERSHIP':img_membership,
                'GO_SOURCE':self.get('S_go_category'),
                'MIN_COUNT':self.get('min_overlap'),
                'P_CUTOFF':self.get('p_cutoff'),
                'ENRICH_CUTOFF':self.get('min_enrichment'),
                'IMG_GO_HEATMAP':img_go_heatmap,
                'IMG_GO_NETWORK_CLUSTER':img_go_nw_clstr,
                'IMG_GO_NETWORK_PVALUE':img_go_nw_pval,
                'IMG_GO_NETWORK_COUNT':img_go_nw_cnt,
                'TARGET_SPECIES': ez.Cache.C_TAX_NAME.get(self.get('tax_id'), 'Unkown'),
                'GPEC':l_GPEC,
                'GPEC_BENCHMARK':bm,
                'BACKGROUND_GENE_COUNT':bg,
                'engine':self,
                #'l_WEB_MODE':l_WEB_MODE,
                'PPI':{'disablePPI': not l_PPI, 'NETWORK_COUNT':1, 'DATASOURCE': ppi_ds, 'minSize': self.get('min_ppi_size'), 'maxSize': self.get('max_ppi_size'), 'img_fn':ppi_img_fn, 'GO_MCODE_Top3_fn':s_out+'/Enrichment_PPI/GO_MCODE_Top3.csv', 'MCODE_CSV_fn':s_out+'/Enrichment_PPI/MCODE.csv', '_FINAL_MCODE_CSV_fn':s_out+'/Enrichment_PPI/_FINAL_MCODE.csv'}
            }
            s_rpt, s_rpt_offline=anarpt.report(myopt, dbcon=dbcon)
            if s_rpt_offline is not None and s_rpt_offline!='':
                util.save_string(s_out+'/AnalysisReport.html', s_rpt_offline)
            s_out2=self.get_subfolder('icon')
            import shutil
            for x in ['PDF48.png','SVG48.png','CYS48.png']:
                shutil.copyfile(os.path.join(os.path.dirname(__file__), "report", "icon", x), os.path.join(s_out2, x))
        # produce zip package
        report.make_analysis_report(s_src=s_out, s_out= s_out + '/all.zip')
        dur_sum=sum([ x[1] for x in data])
        data.append(['Total', dur_sum])
        return pd.DataFrame(data, columns=['Step','Duration'])

class MetaAnalysis:

    def __init__(self, s_output, C_session=None, l_WEB_MODE=False):
        self.engine=MetaEngine(s_output, S_session_files=C_session, l_WEB_MODE=l_WEB_MODE)

    def create_gene_lists(self, S_name, S_hits, type='Entrez', l_skip_convert=False, source_tax_id=9606, target_tax_id=9606, S_color=None):
        if type!='Entrez' or not l_skip_convert:
            S_hits=[ self.id_convert(x, type, target_tax_id=target_tax_id, source_tax_id=source_tax_id) for x in S_hits]
        #print S_hits
        list_benchmark=None
        S_background=None
        lists=[]
        for i,x in enumerate(S_name):
            if x=='_BENCHMARK':
                list_benchmark=bl.GeneList(x, S_hits[i])
            elif x=='_BACKGROUND':
                S_background=S_hits[i]
            elif x in ('_MERGE', '_GPEC', '_FINAL'):
                util.error_msg('Name %s is a reserved gene list name, please avoid starting names with underscore!' % x)
            else:
                lists.append(bl.GeneList(x, S_hits[i]))
        self.engine.put('list_benchmark', list_benchmark)
        self.engine.put('S_background', S_background)
        self.engine.put('tax_id', target_tax_id)
        self.engine.reset_gene_lists(bl.GeneLists(lists), S_color, None)

    @staticmethod
    def hits2csv(S_name, S_hits, s_file="hits.csv"):
        if not s_file.lower().endswith('.csv'):
            util.error_msg('file must ends with .csv, otherwise Metascape cannot load!')
        S_hits=[ util.unique(x) for x in S_hits]
        n=int(np.array([len(x) for x in S_hits]).max())
        S2=[]
        for i,s in enumerate(S_name):
            SS=[""]*n
            n2=len(S_hits[i])
            for j in range(n2):
                SS[j]=S_hits[i][j]
            S2.append(SS)
        util.save_string(s_file, ",".join(S_name)+"\n"+"\n".join([",".join(x) for x in zip(*S2)]))

    @staticmethod
    def hits2txt(S_name, S_hits, s_file="hits.txt"):
        if not s_file.lower().endswith('.txt'):
            util.error_msg('file must ends with .txt, otherwise Metascape cannot load!')
        S_hits=[ util.unique(x) for x in S_hits]
        S_hits=[ " ".join(x) for x in S_hits]
        t=pd.DataFrame({'Name':S_name, 'Hits':S_hits})
        S=[ s+"\t"+S_hits[i] for i,s in enumerate(S_name)]
        util.save_list(s_file, S, s_end="\n")

    def id_convert(self, S_hit, type='Entrez', target_tax_id=9606, source_tax_id=9606):
        myez=self.engine.get('myez')
        t=myez.id_conversion(S_hit, s_source_type=type, target_tax_id=target_tax_id, Source_tax_id=[source_tax_id])
        return util.unique(list(t['gene_id']))

    def load_csv(self, s_csv, type="Entrez", l_skip_convert=True, l_one_list=False, source_tax_id=9606, target_tax_id=9606, S_color=None):
        s_file, s_ext=os.path.splitext(s_csv)
        if s_ext.lower()=='.csv':
            ds=pd.read_csv(s_csv)
        elif l_one_list and s_ext.lower()=='.txt':
            ds=pd.read_table(s_csv)
        else:
            util.error_msg("Unsupport file format: %s" % s_csv)
        S_name=ds.header()
        if l_one_list:
            S_name=S_name[:1]
        if type=='Entrez':
            S_lists=[ util.unique([util.r2i2s(y) for y in list(ds[x]) if not pd.isnull(y)]) for x in S_name ]
        else:
            S_lists=[ util.unique([y for y in list(ds[x]) if not pd.isnull(y)]) for x in S_name ]
        self.create_gene_lists(S_name, S_lists, type=type, l_skip_convert=l_skip_convert, source_tax_id=source_tax_id, target_tax_id=target_tax_id, S_color=S_color)

    def load_txt(self, s_txt, type="Entrez", l_skip_convert=True, l_one_list=False, source_tax_id=9606, target_tax_id=9606, S_color=None):
        if l_one_list:
            return self.load_csv(s_txt, type=type, l_skip_convert=l_skip_convert, l_one_list=l_one_list, source_tax_id=source_tax_id, target_tax_id=target_tax_id, S_color=None)
        S=util.read_list(s_txt)
        S_name=[]
        S_lists=[]
        for i,s in enumerate(S):
            if s.startswith('#'): continue
            SS=s.split("\t")
            if i==0 and SS[0].lower()=='name': continue
            S_name.append(SS[0])
            S_lists.append(util.unique(re.split(r'[ ,;]', SS[1])))
            S_lists=[ [y for y in x if y!='' ] for x in S_lists]
        self.create_gene_lists(S_name, S_lists, type=type, l_skip_convert=l_skip_convert, source_tax_id=source_tax_id, target_tax_id=target_tax_id)

    def settings_from_file(self, s_file, l_one_list, type='Entrez', l_skip_convert=True, source_tax_id=9606, target_tax_id=9606, S_color=None):
        if s_file.lower().endswith('.csv'):
            self.load_csv(s_file, l_one_list=l_one_list, type=type, l_skip_convert=l_skip_convert, source_tax_id=source_tax_id, target_tax_id=target_tax_id, S_color=S_color)
        else:
            self.load_txt(s_file, l_one_list=l_one_list, type=type, l_skip_convert=l_skip_convert, source_tax_id=source_tax_id, target_tax_id=target_tax_id, S_color=S_color)


if __name__=="__main__":
    import argparse as arg
    opt=arg.ArgumentParser(description='Metascape Analysis of Gene Lists')
    #opt.add_argument('-u','--ub', type=float, default=1, help='upper bound, defaults to 1')
    opt.add_argument('-o', '--output', help='output file directory')
    opt.add_argument('-t','--type', type=str, default='Entrez', help='source id type. Entrez,RefSeq,Symbol,dbxref')
    opt.add_argument('-s','--skip_convert', default=False, action='store_true', help='Skip ID conversion, if we are sure ID are updated.')
    opt.add_argument('-u','--one_list', default=False, action='store_true', help='Only one gene list in input.')
    #opt.add_argument('-b','--background', type=str, default='', help='a text file containing a list of background Gene IDs, Entrez ID format, space/comma/semicolon separated')
    opt.add_argument('-c','--CPU', type=int, default=12, help='number of CPUs to be used')
    opt.add_argument('-r','--color', type=str, default=None, help='colors used for lists, hex comma separated string')
    opt.add_argument('-w','--weight', type=str, default=None, help='weight used for lists, hex comma separated string, (+,0,-)')
    opt.add_argument('-p','--PPI', default=False, action='store_true', help='Perform PPI analysis')
    opt.add_argument('-G','--skip_go', default=False, action='store_true', help='Skip GO analysis')
    opt.add_argument('-S','--source_tax_id', type=int, default=9606, help='Source tax_id')
    opt.add_argument('-T','--target_tax_id', type=int, default=9606, help='Target tax_id')
    opt.add_argument('-l','--circos_symbol', default=False, action='store_true', help='Target tax_id')
    opt.add_argument('-E','--skip_GPEC', default=False, action='store_true', help='Skip GPEC prioritization')
    opt.add_argument('-W','--web_mode', default=False, action='store_true', help='Mimic Web Mode, use Cache')

    opt.add_argument('input', nargs=1, help='input file must a in CSV format')

    S_known_tax_id=set(ez.Cache.C_TAX_ID.values())

    args=opt.parse_args()
    if not args.input:
        util.error_msg('Input .csv file is missing!')
    else:
        if not os.path.exists(args.input[0]):
            util.error_msg('File not found: %s' % args.input[0])
    if not args.output:
        util.error_msg('Output folder name is missing, use -o')
    s_out=args.output

    if args.source_tax_id not in S_known_tax_id:
        util.error_msg('Source tax_id: %d is not supported!' % args.source_tax_id)
    if args.target_tax_id not in S_known_tax_id:
        util.error_msg('Target tax_id: %d is not supported!' % args.target_tax_id)
    if args.source_tax_id!=args.target_tax_id and args.skip_convert:
        util.warn_msg('ID conversion is required when source and target tax_id are differnet!')
        args.skip_convert=False

    S_color=args.color
    if S_color is not None:
        S_color=S_color.split(",")
        if len(S_color)==0: S_color=None
    R_weight=args.weight
    if R_weight is not None:
        #c_sign={'+':1, '-':-1, '0':0}
        R_weight=np.array([float(x) for x in R_weight.split(",")])
        if np.isnan(R_weight).any():
            util.error_msg('Bad symbols in weight: only take floating numbers') #+/-/0')
    #S_background=None
    #if args.background!='':
    #    s=util.read_string(args.background)
    #    S_background=re.split(r'[ ;,]', s)
    #    S_background=set([x for x in S_background if x!='' and re.search(r'^\d+$', x)])
    l_GO=not args.skip_go
    l_PPI=args.PPI
    l_GPEC=not args.skip_GPEC

    c_settings={'n_CPU':args.CPU, 'S_color':S_color, 'l_PPI':l_PPI, 'l_GO':l_GO, 'l_GPEC':l_GPEC, 'tax_id':args.target_tax_id, 'circos_show_symbol':args.circos_symbol}
    x=MetaAnalysis(s_output=s_out, C_session=c_settings, l_WEB_MODE=args.web_mode)
    x.settings_from_file(args.input[0], l_one_list=args.one_list, type=args.type, l_skip_convert=args.skip_convert, source_tax_id=args.source_tax_id, target_tax_id=args.target_tax_id, S_color=S_color)
    engine=x.engine
    #engine.puts({S_background':S_background, 'S_color':S_color, 'R_weight':R_weight})
    engine.puts({'S_color':S_color, 'R_weight':R_weight})
    engine.is_okay() # check if S_color and R_weight match the lenght of gene lists
    #engine.old_analysis()
    engine.analysis()

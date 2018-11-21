#!/usr/bin/env python
from os import sys, path
sys.path.insert(0, path.join(path.dirname(path.abspath(__file__)),'../'))
import pandas as pd
import util
import db
from core import *
import urllib
from IPython.core.debugger import Tracer
import re


class GeneGo:
    def __init__(self, taxidList):
        self.con=db.get_con('GENEGO')
        self.taxidList=taxidList


    def get_interaction(self, l_filter=True, l_physical=True):
        ### dump interactome
        sw=util.StopWatch()
        # organism column here is obsolete, should not be used, see Reference 3.1 Species information and interactions, consider they are generic network objects
        s_sql_regulation_rels = 'select distinct id1, id2, type as effect, mechanism, trust, link_id,0 org_link from regulation_rels where nvl(trust, -2) <> -1' # -1, not exist
        # calculated interactions, group and complex relationships
        s_sql_reg_r = 'select distinct id1, id2, 0 as type, mechanism, null as trust, null as link_id,org as org_link from reg_r where mechanism in (10,14)' # 10: Group relation, 14: complex subunit
        s_sql_edge = "select * from (%s union all %s) r where id1!=id2" % (s_sql_regulation_rels, s_sql_reg_r)

        t = db.from_sql(self.con, s_sql_edge)
        sw.check('Interaction data loaded')
        if l_filter:
            # 1: NLP, -1: No Link
            t=t[t.TRUST.apply(lambda x: x not in (-1, 1))].copy()
            #ID Value Meaning Level
            #0 Present Interaction is proven by trusted methods on this organism High
            #8 Approved Interaction is proven for all protein group members (with Present trust) High
            #9 Conflicting data Proven interaction, but different effects in different papers High
            #3 Animal model Proven on animal model High
            #7 Possible common Proven for some protein group members, but not all Medium
            #6 Mix Proven for the protein group as a whole, but not for individual members Medium
            #2 Domain interaction Interaction derived using unreliable methods (yeast2hybrid), only binding site for trans. Factors Low
            #10 Signaling pathway Interaction is made specially for signaling pathway map, may be indirect Low
            #1 NLP Result of data mining, or paper with high-throughput screen (chip on chip, prediction) Low
            #-1 No link Means that this interaction is absent for the particular species No link
            sw.check('Weak link filtered')

        t_tax=db.from_sql(self.con, 'select orgid,taxonomyid from orgs where taxonomyid is not NULL')
        c_tax={t_tax.ix[i,'ORGID']:t_tax.ix[i,'TAXONOMYID'] for i in t_tax.index}

        # filter out undesirable mechanisms
        t_m=self.get_mechanism(l_filter=l_filter)
        t=t.merge(t_m, left_on='MECHANISM', right_on='ID')
        sw.check('Extract Mechanism')
        t_e=self.get_effect()
        t=t.merge(t_e, left_on='EFFECT', right_on='ID')
        sw.check('Extract Effect')
        t_no=self.get_network_object()
        print ">> t", len(t)
        print ">> t_no", util.unique_count(t_no.ORG)
        t=t.merge(t_no, left_on='ID1', right_on='NETW_OBJ_ID')
        t.rename2({'GENE_ID':'GENE_A','ORG':'ORG_A'})
        t=t.merge(t_no, left_on='ID2', right_on='NETW_OBJ_ID')
        t.rename2({'GENE_ID':'GENE_B','ORG':'ORG_B'})
        print ">> A", len(t)
        t=t[(t.GENE_A!=t.GENE_B)&(t.ORG_A==t.ORG_B)]
        print ">> B", len(t)
        t1=t[t.ORG_LINK==0].copy() #LINK_ID is not NULL
        print ">> t1", len(t1)
        t2=t[t.ORG_LINK!=0].copy() #LINK_ID is NULL
        print ">> t2", len(t2)
        sw.check('Add Entrez Gene ID')
        t_p=self.get_pubmed()
        t_p.rename2({'ORG':'ORG_PUBMED'})
        sw.check('Extract PubMed, merging ...')
        t1=t1.merge(t_p, left_on=['LINK_ID'], right_on=['LINK_ID'], how='left')
        t1_1=t1[t1.ORG_PUBMED.isnull()].copy()
        t1_2=t1[~t1.ORG_PUBMED.isnull()].copy()
        t1_2=t1_2[t1_2.ORG_A==t1_2.ORG_PUBMED].copy()
        print ">> t1+pubmed", len(t1_1), len(t1_2)
        t2=t2[t2.ORG_A==t2.ORG_LINK].copy()
        print ">> t2, ORG_LINK", len(t2)
        t=pd.concat([t1_1,t1_2,t2], ignore_index=True)
        t['TRUST']=t['TRUST'].fillna(-2)
        t['TRUST']=t.TRUST.astype(int)
        t=t[['GENE_A','GENE_B','EFFECT_NAME','MECHANISM_NAME','TRUST','PUBMED','ORG_A']] #,'ORG_B','ORG_PUBMED','ORG_LINK']]
        t.rename2({'ORG_A':'ORG'})
        t['tax_id_A']=t.ORG.apply(lambda x: c_tax.get(x, 0))
        t['tax_id_B'] = t['tax_id_A']
        t=t[['GENE_A','GENE_B', 'tax_id_A', 'tax_id_B', 'EFFECT_NAME','MECHANISM_NAME','TRUST','PUBMED']]
        t = t.query('tax_id_A in [%s]' % ','.join(self.taxidList))
        print "DONE", len(t), util.unique_count(t.tax_id_A)
        return t

    def get_mechanism(self, l_filter=True):
        s_sql = "select id, abbr as mechanism_label, tmp as mechanism_name, direct from regulation_mechanisms order by tmp"
        t=db.from_sql(self.con, s_sql)
        if l_filter:
            t=t[ t.DIRECT>0 ]
            t=t[ t.ID.apply(lambda x: x not in (10,14,31,11)) ].copy()
        return t
        #    ID MECHANISM_LABEL                  MECHANISM_NAME  DIRECT
        #0   35              CM                ADP-ribosylation     100
        #1   33              CM                     Acetylation     100
        #2    5               B                         Binding     100
        #3   15               Z                       Catalysis     100
        #4    8               C                        Cleavage     100
        #5   32               C              Cleavage catalysis     100
        #6    6              Cn                     Competition       0
        #7   23               B               Complex formation     100
        #8   11              CR          Complex+Group relation     100
        #9    2              CM           Covalent modification     100
        #10  34              CM                   Deacetylation     100
        #11  39              CM                   Demethylation     100
        #12  30              CM                   Deneddylation     100
        #13   4              -P               Dephosphorylation     100
        #14  28              CM                   Desumoylation     100
        #15  26              CM                Deubiquitination     100
        #16  42              CM                      GPI-anchor     100
        #17  36              CM                   Glycosylation     100
        #18  10              GR                  Group relation     100
        #19  37              CM                   Hydroxylation     100
        #20  12              IE         Influence on expression       0
        #21  38              CM                     Methylation     100
        #22  29              CM                     Neddylation     100
        #23  41              CM                       Oxidation     100
        #24  21              PE          Pharmacological effect       0
        #25   3              +P                 Phosphorylation     100
        #26  44              Rg                      Regulation     100
        #27  40              CM                 S-nitrosylation     100
        #28  31              SR             Similarity relation     100
        #29  43              CM                       Sulfation     100
        #30  27              CM                     Sumoylation     100
        #31  -1            None                       Technical     100
        #32  22              TE                    Toxic effect       0
        #33   9              TR        Transcription regulation      10
        #34   7               T                  Transformation     100
        #35  16              Tn                       Transport     100
        #36  25              CM                  Ubiquitination     100
        #37   0               ?                     Unspecified       0
        #38  20             cRT  co-regulation of transcription     100
        #39  14              CS                 complex subunit     100
        #40  24               M                   miRNA binding     100
        #41  17               B                receptor binding     100
        #42  18              Tn             transport catalysis     100

    def get_effect(self):
        s_sql = "select id, desc_ as effect_name from regulation_types"
        t=db.from_sql(self.con, s_sql)
        return t

    def get_pubmed(self):
        # stragg probably is too slow
        #s_sql="select link_id,stragg(distinct rpo.pubmed) pubmeds from regulation_rels_org rro, regulation_pubmed_org rpo where rro.koid=rpo.koid and rro.org=1 group by link_id"
        s_sql="select link_id,rpo.pubmed,org from regulation_rels_org rro join regulation_pubmed_org rpo on rro.koid=rpo.koid order by link_id,org"
        t=db.from_sql(self.con, s_sql)
        t['PUBMED']=t.PUBMED.astype(str)
        n=len(t)
        iB=iE=0
        I=[]
        S=[]
        O=[]
        for i in xrange(1,n+1):
            #if (i+1)%10000==0:
            #    print "> %d of %d" % (i+1, n)
            if i==n or t.ix[i,'LINK_ID']!=t.ix[i-1,'LINK_ID'] or t.ix[i,'ORG']!=t.ix[i-1,'ORG']:
                iE=i-1
                I.append(t.ix[iB,'LINK_ID'])
                O.append(t.ix[iB,'ORG'])
                if iE>iB:
                    S.append("|".join(set(t.ix[iB:iE, 'PUBMED'])))
                else:
                    S.append(t.ix[iB, 'PUBMED'])
                iB=i
        t=pd.DataFrame(data={'LINK_ID':I, 'ORG':O, 'PUBMED':S})
        return t

    def get_network_object(self):
        s_sql="select g.id netw_obj_id,d.ref gene_id,o.org from gene_netw g,genedbs d,geneorgs o,genes gs where g.gene=d.gene and o.gene=g.gene and g.gene=gs.geneid and gs.type=1 and d.db=17"
        t=db.from_sql(self.con, s_sql)
        return t




'''

    def get_interaction(self, l_filter=True, l_physical=True):
        ### dump interactome
        sw=util.StopWatch()
        s_sql_regulation_rels = 'select distinct id1, id2, type as effect, mechanism, trust, link_id from regulation_rels where nvl(trust, -2) <> -1'
        s_sql_reg_r = 'select distinct id1, id2, 0 as effect, mechanism, null as trust, null as link_id from reg_r where mechanism in (10,14)'

        s_sql_edge = "select * from (%s union all %s) r where id1!=id2" % (s_sql_regulation_rels, s_sql_reg_r)
        t = db.from_sql(self.con, s_sql_edge)#, l_select=True)
        sw.check('Interaction data loaded')
        if l_filter:
            # 1: NLP, -1: No Link
            t=t[t.TRUST.apply(lambda x: x not in (-1, 1))].copy()
            #ID Value Meaning Level
            #0 Present Interaction is proven by trusted methods on this organism High
            #8 Approved Interaction is proven for all protein group members (with Present trust) High
            #9 Conflicting data Proven interaction, but different effects in different papers High
            #3 Animal model Proven on animal model High
            #7 Possible common Proven for some protein group members, but not all Medium
            #6 Mix Proven for the protein group as a whole, but not for individual members Medium
            #2 Domain interaction Interaction derived using unreliable methods (yeast2hybrid), only binding site for trans. Factors Low
            #10 Signaling pathway Interaction is made specially for signaling pathway map, may be indirect Low
            #1 NLP Result of data mining, or paper with high-throughput screen (chip on chip, prediction) Low
            #-1 No link Means that this interaction is absent for the particular species No link
            sw.check('Weak link filtered')

        # filter out undesirable mechanisms
        t_m=self.get_mechanism(l_filter=l_filter)
        t=t.merge(t_m, left_on='MECHANISM', right_on='ID')
        sw.check('Extract Mechanism')
        t_e=self.get_effect()
        t=t.merge(t_e, left_on='EFFECT', right_on='ID')
        sw.check('Extract Effect')
        t_no=self.get_network_object()
        #Tracer()()
        t=t.merge(t_no, left_on='ID1', right_on='NETW_OBJ_ID')
        t.rename2({'GENE_ID':'GENE_A'})
        t.rename2({'TAX_ID': 'tax_id_A'})
        t=t.merge(t_no, left_on='ID2', right_on='NETW_OBJ_ID')
        t.rename2({'GENE_ID':'GENE_B'})
        t.rename2({'TAX_ID': 'tax_id_B'})
        t = t[(t.GENE_A!=t.GENE_B) & (t.tax_id_A == t.tax_id_B)]
        t = t.query('tax_id_A in [%s]'%','.join(self.taxidList))
        sw.check('Add Entrez Gene ID')
        t_p=self.get_pubmed()
        sw.check('Extract PubMed, merging ...')
        t=t.merge(t_p, left_on='LINK_ID', right_on='LINK_ID', how='left')
        #t['TRUST']=t.TRUST.astype(int)
        t=t[['GENE_A','GENE_B', 'tax_id_A', 'tax_id_B', 'EFFECT_NAME','MECHANISM_NAME','TRUST','PUBMED']]
        return t

    def get_mechanism(self, l_filter=True):
        s_sql = "select id, abbr as mechanism_label, tmp as mechanism_name, direct from regulation_mechanisms order by tmp"
        t=db.from_sql(self.con, s_sql)
        if l_filter:
            t=t[ t.DIRECT>0 ]
            t=t[ t.ID.apply(lambda x: x not in (10,14,31,11)) ].copy()
        return t
        #    ID MECHANISM_LABEL                  MECHANISM_NAME  DIRECT
        #0   35              CM                ADP-ribosylation     100
        #1   33              CM                     Acetylation     100
        #2    5               B                         Binding     100
        #3   15               Z                       Catalysis     100
        #4    8               C                        Cleavage     100
        #5   32               C              Cleavage catalysis     100
        #6    6              Cn                     Competition       0
        #7   23               B               Complex formation     100
        #8   11              CR          Complex+Group relation     100
        #9    2              CM           Covalent modification     100
        #10  34              CM                   Deacetylation     100
        #11  39              CM                   Demethylation     100
        #12  30              CM                   Deneddylation     100
        #13   4              -P               Dephosphorylation     100
        #14  28              CM                   Desumoylation     100
        #15  26              CM                Deubiquitination     100
        #16  42              CM                      GPI-anchor     100
        #17  36              CM                   Glycosylation     100
        #18  10              GR                  Group relation     100
        #19  37              CM                   Hydroxylation     100
        #20  12              IE         Influence on expression       0
        #21  38              CM                     Methylation     100
        #22  29              CM                     Neddylation     100
        #23  41              CM                       Oxidation     100
        #24  21              PE          Pharmacological effect       0
        #25   3              +P                 Phosphorylation     100
        #26  44              Rg                      Regulation     100
        #27  40              CM                 S-nitrosylation     100
        #28  31              SR             Similarity relation     100
        #29  43              CM                       Sulfation     100
        #30  27              CM                     Sumoylation     100
        #31  -1            None                       Technical     100
        #32  22              TE                    Toxic effect       0
        #33   9              TR        Transcription regulation      10
        #34   7               T                  Transformation     100
        #35  16              Tn                       Transport     100
        #36  25              CM                  Ubiquitination     100
        #37   0               ?                     Unspecified       0
        #38  20             cRT  co-regulation of transcription     100
        #39  14              CS                 complex subunit     100
        #40  24               M                   miRNA binding     100
        #41  17               B                receptor binding     100
        #42  18              Tn             transport catalysis     100

    def get_effect(self):
        s_sql = "select id, desc_ as effect_name from regulation_types"
        t=db.from_sql(self.con, s_sql)
        return t

    def get_pubmed(self):
        # stragg probably is too slow
        #s_sql="select link_id,stragg(distinct rpo.pubmed) pubmeds from regulation_rels_org rro, regulation_pubmed_org rpo where rro.koid=rpo.koid and rro.org=1 group by link_id"
        s_sql="select link_id,rpo.pubmed from regulation_rels_org rro join regulation_pubmed_org rpo on rro.koid=rpo.koid order by link_id"

        t=db.from_sql(self.con, s_sql)
        t['PUBMED']=t.PUBMED.astype(str)
        n=len(t)
        iB=iE=0
        I=[]
        S=[]
        for i in xrange(1,n+1):
            if (i+1)%10000==0:
                print "> %d of %d" % (i+1, n)
            if i==n or t.ix[i,'LINK_ID']!=t.ix[i-1,'LINK_ID']:
                iE=i-1
                I.append(t.ix[iB,'LINK_ID'])
                if iE>iB:
                    S.append("|".join(set(t.ix[iB:iE, 'PUBMED'])))
                else:
                    S.append(t.ix[iB, 'PUBMED'])
                iB=i
        t=pd.DataFrame(data={'LINK_ID':I, 'PUBMED':S})
        return t

    def get_network_object(self):
        s_sql="select g.id netw_obj_id,d.ref gene_id, orgs.taxonomyid tax_id from gene_netw g,genedbs d,genes gs, orgs, geneorgs o where g.gene=d.gene and o.gene=g.gene and g.gene=gs.geneid and gs.type=1 and d.db=17 and orgs.orgid=o.org "
        t=db.from_sql(self.con, s_sql)
        return t
'''
# class BioGrid:
#     def __init__(self, taxidList):
#         self.taxidList=taxidList
#
#     def get_data (self):
#         if not os.path.exists(os.path.join(SyncDB.DOWNLOAD_DIR(),"BIOGRID-ALL-3.4.134.tab2.txt")):
#             urllib.urlretrieve("http://thebiogrid.org/downloads/archives/Release%20Archive/BIOGRID-3.4.134/BIOGRID-ALL-3.4.134.tab2.zip", os.path.join(SyncDB.DOWNLOAD_DIR(),"BIOGRID-ALL-3.4.134.tab2.zip"))
#             cmd = "unzip " + os.path.join(SyncDB.DOWNLOAD_DIR(),"BIOGRID-ALL-3.4.134.tab2.zip") + " -d " + SyncDB.DOWNLOAD_DIR();
#
#             print cmd;
#             util.unix(cmd);
#
#         t=pd.read_table(os.path.join(SyncDB.DOWNLOAD_DIR(),"BIOGRID-ALL-3.4.134.tab2.txt") , dtype=str)
#         #print t.header()
#         #['#BioGRID Interaction ID', 'Entrez Gene Interactor A', 'Entrez Gene Interactor B', 'BioGRID ID Interactor A', 'BioGRID ID Interactor B', 'Systematic Name Interactor A', 'Systematic Name Interactor B', 'Official Symbol Interactor A', 'Official Symbol Interactor B', 'Synonyms Interactor A', 'Synonyms Interactor B', 'Experimental System', 'Experimental System Type', 'Author', 'Pubmed ID', 'Organism Interactor A', 'Organism Interactor B', 'Throughput', 'Score', 'Modification', 'Phenotypes', 'Qualifications', 'Tags', 'Source Database']
#         #print util.unique(t['Experimental System Type'])
#         t=t[['Entrez Gene Interactor A', 'Entrez Gene Interactor B', 'Organism Interactor A', 'Organism Interactor B', 'Experimental System', 'Experimental System Type', 'Score', 'Pubmed ID']].copy()
#         t.rename2({'Entrez Gene Interactor A':'gid_A', 'Entrez Gene Interactor B':'gid_B', 'Experimental System Type':'interaction_category', 'Experimental System':'interaction_type',  'Pubmed ID':'support','Organism Interactor A':'tax_id_A', 'Organism Interactor B':'tax_id_B', 'Score':'score' })
#         if len(self.taxidList) !=0:
#             filter = ','.join(['"%s"'%x for x in self.taxidList])
#             t = t.query('tax_id_A in [%s] and tax_id_B in [%s] and tax_id_B==tax_id_A'%(filter, filter))
#
#         t.drop_duplicates(['gid_A', 'gid_B'], inplace=True)
#
#         return t;

class PPI(XmlClass):

    def __init__(self, xe=None):
        self.c_gid=None
        self.c_tax=None
        XmlClass.__init__(self,xe=xe)
        self.tag = "PPI"
        self.fn_dest=os.path.join(SyncDB.DOWNLOAD_DIR(),xe.attrib['dest'])
        self.taxidList = SyncDB.SUPPORTED_SPECIES;
        for child in self.children:
            if type(child).__name__ == "Species":
                self.taxidList = child.supported_species
                break;
        self.inputs = [os.path.join(SyncDB.UPLOAD_DIR(),"gid2source_id.csv")]


    def uniprot(self):
        if self.c_gid is None or self.c_tax is None:
            fn = os.path.join(SyncDB.UPLOAD_DIR(),"gid2source_id.csv")
            t_uniprot = pd.read_csv(fn, header=None,names=['gid','id_type_id','source_id', 'col1','col2','tax_id'])
            t_uniprot = t_uniprot[t_uniprot['id_type_id']==7]
            print(len(t_uniprot))
            print(t_uniprot[:5])
            self.c_gid=dict(t_uniprot[['source_id','gid']].values)
            self.c_tax=dict(t_uniprot[['source_id','tax_id']].values)
        return (self.c_gid, self.c_tax)

    @staticmethod
    def df2dict(t_ppi):
        c_ppi={}
        for i in t_ppi.index:
            #if i%1000==0: print i
            g1=t_ppi.ix[i,'gid_A']
            g2=t_ppi.ix[i,'gid_B']
            if g1 not in c_ppi: c_ppi[g1]={}
            if g2 not in c_ppi: c_ppi[g2]={}
            c_ppi[g1][g2]=1
            c_ppi[g2][g1]=1
        return c_ppi

    @staticmethod
    def overlap(c1, c2):
        N=n_a=n_b=n_ab=0
        for k,v in c1.items():
            n_a+=len(v)
            if k not in c2: continue
            for k2 in v.keys():
                if k2 in c2[k]: n_ab+=1
        N=n_a
        for k,v in c2.items():
            n_b+=len(v)
            if k not in c1:
                N+=len(v)
                #print k
                continue
            for k2 in v.keys():
                if k2 not in c1[k]: N+=1
        return (N/2, n_a/2, n_b/2, n_ab/2, n_ab*1.0/N)

    @staticmethod
    def share(c1, c2):
        c={}
        for k,v in c1.items():
            if k not in c2: continue
            for k2 in v.keys():
                if k2 not in c2[k]: continue
                if k not in c: c[k]={}
                c[k][k2]=1
        return c

    def biogrid(self, l_human_only=False):
        fn_source = os.path.join(SyncDB.DOWNLOAD_DIR(),"BIOGRID-ALL-3.4.134.tab2.txt")
        if not os.path.exists(fn_source):
            urllib.urlretrieve("http://thebiogrid.org/downloads/archives/Release%20Archive/BIOGRID-3.4.134/BIOGRID-ALL-3.4.134.tab2.zip",
                               os.path.join(SyncDB.DOWNLOAD_DIR(),"BIOGRID-ALL-3.4.134.tab2.zip"))
            cmd = "unzip " + os.path.join(SyncDB.DOWNLOAD_DIR(),"BIOGRID-ALL-3.4.134.tab2.zip") + " -d " + SyncDB.DOWNLOAD_DIR();
            print cmd;
            util.unix(cmd);

        t=pd.read_table(fn_source , dtype=str)
        #print t.header()
        #['#BioGRID Interaction ID', 'Entrez Gene Interactor A', 'Entrez Gene Interactor B', 'BioGRID ID Interactor A', 'BioGRID ID Interactor B', 'Systematic Name Interactor A', 'Systematic Name Interactor B', 'Official Symbol Interactor A', 'Official Symbol Interactor B', 'Synonyms Interactor A', 'Synonyms Interactor B', 'Experimental System', 'Experimental System Type', 'Author', 'Pubmed ID', 'Organism Interactor A', 'Organism Interactor B', 'Throughput', 'Score', 'Modification', 'Phenotypes', 'Qualifications', 'Tags', 'Source Database']
        print util.unique(t['Experimental System Type'])
        #t=t[(t['Organism Interactor A']=='9606') & (t['Organism Interactor B']=='9606')]
        t.rename2({'Entrez Gene Interactor A':'gid_A', 'Entrez Gene Interactor B':'gid_B', 'Experimental System Type':'interaction_category', 'Experimental System':'interaction_type',  'Pubmed ID':'support', 'Source Database':'ds', 'Organism Interactor A':'tax_id_A', 'Organism Interactor B':'tax_id_B', 'Score':'score'})
        #print t.header()
        t['interaction_type_id']=2
        t=t[['gid_A','gid_B','tax_id_A','tax_id_B','interaction_type_id','interaction_category','interaction_type','score','support','ds']]
        t=t[(t.gid_A!='-') & (t.gid_B!='-')]
        t['gid_A']=t.gid_A.astype(int)
        t['gid_B']=t.gid_B.astype(int)
        t['tax_id_A']=t.tax_id_A.astype(int)
        t['tax_id_B']=t.tax_id_B.astype(int)
        t['ds']='BioGrid'
        t=t[(t.gid_A!=t.gid_B) & (t.tax_id_A==t.tax_id_B) & (t.gid_A>0) & (t.gid_B>0)].copy()
        self.bio=t
        return self.bio


    def InWebIM(self):
        t=pd.read_table(os.path.join(SyncDB.DOWNLOAD_DIR(), 'core.psimitab'),
                        names=['source','target','source_alt','target_alt','source_alias','target_alias','interaction_method','first_author','publication','tax_id_A','tax_id_B','interaction_type','ds','interaction_id','score','other']) #, nrows=10)
        t['source']=t.source.apply(lambda x: x.replace('uniprotkb:',''))
        t['target']=t.target.apply(lambda x: x.replace('uniprotkb:',''))
        #taxid:9606(Homo sapiens)
        t['tax_id_A']=t['tax_id_A'].apply(lambda x: re.sub(r'\(.*$', '', x.replace('taxid:','')))
        t['tax_id_B']=t['tax_id_B'].apply(lambda x: re.sub(r'\(.*$', '', x.replace('taxid:','')))

        #print t[:2]

        #source  target  is_directed     is_stimulation  is_inhibition   dip_url
        #Q13873  Q969L4  1       0       0
        c_gid, c_tax=self.uniprot()

        t['source_gene_id']=t.source.apply(lambda x: c_gid.get(x, 0))
        t['tax_id_A']=t.source.apply(lambda x: c_tax.get(x, 0))
        t['target_gene_id']=t.target.apply(lambda x: c_gid.get(x, 0))
        t['tax_id_B']=t.target.apply(lambda x: c_tax.get(x, 0))
        def f_meth(x):
            m=re.search(r'\((?P<meth>.+?)\)', x)
            if m:
                return m.group('meth')
            return ''
        t['interaction_method']=t.interaction_method.apply(lambda x: f_meth(x))
        t['ds']=t.ds.apply(lambda x: f_meth(x))

        #R1=t.score.apply(lambda x: x.split("|")[0]).values
        #R2=t.score.apply(lambda x: x.split("|")[1]).values
        #t=pd.DataFrame({'X':R1, 'Y':R2})
        #t.to_csv('t.csv', index=False)

        # first score is the final score, 2nd is the initial, see p4 of the suppl. method
        # The initial scores were then transformed into the final scores using this function and by artificially assigning a final score of 1 to the inferred interactions with evidence from at least one pathway database.
        # that means some pairs were set to final score of 1, despite what its initial score is.
        # this observation make me believe the first score is the final score

        t['score']=t.score.apply(lambda x: x.split("|")[0])
        #R=util.sarray2rarray(t.score)
        #print len(t)
        #tmp=t[(t.source=='Q9BTZ2')&(t.target=='Q9NY65')]
        t=t[(t.source_gene_id > 0) & (t.target_gene_id>0) & (t.tax_id_A==t.tax_id_B) & (t.source_gene_id!=t.target_gene_id)]
        #print len(t)
        c_cnt=util.unique_count(t.tax_id_A)
        print c_cnt
        t=t[t.tax_id_A==9606].copy()
        t['interaction_type_id']=4
        t['score']='-'
        t['ds']='InWeb_IM'
        t.rename2({'source_gene_id':'gid_A', 'target_gene_id':'gid_B', 'publication':'support','interaction_method':'interaction_category'})
        t=t[['gid_A','gid_B','tax_id_A','tax_id_B','interaction_type_id','interaction_category','interaction_type','score','support','ds']].copy()
        self.inweb=t
        return self.inweb

    # table can be downloaded using
    def OmniPath(self):
        fn_soure = os.path.join(SyncDB.DOWNLOAD_DIR(),"interactions")
        if not os.path.exists(fn_soure):
            urllib.urlretrieve("http://omnipathdb.org/interactions", fn_soure)
        t=pd.read_table(fn_soure)
        #source  target  is_directed     is_stimulation  is_inhibition   dip_url
        #Q13873  Q969L4  1       0       0
        c_gid, c_tax=self.uniprot()
        t['source_gene_id']=t.source.apply(lambda x: c_gid.get(x, 0))
        t['tax_id_A']=t.source.apply(lambda x: c_tax.get(x, 0))
        t['target_gene_id']=t.target.apply(lambda x: c_gid.get(x, 0))
        t['tax_id_B']=t.target.apply(lambda x: c_tax.get(x, 0))

        #print len(t)
        t=t[(t.source_gene_id > 0) & (t.target_gene_id>0) & (t.tax_id_A==t.tax_id_B) & (t.source_gene_id!=t.target_gene_id)].copy()
        c_cnt=util.unique_count(t.tax_id_A)
        print c_cnt
        if len(c_cnt)>1:
            util.warn_msg('Found non-human genes')
        t=t[t.tax_id_A==9606]
        t['interaction_type_id']=3
        t['interaction_categy']='Signaling Pathway'
        t['interaction_type']='Literature'
        t['score']='-'
        t['ds']='OmniPath'
        t.rename2({'source_gene_id':'gid_A', 'target_gene_id':'gid_B', 'dip_url':'support'})
        t=t[['gid_A','gid_B','tax_id_A','tax_id_B','interaction_type_id','interaction_categy','interaction_type','score','support','ds']].copy()
        self.omni=t
        return self.omni


    def populate_more(self,root):
        self.outputs = [self.fn_dest]
   
    def do_update(self):
        ts = []
        if SyncDB.GET_APP_SETTIGNS("server") != 'METASCAPE':
            gg=GeneGo(self.taxidList)
            gg_res=gg.get_interaction(l_filter=False)
            gg_res.rename2({"GENE_A":"gid_A",
                            "GENE_B":"gid_B",
                            "EFFECT_NAME":"interaction_category",
                            "MECHANISM_NAME":"interaction_type",
                            "TRUST":"score",
                            "PUBMED":"support",
                            "source":"ds"
                            })
            gg_res["ds"]="GeneGO"
            gg_res["interaction_type_id"]=1
            ts.append(gg_res)

        ts.append(self.InWebIM())
        ts.append(self.OmniPath())
        ts.append(self.biogrid())
        t_gp = pd.concat(ts, ignore_index=True)
        t_gp.to_csv(self.fn_dest, index=False)


if __name__=="__main__":

    #gg=GeneGo(["9606,4932"])
    #t=gg.get_interaction(l_filter=False)
    #t.to_csv('GeneGo_PPI.csv', index=False)

    bg = GeneGo(["4932"])
    t = bg.get_data();
    t.to_csv('BioGrid_PPI.csv', index=False)

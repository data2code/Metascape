#!/usr/bin/env python
from __future__ import absolute_import
from __future__ import print_function
#from __future__ import unicode_literals
import db
import pandas as pd
import util
import eutils
import  Entity.PubMed as pm
import re
from six.moves import range

#changed
#import MySQLdb as mysql
import json

ICON_PATH='Content/Images/'
VIS_ICONS= {
    "Circos Plot": "circos.png",
    "Enrichment Bar Graph": "bar_graph.png",
    "Enrichment Heatmap": "heatmap.png",
    "Enrichment Network": "network.png",
    "Protein Complex Network": "protein.png",
    "Membership Pie Chart": "pie_chart.png",
    "Spreadsheet": "data.png"
}

class Citation:

    def __init__(self, db_server='METASCAPE'):
        self.db=db.DB(db_server)
        #self.con = mysql.connect(host='localhost', user='ivy', passwd='ivy_2018', db='ivy_test')

    def get_id(self, PMID):
        t=self.db.from_sql('select id from gp.citation where PMID=?', params=[PMID])
        if len(t)==0:
            return 0
        else:
            return t.ix[0, 'id']

    def fetch_pmid(self, PMID):
        eu=eutils.EUtils()
        args={}
        args['db']='pubmed'
        args['id']=[PMID]
        out=eu.efetch(args)
        X=pm.FetchList(out)
        X=X.to_list(['pubmed_id', 'journal', 'title', 'author'])
        if len(X)==0:
            return {}
        #{'journal.title': 'Pharmaceutical medicine', 'journal.day': None, 'journal.year': '2017', 'title': 'Measuring and Improving Physician Knowledge of Safety Risks Using Traditional and Online Methods in Pharmacovigilance.', 'journal.month': None, 'journal.volume': '31', 'journal.issue': '4', 'author': ['Liede A', 'Amelio J', 'Bennett J', 'Goodman H', 'Peters PM', 'Barber R', 'Kehler E', 'Michael Sprafka J'], 'pubmed_id': '28824275', 'journal.page': '257-266'}
        X=X[0]
        data={'PMID': X['pubmed_id'], 'Title': X['title'], 'Year': X['journal.year'], 'Journal':X['journal.title'], 'Month':X['journal.month'], 'Volume':X['journal.volume'], 'Issue':X['journal.issue'], 'Page':X['journal.page'], 'Authors':", ".join(X['author'])}
        return data


    def insert(self, PMID, PStatus='N', Graph='N'):
        id=self.get_id(PMID)
        if id!=0:
            util.warn_msg('PMID: %s already in collection, id=%s!' % (PMID, id))
            data={}
        else:
            data=self.fetch_pmid(PMID)
            if len(data)==0:
                util.warn_msg('PMID: %s cannot be fetched!' % (PMID))
            else:
                #print data
                self.insert_by_data(data, PStatus, Graph)
            return data

    def insert_by_data(self, data, PStatus='N', Graph='N'):
        data['Title']=data['Title'].encode('utf-8','ignore')
        data['Authors']=data['Authors'].encode('utf-8','ignore')
        self.db.from_sql("insert into gp.citation (id,PMID,Title,Authors,Journal,PYear,Volume,Issue,Page,PStatus,Graph,PMonth,URL,Created) values (null,?,?,?,?,?,?,?,?,?,?,?,?,now())", params=[data['PMID'], data['Title'], data['Authors'], data['Journal'], data['Year'], data['Volume'], data['Issue'], data['Page'], PStatus, Graph, data['Month'],data.get('URL',None)])

    def insert_by_web(self, PMID, Code="", Graph="N"):
        PStatus = 'N'
        s_msg="This is an internal tool!"
        code=Code.strip()
        if code in ("92037","92038"):
            if code=="92038": PStatus="Y"
            data = self.insert(PMID, PStatus, Graph)
            t=self.db.from_sql("select * from gp.citation where PMID=?", params=[PMID])
            s_msg="Invalid PubMed ID: %s!" % PMID
            if len(t):
                s_msg=self.format_entry(t.ix[0])
        s="""
<div class="container">
  <h2>Metascape Citations</h2>
  <p/>
  <div class="panel-group">
    <div class="panel panel-info">
      <div class="panel-heading">New Citation Entry</div>
      <div class="panel-body">"""
        s+=s_msg
        s+="""
    </div>
    <p/>
  </div>
</div>"""
        return s

    def approve(self, PMID):
        id=self.get_id(PMID)
        if id==0:
            util.warn_msg('PMID: %s is not in collection!' % PMID)
        else:
            self.db.from_sql("update gp.citation set PStatus='Y' where PMID=?", params=[PMID])

    def disapprove(self, PMID):
        id=self.get_id(PMID)
        if id==0:
            util.warn_msg('PMID: %s is not in collection!' % PMID)
        else:
            self.db.from_sql("update gp.citation set PStatus='N' where PMID=?", params=[PMID])

    def graph(self, PMID, has_graph='Y'):
        id=self.get_id(PMID)
        if id==0:
            util.warn_msg('PMID: %s is not in collection!' % PMID)
        else:
            self.db.from_sql("update gp.citation set Graph=? where PMID=?", params=[has_graph])

    def remove(PMID):
        id=self.get_id(PMID)
        if id==0:
            util.warn_msg('PMID: %s is not in collection!' % PMID)
        else:
            self.db.from_sql("delete from gp.citation where PMID=?", params=[PMID])

    def get_all(self, l_full_page=False):
        #sw=util.StopWatch()
        t=self.db.from_sql("select * from gp.citation where PStatus='Y' order by PYear desc, Abs(PMID) desc")
        #t=db.from_sql(self.con, "select * from citation where PStatus='Y' order by PYear desc, Abs(PMID) desc")
        #print ">>>>>>>>>", len(t)
        #print t[:]
        #sw.check('Fetched')
        S=[]
        if l_full_page:
            S.append("""
<!DOCTYPE html>
<html lang="en">
  <head>
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css">
    <script src="https://ajax.googleapis.com/ajax/libs/jquery/3.2.1/jquery.min.js"></script>
    <script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/js/bootstrap.min.js"></script>
  </head>
<body ng-app="app">
""")
        S.append("""
<div class="container">
  <h2>Metascape Citations</h2>
  <p/>
  <div class="panel-group">
    <div class="panel panel-info">
      <div class="panel-heading">Citation Info</div>
      <div class="panel-body">Please cite the following paper: Tripathi et al., Cell Host & Microbe (2015), 18: 723-735.; it is our first application. We will really appreciate if you show our URL [http://metascape.org] in your text.  If you like Metascape, help us spread the words.<p/><p/>
The different icons represent different types of visualizations used by the publications, as explained below:<ul>
<!-- <span class="glyphicon glyphicon-eye-open"></span>: graphics<br/> -->
<img src="%s/circos.png" style="width:32px;height:32px"> circos plot;
<img src="%s/bar_graph.png" style="width:32px;height:32px"> enrichment bar graph;
<img src="%s/heatmap.png" style="width:32px;height:32px"> enrichment heatmap;
<img src="%s/network.png" style="width:32px;height:32px"> enrichment network;<br>
<img src="%s/protein.png" style="width:32px;height:32px"> protein complex network;
<img src="%s/pie_chart.png" style="width:32px;height:32px"> membership pie chart;
<img src="%s/data.png" style="width:32px;height:32px"> spreadsheet.</ul>
    </div>
    <p/>
""" % (ICON_PATH, ICON_PATH, ICON_PATH, ICON_PATH, ICON_PATH, ICON_PATH, ICON_PATH))
        iB=iE=0
        n=len(t)
        i_cnt=n
        for i in range(1, n+1):
            if i>=n or t.ix[i, 'PYear']!=t.ix[i-1, 'PYear']:
                iE=i
                k=t.ix[iB, 'PYear']
                #print ">>>>>>>>>> ", k, iE-iB

                S.append("""
    <div class="panel panel-info">
      <div class="panel-heading">Year %s</div>""" % k)
                S.append("""
      <div class="panel-body">
""")
                for j in range(iB, iE):
                    S.append(("%d. " % i_cnt)+self.format_entry(t.ix[j]))
                    i_cnt-=1
                S.append("""
      </div>
    </div>
    <p/>
""")
                iB=iE
        S.append("""
  </div>
</div>""")
        if l_full_page:
            S.append("\n</body></html>")
        #sw.check('DONE')
        return "".join(S)

    def format_entry(self, data):
        data=dict(data)
        id=data['PMID']
        if id>0:
            s_url='http://www.ncbi.nlm.nih.gov/pubmed/%s' % id
        else:
            s_url=data['URL']
        S=[]
        if data['Journal'] is not None:
            S.append(data['Journal'])
        if data['PYear'] is not None:
            s_month=data['PMonth']
            if s_month is not None:
                s_month=re.sub(r'^\d\d_', '', s_month)
                S.append(str(data['PYear']))
                S.append(s_month+";")
            else:
                S.append(str(data['PYear'])+";")
        if data['Volume'] is not None:
            s_page=""
            if data['Page'] is not None:
                s_page=":"
            if data['Issue'] is not None:
                S.append(str(data['Volume'])+"("+str(data['Issue'])+")"+s_page)
            else:
                S.append(str(data['Volume'])+s_page)
        if data['Page'] is not None:
                S.append(str(data['Page'])+".")
        #data['Title']=data['Title'].decode('utf-8', 'ignore')
        #data['Authors']=data['Authors'].decode('utf-8', 'ignore')
        #print "----------------------"
        #print S
        s_journal=" ".join(S)
        s_graph='' # <span class="glyphicon glyphicon-eye-open"></span>' if (data['Graph'] is not None and data['Graph']=='Y') else ''
        s_vis = data['visualization']
        if not pd.isnull(s_vis) and s_vis!='':
            vis_data = json.loads(data['visualization'])
            s_graph="".join(['  <img src="%s" style="width:16px;height:16px" title="%s">' % (ICON_PATH+VIS_ICONS[k], v) for k,v in vis_data.items() ])
        s='<a href="%s">%s</a><br>%s.<br>%s%s<p/>\n' % (s_url, data['Title'], data['Authors'], s_journal, s_graph)
        return s

    def export_visualization(self):
        t = self.db.from_sql('select PMID, Journal, Title, PYear, visualization from gp.citation where visualization is not null order by abs(PMID) desc')
        t2 = pd.DataFrame(t, columns=['Circos Plot', 'Enrichment Bar Graph', 'Enrichment Heatmap', 'Enrichment Network', 'Protein Complex Network', 'Membership Pie Chart', 'Spreadsheet'])
        t2.insert(0, 'Title', t.Title)
        t2.insert(0, 'Year', t.PYear)
        t2.insert(0, 'Journal', t.Journal)
        t2.insert(0, 'PMID', t.PMID)
        for i in t.index:
            vis_data = json.loads(t.ix[i, 'visualization'])
            for k,v in vis_data.items():
                t2.ix[i, k] = v
        t2.fillna('', inplace=True)
        return t2

if __name__=="__main__":
    ids_graph=[26651948, 28683259,28272987,28558796,28445719,28416637,28625551,28248290,28065413,28514688,28157506,27980230,28476932,28801683,28428242,28529766,28556776,28630416,28754678,28679947,28794002,28296633,28322580,28193876]
    ids_nograph=[28742165,28719624,27489048,28765524,28167044,27633323,28537707,27376549,28757253,28602351,28111071,28397821,28306507,26657882,28394251,27713570,27980103,28698297,27656401,27566152]
    # http://www.biorxiv.org/content/biorxiv/early/2017/06/20/152652.full.pdf
    # http://www.sciencedirect.com/science/article/pii/S0168170215301337
    # http://www.biorxiv.org/content/early/2017/06/20/152769
    # http://www.biorxiv.org/content/biorxiv/early/2017/04/14/127522.full.pdf
    # http://www.biorxiv.org/content/biorxiv/early/2017/07/05/159384.full.pdf
    # http://archiv.ub.uni-heidelberg.de/volltextserver/21963/
    # http://archiv.ub.uni-heidelberg.de/volltextserver/21963/1/PakKin_LOU_thesis_merged_e.pdf
    # http://www.biorxiv.org/content/early/2017/05/24/141952.full.pdf
    x=Citation()
    #print x.insert_by_web(28824275, '92130', 'Y')
    #exit()
    #x.fetch_pmid(28824275)
    #x.insert(28824275)
    #x.insert(28381584)
    #x.insert(24466558)
    #x.insert(27811064)
    #x.insert(26657882)
    l_insert=True
    if l_insert:
        for id in ids_nograph:
            print(">>%d" % id)
            x.insert(id, PStatus='Y', Graph='N')
        for id in ids_graph:
            x.insert(id, PStatus='Y', Graph='Y')
    print(x.get_all(l_full_page=True))

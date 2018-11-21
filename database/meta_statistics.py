#!/usr/bin/env python
import numpy as np
import pandas as pd
import util
import db
import matplotlib
import scipy.misc

from core import SyncDB

matplotlib.use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from textwrap import wrap
import os
import glob

class Statistics:

    def __init__(self):
        self.db=db.DB(name=SyncDB.CONNECTION_ID, db=SyncDB.DATABASE)
        print SyncDB.CONNECTION_ID
        self.db.from_sql('use {0}'.format(SyncDB.DATABASE))
        self.data=None
        self.get_date()
        #self.d1='2016-10-03'
        #self.d2='2016-06-14'
        S_tax_id=self.db.from_sql('select distinct tax_id from gid2source_id')['tax_id'].tolist()
        self.t_tax=self.db.sql_in('select tax_id,display_name tax_name from taxid2name where tax_id in (', ') order by display_name', S_tax_id)

    def get_date(self):
        t=self.db.from_sql('select distinct history from statistics order by history desc limit 0,2')
        print t
        self.d1=self.d2=None
        if len(t)==0:
            util.error_msg('No history data is found')
        elif len(t)==1:
            util.warn_msg('Only one history entry')
            self.d1=str(t.ix[0, 'history'])
        else:
            self.d1=str(t.ix[0, 'history'])
            self.d2=str(t.ix[1, 'history'])
        print self.d1
        print self.d2

    def get_data(self):
        t=pd.DataFrame()
        if self.d1 is None: return t
        if self.d2 is None:
            s_date="('%s')" % self.d1
        else:
            s_date="('%s','%s')" % (self.d1,self.d2)
        self.data=self.db.from_sql('select * from statistics where history in %s order by history desc' % s_date)
        self.data['tax_id'].fillna(9606, inplace=True) #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        self.data['tax_id']=self.data['tax_id'].astype(int)
        self.data['history']=self.data['history'].astype(str)
        self.data=self.data.merge(self.t_tax, left_on='tax_id', right_on='tax_id', how='left')
        self.data['tax_name'].fillna('Any', inplace=True)

    @staticmethod
    def split_value(R):
        if len(R)==1:
            return (R[0], 0, 0, R[0])
        else:
            return (min(R[0], R[1]), max(R[1]-R[0], 0), max(R[0]-R[1], 0), max(R[0], R[1], 1))

    def plot_by_type(self, s_png, t_data=None, n_per_img=4, n_per_row=2, l_pct=True, ylabel='Genes'):
        """l_pct: controls whether to plot the normalized percentage data"""
        if t_data is None:
            if self.data is None:
                self.get_data()
            t_data=self.data
        X_label=['Any']+self.t_tax.tax_name.tolist()
        X_label = [ '\n'.join(wrap(l, 10)) for l in X_label]
        X_id=[0]+self.t_tax.tax_id.tolist()

        n=len(X_id)
        c_idx={tax_id: i for i,tax_id in enumerate(X_id) }
        X=np.array([i+1 for i in range(n)])
        Y_cnt=[]
        Y_pct=[]
        S_title=[]
        S_out=[]

        if l_pct:
            factor=2
        else:
            factor=1

        for k,t_v in t_data.groupby(['table_name','type_id','type_name']):
            t_v=t_v.copy()
            S_title.append(k[2]+"."+str(k[1]))
            Y_blue=np.zeros(n) # base=min(new_val, old_val)
            Y_red=np.zeros(n)  # red=max(old_val-new_val, 0)
            Y_green=np.zeros(n)# green=max(new_val-old_val, 0)
            # normalized version, percentage
            Y_blue2=np.zeros(n) # base=min(new_val, old_val)
            Y_red2=np.zeros(n)  # red=max(old_val-new_val, 0)
            Y_green2=np.zeros(n)# green=max(new_val-old_val, 0)

            for k2,t_v2 in t_v.groupby(['tax_id']):
                val=np.ones(2)
                for i in t_v2.index:
                    val[ 0 if t_v2.ix[i, 'history']==self.d1 else 1 ]=t_v2.ix[i, 'total']
                b,r,g,m=Statistics.split_value(val)
                idx=c_idx.get(k2)
                Y_blue[idx]=b
                Y_red[idx]=r
                Y_green[idx]=g
                Y_blue2[idx]=b*100/m
                Y_red2[idx]=r*100/m
                Y_green2[idx]=g*100/m
            Y_cnt.append([Y_blue, Y_red, Y_green])
            Y_pct.append([Y_blue2, Y_red2, Y_green2])

        n_type=len(Y_cnt)
        width=0.8
        cnt=0
        font = {'family' : 'normal', 'size'   : 6}
        matplotlib.rc('font', **font)
        for i in range(n_type):
            i_offset=i % n_per_img
            #print ">>> ", n_per_img/n_per_row*2, n_per_row, i_offset+1
            if i_offset==0:
                plt.figure(figsize=(7,6))
            plt.subplot(n_per_img/n_per_row*factor, n_per_row, i_offset+1+ (i_offset/n_per_row*n_per_row if l_pct else 0))
            plt.bar(X, Y_cnt[i][0], width, color='#3182bd', edgecolor = "none")
            plt.bar(X, Y_cnt[i][1], width, bottom=Y_cnt[i][0], color='#fdae6b', edgecolor = "none")
            plt.bar(X, Y_cnt[i][2], width, bottom=Y_cnt[i][0]+Y_cnt[i][1], color='#a1d99b', edgecolor = "none")
            if i_offset % n_per_row==0:
                plt.ylabel(ylabel, fontsize=8)

            if l_pct:
                plt.gca().axes.get_xaxis().set_visible(False)
                plt.title(S_title[i])
                plt.subplot(n_per_img/n_per_row*factor, n_per_row, i_offset+1+n_per_row+i_offset/n_per_row*n_per_row)
                plt.bar(X, Y_pct[i][0], width, color='#3182bd', edgecolor = "none")
                plt.bar(X, Y_pct[i][1], width, bottom=Y_pct[i][0], color='#fdae6b', edgecolor = "none")
                plt.bar(X, Y_pct[i][2], width, bottom=Y_pct[i][0]+Y_pct[i][1], color='#a1d99b', edgecolor = "none")
                if i_offset % n_per_row==0:
                    plt.ylabel('%', fontsize=8)
                if i_offset>=n_per_img-n_per_row:
                    plt.subplots_adjust(bottom=0.5)
                    plt.xticks(X+width/2.0, X_label, rotation=90, fontsize=6)
                else:
                    plt.gca().axes.get_xaxis().set_visible(False)
            else:
                plt.subplots_adjust(bottom=0.5)
                plt.xticks(X+width/2.0, X_label, rotation=90, fontsize=6)

            plt.title(S_title[i])
            if i_offset==n_per_img-1 or i==n_type-1:
                #plt.tight_layout()
                for j in range(max(n_type-1, n_per_row), n_per_img):
                    plt.subplot(n_per_img/n_per_row*2, n_per_row, j+1+(n_per_row+j/n_per_row*n_per_row if l_pct else 0))
                    plt.bar(X, np.zeros(n), width, color='#3182bd', edgecolor = "none")
                    plt.subplots_adjust(bottom=0.5)
                    plt.xticks(X+width/2.0, X_label, rotation=90, fontsize=6)
                    plt.gca().axes.get_yaxis().set_visible(False)
                cnt+=1
                s_file=s_png+"_%d.png" % cnt
                plt.savefig(s_file, dpi=100)
                Statistics.crop(s_file)
                plt.close()
                src_url = "Content/MenuPages/out/{0}".format(os.path.basename(s_file))
                S_out.append('<img src="%s">\n' % src_url)

        return "".join(S_out)

    def plot_by_organism(self, s_png, t_data=None, l_pct=True, ylabel='Genes'):
        """l_pct: controls whether to plot the normalized percentage data"""
        if t_data is None:
            if self.data is None:
                self.get_data()
            t_data=self.data

        S_out=[]
        for k,t_v in t_data.groupby(['tax_id','tax_name']):
            if k[0] == 4932:
                a =5
            t_v=t_v.copy()
            X_label=[]
            S_title=[]
            Y_cnt=[]
            Y_pct=[]
            tax_id=k[0]

            Y_blue=[]
            Y_red=[]
            Y_green=[]
            Y_blue2=[]
            Y_red2=[]
            Y_green2=[]
            X=[]

            for k2,t_v2 in t_v.groupby(['table_name','type_id','type_name']):
                t_v2=t_v2.copy()

                X.append(k2[0]+"."+k2[2])
                val=np.ones(2)
                for i in t_v2.index:
                    val[ 0 if t_v2.ix[i, 'history']==self.d1 else 1 ]=t_v2.ix[i, 'total']
                b,r,g,m=Statistics.split_value(val)
                Y_blue.append(b)
                Y_red.append(r)
                Y_green.append(g)
                Y_blue2.append(b*100/m)
                Y_red2.append(r*100/m)
                Y_green2.append(g*100/m)

            n_chunk=max(round(len(X)/20.0), 1)
            if n_chunk==1:
                S_title.append(k[1])
                X_label.append(X)
                Y_cnt.append([np.array(Y_blue), np.array(Y_red), np.array(Y_green)])
                Y_pct.append([np.array(Y_blue2), np.array(Y_red2), np.array(Y_green2)])
            else:
                I=util.split(range(len(X)), n_chunk=n_chunk)
                for i,I_ in enumerate(I):
                    S_title.append(k[1]+'.'+str(i+1))
                    X_label.append(X[I_[0]:I_[-1]])
                    Y_cnt.append([np.array(Y_blue[I_[0]:I_[-1]]), np.array(Y_red[I_[0]:I_[-1]]), np.array(Y_green[I_[0]:I_[-1]])])
                    Y_pct.append([np.array(Y_blue2[I_[0]:I_[-1]]), np.array(Y_red2[I_[0]:I_[-1]]), np.array(Y_green2[I_[0]:I_[-1]])])

            n=len(X_label)
            width=0.8
            cnt=0
            if l_pct:
                n_per_img=1
                factor=2
            else:
                n_per_img=2
                factor=1
            font = {'family' : 'normal', 'size'   : 6}
            matplotlib.rc('font', **font)
            factor=2 if l_pct else 1
            for i in range(n):
                i_offset=i % n_per_img
                if i_offset==0:
                    plt.figure(figsize=(7,6))
                plt.subplot(n_per_img*factor, 1, i_offset+1)
                m=len(X_label[i])
                X=np.array([j+1 for j in range(m)])
                plt.bar(X, Y_cnt[i][0], width, color='#3182bd', edgecolor = "none")
                plt.bar(X, Y_cnt[i][1], width, bottom=Y_cnt[i][0], color='#fdae6b', edgecolor = "none")
                plt.bar(X, Y_cnt[i][2], width, bottom=Y_cnt[i][0]+Y_cnt[i][1], color='#a1d99b', edgecolor = "none")
                plt.ylabel(ylabel, fontsize=8)
                plt.title(S_title[i])

                labels = [ '\n'.join(wrap(l, 20)) for l in X_label[i]]
                if l_pct:
                    plt.gca().axes.get_xaxis().set_visible(False)
                    plt.subplot(n_per_img*factor, 1, i_offset+2)
                    plt.bar(X, Y_pct[i][0], width, color='#3182bd', edgecolor = "none")
                    plt.bar(X, Y_pct[i][1], width, bottom=Y_pct[i][0], color='#fdae6b', edgecolor = "none")
                    plt.bar(X, Y_pct[i][2], width, bottom=Y_pct[i][0]+Y_pct[i][1], color='#a1d99b', edgecolor = "none")
                    plt.ylabel('%', fontsize=8)
                plt.subplots_adjust(bottom=0.5)
                plt.xticks(X+width/2.0, labels, rotation=90, fontsize=6)
                if i_offset==n_per_img-1 or i==n-1:
                    plt.tight_layout()
                    cnt+=1
                    s_file=s_png+"_%d_%d.png" % (tax_id, cnt)
                    plt.savefig(s_file, dpi=100)
                    Statistics.crop(s_file)
                    plt.close()
                    src_url = "Content/MenuPages/out/{0}".format(os.path.basename(s_file))
                    S_out.append('<img src="%s">\n' % src_url)
        return "".join(S_out)

    def plot_all(self, s_dir='out', t_data=None, l_pct=True):

        def scale(t_data, table_name=None, ds=None, type_name=None, factor=1):
            tmp=t_data
            if table_name is not None:
                tmp=tmp[tmp.table_name==table_name]
            if ds is not None:
                tmp=tmp[tmp.ds==ds]
            if type_name is not None:
                tmp=tmp[tmp.type_name==type_name]

            for i in tmp.index:
                t_data.ix[i, 'total']/=factor
                t_data.ix[i, 'table_name']+=' (x%d)' % factor

        if t_data is None:
            if self.data is None:
                self.get_data()
            t_data=self.data

        I=np.array(t_data.ds.apply(lambda x: x in ('GeneGo')))
        t_genego=t_data[I].copy()
        t_data=t_data[~I].copy()
        scale(t_genego, table_name='interaction', factor=25)

        print ">>> GeneGo", len(t_genego)

        I=np.array(t_data.table_name.apply(lambda x: x in ('term','term2term')))
        t_term=t_data[I].copy()
        t_data=t_data[~I].copy()
        print ">>> Term", len(t_term)

        I=np.array(t_data.table_name.apply(lambda x: x in ('interaction')))
        t_ppi=t_data[I].copy()
        t_data=t_data[~I].copy()
        print ">>> PPI", len(t_ppi)

        # those human-only annotations
        t_type_id=self.db.from_sql('''SELECT type_id,type_name, ds
                                      FROM statistics
                                      WHERE table_name='annotation'
                                      AND history IN ('{0}')
                                      GROUP BY type_id,type_name
                                      HAVING COUNT(DISTINCT tax_id)=1'''.format(self.d2))
        S_type_id=set(t_type_id.type_id)
        I=np.array(t_data.table_name=='annotation') & np.array(t_data.type_id.apply(lambda x: x in S_type_id))
        t_human=t_data[I].copy()
        print ">>> Human-only", len(t_human)

        t_by_org=t_data[~I].copy()
        print ">>> Remaining", len(t_by_org)

        #print len(t_by_type), len(t_by_org), len(t_genego), len(t_data)

        S_out=[]
        if len(t_genego) and False:
            S_out.append("<h3>Non-public Data Sources</h3>\n")
            S_out.append(self.plot_by_organism(s_dir+'/NonPublic', t_data=t_genego, l_pct=l_pct))
        scale(t_term, table_name='term2term', factor=100)
        scale(t_term, table_name='term', type_name='GO Biological Processes', factor=2)
        scale(t_term, table_name='term', type_name='L1000 shRNA', factor=5)
        scale(t_term, table_name='term', type_name='L1000 Compound', factor=20)
        scale(t_term, table_name='term', type_name='L1000 cDNA', factor=5)

        #print(util.unique_count(t_term.table_name), util.unique_count(t_term.type_name))
        scale(t_by_org, table_name='term2gids', type_name='GO Biological Processes', factor=2)
        scale(t_by_org, table_name='term2gids', type_name='L1000 shRNA', factor=10)
        scale(t_by_org, table_name='term2gids', type_name='L1000 Compound', factor=40)
        scale(t_by_org, table_name='term2gids', type_name='L1000 cDNA', factor=10)
        scale(t_by_org, table_name='term2gids', type_name='L1000 Ligand', factor=2)
        scale(t_by_org, table_name='term2gids', type_name='DisGeNET', factor=2)
        #print(util.unique_count(t_data.table_name), util.unique_count(t_data.type_name))
        panel_str = '''
            <div class="panel panel-info">
                <div class="panel-heading">{0}</div>
                <div class="panel-body">{1}</div>
            </div>
        '''
        S_out.append('<div class="container"><div class="panel-group">')
        import datetime
        now = datetime.datetime.now()
        S_out.append('<div style="font-size: 150%;">Database Last Update Date: {0}</div>'.format(now.strftime("%Y-%m-%d")))
        S_out.append('<img src="Images/UpdateLegend.png"/>')
        S_out.append(panel_str.format("Ontology Terms and Term-Term Relationships",
                    self.plot_by_organism(s_dir+'/Term', t_data=t_term, l_pct=l_pct, ylabel='Counts')))
        S_out.append(panel_str.format("Protein-protein Interaction",
                    self.plot_by_type(s_dir+'/PPI', t_data=t_ppi, n_per_img=1, n_per_row=1, l_pct=l_pct, ylabel='Pairs')))
        S_out.append(panel_str.format("Human-specific Gene Annotations",
                    self.plot_by_organism(s_dir+'/HumanOnly', t_data=t_human, l_pct=l_pct)))
        S_out.append(panel_str.format("Remaining Data Sources by Organism",
                    self.plot_by_organism(s_dir+'/Organism', t_data=t_by_org, l_pct=l_pct, ylabel='Counts')))
        S_out.append('</div></div>')
        return "".join(S_out)

    @staticmethod
    def crop(s_png, margin=10, left=0, right=0, top=0, bottom=0):
        M=scipy.misc.imread(s_png)
        h,w,z=M.shape
        r=M[:,:,0]
        g=M[:,:,1]
        b=M[:,:,2]
        bg=np.logical_and(r==255, g==255, b==255)
        xa=xb=ya=yb=0
        X=np.sum(bg, axis=0)
        #print h,w,X, len(X)
        for i in xrange(w):
            if X[i]<h:
                xa=i
                break
        for i in xrange(w-1,-1,-1):
            if X[i]<h:
                xb=i
                break
        X=np.sum(bg, axis=1)
        for i in xrange(h):
            if X[i]<w:
                ya=i
                break
        for i in xrange(h-1,-1,-1):
            if X[i]<w:
                yb=i
                break
        if margin>0 or left>0: xa=max(0,xa-max(margin, left))
        if margin>0 or right>0: xb=min(w-1,xb+max(margin, right)+1)
        if margin>0 or top>0: ya=max(0,ya-max(margin,top))
        if margin>0 or bottom>0: yb=min(h-1,yb+max(margin, bottom)+1)
        if xa==0 and xb==w-1 and ya==0 and yb==h-1:
            print "No cropping"
            return # no cropping
        else:
            print "Crop: width %d:%d, height: %d:%d" % (ya, yb, xa, xb)
            M=M[ya:yb+1,xa:xb+1,:]
            scipy.misc.imsave(s_png, M)

    def generate_statistic_html(self):
        s_dir = '/home/meta_test/metascape_keep/Content/MenuPages'


        s_out_dir = os.path.join(s_dir,'out')
        if not os.path.exists(s_out_dir):
            os.mkdir(s_out_dir)
        else:
            for x in glob.glob(s_out_dir+"/*.png"):
                os.remove(x)
        self.get_data()
        s_html=self.plot_all(l_pct=False, s_dir=s_out_dir)
        util.save_string(os.path.join(s_dir,'statistics.html'), s_html)



if __name__=="__main__":
    import os
    import shutil
    import glob
    import argparse as arg
    opt=arg.ArgumentParser(description='Visualize Metascape Database Statistics')
    opt.add_argument('-d', '--output', default='out', help='output file directory')
    opt.add_argument('-o', '--html', default='statistics.html', help='output html summary file')

    args=opt.parse_args()
    if not args.output:
        util.error_msg('Output folder is missing!')
    else:
        s_dir=args.output
        s_stat=args.html

    if not os.path.exists(s_dir):
        os.mkdir(s_dir)
    else:
        for x in glob.glob(s_dir+"/*.png"):
            os.remove(x)

    x=Statistics()
    x.get_data()
    s_html=x.plot_all(l_pct=False)
    util.save_string(s_stat, s_html)


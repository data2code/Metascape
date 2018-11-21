#!/usr/bin/env python -u
import pandas as pd
import numpy as np
import os

import socket
import util
import xml.etree.ElementTree as ET
from pprint import pprint
import sys
import db
import itertools
import pprint as pprint
import pdb
from IPython.core.debugger import Tracer

class XmlClass(object):
    def __init__(self,xe):
        self.options=xe.attrib
        self.tag=xe.tag
        self.id=xe.tag
        self.inputs=[]
        self.outputs=[]
        self.ds = None
        self.fn_dest = None #initialed in populate_more,  which call subclass's get_fn_dest
        if 'ds' in xe.attrib:
            self.ds = xe.attrib['ds']
        self.children = [ SyncDB.task_factory(c) for c in xe ]
        if 'dependence' in xe.attrib:
            self.dependence = xe.attrib['dependence']
        self.has_calculated = None

    def do_work(self):
        #get input objects
        print '***enter do_work', self.id
        #Tracer()();
        if self.has_calculated is not None:
            print 'has calculated', self.id, 'return:', self.has_calculated
            return self.has_calculated
        print self.print_me()
        input_dirty = self.prepare_inputs()
        #check all the output is existing , len is not 0
        output_dirty = any([self.is_dirty(f) for f in self.outputs])
        print 'input, output dirty checking: ',input_dirty,output_dirty,'in', self.id

        update = True
        if ((not input_dirty)\
           and  (not output_dirty)\
           and not SyncDB.IS_REBUILD):
            update = False
        else:
            update = True
        if len(SyncDB.DEBUG_NODE_IDS):
            if self.id in SyncDB.DEBUG_NODE_IDS:
                update = True
                print 'input, output dirty checking: turn on by debug nodes'
            elif SyncDB.DEBUG_NODE_IDS_ONLY:
                update = False
        if not update:
            print '***do nothing,return from:',self.id
            r = False
        else:
            print '***calling do_update in ',self.id
            r = True
            if SyncDB.DONT_RECOVER:
                r = self.do_update()
            else:
                try:
                    r = self.do_update()
                except Exception as exp:
                    not_recovered = True
                    import traceback
                    err_msg = '!!!!!!!!!!!Error occurred while running %s (%s) \n (%s).'% (self.id,exp, traceback.format_exc())
                    print err_msg
                    SyncDB.ERROR_LOG.write(err_msg)
                    if hasattr(SyncDB,'RECOVERY_FOLDER') and os.path.exists(SyncDB.RECOVERY_FOLDER):
                        #recover files from this folder
                        for f in self.outputs:
                            file = os.path.basename(f)
                            if os.path.exists(SyncDB.RECOVERY_FOLDER+"/"+file):
                                SyncDB.ERROR_LOG.write('\tRecovered %s from folder %s\n'%(file,SyncDB.RECOVERY_FOLDER));
                                print '\tRecovered %s from folder %s\n'%(file,SyncDB.RECOVERY_FOLDER);
                                util.unix('cp ' + SyncDB.RECOVERY_FOLDER+"/"+file + ' ' + SyncDB.DOWNLOAD_DIR());
                            else:
                                SyncDB.ERROR_LOG.write('\tDid not recover %s from folder %s (not such file)\n'%(file,SyncDB.RECOVERY_FOLDER));
                                print '\tDid not recover %s from folder %s (not such file)\n'%(file,SyncDB.RECOVERY_FOLDER);
                    else:
                        print '\tRecovery folder does not exist for %s'%(str(self.outputs));
                        SyncDB.ERROR_LOG.write('\tRecovery folder does not exist\n');
                        SyncDB.ERROR_LOG.write(str(self.outputs));
                        
            if r is None:
                r = True
        self.has_calculated = r
        return r

    def do_update(self):
        if not isinstance(self,Main):
            print 'Error:we must overide this method' 
        #only Main.py is using the belows
        for c in self.children:
             c.do_work()
    
    def is_dirty(self,fn):
        dirty = True
        import os.path, time,datetime
        if hasattr(SyncDB,'DEBUG_NODE_IDS'):
            if (self.id in SyncDB.DEBUG_NODE_IDS):
                return True
            if os.path.isfile(fn) or fn.startswith('CsvInMem_'):
                return False
            if SyncDB.DEBUG_NODE_IDS_ONLY:
                return  False
            return True

        if fn.startswith('CsvInMem_'):
            return False

        if os.path.isfile(fn): 
            statinfo = os.stat(fn)
            if statinfo.st_size > 0:
                file_time = os.path.getmtime(fn)
                current_time = datetime.datetime.now()
                time_delta_hour = (current_time -datetime.datetime.fromtimestamp(file_time)).total_seconds()/3600
                print 'time_delta_hour', time_delta_hour
                if time_delta_hour <= SyncDB.IGNORE_HOURS_DOWNLOAD:
                    dirty  = False
        return dirty


    def prepare_inputs(self):
        #Tracer()();
        dirty = False
        for fn in self.inputs:
            if fn.startswith('ds:'):
                ds = fn.replace('ds:','')
                if (SyncDB.IS_REBUILD and  ((len(SyncDB.TO_SYNC_DS)==0) or (ds in SyncDB.TO_SYNC_DS))):
                    dirty = True
                continue
            o = SyncDB.ALL_TASKS.get_task_by_output(fn)
            dirty = dirty | o.do_work()
        return dirty

    def get_id(self):
        return self.tag

    def populate_more(self,root):
        if hasattr(self,'get_fn_dest'): #termgenepair doesn't have get_fn_dest
            self.fn_dest = self.get_fn_dest()
        self.id = self.get_id()
        self.outputs.append(self.fn_dest)
        
        if hasattr(self,'dependence'):
            for f in self.dependence.split(','):
                f = f.strip()
                self.inputs.append(self.get_dependence_fn(f))

    def get_dependence_fn(self,f):
        return os.path.join(SyncDB.DOWNLOAD_DIR(),f)

    #def __str__(self):
    def print_me(self):
        return '-------------\nid:%s\ninputs:%s\noutputs:%s\nds:%s' % (self.id,self.inputs,self.outputs,self.ds)

    def get_xmlclass_instances(self):
        return [self]


class ConfigRoot(XmlClass):
    def __init__(self, xe):
         XmlClass.__init__(self, xe = xe)
         self.all_children = []

    def get_task_by_output(self, input):
        for c in self.all_children:
            if input in c.outputs:
                return c
        return None

    def get_ds(self,child):
        if child.ds is not None:
            return child.ds
        ds=set()
        for fn in child.inputs:
            o = SyncDB.ALL_TASKS.get_task_by_output(fn)
            if o is None:
                print '%s has input: %s and cannot find any class which generates this file' % (child.id,fn)
            print '%s has input: %s and %s generates this file' % (child.id,fn,o.id)
            o.ds = self.get_ds(o)
            if o.ds is not None:
                 ds_c_set =  set(o.ds.split(','))
                 ds.update(ds_c_set)
        if len(ds) > 0:
            return ','.join(ds)
        return None

    def process_settings(self, settings):
        for c in settings.children:
            pass
                    
                
    def do_work(self):
        #populate inputs, outputs
        #get all children
        #get all upload children
        all_upload_children = []

        print '---------inputs, outputs  begin----------------------------------------------------'
        for c in self.children:
            if 'settings.Settings' in str(c):
                self.process_settings(c)
            elif 'main.Main' in str(c):
                self.all_children.extend(c.children)
                all_upload_children.extend(c.children)
            else:
                xs = c.get_xmlclass_instances()
                self.all_children.extend(xs)
        for c in self.all_children:
            c.populate_more(self)
            print c.print_me()
        for c in self.all_children:
            print 'to do get_ds for:',c.id
            c.ds = self.get_ds(c)
            print c.print_me()
        print '---------inputs, outputs  end-----------------------------------------------------------'

        for c in all_upload_children:
            c.do_work()

    def check_inputs(self):
        #check input files are available
        passed = True
        from gputil import GPUtils
        if not GPUtils.check_inputs():
            passed=False
            print 'Inputs are missing in GPUtils'            
            
        all_children = []
        print '---------inputs, outputs  begin----------------------------------------------------'
        for c in self.children:
            if 'main.Main' in str(c):
                all_children.extend(c.children)
            else:
                xs = c.get_xmlclass_instances()
                all_children.extend(xs)

        for c in all_children:
            #print 'check inputs for:',c.id
            if hasattr(c, 'check_inputs'):
                if not c.check_inputs():
                    print 'Inputs are missing for', c.id
                    passed=False;            
        return passed;


class CsvConvert(XmlClass):
    def __init__(self, xe):
        XmlClass.__init__(self,xe=xe)
        self.fn_source = os.path.join(SyncDB.DOWNLOAD_DIR(),xe.attrib['source'])
        if 'newCols' in xe.attrib:
            self.new_cols = map(str.strip, xe.attrib['newCols'].split(','))
        if 'dest' in xe.attrib:
            self.dest = xe.attrib['dest']
        self.add_col_element(xe)
   
    def add_col_element(self,xe):
        if 'newCols' in xe.attrib: # term gene pair doesn't have 'newCols'
            defined_cols = []
            for x in xe:
                if 'colName' in x.attrib:
                    defined_cols.append(x.attrib['colName'])
            todo_names = {}
            for x in util.minus(self.new_cols,defined_cols):
                x = x.strip()
                xml = '<Col colName="%s" sourceCol="%s"/>' % (x,x)
                root = ET.fromstring(xml)
                c = SyncDB.task_factory(root)
                self.children.append(c)

    def get_fn_dest(self):
        return os.path.join(SyncDB.DOWNLOAD_DIR(),self.dest)

    def generate_new_row(self,row):
        result={}
        '''
        #when taxidList tag is missing the default species list defined in SyncDB.SUPPORTED_SPECIES is used.
        sp_xml= """<Species 
                        taxidCol='tax_id' />"""
        species = SyncDB.task_factory(ET.fromstring(sp_xml))
        
        for child in self.children:
            if type(child).__name__ == "Species":
                species = child;
                break;
                
        if not species.check(row=row):
            return None;            
        '''
        
        for col in self.children:
            if type(col).__name__ != "Species":
                r = col.do_work(row=row)
                result[col.col_name] = r
                if col.required  and ((r is None) | (len(str(r)) == 0)):
                    return None
        return result

    def get_chunksize(self):
        return 10000

    def do_update(self):
        #read csv by chunk
        kwargs={}
        taxidList = SyncDB.SUPPORTED_SPECIES;
        taxidCol = 'tax_id'
        for child in self.children:
            if type(child).__name__ == "Species":
                taxidList = child.supported_species
                taxidCol = child.taxid_col
                break;

        if 'oldCols' in self.options:
            kwargs['names']= self.options['oldCols'].split(',')
        if 'read_csv' in self.options:
            for kv_str in  self.options['read_csv'].split(','):
                kv = kv_str.split('=')
                kwargs[kv[0]] = kv[1]
                if kv[1] == 'None':
                    kwargs[kv[0]] = None
                if kv[0].lower() == 'skiprows':
                    kwargs[kv[0]] = int(kv[1])
                if kv[0].lower() == 'quoting':
                    kwargs[kv[0]] = int(kv[1])                    
                if kv[0].lower() == 'sep':
                    kwargs[kv[0]] = kv[1].replace('\\t', '\t')

        if 'read_csv' not in self.options or 'dtype' not in self.options['read_csv']:
            kwargs['dtype'] = str
                    
        #need to check if tax_id is in the file.
        kwargs['nrows']=1
        #Tracer()()
        check=util.read_csv(self.fn_source,**kwargs)

        if not taxidCol in check.columns:
            taxidCol = None
        del kwargs["nrows"]
        
        self.filter_str = None
        if 'filter_str' in self.options:
            self.filter_str = self.options['filter_str'];

        #['0'] as tax_id means everything.
        if taxidCol is not None and len(taxidList) != 0:
            if self.filter_str is None:
                self.filter_str = taxidCol + ' in [' + ','.join(['"%s"'%t for t in (taxidList + ['0'])]) + ']';
            else:
                self.filter_str = '(' + taxidCol + ' in [' + ','.join(['"%s"'%t for t in (taxidList + ['0'])]) + ']) and (' + self.filter_str + ')';
        
        print 'read file: %s with options %s' % (self.fn_source,kwargs)
        
        iter_csv=None
        if 'CsvInMem' in self.fn_source:
            iter_csv = SyncDB.get_csvinmem(self.fn_source)
        else:
            print kwargs
            iter_csv = util.read_csv(self.fn_source, iterator=True, chunksize=self.get_chunksize(),**kwargs)
            
        with open(self.fn_dest, "w") as myfile:
            has_output_header  = False
            count=0
            for chunk in iter_csv:
                count = count + 1
                print '%s %d' % (self.fn_source,count)
                #populate new data by chunk
                if not hasattr(self,'new_cols'):
                    self.new_cols= chunk.header()
                if self.filter_str is not None:
                    chunk= chunk.query(self.filter_str)                   
                chunk = chunk.fillna('');    
                output_df = pd.DataFrame(columns=self.new_cols)
                newRows = self.do_one_chunk(chunk)
                if len(newRows)>0 :
                    output_df = output_df.append(newRows,ignore_index=True)
                    output_df = output_df.reindex(columns=self.new_cols)

                    header = not has_output_header
                    if not has_output_header:
                        has_output_header = True
                    myfile.write(util.to_csv_string(output_df,index=False,header=header))

        
        if 'group_by' in self.options:
            #group_by='group_cols:col1,col2;aggregate_method:col1|count,col1|count_unique,col2|max,col3|min,col4|first,col4|last,col5|concat,col6|sum;' default aggregate_method is first.
            import re;
            m = re.search("group_cols:([^;]*;)",self.options['group_by'])
            if m:
                group_cols=m.groups()[0][:-1].split(',')
                m = re.search("aggregate_method:([^;]*;)",self.options['group_by'])
                aggr_methods = {};
                if m:
                    aggrs=m.groups()[0][:-1].split(',')
                    for aggr in aggrs:
                        aggr_methods[aggr.split('|')[0]] = aggr.split('|')[1]
                        
                df = util.read_csv(self.fn_dest);
                #Tracer()()
                for c in group_cols:
                    if c not in df.columns:
                        print 'Error in group_by options. Column %s does not exist'%c;
                        exist();
                        
                data=[];
                for k, row in df.groupby(group_cols):
                    cur_row={}
                    for i in range(len(group_cols)):
                        cur_row[group_cols[i]] = k[i]
                    #Tracer()()    
                    for c in df.columns:
                        if c not in group_cols:
                            if (c not in aggr_methods) or aggr_methods[c]=='first': #default is first
                                cur_row[c] = row[c].values[0];
                            elif aggr_methods[c]=='last':
                                cur_row[c] = row[c].values[len(row)-1];
                            elif aggr_methods[c]=='max':
                                cur_row[c] = np.max(row[c]);
                            elif aggr_methods[c]=='min':
                                cur_row[c] = np.min(row[c]);
                            elif aggr_methods[c]=='sum':
                                cur_row[c] = np.sum(row[c]);
                            elif aggr_methods[c]=='count':
                                cur_row[c] = len(row[c]);
                            elif aggr_methods[c]=='count_unique':
                                cur_row[c] = len(np.unique(row[c]));
                            elif aggr_methods[c]=='concat':
                                cur_row[c] = ','.join(pd.DataFrame(np.unique(row[c]))[0].map(str));
                            else: #default is first                                
                                print 'Warning: %s is not supported as aggregate_method. Default is used'%aggr_methods[c]
                                cur_row[c] = row[c].values[0];
                                
                    data.append(cur_row);    
                pd.DataFrame(data).to_csv(self.fn_dest, index=False);
            else:        
                print 'Error in group_by options. Group columns is not specified';
                exit();
            
        self.do_post_update()                    

    def do_post_update(self):
        pass

    def do_post_chunk(self):
        pass


    def do_one_chunk(self,chunk):
        newRows = []
        for index, row in chunk.iterrows():
            r = self.generate_new_row(row)
            if r is not None:
                newRows.append(r)
        return newRows

    def populate_more(self,root):
        XmlClass.populate_more(self,root)
        self.inputs.append(self.fn_source)

    def get_id(self):
        r = self.tag+":"
        if hasattr(self,'type_name') and  self.type_name is not None:
            r = r+ self.type_name
        elif hasattr(self,'dest') and  self.dest is not None:
            r = r + self.dest
        return r


class DownloadCsvConvert(CsvConvert):
    def __init__(self, xe):
        CsvConvert.__init__(self,xe=xe)
        if not 'newCols' in xe.attrib:
            self.new_cols = self.options['oldCols']

class UploadCsvConvert(CsvConvert):
    def __init__(self, xe,dest): # only be called by subclass
        CsvConvert.__init__(self,xe=xe)
        if 'typeName' in xe.attrib:
            self.type_name = xe.attrib['typeName']
        else:
            self.type_name = None # homologene doens't have typeName
        self.con = None
        self.type_col_value = None
        self.dest = dest

    def get_id(self):
        r = self.tag+":"
        if self.type_name is not None:
            r = r+ self.type_name
        else:
            r = r + self.dest
        return r

    def get_connection(self):
        if self.con is None:
            self.con = db.get_con(SyncDB.CONNECTION_ID, db=SyncDB.DATABASE)

        return self.con

    def get_type_col_value_sql(self):
        return 'SELECT 2 FROM term_category t WHERE t.category_name = ?'

    def generate_new_row(self,row):
        r = super(UploadCsvConvert,self).generate_new_row(row)
        if r is not None:
            r['ds'] = self.ds
            if self.type_name is not None:
                r[self.type_col] = self.get_type_col_value()
        return r

    def get_type_col_value(self):
        if self.type_col_value is None:
            con = self.get_connection()
            t = db.from_sql(con,self.get_type_col_value_sql(),
                        params=[self.type_name])
            self.type_col_value = t.ix[0,0].astype(str)
        return self.type_col_value

    def get_fn_dest(self):
        dir = os.path.join(SyncDB.UPLOAD_DIR(), self.dest)
        if not os.path.exists(dir):
            os.makedirs(dir)
        
        if hasattr(self, 'use_dest'):
            fn = os.path.join(dir, "%s.csv" % (self.use_dest))
        elif self.type_name is None: # for homologene, only one source, not many types
            fn = os.path.join(dir, "%s.csv" % (self.dest))
        else:
            fn = os.path.join(dir, "%s.csv" % (self.type_name))
        return fn
    
    def do_post_update(self):
        #return;
        if util.lines_in_file(self.fn_dest) > 0:
            lines_seen = set() # holds lines already seen
            unique_lines = [];
            for line in open(self.fn_dest, "r"):
                if line not in lines_seen: # not a duplicate
                    unique_lines.append(line)
                    lines_seen.add(line)
                    
            outfile = open(self.fn_dest, 'w')
            outfile.writelines(unique_lines)
            outfile.close()
            
            #df = util.read_csv(self.fn_dest,dtype=str)
            #df.drop_duplicates(inplace=True)
            #df.to_csv(self.fn_dest, index=False)

class SyncDB(object):
    #def __init__(self):
        #self.tasks=[]
    DEBUG=True
    __HOME_DIR=os.path.dirname(os.path.realpath(__file__)) 
    __TAX_ID_LIST = {9606,10090,10116}
    IGNORE_HOURS_DOWNLOAD = 48000 # < 24 hour will not downaload again
    ALL_TASKS = None
    CSV_IN_MEMS = {}
    TO_SYNC_DS = {}
    IS_REBUILD = True # flase means sync , if file existing then file is good.
    DEBUG_NODE_IDS=[]
    DEBUG_NODE_IDS_ONLY=False
    DATABASE='gp'
    APP_SETTINGS=None
    
    
    @staticmethod
    def GET_APP_SETTIGNS(setting_name):
        if SyncDB.APP_SETTINGS is None:
            df=util.read_csv(os.path.join(SyncDB.HOME_DIR(),"application_settings.csv"))
            SyncDB.APP_SETTINGS = {df.at[i, "setting_name"]:df.at[i, "setting_value"] for i in df.index}
            SyncDB.APP_SETTINGS['server'] = 'METASCAPE' 

        return SyncDB.APP_SETTINGS[setting_name] if setting_name in SyncDB.APP_SETTINGS else None;
            

    @staticmethod
    def TAX_ID_LIST():
        return SyncDB.__TAX_ID_LIST

    @staticmethod
    def HOME_DIR():
        return SyncDB.__HOME_DIR

    @staticmethod
    def set_HOME_DIR(home_dir):
        if home_dir:
            SyncDB.__HOME_DIR=home_dir

    @staticmethod
    def LIB_DIR():
        return os.path.join(SyncDB.HOME_DIR(),"task_class")

    @staticmethod
    def UPLOAD_DIR():
        dir = os.path.join(SyncDB.HOME_DIR(),"upload")
        if not os.path.exists(dir):
                os.makedirs(dir)
        return dir

    @staticmethod
    def DOWNLOAD_DIR():
        dir = os.path.join(SyncDB.HOME_DIR(),"download")
        if not os.path.exists(dir):
                os.makedirs(dir)
        return dir

    @staticmethod
    def load_xml_file(s_file=None):
        if SyncDB.DEBUG:
            util.debug_msg('Loading config file: '+s_file)
        if s_file is None:
            s_file=os.path.join(SyncDB.HOME_DIR(),"config.xml")
        tree=ET.parse(s_file)
        root=tree.getroot()
        cr = SyncDB.task_factory(root)
        SyncDB.ALL_TASKS = cr
        return cr
            
    @staticmethod
    def task_factory(xe):
        class_name = xe.tag
        if class_name in globals():
            task_class = globals()[class_name]
        else:
            module_name =  class_name.lower()
            if os.path.exists(os.path.join(SyncDB.LIB_DIR(),module_name+".py")):
                m=__import__('task_class.'+module_name, globals(), locals(), ['*'], -1)
                task_class=m.__dict__[class_name]
                globals()[class_name]=task_class
            else:
                print util.error_msg('Module not found: '+class_name)
        return task_class(xe=xe)
    @staticmethod
    def get_csvinmem(fn):
        if fn not in SyncDB.CSV_IN_MEMS:
            df = util.read_csv(fn,dtype=str)
            SyncDB.CSV_IN_MEMS[fn] = [df]
        return SyncDB.CSV_IN_MEMS.get(fn)

    @staticmethod
    def history():
        con = db.get_con(SyncDB.CONNECTION_ID, db=SyncDB.DATABASE)
        query = "DELETE FROM {0}.statistics where history = CURDATE()"
        query = query.format(SyncDB.DATABASE)
        db.from_sql(con,query)
        
        query = """
            INSERT INTO {0}.statistics
            SELECT a.*, CURDATE() AS history
            FROM(
            SELECT 'gid2source_id' AS table_name,  it.display_name  AS type_name, it.id_type_id as type_id, it.ds_name as ds, gs.tax_id, COUNT(*) AS total FROM {0}.gid2source_id gs , {0}.id_type it
            WHERE gs.id_type_id = it.id_type_id
            GROUP BY gs.id_type_id, gs.tax_id
            UNION ALL
            SELECT 'annotation' AS table_name, a_t.display_name AS type_name, a_t.annotation_type_id as type_id, a_t.ds_name as ds,  a.tax_id, COUNT(*) AS total
            FROM {0}.annotation a, {0}.annotation_type a_t
            WHERE a.annotation_type_id = a_t.annotation_type_id
            GROUP BY a.annotation_type_id, a.tax_id
            UNION ALL
            SELECT 'gid2terms' AS table_name, tc.category_name AS type_name, tc.term_category_id as type_id, tc.ds_name as ds, gt.tax_id , COUNT(*)
            FROM {0}.gid2terms gt, {0}.term_category tc
            where tc.term_category_id = gt.term_category_id 
            GROUP BY tc.term_category_id, gt.tax_id
            UNION ALL
            SELECT 'homologene' AS table_name, 'Homologene' AS type_name, 1 as type_id, 'NCBI' as ds, hg.tax_id,  COUNT(*) AS total FROM {0}.homologene hg group by hg.tax_id
            UNION ALL
            SELECT 'term' AS table_name, tc.category_name AS type_name, tc.term_category_id as type_id, tc.ds_name as ds, NULL as tax_id, COUNT(*) AS total FROM {0}.term t, {0}.term_category tc
            WHERE t.term_category_id = tc.term_category_id
            GROUP BY t.term_category_id
            UNION ALL
            SELECT 'term2gids' AS table_name, tc.category_name AS type_name, tc.term_category_id as type_id, tc.ds_name as ds, gt.tax_id, COUNT(*) as total FROM {0}.term2gids gt, {0}.term_category tc
            where tc.term_category_id = gt.term_category_id 
            GROUP BY tc.term_category_id, gt.tax_id
            UNION ALL
            SELECT 'term2term' AS table_name, 'Term relations' AS type_name, tc.term_category_id as type_id, tc.ds_name as ds, NULL as tax_id, COUNT(*) AS total FROM {0}.term2term tt, {0}.term_category tc
            WHERE tt.term_category_id = tc.term_category_id
            GROUP BY tt.term_category_id
            UNION ALL
            SELECT 'interaction' AS table_name, i_t.interaction_type_name AS type_name, i_t.interaction_type_id as type_id, i_t.ds_name as ds, i.tax_id_A as tax_id, COUNT(*) AS total FROM {0}.interaction i, {0}.interaction_type i_t
            WHERE i.interaction_type_id = i_t.interaction_type_id
            GROUP BY i.interaction_type_id, i.tax_id_A
            )a order by a.table_name, a.ds, a.type_name, a.tax_id
        """
        query = query.format(SyncDB.DATABASE)
        db.from_sql(con,query)

    @staticmethod
    def report_html(rpath):
        import re
        con = db.get_con(SyncDB.CONNECTION_ID, db=SyncDB.DATABASE)
        #get last update date;
        dt=db.from_sql(con,"Select distinct history FROM {0}.statistics order by history desc".format(SyncDB.DATABASE))
        if (len(dt['history']) < 2):
            print "No previous build statistics found."
            return;

        cur_date=str(dt['history'][0])            
        last_date=str(dt['history'][1])
        query = """
            select t.* from (
            SELECT new.ds as data_source, new.table_name, new.type_name, new.total as size, (new.total - old.total) as growth, (new.total - old.total)/old.total *100 AS delta, 0 as new_missing
            FROM {2}.statistics old, {2}.statistics new
            WHERE old.history = '{0}'
            AND new.history = '{1}'
            AND old.table_name = new.table_name
            AND old.type_name = new.type_name
            UNION
            SELECT new.ds as data_source, new.table_name, new.type_name, new.total as size, new.total as growth, 101 AS delta, 1 as new_missing
            FROM {2}.statistics new left Join (select * from {2}.statistics where history = '{0}') old on old.table_name = new.table_name and old.type_name = new.type_name
            WHERE new.history = '{1}' 
            AND old.table_name is NULL
            UNION
            SELECT old.ds as data_source, old.table_name, old.type_name, old.total as size, -old.total as growth, -101 AS delta, -1 as new_missing
            FROM {2}.statistics old left Join (select * from {2}.statistics where history = '{1}') new on old.table_name = new.table_name and old.type_name = new.type_name
            WHERE old.history = '{0}'
            AND old.table_name is NULL
            ) t
            ORDER BY data_source, table_name, delta DESC;
        """
        query = query.format(last_date, cur_date, SyncDB.DATABASE)
        dt=db.from_sql(con,query)
        #Tracer()()
        missing_data = dt.query('new_missing < 0')
        new_data = dt.query('new_missing > 0')
        expanded_data = dt.query('growth > 0')
        reduced_data = dt.query('growth < 0')
        unchanged_data = dt.query('growth == 0 and new_missing==0')
        
       

        reports = [r for r in sorted(os.listdir(rpath)) if "Report_" in r and r != "Report_"+cur_date+".html"]
        last_report = reports[len(reports)-1] if len(reports) > 0 else None
      
        #import datetime
        #time_now = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        html='<!DOCTYPE html><html><head><title> GP Build Report ' + cur_date + '</title></head><body><h1> GP Build Report ' + cur_date + '</h1>'
        if last_report:
            html += '<span id="lastreportlink"><a href="' + last_report + '"><<</a></span>\n'
        html += '<span id="nextreportlink"></span>\n'            
        html += '<div>'
        if (len(missing_data)>0):
            html += '<h4> Missing Data: </h4><table><tr><th>Data Source</th><th>Table Name</th><th>Type Name</th><th>Size</th></tr>'
        for i in range (len(missing_data)):
            html += '<tr><td>' + missing_data.irow(i)['data_source'] + '</td><td>' + missing_data.irow(i)['table_name'] + '</td><td>' + missing_data.irow(i)['type_name'] + '</td><td>' + str(missing_data.irow(i)['growth']) + '</td></tr>';
        if (len(missing_data)>0):
            html += '</table>'
        html += '</div><div>'
        
        if (len(new_data)>0):
            html += '<h4> New Data: </h4><table><tr><th>Data Source</th><th>Table Name</th><th>Type Name</th><th>Size</th></tr>'
        for i in range (len(new_data)):
            html += '<tr><td>' + new_data.irow(i)['data_source'] +'</td><td>' + new_data.irow(i)['table_name'] + '</td><td>' + new_data.irow(i)['type_name'] + '</td><td>' + str(new_data.irow(i)['growth']) + '</td></tr>';
        if (len(new_data)>0):
            html += '</table>'
        html += '</div><div>'
        if (len(expanded_data)>0):
            html += '<h4> Expanded Data: </h4><table><tr><th>Data Source</th><th>Table Name</th><th>Type Name</th><th>Size Diff</th> <th>Ratio</th></tr>'
        for i in range (len(expanded_data)):
            html += '<tr><td>' + expanded_data.irow(i)['data_source'] +  '</td><td>' + expanded_data.irow(i)['table_name'] + '</td><td>' + expanded_data.irow(i)['type_name'] + '</td><td>' + str(expanded_data.irow(i)['growth']) + '</td><td>' + str(expanded_data.irow(i)['delta']) + '</td></tr>';
        if (len(expanded_data)>0):
            html += '</table>'
        html += '</div><div>'
        if (len(reduced_data)>0):
            html += '<h4> Reduced Data: </h4><table><tr><th>Data Source</th><th>Table Name</th><th>Type Name</th><th>Size Diff</th> <th>Ratio</th></tr>'
        for i in range (len(reduced_data)):
            html += '<tr><td>' + reduced_data.irow(i)['data_source']  + '</td><td>' + reduced_data.irow(i)['table_name'] + '</td><td>' + reduced_data.irow(i)['type_name'] + '</td><td>' + str(reduced_data.irow(i)['growth']) + '</td><td>' + str(reduced_data.irow(i)['delta']) + '</td></tr>';
        if (len(reduced_data)>0):
            html += '</table>'
        html += '</div><div>'
        if (len(unchanged_data)>0):
            html += '<h4> Unchanged Data: </h4><table><tr><th>Data Source</th><th>Table Name</th><th>Type Name</th><th>Size</th> </tr>'
        for i in range (len(unchanged_data)):
            html += '<tr><td>' + unchanged_data.irow(i)['data_source']  + '</td><td>' + unchanged_data.irow(i)['table_name'] + '</td><td>' +  unchanged_data.irow(i)['type_name'] + '</td><td>' + str(unchanged_data.irow(i)['size']) + '</td></tr>';
        if (len(unchanged_data)>0):
            html += '</table>'
            
        html += '</div>'
        html += '</body></html>'

        with open(rpath+"Report.html", "w") as report_file:
            report_file.write(html)
            
        files=util.unix('cp ' + rpath + 'Report.html ' + rpath + 'Report_' + cur_date + '.html')
        
        if last_report:
            with open(rpath+last_report, "r") as last_report_file:
                last_html=last_report_file.read()
                last_html = re.sub(r"<span id=\"nextreportlink\">.*</span>", '<span id="nextreportlink"><a href="' + 'Report_' + cur_date + '.html">>></a></span>', last_html)

            with open(rpath+last_report, "w") as last_report_file:
                last_report_file.write(last_html)
     
    @staticmethod
    def report_js(rpath):
        import re
        con = db.get_con(SyncDB.CONNECTION_ID, db=SyncDB.DATABASE)
        #get last update date;
        dt=db.from_sql(con,"Select * FROM {0}.statistics order by history desc".format(SyncDB.DATABASE))
        dt['history'] = dt['history'].apply(str);
        cur_date = max(dt['history'])
        json_arr = []
        for c in dt.columns:
            json_arr.append('"' + c + '":' + dt[c].to_json(orient='values'))
                
        with open(rpath+"gp_stats.js", "w") as report_file:
            report_file.write("window.buildLogStatistics=")
            report_file.write(dt.to_json(orient='records'))
            
        files=util.unix('cp ' + rpath + 'gp_stats.js ' + rpath + 'gp_stats_' + cur_date + '.js')

    @staticmethod
    def stat_js(rpath):
        import re
        con = db.get_con(SyncDB.CONNECTION_ID, db=SyncDB.DATABASE)
        #get last update date;
        dt=db.from_sql(con,"Select * FROM {0}.statistics order by history desc".format(SyncDB.DATABASE))
        dt['history'] = dt['history'].apply(str);
        latest_build_date = max(dt['history'])
        #Tracer()()
        dt = dt.query("history=='" + latest_build_date + "'")
        jquery_str = "window.descriptionStatistics ={'latest_build':'" + latest_build_date + "','ds_counts':{\n"

        totals = [];
        for table in [['i_ann_','annotation'], ['i_gene_','gid2terms'], ['i_term_', 'term2gids'], ['i_interaction_', 'interaction']]:
            dt_type = dt.query("table_name=='" + table[1] + "'")
            for k, g in dt_type.groupby(by=['type_id']):
                total_in_one_type = g['total'].sum()
                totals.append("'" + table[0] + str(k) + "':" + str(total_in_one_type))
        
        jquery_str += ',\n'.join(totals) + "}}"
                       
        with open(rpath+"gp_description_stats.js", "w") as report_file:
            report_file.write(jquery_str)
           
        files=util.unix('cp ' + rpath + 'gp_description_stats.js ' + rpath + 'gp_description_stats_' + latest_build_date + '.js')
        
    @staticmethod
    def report(rpath):
        SyncDB.stat_js(rpath);    
        SyncDB.report_js(rpath);
        SyncDB.report_html(rpath);
        

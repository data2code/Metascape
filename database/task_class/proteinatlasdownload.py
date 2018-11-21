#!\usr\bin\env python

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
import csv
from gputil import GPUtils
from lxml import etree

from IPython.core.debugger import Tracer

class ProteinAtlasDownload(XmlClass):

    def __init__(self, xe=None):
        XmlClass.__init__(self,xe=xe)
        self.tag = "ProteinAtlasDownload"
        self.fn_dest=os.path.join(SyncDB.DOWNLOAD_DIR(),xe.attrib['dest'])
        self.protein_function_classes = [
            "Enzymes",
            
            "Enzymes\ENZYME proteins",
            
            "Enzymes\ENZYME proteins\Oxidoreductases",
            "Enzymes\ENZYME proteins\Transferases",
            "Enzymes\ENZYME proteins\Hydrolases",
            "Enzymes\ENZYME proteins\Lyases",
            "Enzymes\ENZYME proteins\Isomerase",
            "Enzymes\ENZYME proteins\Ligase",
            
            "Enzymes\Peptidases",
            
            "Enzymes\Peptidases\Aspartic-type peptidases",
            "Enzymes\Peptidases\Cysteine-type peptidases",            
            "Enzymes\Peptidases\Metallopeptidases",
            "Enzymes\Peptidases\Serine-type peptidases",
            "Enzymes\Peptidases\Threonine-type peptidases",            
            
            "Enzymes\Kinases",
            
            "Enzymes\Kinases\AGC Ser/Thr protein kinases",            
            "Enzymes\Kinases\Tyr protein kinases",            
            "Enzymes\Kinases\TKL Ser/Thr protein kinases",                        
            "Enzymes\Kinases\STE Ser/Thr protein kinases",            
            "Enzymes\Kinases\RGC receptor guanylate cyclase kinases",            
            "Enzymes\Kinases\NEK Ser/Thr protein kinases",                        
            "Enzymes\Kinases\CMGC Ser/Thr protein kinases",            
            "Enzymes\Kinases\CK1 Ser/Thr protein kinases",            
            "Enzymes\Kinases\CAMK Ser/Thr protein kinases",                        
            "Enzymes\Kinases\Atypical kinases",        

            "CD markers",
            "Blood group antigen proteins",
            "Nuclear receptors",
            "Transporters",
            
            "Transporters\Transporter channels and pores",
            "Transporters\Electrochemical Potential-driven transporters",
            "Transporters\Primary Active Transporters",
            "Transporters\Transport Electron Carriers",
            "Transporters\Accessory Factors Involved in Transport",
            
            "Ribosomal proteins",
            
            "G-protein coupled receptors",
            
            "G-protein coupled receptors\GPCRs excl olfactory receptors",
            "G-protein coupled receptors\Adenosine and adenine nucleotide receptors",
            "G-protein coupled receptors\Chemokines and chemotactic factors receptors",            
            "G-protein coupled receptors\Lysolipids receptors",
            "G-protein coupled receptors\Odorant/olfactory and gustatory receptors",
            "G-protein coupled receptors\Opsins",            
            "G-protein coupled receptors\Serotonin receptors",
            "G-protein coupled receptors\Family 2 (B) receptors",
            "G-protein coupled receptors\Family T2R receptors (taste receptor GPCRs)",            
            "G-protein coupled receptors\Family fz/smo receptors",
            
            "Voltage-gated ion channels",

            "Voltage-gated ion channels\Calcium-Activated Potassium Channels",
            "Voltage-gated ion channels\CatSper and Two-Pore Channels",
            "Voltage-gated ion channels\Cyclic Nucleotide-Regulated Channels",
            "Voltage-gated ion channels\Inwardly Rectifying Potassium Channels",
            "Voltage-gated ion channels\Transient Receptor Potential Channels",
            "Voltage-gated ion channels\Two-P Potassium Channels",
            "Voltage-gated ion channels\Voltage-Gated Calcium Channels",
            "Voltage-gated ion channels\Voltage-Gated Potassium Channels",
            "Voltage-gated ion channels\Voltage-Gated Sodium Channels",

            "Transcription factors",
            
            "Transcription factors\Yet undefined DNA-binding domains",
            "Transcription factors\Basic domains",
            "Transcription factors\Zinc-coordinating DNA-binding domains",            
            "Transcription factors\Helix-turn-helix domains",
            "Transcription factors\Other all-alpha-helical DNA-binding domains",
            "Transcription factors\alpha-Helices exposed by beta-structures",            
            "Transcription factors\Immunoglobulin fold",
            "Transcription factors\beta-Hairpin exposed by an alpha/beta-scaffold",
            "Transcription factors\beta-Sheet binding to DNA",            
            "Transcription factors\beta-Barrel DNA-binding domains",

            "Mitochondrial proteins",
            "RNA polymerase related proteins",
            "RAS pathway related proteins",
            "Citric acid cycle related proteins",
            "Cytoskeleton related proteins"
        ]
           
        self.protein_location_classes = [ 
            "Predicted membrane proteins",
            "Predicted secreted proteins",
            "Plasma proteins",
        ]
        
    def populate_more(self,root):
        self.outputs = [self.fn_dest]
       
    def prepare (self):
        pass;

    def longest_common_prefix (self, s1, s2):
        lcp = "";
        for k in range(min(len(s1), len(s2))):
            if (s1[k] != s2[k]):
                break;
            else:
                lcp += s1[k];

        import re;
        if not lcp.endswith('/') or len(re.findall('/',lcp)) < 1:
            return "";
        
        return lcp
        
    def colapse_protein_class (self, protein_classes):
        protein_classes.sort()
        out = []
        for i in range(len(protein_classes)):
            if (i==0):
                continue;
            if (not protein_classes[i].startswith(protein_classes[i-1])):
                out.append(protein_classes[i-1]);
        
        out.append(protein_classes[-1]);
        out2 = []
        i=0;
        while (i<len(out)):
            lcp = out[i];
            j = i + 1;
            while (j < len(out)):
                next_lcp = self.longest_common_prefix (lcp, out[j]);
                if (len(next_lcp)==0):
                    break;
                lcp = next_lcp;    
                j += 1;            
            
            #collapse between i and j;
            collapsed = [];            
            for k in range(i, j):
                collapsed.append(out[k][len(lcp):len(out[k])])
                
            if len(collapsed)==1:
                out2.append(lcp)
            else:
                out2.append(lcp + "{" + ",".join(collapsed) + "}");
                
            i = j;
            
        return out2;

    def get_protein_expression(self, type_name):
        self.get_ensembl2gid_map=GPUtils.get_ensembl2gid_map();  
        #"Gene","Tissue","Cell type","Level","Expression type","Reliability"        
        if type_name == 'Protein_Atlas_Normal_Tissue':
            data_file = os.path.join(SyncDB.DOWNLOAD_DIR(),"protein_atlas/normal_tissue.csv.zip")
            urllib.urlretrieve("http://www.proteinatlas.org/download/normal_tissue.csv.zip", data_file)
        else:
            #"Gene","Tumor","Level","Count patients","Total patients","Expression type"        
            data_file = os.path.join(SyncDB.DOWNLOAD_DIR(),"protein_atlas/cancer.csv.zip")
            urllib.urlretrieve("http://www.proteinatlas.org/download/cancer.csv.zip", data_file)
            
        data = util.read_csv(data_file)
        allgenes = set();
        notfoundgenes =  set();

        out = []
        for index, row in data.iterrows():
            if row["Level"] == "Not detected":
                continue;              
            allgenes.add(row['Gene'])
        
            if row['Gene'] in self.get_ensembl2gid_map:
                content = None;                    
                if type_name == 'Protein_Atlas_Normal_Tissue' :
                    content = row["Tissue"] + "/" + row["Cell type"] + "(" + row["Level"] + ")"
                    out.append({'gid':self.get_ensembl2gid_map[row['Gene']], 'content':content, 'annotation_field1':row['Gene']})                 
                elif type_name == 'Protein_Atlas_Cancer_Tissue':
                    content = row["Tumor"] + "(" + row["Level"] + ")"
                    out.append({'gid':self.get_ensembl2gid_map[row['Gene']], 'tumor': row["Tumor"], 'level':row["Level"], 'count_patients':row["Count patients"], 'total_patients':row["Total patients"],'annotation_field1':row['Gene']})                 
                                
            else:
                notfoundgenes.add(row['Gene'])            
                        
        print len(notfoundgenes), '/', len(allgenes), ' gene symbols (', float(len(notfoundgenes))/len(allgenes) if len(allgenes) != 0 else 1, ') cannot be converted to gene ids in ' + data_file;
        
        data=[]

             
        if type_name == 'Protein_Atlas_Cancer_Tissue':
            data0=[]
            for k, g in pd.DataFrame(out).groupby(['gid','tumor']):   
                #Tracer()()            
                if sum([g["count_patients"][i] for i in g.index]) == 0:
                    continue;
                    
                S = [g["level"][i] + ":" + str(g["count_patients"][i]) for i in g.index]
                S.append("ALL:" + str(g["total_patients"].values[0]))
                data0.append({'gid':k[0], 'content':k[1] + "(" + "|".join(S) + ")" , 'type_name': type_name, 'annotation_field1':g['annotation_field1'].tolist()[0]})                
            out = data0;

        for k, g in pd.DataFrame(out).groupby(['gid']):   
            data.append({'gid':k, 'content':";".join(g['content'].tolist()), 'type_name': type_name, 'annotation_field1':g['annotation_field1'].tolist()[0]})
                
        return data;

    def is_cellline (self, sample):
        for i, c in enumerate (sample):
            if c.isupper():
                return True;
        return False;
        
    def get_rna_expression(self):
        self.get_ensembl2gid_map=GPUtils.get_ensembl2gid_map();  
        #"Gene","Tissue","Cell type","Level","Expression type","Reliability"        
        data_file = os.path.join(SyncDB.DOWNLOAD_DIR(),"protein_atlas/rna.csv.zip")
        urllib.urlretrieve("http://www.proteinatlas.org/download/rna.csv.zip", data_file)
            
        data = util.read_csv(data_file)
        allgenes = set();
        notfoundgenes =  set();

        out = []
        for index, row in data.iterrows():
            if row["Abundance"] == "Not detected":
                continue;
            allgenes.add(row['Gene'])
            if row['Gene'] in self.get_ensembl2gid_map:                
                if self.is_cellline(row["Sample"]):
                    continue;
                    
                content = None;                    
                content = row["Sample"] + "(" + row["Abundance"] + "|" + row["Unit"] + ":" + str(row["Value"]) + ")"
                for gid in self.get_ensembl2gid_map[row['Gene']]:
                    out.append({'gid':gid, 'content':content, 'annotation_field1':row['Gene']})                                                
            else:
                #print row['Gene']
                notfoundgenes.add(row['Gene'])
        
        print len(notfoundgenes), '/', len(allgenes), ' gene symbols (', float(len(notfoundgenes))/len(allgenes) if len(allgenes) != 0 else 1, ') cannot be converted to gene ids in ' + data_file;

        data=[]
           
        for k, g in pd.DataFrame(out).groupby(['gid']):   
            data.append({'gid':k, 'content':";".join(g['content'].tolist()), 'type_name': 'Protein_Atlas_RNA', 'annotation_field1':g['annotation_field1'].tolist()[0]})
                
        return data;
        
    def get_subcellular_location(self):
        self.get_ensembl2gid_map=GPUtils.get_ensembl2gid_map();       
        data_file = os.path.join(SyncDB.DOWNLOAD_DIR(),"protein_atlas/subcellular_location.tsv.zip")
        urllib.urlretrieve("http://www.proteinatlas.org/download/subcellular_location.tsv.zip", data_file)
        data = util.read_csv(data_file,  sep="\t")
        
        allgenes = set();
        notfoundgenes =  set();
        out = []
        for index, row in data.iterrows():
            allgenes.add(row['Gene'])
            if row['Gene'] == 'ENSG00000090861':
                d = 9
            if row['Gene'] in self.get_ensembl2gid_map:
                ml = ''
                if not pd.isnull(row['Enhanced']):
                    ml += row['Enhanced'] + '(Enhanced)'
                if not pd.isnull(row['Supported']):
                    ml += row['Supported'] + '(Supported)'
                if not pd.isnull(row['Approved']):
                    ml += row['Approved'] + '(Approved)'
                if not pd.isnull(row['Uncertain']):
                    ml += row['Uncertain'] + '(Uncertain)'
                if len(ml) != 0:
                    for gid in self.get_ensembl2gid_map[row['Gene']]:
                        out.append({'gid':gid, 'content':ml, 'annotation_field1':row['Gene']})
            else:
                notfoundgenes.add(row['Gene'])
                
        print len(notfoundgenes), '/', len(allgenes), ' gene symbols (', float(len(notfoundgenes))/len(allgenes) if len(allgenes) != 0 else 1, ') cannot be converted to gene ids in ' + data_file;
        #
        data=[]
        for k, g in pd.DataFrame(out).groupby(['gid']):   
            #Tracer()()        
            data.append({'gid':k, 'content':";".join(g['content'].tolist()), 'type_name': "Protein_Atlas_Subcellular_Location", 'annotation_field1':g['annotation_field1'].tolist()[0]})               
       
        return data;

    def get_protein_class_data (self, id_map, id_type):
        out = []
        i = 0;
        for prot_class in self.protein_function_classes + self.protein_location_classes:
            pclass = prot_class.split("\\")[-1];
            cur_data_file = os.path.join(SyncDB.DOWNLOAD_DIR(), "protein_atlas/" + pclass.replace('/', '_')+".tab.gz");
            urllib.urlretrieve("http://www.proteinatlas.org/search/protein_class:" + pclass + "?format=tab", cur_data_file)
            data = util.read_csv(cur_data_file, sep='\t')
            allgenes = set()
            notfoundgenes =  set();

            for index, row in data.iterrows():
                allgenes.add(row[id_type])
                if row[id_type] in id_map:                                       
                    if prot_class in self.protein_function_classes:
                        type_name = 'Protein_Atlas_Function'
                    elif prot_class == "Predicted membrane proteins":
                        type_name = 'Protein_Atlas_Membrane'
                    elif prot_class == "Predicted secreted proteins":
                        type_name = 'Protein_Atlas_Secreted'
                    elif prot_class == "Plasma proteins":
                        type_name = 'Protein_Atlas_Plasma'
                    else:
                        print 'unknown protein class in protein atlas: ', prot_class
                        continue;
                    for gid in id_map[row[id_type]]:
                        out.append({'gid':gid, 'content':prot_class.replace('\\', '/'), 'type_name': type_name , 'annotation_field1':row['Gene']})                
                else:
                    notfoundgenes.add(row[id_type])
                    
            print len(notfoundgenes), '/', len(allgenes), ' gene symbols (', float(len(notfoundgenes))/len(allgenes) if len(allgenes) != 0 else 1, ') cannot be converted to gene ids in ' + cur_data_file;
            i +=1;
            #if (i > 4):
            #    break;          
        data=[]
        for k, g in pd.DataFrame(out).groupby(['gid','type_name']):   
            #Tracer()()        
            S=self.colapse_protein_class(util.unique([x for x in g['content'] if not pd.isnull(x)]))
            data.append({'gid':k[0], 'content':";".join(S), 'type_name':k[1], 'annotation_field1':g['annotation_field1'].tolist()[0]})               
        
        return data;
            

    def process_gene_entry(self, entry):
        row = {"protein_expression_tissue":"", "protein_expression_cancer":""};
        row ["Gene"] = entry.find('name').text;
        iden = entry.find('identifier')
        
        if iden is not None and iden.attrib['db']=="Ensembl":
            row ["Ensembl"] = iden.attrib["id"]
        else:
            row ["Ensembl"] = None
            
        TE=entry.findall('tissueExpression')        
        for te in TE:
            if te.find('summary') is None:
                continue;
            if "assayType" in te.attrib and te.attrib["assayType"] == 'tissue':
                row["protein_expression_tissue"] = te.find('summary').text;
            elif "assayType" in te.attrib and te.attrib["assayType"] == 'cancer':
                row["protein_expression_cancer"] = te.find('summary').text;

        antibody=entry.find('antibody')   
        
        if antibody is not None:
            TE=antibody.findall('tissueExpression') 
            for te in TE:
                if te.find('summary') is None:
                    continue;
                if "assayType" in te.attrib and te.attrib["assayType"] == 'tissue':
                    row["protein_expression_tissue"] = te.find('summary').text;
                elif "assayType" in te.attrib and te.attrib["assayType"] == 'cancer':
                    row["protein_expression_cancer"] = te.find('summary').text;
            
        return row;
        
    def get_summary_expression_data (self, id_map, id_type):
        #Tracer()()
        if not os.path.exists(os.path.join(SyncDB.DOWNLOAD_DIR(),"proteinatlas.xml")):
            urllib.urlretrieve("http://www.proteinatlas.org/download/proteinatlas.xml.gz", os.path.join(SyncDB.DOWNLOAD_DIR(),"proteinatlas.xml.gz"))
            util.unix("gunzip " + os.path.join(SyncDB.DOWNLOAD_DIR(),"proteinatlas.xml.gz"))
            
        context = etree.iterparse(os.path.join(SyncDB.DOWNLOAD_DIR(),"proteinatlas.xml"), tag='entry')
        out = [];
        i = 0;
        for action, elem in context:
            i = i + 1;
            if elem.tag=='entry':
                row = self.process_gene_entry (elem);
                if row[id_type] in id_map:
                    content = str(unicode(row["protein_expression_tissue"]).encode('utf-8'))
                    content = content if  content!='nan' else '';
                    for gid in id_map[row[id_type]]:
                        out.append({'gid':gid, 'content':content, 'type_name': "Protein_Atlas_PROTEIN_EXPR_TISSUE" , 'annotation_field1':row['Gene']});      
                    content = str(unicode(row["protein_expression_cancer"]).encode('utf-8'))
                    content = content if  content!='nan' else '';  
                    for gid in id_map[row[id_type]]:
                        out.append({'gid':gid, 'content':content, 'type_name': "Protein_Atlas_PROTEIN_EXPR_CANCER" , 'annotation_field1':row['Gene']});               
                # processing goes here
                pass

            elem.clear()

            while elem.getprevious() is not None:
                del elem.getparent()[0]
            if i%1000 == 0: print "processed " + str(i) + " genes";
        return out;
    
    def get_rna_tissue_category(self, id_map, id_type):
        if not os.path.exists(os.path.join(SyncDB.DOWNLOAD_DIR(),"proteinatlas.tsv.zip")):
            urllib.urlretrieve("http://www.proteinatlas.org/download/proteinatlas.tsv.zip", os.path.join(SyncDB.DOWNLOAD_DIR(),"proteinatlas.tsv.zip"))

        data = util.read_csv(os.path.join(SyncDB.DOWNLOAD_DIR(),"proteinatlas.tsv.zip"), sep='\t')
        allgenes = set();
        notfoundgenes =  set();
        out=[];
        
        for index, row in data.iterrows():
            allgenes.add(row[id_type])
            if row[id_type] in id_map:
                for gid in id_map[row[id_type]]:
                    out.append({'gid':gid, 'content':row['RNA tissue category'], 'type_name': 'Protein_Atlas_RNA_tissue_category', 'annotation_field1':row['Gene']})                
            else:
                notfoundgenes.add(row[id_type])

        print len(notfoundgenes), '/', len(allgenes), ' gene symbols (', float(len(notfoundgenes))/len(allgenes) if len(allgenes) != 0 else 1, ') cannot be converted to gene ids in proteinatlas.tab';
        return out;
        
    def do_update(self):
        self.prepare()
        if not os.path.exists(os.path.join(SyncDB.DOWNLOAD_DIR(),"protein_atlas")):            
            util.unix("mkdir " + os.path.join(SyncDB.DOWNLOAD_DIR(),"protein_atlas"))
        data = self.get_summary_expression_data(GPUtils.get_ensembl2gid_map(), "Ensembl") \
               + self.get_protein_class_data(GPUtils.get_ensembl2gid_map(), "Ensembl") \
               + self.get_subcellular_location() \
               + self.get_rna_tissue_category(GPUtils.get_ensembl2gid_map(), "Ensembl")
        
        dff=pd.DataFrame(data);
        reload(sys)  
        sys.setdefaultencoding('utf8')
        dff['content']=dff['content'].map(lambda x: unicode(x).encode('utf-8'))
        dff['content']=dff['content'].map(lambda x: x.replace('\r','').replace('\n','') if x!='nan' else '')
        dff['tax_id']='9606'
        dff.to_csv(self.fn_dest, index=False, sep=',');                        

    def check_inputs(self):
        passed = True
        urls=[]
        print "Checking Urls for proteinatlas.org ..."
        for prot_class in self.protein_function_classes + self.protein_location_classes:
            pclass = prot_class.split("\\")[-1];
            urls.append("http://www.proteinatlas.org/search/protein_class:" + pclass + "?format=tab")
            
        urls.append("http://www.proteinatlas.org/download/subcellular_location.tsv.zip")
        urls.append("http://www.proteinatlas.org/download/proteinatlas.xml.gz")
        urls.append("http://www.proteinatlas.org/download/proteinatlas.tsv.zip")
        
        #Tracer()()
        for url in urls:
            if not GPUtils.check_url(url):
                passed = False;                
                
        return passed
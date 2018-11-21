#!/usr/bin/env python
from __future__ import absolute_import
from __future__ import print_function
from . import Entity as el
import xml.etree.ElementTree as ET
import re

# NCBI's efetch does not support pcassay and pccompound database
# so XML are obtained from other means, e.g.,
# https://pubchem.ncbi.nlm.nih.gov/assay/assay.cgi?aid=651919&version=1.2&q=expdesc_xmldisplay
# However, Eutil.efetch() method has wrapp this details inside, so we don't need to worry about it

class FetchList(el.EntityList):

    def __init__(self,xml):
        el.EntityList.__init__(self, 'Record', xml)

    def attributes(self):
        return ['pccompound_id', 'smiles', 'InChIKey', 'name', 'date', 'synonym', 'mesh_synonym', 'mesh_moa']

    def to_list(self, S_attr=None):
        S_attr=['name','smiles','InChiKey'] if S_attr is None else S_attr
        if 'pccompound_id' not in S_attr:
            S_attr.append('pccompound_id')
        return el.EntityList.to_list(self, S_attr)

    def pccompound_id(self):
        return self._get_value('./RecordNumber')

    def name(self):
        out=self._get_nodes('./Section/Section')
        for i,x in enumerate(out):
            if x is not None:
                l_found=False
                for y in x:
                    if y.find('./TOCHeading').text=='Record Title':
                        out[i]=y.find('./Information/StringValue').text
                        l_found=True
                        break
                if not l_found: out[i]=None
        return out

    def smiles(self):
        out=self._get_nodes('./Section/Section/Section')
        for i,x in enumerate(out):
            if x is not None:
                l_found=False
                for y in x:
                    if y.find('./TOCHeading').text=='Canonical SMILES':
                        out[i]=y.find('./Information/StringValue').text
                        l_found=True
                        break
                if not l_found: out[i]=None
        return out

    def InChIKey(self):
        out=self._get_nodes('./Section/Section/Section')
        for i,x in enumerate(out):
            if x is not None:
                l_found=False
                for y in x:
                    if y.find('./TOCHeading').text=='InChI Key':
                        out[i]=y.find('./Information/StringValue').text
                        l_found=True
                        break
                if not l_found: out[i]=None
        return out

    def synonym(self):
        out=self._get_nodes('./Section/Section/Section')
        for i,x in enumerate(out):
            if x is not None:
                l_found=False
                for y in x:
                    if y.find('./TOCHeading').text=='Depositor-Supplied Synonyms':
                        out[i]=[ syn.text for syn in y.findall('./Information/StringValueList')]
                        l_found=True
                        break
                if not l_found: out[i]=None
        return out

    def mesh_synonym(self):
        out=self._get_nodes('./Section/Section/Section')
        for i,x in enumerate(out):
            if x is not None:
                l_found=False
                for y in x:
                    if y.find('./TOCHeading').text=='MeSH Synonyms':
                        out[i]=[ syn.text for syn in y.findall('./Information/StringValueList')]
                        l_found=True
                        break
                if not l_found: out[i]=None
        return out

    def mesh_moa(self):
        out=self._get_nodes('./Section/Section')
        for i,x in enumerate(out):
            if x is not None:
                l_found=False
                for y in x:
                    if y.find('./TOCHeading').text=='MeSH Pharmacological Classification':
                        out[i]=[z.text for z in y.findall('./Information/Name')]
                        l_found=True
                        break
                if not l_found: out[i]=None
        return out

    def date(self):
        S_year=[None]*len(self.data)
        S_month=[None]*len(self.data)
        S_day=[None]*len(self.data)
        out=self._get_nodes('./Section/Section')
        for i,x in enumerate(out):
            if x is not None:
                for y in x:
                    if y.find('./TOCHeading').text=='Create Date':
                        z=y.find('./Information/DateValue').text
                        s=re.match(r'(?P<year>\d{4})(-0?(?P<month>\d+)(-0?(?P<day>\d{1,2}))?)?', z)
                        S_year[i]=s.group('year')
                        S_month[i]=s.group('month')
                        S_day[i]=s.group('day')
        return {'year':S_year, 'month':S_month, 'day':S_day}

class SummaryList(el.EntityList):

    def __init__(self,xml):
        el.EntityList.__init__(self, 'DocSum', xml)

    def attributes(self):
        return ['pccompound_id', 'smiles', 'InChIKey', 'name', 'date', 'synonym', 'mesh_synonym', 'mesh_moa']

    def to_list(self, S_attr=None):
        S_attr=['name','smiles','InChiKey'] if S_attr is None else S_attr
        if 'pccompound_id' not in S_attr:
            S_attr.append('pccompound_id')
        return el.EntityList.to_list(self, S_attr)

    def pccompound_id(self):
        return self._get_value('./Id')

    def name(self):
        out=self._get_node("./Item[@Name='MeSHHeadingList']")
        for i,x in enumerate(out):
            if x is not None and len(x)>0:
                out[i]=x[0].text
        return out

    def date(self):
        S_date=self._get_value("./Item[@Name='CreateDate']")
        S_year=[None]*len(S_date)
        S_month=[None]*len(S_date)
        S_day=[None]*len(S_date)
        for i,x in enumerate(S_date):
            if x is None: continue
            y=re.match(r'(?P<year>\d{4})(/0?(?P<month>\d+)(/0?(?P<day>\d{1,2}))?)?', x)
            S_year[i]=y.group('year')
            S_month[i]=y.group('month')
            S_day[i]=y.group('day')
        return {'year':S_year, 'month':S_month, 'day':S_day}

    def synonym(self):
        return self._get_values("./Item[@Name='SynonymList']/Item[@Name='string']")

    def mesh_synonym(self):
        return self._get_values("./Item[@Name='MeSHTermList']/Item[@Name='string']")

    def mesh_moa(self):
        return self._get_values("./Item[@Name='PharmActionList']/Item[@Name='string']")

    def smiles(self):
        return self._get_value("./Item[@Name='CanonicalSmiles']")

    def InChIKey(self):
        return self._get_value("./Item[@Name='InChIKey']")

if __name__=="__main__":
    import eutils
    eu=eutils.EUtils()
    out=eu.efetch({'db':'pccompound', 'id':'123596,123540'})
    x=FetchList(out)
    print(">>>Fetch")
    print(x.to_list(x.attributes()))

    out=eu.esummary({'db':'pccompound', 'id':'123596,123540'})
    x=SummaryList(out)
    print(">>>Summary")
    print(x.to_list(x.attributes()))
    exit()


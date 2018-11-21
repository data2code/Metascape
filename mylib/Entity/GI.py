#!/usr/bin/env python
from __future__ import absolute_import
from __future__ import print_function
from . import Entity as el
import xml.etree.ElementTree as ET
import re

class FetchList(el.EntityList):

    def __init__(self,xml):
        el.EntityList.__init__(self, 'GBSeq', xml)

    def attributes(self):
        return ['gi', 'tax_id', 'type', 'gene', 'CDS', 'accession', 'description', 'sequence']

    def to_list(self, S_attr=None):
        S_attr=['tax_id','accession','gene'] if S_attr is None else S_attr
        if 'gi' not in S_attr:
            S_attr.append('gi')
        return el.EntityList.to_list(self, S_attr)

    def gi(self):
        out=self._get_values('./GBSeq_other-seqids/GBSeqid')
        for i,x in enumerate(out):
            if x is not None:
                S=[y.replace('gi|','') for y in x if y.startswith('gi|')]
                out[i]=S[0] if len(S)>0 else None
        return out

    def type(self):
        return self._get_value('./GBSeq_moltype')

    def tax_id(self):
        out=self._get_nodes('./GBSeq_feature-table/GBFeature')
        for i,x in enumerate(out):
            if x is not None:
                l_found=False
                for y in x:
                    if l_found: break
                    if y.find('./GBFeature_key').text=='source':
                        for z in y.findall('./GBFeature_quals/GBQualifier'):
                            if z.find('GBQualifier_name').text=='db_xref':
                                if z.find('GBQualifier_value').text.startswith('taxon'):
                                    out[i]=z.find('GBQualifier_value').text.replace('taxon:', '')
                                    l_found=True
                                    break
                if not l_found: out[i]=None
        return out

    def gene(self):
        out=self._get_nodes('./GBSeq_feature-table/GBFeature')
        for i,x in enumerate(out):
            if x is not None:
                l_found=False
                for y in x:
                    if l_found: break
                    if y.find('./GBFeature_key').text=='gene':
                        for z in y.findall('./GBFeature_quals/GBQualifier'):
                            if z.find('GBQualifier_name').text=='db_xref':
                                if z.find('GBQualifier_value').text.startswith('GeneID'):
                                    out[i]=z.find('GBQualifier_value').text.replace('GeneID:', '')
                                    l_found=True
                                    break
                if not l_found: out[i]=None
        return out

    def CDS(self):
        out=self._get_nodes('./GBSeq_feature-table/GBFeature')
        for i,x in enumerate(out):
            if x is not None:
                S=[]
                for y in x:
                    if y.find('./GBFeature_key').text=='CDS':
                        S.append(y.find('./GBFeature_location').text)
                out[i]=S if len(S)>0 else None
        return out

    def accession(self):
        return self._get_value('./GBSeq_accession-version')

    def description(self):
        return self._get_value('./GBSeq_definition')

    def sequence(self):
        return self._get_value('./GBSeq_sequence')

class SummaryList(el.EntityList):

    def __init__(self,xml):
        el.EntityList.__init__(self, 'DocSum', xml)

    def attributes(self):
        return ['gi', 'tax_id', 'extra', 'title', 'status', 'replace_by']

    def to_list(self, S_attr=None):
        S_attr=['tax_id','extra','title'] if S_attr is None else S_attr
        if 'gi' not in S_attr:
            S_attr.append('gi')
        return el.EntityList.to_list(self, S_attr)

    def gi(self):
        return self._get_value('./Id')

    def tax_id(self):
        return self._get_value("./Item[@Name='TaxId']")

    def extra(self):
        return self._get_value("./Item[@Name='Extra']")

    def title(self):
        return self._get_value("./Item[@Name='Title']")

    def status(self):
        return self._get_value("./Item[@Name='Status']")

    def replace_by(self):
        return self._get_value("./Item[@Name='ReplacedBy']")

if __name__=="__main__":
    import eutils
    eu=eutils.EUtils()
    out=eu.efetch({'db':'nucleotide', 'id':'67944638,578820969'})
    x=FetchList(out)
    print(">>>Fetch")
    print(x.to_list(x.attributes()))

    out=eu.esummary({'db':'nucleotide', 'id':'67944638,578820969'})
    x=SummaryList(out)
    print(">>>Summary")
    print(x.to_list(x.attributes()))


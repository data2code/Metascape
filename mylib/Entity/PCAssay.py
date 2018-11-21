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
        el.EntityList.__init__(self, 'PC-AssayDescription', xml)

    def attributes(self):
        return ['pcassay_id', 'source', 'date', 'title', 'description', 'gene_target', 'target']

    def to_list(self, S_attr=None):
        S_attr=['title','source','target'] if S_attr is None else S_attr
        if 'pcassay_id' not in S_attr:
            S_attr.append('pcassay_id')
        return el.EntityList.to_list(self, S_attr)

    def pcassay_id(self):
        return self._get_value('./PC-AssayDescription_aid/PC-ID/PC-ID_id')

    def source(self):
        return self._get_value('./PC-AssayDescription_aid-source/PC-Source/PC-Source_db/PC-DBTracking/PC-DBTracking_name')

    def date(self):
        S_year=self._get_value('./PC-AssayDescription_aid-source/PC-Source/PC-Source_db/PC-DBTracking/PC-DBTracking_date/Date/Date_std/Date-std/Date-std_year')
        S_month=self._get_value('./PC-AssayDescription_aid-source/PC-Source/PC-Source_db/PC-DBTracking/PC-DBTracking_date/Date/Date_std/Date-std/Date-std_month')
        S_day=self._get_value('./PC-AssayDescription_aid-source/PC-Source/PC-Source_db/PC-DBTracking/PC-DBTracking_date/Date/Date_std/Date-std/Date-std_day')
        return {'year':S_year, 'month':S_month, 'day':S_day}

    def title(self):
        return self._get_value('./PC-AssayDescription_name')

    def description(self):
        """The XML is problematic, as it does not distinguish different fields, that's probably why efetch does nto work, as the XML is not really well defined."""
        out=self._get_values('./PC-AssayDescription_description/PC-AssayDescription_description_E')
        out=[" ".join([y for y in X if y]) for X in out ]
        return out

    def gene_target(self):
        return self._get_value('./PC-AssayDescription_xref/PC-AnnotatedXRef/PC-AnnotatedXRef_xref/PC-XRefData/PC-XRefData_gene')

    def target(self):
        S_id=self._get_value('./PC-AssayDescription_target/PC-AssayTargetInfo/PC-AssayTargetInfo_mol-id')
        S_type=self._get_attrib('./PC-AssayDescription_target/PC-AssayTargetInfo/PC-AssayTargetInfo_molecule-type', 'value')
        return {'id':S_id, 'type':S_type}

class SummaryList(el.EntityList):

    def __init__(self,xml):
        el.EntityList.__init__(self, 'DocSum', xml)

    def get_uid(self):
        return 'pcassay_id';

    def attributes(self):
        return ['pcassay_id', 'source', 'date', 'title', 'description', 'target']

    def to_list(self, S_attr=None):
        S_attr=['title','source','target'] if S_attr is None else S_attr
        if 'pcassay_id' not in S_attr:
            S_attr.append('pcassay_id')
        return el.EntityList.to_list(self, S_attr)

    def pcassay_id(self):
        return self._get_value('./Id')

    def source(self):
        out=self._get_node("./Item[@Name='SourceNameList']")
        for i,x in enumerate(out):
            if x is not None:
                out[i]=[ y.text for y in x]
        return out

    def date(self):
        """Deposit Date is different from the published date obtained from esummary"""
        S_date=self._get_value("./Item[@Name='DepositDate']")
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

    def title(self):
        return self._get_value("./Item[@Name='AssayName']")

    def description(self):
        return self._get_value("./Item[@Name='AssayDescription']")

    def target(self):
        """May need to look for other target tag, such as GeneTargetList?"""
        out=self._get_nodes("./Item[@Name='ProteinTargetList']/Item[@Name='ProteinTarget']")
        S_type=[None]*len(out)
        S_id=[None]*len(out)
        for i,x in enumerate(out):
            if x is not None:
                S_id[i]=[y.find("./Item[@Name='GI']").text for y in x if y.find("./Item[@Name='GI']") is not None]
                S_type[i]=['protein']*len(S_id)
        return {'id':S_id, 'type':S_type}


if __name__=="__main__":
    import eutils
    eu=eutils.EUtils()
    out=eu.efetch({'db':'pcassay', 'id':'651919,651920'})
    x=FetchList(out)
    print(">>>Fetch")
    print(x.to_list(x.attributes()))

    out=eu.esummary({'db':'pcassay', 'id':'651919,651920'})
    x=SummaryList(out)
    print(">>>Summary")
    print(x.to_list(x.attributes()))


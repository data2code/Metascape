#!/usr/bin/env python
from __future__ import absolute_import
from __future__ import print_function
from . import Entity as el
import xml.etree.ElementTree as ET
import re

class FetchList(el.EntityList):

    def __init__(self,xml):
        el.EntityList.__init__(self, 'Entrezgene', xml)

    def attributes(self):
        return ['gene_id', 'type', 'tax_id', 'tax_name', 'tax_common', 'summary', 'symbol', 'description', 'chromosome', 'mim', 'ensembl']

    def to_list(self, S_attr=None):
        S_attr=['tax_id','symbol','description'] if S_attr is None else S_attr
        if 'gene_id' not in S_attr:
            S_attr.append('gene_id')
        return el.EntityList.to_list(self, S_attr)

    def gene_id(self):
        return self._get_value('./Entrezgene_track-info/Gene-track/Gene-track_geneid')

    def type(self):
        return self._get_attrib('./Entrezgene_type', 'value')

    def tax_name(self):
        return self._get_value('./Entrezgene_source/BioSource/BioSource_org/Org-ref/Org-ref_taxname')

    def tax_common(self):
        return self._get_value('./Entrezgene_source/BioSource/BioSource_org/Org-ref/Org-ref_common')

    def tax_id(self):
        return self._get_value('./Entrezgene_source/BioSource/BioSource_org/Org-ref/Org-ref_db/Dbtag/Dbtag_tag/Object-id/Object-id_id')

    def summary(self):
        return self._get_value('./Entrezgene_summary')

    def symbol(self):
        return self._get_value('./Entrezgene_gene/Gene-ref/Gene-ref_locus')

    def description(self):
        return self._get_value('./Entrezgene_gene/Gene-ref/Gene-ref_desc')

    def chromosome(self):
        return self._get_value('./Entrezgene_gene/Gene-ref/Gene-ref_maploc')

    def mim(self):
        out=self._get_nodes('./Entrezgene_gene/Gene-ref/Gene-ref_db/Dbtag')
        for i,x in enumerate(out):
            l_found=False
            for y in x:
                if y.find('./Dbtag_db').text == 'MIM':
                    out[i]=y.find('./Dbtag_tag/Object-id/Object-id_id').text
                    l_found=True
                    break
            if not l_found:
                out[i]=None
        return out

    def ensembl(self):
        out=self._get_nodes('./Entrezgene_gene/Gene-ref/Gene-ref_db/Dbtag')
        for i,x in enumerate(out):
            l_found=False
            for y in x:
                if y.find('./Dbtag_db').text == 'Ensembl':
                    out[i]=y.find('./Dbtag_tag/Object-id/Object-id_str').text
                    l_found=True
                    break
            if not l_found:
                out[i]=None
        return out

class SummaryList(el.EntityList):

    def __init__(self,xml):
        el.EntityList.__init__(self, './*/DocumentSummary', xml)

    def get_uid(self):
        return 'gene_id';
    def attributes(self):
        return ['gene_id', 'tax_id', 'summary', 'symbol', 'description', 'chromosome', 'mim']

    def to_list(self, S_attr=None):
        S_attr=['tax_id','symbol','description'] if S_attr is None else S_attr
        if 'gene_id' not in S_attr:
            S_attr.append('gene_id')
        return el.EntityList.to_list(self, S_attr)

    def gene_id(self):
        return self._get_attrib('.', 'uid')

    #def type(self):
    #    #Type seems missing on 2/13/2015
    #    #https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id=51284&retmode=xml
    #    return self._get_attrib('./Entrezgene_type', 'value')

    def tax_id(self):
        return self._get_value("./Organism/TaxID")

    def summary(self):
        return self._get_value("./Summary")

    def symbol(self):
        return self._get_value("./Name")

    def description(self):
        return self._get_value("./Description")

    def chromosome(self):
        return self._get_value("./Chromosome")

    def mim(self):
        out=self._get_nodes("./Mim/int")
        for i,x in enumerate(out):
            if x is None: continue
            out[i]=" ".join([y.text for y in x])
        return out

if __name__=="__main__":
    import eutils
    eu=eutils.EUtils()
    out=eu.efetch({'db':'gene', 'id':'170743,51284'})
    x=FetchList(out)
    print(">>>Fetch")
    print(x.to_list(x.attributes()+["summary","type"]))

    out=eu.esummary({'db':'gene', 'id':'170743,51284,466142'})

    x=SummaryList(out)
    print(">>>Summary")
    print(x.to_list(x.attributes()))


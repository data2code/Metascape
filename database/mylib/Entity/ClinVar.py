#!/usr/bin/env python
import Entity as el
import xml.etree.ElementTree as ET

class SummaryList(el.EntityList):

    def __init__(self,xml):
        el.EntityList.__init__(self, './*/DocumentSummary', xml)
        #name = 'DocumentSummary'
        #if type(xml) is str:
        #    out = [ET.fromstring(xml)]
        #else:
        #    out =[ET.fromstring(x) for x in xml]
        #
        #self.data=[]
        #for x in out:
        #    for y in x.iter(name):
        #        self.data.append(y)

    def to_list(self, S_attr=None):
        S_attr=['variant_id', 'genes', 'trait', 'variant'] if S_attr is None else S_attr
        return el.EntityList.to_list(self, S_attr)

    def attributes(self):
        return ['variant_id','variant','genes','trait']

    def genes(self):
        out = []
        for x in self.data:
            y = x.find('genes')
            datas = y.findall('gene')
            genes = [[entry.find('GeneID').text, entry.find('symbol').text] for entry in datas]
            # remove duplicates
            genes_dict = dict(genes)
            genes = [[symbol,id] for id,symbol in genes_dict.items()]
            out.append(genes)
        return out

    def variant_id(self):
        return [x.attrib['uid'] for x in self.data]

    def variant(self):
        out = []
        for x in self.data:
            variant = [next(x.iter('variant_type')).text, next(x.iter('variation_name')).text]
            variant.append(next(x.iter('clinical_significance')).find('description').text)
            out.append(variant)
        return out

    def trait(self):
        out = []
        for x in self.data:
            trait = {}
            for y in x.iter('trait'):
                 datas = y.find('trait_xrefs').findall('trait_xref')
                 sources = [entry.find('db_source').text + ":" + entry.find('db_id').text for entry in datas]
                 trait[y.find('trait_name').text] = sources
            out.append(trait)
        return out

if __name__=="__main__":
    import eutils
    eu=eutils.EUtils()

    out=eu.esummary({'db':'clinvar', 'id':'187763'})
    x=SummaryList(out)
    print ">>>Summary"
    print x.to_list(x.attributes())


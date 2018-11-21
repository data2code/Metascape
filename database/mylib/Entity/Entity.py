import xml.etree.ElementTree as ET
import re

class EntityList(object):
    def __init__(self, name, xml):
        '''name is the xml tag for the entity, e.g., Entrezgene for efetch gene db entries
        It takes xml text (or a list of xml, it is resulted from multiple efetch)
        and convert that into a list of XML nodes, each node represents one entity.'''
        if type(xml) is str:
            out = [ ET.fromstring(xml)]
        else: # a list of xml strings
            out =[ET.fromstring(x) for x in xml]
        self.data=[]
        for x in out:
            for y in x.findall(name):
                self.data.append(y)
        #self.data is a list of gene node objects

    def attributes(self):
        """Return the list of attributes, that can be used in to_list()"""
        return []

    def to_list(self, S_attr=None):
        """Return a list of dict, each dict represents an object, each object contains keys specified by S_attr.
        Generally, this is the only method we need to call to convert XML into a list of dict."""
        out=[{} for x in self.data]
        for x in S_attr:
            if getattr(self, x) is not None:
                S_value=getattr(self, x)()
                if type(S_value) is list:
                    # for gene description, it's a list
                    self._add_attr(out, x, S_value)
                elif type(S_value) is dict: # a dict where each value is a list
                    # for pubmed journal, it's a list of dict
                    for k,v in S_value.items():
                        self._add_attr(out, x+'.'+k, v)
            else:
                util.error_msg('Attribute %s not implemented!' % x)
        return out

    def _get_value(self, xpath):
        """Each object returns one value str, so it returns list(str)"""
        out=[]
        for x in self.data:
            X=x.find(xpath)
            if X is not None:
                out.append(X.text)
            else:
                out.append(None)
        return out

    def _get_values(self, xpath):
        """Each object returns one list(str), so it returns list(list(str))"""
        # e.g., PCassay description_E
        out=[]
        for x in self.data:
            X=x.findall(xpath)
            if X is not None:
                out.append([y.text for y in X])
            else:
                out.append(None)
        return out

    def _get_attrib(self, xpath, attr):
        """Each object returns one attr str, so it returns list(str)"""
        out=[]
        for x in self.data:
            X=x.find(xpath)
            if X is not None:
                out.append(X.attrib[attr])
            else:
                out.append(None)
        return out

    def _get_node(self, xpath):
        """Each object returns one node, returns list(node)"""
        out=[]
        for x in self.data:
            X=x.find(xpath)
            if X is not None:
                out.append(X)
            else:
                out.append(None)
        return out

    def _get_nodes(self, xpath):
        """Each object returns a list of nodes, returns list(list(node))"""
        out=[]
        for x in self.data:
            X=x.findall(xpath)
            if X is not None:
                out.append(X)
            else:
                out.append(None)
        return out

    def _add_attr(self, objs, property_name, S_value):
        """Used by to_list()"""
        if len(objs)!=len(S_value):
            util.error_msg('Attribute: %s, length of Objects and S_value do not match: %d vs %d!' %
                (property_name, len(objs), len(S_value)))
        for i,x in enumerate(objs):
            objs[i][property_name]=S_value[i]

    @staticmethod
    def factory(db,xml):
        if db=='pcassay':
            import PCAssay
            return PCAssay.SummaryList(xml)
        if db=='gene':
            import Gene
            return Gene.SummaryList(xml)
        if db=='omim':
            import Omim
            return Omim.SummaryList(xml)
        if db=='pubmed':
            import PubMed
            return PubMed.SummaryList(xml)
        if db=='mesh':
            import MeSH
            return MeSH.SummaryList(xml)



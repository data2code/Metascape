from core import *

class Settings(XmlClass):
    def __init__(self, xe):
        #add default instance to children
        XmlClass.__init__(self,xe=xe)

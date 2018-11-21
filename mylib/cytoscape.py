#!/usr/bin/env python
from __future__ import absolute_import
from __future__ import print_function
from py2cytoscape.data.cyrest_client import CyRestClient
from py2cytoscape.data.cynetwork import CyNetwork
import xgmml
import pandas as pd
import util
import traceback
from six.moves import range

### READ ME ###
# This package required running Cytoscape (with CyREST API installed, required Cytoscape 3.2 and Java 1.8)
# Cytoscape cannot be launched in headless mode
# So we first need to install Xvbf, Java1.8, Cytoscape 3.2.1
#
# Run: "Xvbf :99" open a display port at say 99
#      It seems we don't have to specify resolution, if need, use Xvfb :25 -screen 0 1024x768x16
###########################
# Cytoscape Installation:
###########################
# We cannot change Port from 1234 to something else via -R, need to be changed via GUI or
#            use configure file.
# To be checked: how to install cyREST app in Cytoscape without GUI??? (find out which folder)
#      Probably under
#        ~/CytoscapeConfiguration/3/apps/installed> ls -lt
# Port is set at
#    ~/CytoscapeConfiguration/cytoscape3.props
#    ### rest.port=1234
#
# Run: export DISPLAY=:99; nohup ./cytoscape.sh &
#

class CytoscapeClient:
    """This is a wrapper to the py2cytoscape wrapper. The local installation of Py2cytoscape has been modified for it to work"""
    def __init__(self, host="localhost", port="1234", version="v1"):
        self.cy = CyRestClient(ip=host, port=port, version=version)

    def gc(self):
        self.cy.gc()

    def free_memory(self):
        self.cy.gc()
        X=self.cy.status()
        return max(X['memoryStatus']['maxMemory']-X['memoryStatus']['usedMemory'], X['memoryStatus']['freeMemory'])

    def random_id(self, prefix="MY_"):
        import random
        s="ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"
        return prefix+''.join([s[random.randint(0, len(s)-1)] for i in range(10)])

    def cynet_from_xgmml(self, s_file, name=None, s_layout="force-directed", l_bundle=True, t_xy=None, l_digraph=False):
        """If t_xy is provided, s_layout is ignored
        return CyNetwork"""
        import os
        if not os.path.exists(s_file):
            import util
            util.warn_msg('Missing file: %s' % s_file)
            return None
        if not l_digraph:
            net=xgmml.Network(s_file, name=name)
        else:
            import digraph
            net=digraph.Digraph(s_file, name=name)
        if t_xy is not None:
            net.set_XY(t_xy)
            s_layout=None
        S=net.T_node.header()
        #print "---------------->", S
        #S=[x for x in S if x not in ('Gene','Name','Symbol')]
        #net.T_node.drop(S, axis=1, inplace=True)
        ##############################
        return self.cynet_from_network(net, name=name, s_layout=s_layout, l_bundle=l_bundle)

    def cynet_from_network(self, net, name=None, s_layout="force-directed", l_bundle=True):
        """Input net: xgmml.Network object
        return CyNetwork"""
        #print "====>>>>>>>>", net.T_node.header()
        data=net.to_json()
        if name is None:
            name=data['data']['name']
        #cynet=self.cy.create(suid=net)
        cynet=None
        for trial in range(3): # something Cytoscape errors, try up to three times
            try:
                if cynet is None:
                    cynet=self.cy.network.create(name=name, data=data)
            except:
                cynet=None
                print("Fail at trail %d" % (trial+1))
                print("JSON data>", data)
        if cynet is None:
            import sys
            util.warn_msg('Error during plot_go_network(): %s' % sys.exc_info()[0])
            print(traceback.format_exc())

        id=cynet.get_id()

        if s_layout is not None:
            #print "apply layout %s" % s_layout
            self.cy.layout.apply(name=s_layout, network=cynet)
        if l_bundle:
            self.cy.edgebundling.apply(cynet)
            #print "apply bundle %s" % s_layout
        self.cy.layout.fit(cynet)
        return cynet

    def delete_cynet(self, cynet=None):
        """None means delete all networks"""
        if cynet is None:
            self.cy.network.delete_all()
        else:
            self.cy.network.delete(cynet)

    def get_network(self, suid):
        """Get xgmml.Network with X,Y coordinates"""
        cynet=self.cy.network.create(suid=suid)
        net=xgmml.Network.from_json(cynet.to_json())
        net.T_node.drop([x for x in ['graphics_x','graphics_y','SUID','id_original'] if x in net.T_node.header() ], axis=1, inplace=True)
        t_xy=self.cynet_get_XY(cynet)
        if t_xy is not None:
            t_xy.drop('Gene', axis=1, inplace=True)
            t_xy.rename2({'x':'graphics_x', 'y':'graphics_y'})
            net.T_node=net.T_node.merge(t_xy, left_on='id', right_on='id', how='left')
        return net

    def get_cynet(self, cynet):
        """make sure we return a cynet object"""
        if type(cynet) in (str, int):
            return self.cy.network.create(suid=int(cynet))
        return cynet

    def cynet_get_XY(self, cynet):
        cynet=self.get_cynet(cynet)
        data=cynet.get_first_view()
        if data is None: return None # view does not exist
        nodes=data['elements']['nodes']
        X=[]
        for node in nodes:
            # id is an internal id, Gene is our key
            X.append({'id': str(node['data']['id']), 'x':node['position']['x'], 'y':node['position']['y'], 'Gene':node['data']['Gene']})
        t_xy=pd.DataFrame(X)
        return t_xy

    def cynet_save(self, cynet, s_file="x.png"):
        def save_image(s_file, data):
            f=open(s_file, 'wb')
            f.write(data)
            f.close()
        cynet=self.get_cynet(cynet)
        s_ext=s_file.upper()
        if s_ext.endswith('.PNG'):
            save_image(s_file, cynet.get_png())
        elif s_ext.endswith('.PDF'):
            save_image(s_file, cynet.get_pdf())

if __name__=="__main__":
    cc=CytoscapeClient(host="localhost", port="1234")
    cynet=cc.cynet_from_xgmml("~/RIGI/Apr2015/GONetwork.xgmml", name="My_TEST")
    suid=cynet.get_id()
    print(suid, type(suid))
    net=cc.get_network(suid=suid)
    cc.cynet_save(cynet, 'x3.pdf')
    cc.cy.session.save('~/tmp/test.cys')
    cc.cy.session.delete()
    exit()
    print(cc.cy.session.name())
    cc.cy.session.save('~/tmp/test.cys')
    cc.cy.session.delete()
    cc.cy.session.load('~/tmp/test.cys')
    suid=cc.cy.network.get_all()[0]
    print(type(suid))
    print(suid) # changes after reload
    cynet=cc.get_cynet(suid)
    net=cc.get_network(suid=suid)
    cc.cynet_save(cynet, 'x2.pdf')
    cc.delete_cynet(cynet)
    net.T_node.to_csv('node.csv', index=False)
    net.T_edge.to_csv('edge.csv', index=False)

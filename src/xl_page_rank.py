"""Author: Michael Schneider
""" 
import os 
import sys 
import random
sys.path.append("../src")
sys.path.append("/scratch/schneider/libs/lib/python2.7/site-packages")

import networkx as nx
import InputOutput
import numpy
from optparse import OptionParser
import pdb
## @var parser
#  Global parser object used to call program options from anywhere in the program.
parser = OptionParser()

def add_options( parser ):
    """Generic options function to specify command line inputs
    """
    parser.add_option("-e", type="string", dest="example", help="An example")
    options, args = parser.parse_args()
    return options, args 

options, args  = add_options( parser )

"""Add edges for all the "double" cross-links AND loops"""

def load_xl_data( xl_file ):
    file  = open(xl_file)
    from_site = 0
    to_site = 0
    score = 0
    gt_data = []
    xls = []
    is_decoy = 0
    decoy_scores = []
    true_scores = []
    decoy_dict = {}
    for line in file: 
        strline = str(line).strip().split(',')
        from_site = int(strline[7])-28
        to_site = int(strline[8])-28
        score = float(strline[9])
        is_decoy = strline[10]

        if from_site > 0 and to_site > 0 and abs(from_site-to_site) >= 12:


            if is_decoy == 'FALSE':
                xls.append((return_sorted_tuple((from_site, to_site)), score/30.0))
                gt_data.append((from_site, 'CA', to_site, 'CA', score))
                true_scores.append((from_site, 'CA', to_site, 'CA', score))
                decoy_dict[return_sorted_tuple((from_site, to_site))] = 'TRUE'
            else:
                decoy_scores.append((from_site, 'CA', to_site, 'CA', score))
                decoy_dict[return_sorted_tuple((from_site, to_site))] = 'DECOY'
    file.close()
    #InputOutput.InputOutput.write_contact_file(gt_data, 'gt', upper_distance = 20)
    InputOutput.InputOutput.write_contact_file(gt_data, 'gt', upper_distance = 20)
    #InputOutput.InputOutput.write_contact_file(decoy_scores, 'decoy_scores.txt', upper_distance = 20)
    #InputOutput.InputOutput.write_contact_file(true_scores, 'true_scores.txt', upper_distance = 20)

    return xls,decoy_dict

def return_sorted_tuple(tuple):
    """Return the sorted tuple, such as the lower number is always first"""
    list_tup = list(tuple)
    list_tup.sort()
    tuple = (list_tup[0], list_tup[1])
    return tuple

# Version that already worked!
def has_loop(node_1, node_2, graph, loop_len = 1):
    loop = False
    for e in graph.neighbors(node_1[0]):    
        for f in graph.neighbors(node_2[0]):
            if e == f:
            #if is_neighbourhood(graph.node[e]['xl'], graph.node[f]['xl'], delta=2):
		
                loop = True
    return loop
 
def is_neighbourhood(tuple_1, tuple_2, delta = 1, double=True):
    is_nei = False
    t_1 = return_sorted_tuple(tuple_1)
    t_2 = return_sorted_tuple(tuple_2)
    t_1_lower = t_1[0]
    t_2_lower =	t_2[0]
    t_1_upper =	t_1[1]
    t_2_upper =	t_2[1]
    #pdb.set_trace()
    t_1_low_nei = [ t_1_lower + i for i in xrange(-1*delta, 1*delta+1)]    
    t_2_low_nei = [ t_2_lower + i for i in xrange(-1*delta, 1*delta+1)]
    t_1_up_nei =  [ t_1_upper + i for i in xrange(-1*delta, 1*delta+1)]
    t_2_up_nei	= [ t_2_upper + i for i in xrange(-1*delta, 1*delta+1)]
    m_1 = 0
    for i in t_1_low_nei:
        for j in t_2_low_nei:
            if i == j:
                m_1 = 1
    m_2	= 0
    for	i in t_1_up_nei:
        for j in t_2_up_nei:
            if i == j:
                m_2 = 1
    if double:
        if m_1 == 1 and m_2 == 1:
            is_nei = True
    else:
        if m_1 == 1 or m_2 == 1:
            is_nei = True
    return is_nei

   
def do_page_rank (xl_graph,pers,decoy_dict):

    ranked_nodes = nx.pagerank(xl_graph,max_iter=10000, alpha=0.85, tol=1e-06,personalization=pers)#, weight = 'weight')
    for_sorting = [ (score, node) for node, score in ranked_nodes.iteritems()]   
    for_sorting.sort()
    for_sorting.reverse()
    xl_ranked = [] 
    xl_ranked_true = []
    xl_ranked_decoy = []
    for score, n in for_sorting:   
        res_lower = xl_graph.node[n]['xl'][0]
        res_upper = xl_graph.node[n]['xl'][1]
        xl_ranked.append((res_lower,'CA', res_upper,'CA', score))
        if decoy_dict[(res_lower,res_upper)] == 'DECOY':
            xl_ranked_decoy.append((res_lower,'CA', res_upper,'CA', score))
        else:
            xl_ranked_true.append((res_lower,'CA', res_upper,'CA', score))

    InputOutput.InputOutput.write_contact_file(xl_ranked, 'foo', upper_distance = 20)
    #InputOutput.InputOutput.write_contact_file(xl_ranked_true, 'true_pagerank.txt', upper_distance = 20)
    #InputOutput.InputOutput.write_contact_file(xl_ranked_decoy, 'decoy_pagerank.txt', upper_distance = 20)
def add_loops( xl_graph ):
    import itertools
    all_connections = []
    for n in xl_graph.nodes(data=True):
        connecting_nodes = []
	for o in xl_graph.nodes(data=True):
             if o[0] > n[0]:
                 if is_neighbourhood(n[1]['xl'], o[1]['xl'], delta=2, double=False):
                     if len(xl_graph.neighbors(o[0])) > 0:
                         connecting_nodes.append(o)
        if len(connecting_nodes) >= 2:
  	    for o in connecting_nodes:
                all_connections.append((n,o))
    for n,o in all_connections:
        xl_graph.add_edge(n[0],o[0])	
					
def toy_graph():
    y = nx.Graph()
    y.add_node(1)
    y.add_node(2)
    y.add_node(3)
    y.add_node(4)
    y.add_edge(1,2)
    y.add_edge(2,3)
    y.add_edge(3,4)
    return y

def build_xl_graph( xl_data ):
   
    g = nx.Graph()
    index = 1
    pers = {}
    for i, score in xl_data:
        #print i
        g.add_node(index, xl=i, weight = score)
        pers[index] = score
        index += 1    

    for n in g.nodes(data=True):
        for o in g.nodes(data=True):
            if o[0] > n[0]:
                if is_neighbourhood(n[1]['xl'], o[1]['xl'], delta=1):
                    g.add_edge(n[0],o[0], weight=1.0)
                    #pass
            	#elif is_neighbourhood(n[1]['xl'], o[1]['xl'], delta=2):
                #    g.add_edge(n[0],o[0] )

    #add_loops( g )
    
    for i in xrange(0,3):
        to_add = []
        for n in g.nodes(data=True):
            for o in g.nodes(data=True):
                if o[0] > n[0]:
                    if has_loop( o, n,g ):
                        to_add.append((o[0],n[0]))
                    #g.add_edge(n[0],o[0])
        for i,j in to_add:
            g.add_edge(i,j, weight=1.0)

    

    #add_loops_node_graph(g)
    return g, pers

def get_node_map( xl_graph ):
    node_map = {}
    for i in xl_graph.nodes(data=True):
        node_map[i[1]['xl']] = i[0]
    return node_map
def add_cycles_to_graph( cycles, xl_graph, cycle_len=3 ):
    node_map = get_node_map(xl_graph)

    to_sort = [ (len(c), c) for c in cycles ]
    to_sort.sort()
    #print to_sort
    cycles = [t[1] for t in to_sort]
    for c in cycles:
        print c
        if len(c) <= cycle_len:
            link_tuples = []
            """This does full connection!"""
            """
            for i in c:
                for j in c:
                    if j > i:
                        link_tuples.append((i,j))
            """
            for i in xrange(0,len(c)):
                slice = c[i:i+2]
    #print a[i:i+2]
                if len(slice) == 2:
                    link_tuples.append(return_sorted_tuple((slice[0],slice[1])))
                    #c_tuple = return_sorted_tuple((slice[0],slice[1]))
                    #xl_graph.add_edge(node_map[c_tuple[0]])
                else:
                    slice.append(c[0])
                    link_tuples.append(return_sorted_tuple((slice[0],slice[1])))
               # print slice
        #for node in xl_graph.nodes(data=True):

            for i in xrange(0,len(link_tuples)):
                slice = link_tuples[i:i+2]
    #print a[i:i+2]
                if len(slice) == 2:
                    #link_tuples.append(return_sorted_tuple((slice[0],slice[1])))
                    #c_tuple = return_sorted_tuple((slice[0],slice[1]))
                    if node_map.has_key ( slice[0] ) and  node_map.has_key ( slice[1] ):
                        xl_graph.add_edge(node_map[slice[0]], node_map[slice[1]], weight=gauss(float(len(c)),b=4.0,c=0.5))
                else:
                    slice.append(link_tuples[0])
                    #link_tuples.append(return_sorted_tuple((slice[0],slice[1])))
                    if node_map.has_key ( slice[0] ) and  node_map.has_key ( slice[1] ):
                        xl_graph.add_edge(node_map[slice[0]], node_map[slice[1]], weight=gauss(float(len(c)),b=4.0,c=0.5))
                print slice


            """
            for l1 in xrange(0,len(link_tuples)):
                for l2 in xrange(l1+1, len(link_tuples)):

                    if node_map.has_key(link_tuples[l1]) and node_map.has_key( link_tuples[l2]):


                        if xl_graph.has_edge(node_map[link_tuples[l1]],node_map[link_tuples[l2]]):
                            pass
                        else:
                            xl_graph.add_edge(node_map[link_tuples[l1]],node_map[link_tuples[l2]],weight=1.0)

                #print xl_graph.has_edge(node_map[link_tuples[l1]],node_map[link_tuples[l2]])
		    """
		#print xl_graph.has_edge
            #print link_tuples
def gauss(x,a=1.0,b=1.0,c=1.0):
    return a * numpy.exp( -1.0 * ( (x-b)**2/2*c**2) )

def add_loops_node_graph(xl_graph):
    ng = nx.Graph()
    counter = 1
    for i in xl_graph.nodes(data=True):
        r_lower = i[1]['xl'][0]
        r_upper = i[1]['xl'][1]
        if ng.has_node(r_lower):
            pass
        else:
            ng.add_node(r_lower)
        if ng.has_node(r_upper):
            pass
        else:
            ng.add_node(r_upper)
        ng.add_edge(r_lower,r_upper)
    #for i in ng.nodes():
    #    print i
        #print nx.cycle_basis(ng,root=i)
    #for i in ng.edges():
    #    print i
    cycles = nx.cycle_basis(ng)
    add_cycles_to_graph( cycles, xl_graph, cycle_len=4 )
    #print nx.cycle_basis(ng,40)

        #ng.add_node(counter)


def main():        
   """Generic main function. Executes main functionality of program
   """
   xl_data, decoy_dict = load_xl_data( options.example )
   xl_graph, pers = build_xl_graph(xl_data)
   for n in xl_graph.nodes(data=True):
       print n
   for n in xl_graph.edges(data=True):
       print n
   print  xl_graph.number_of_edges()
   do_page_rank(xl_graph,pers,decoy_dict)

if __name__ == '__main__':
    sys.exit(main())
else:
    print("Loaded as a module!")  

"""Author: Michael Schneider
""" 
import os 
import sys 
import random
sys.path.append("../src")
sys.path.append("/scratch/schneider/libs/lib/python2.7/site-packages")
sys.path.append("/scratch/schneider/projects/contact_prediction_git/src/contact_git_code/contact_prediction/features/")

import networkx as nx
import InputOutput
import numpy
from optparse import OptionParser
import pdb
import ResidueFeatureSecStruct
import ResidueFeatureRelSasa
import cPickle
## @var parser
#  Global parser object used to call program options from anywhere in the program.
parser = OptionParser()

def add_options( parser ):
    """Generic options function to specify command line inputs
    """
    parser.add_option("-c", type="string", dest="example", help="An example")
    parser.add_option("-l", type="int", dest="length", help="An example")
    parser.add_option("-p", type="string", dest="pdb_id", help="An example")
    parser.add_option("-f", type="string", dest="pdb_file", help="An example")
    parser.add_option("-s", type="string", dest="psipred_file", help="An example")
    parser.add_option("-t", type="float", dest="top", help="An example")
    parser.add_option("-a", type="float", dest="alpha", help="An example")
    parser.add_option("-o", type="string", dest="out_folder", help="An example", default="../../results/29-08-14/")
    options, args = parser.parse_args()
    return options, args 

options, args  = add_options( parser )

"""Add edges for all the "double" cross-links AND loops"""

def parse_psipred(psipred_file):
    ss = ''
    conf = ''
    for line in open(psipred_file):
        if line.startswith('Conf:'):
            conf+=(line[6:].strip())
        elif line.startswith('Pred:'):
             ss+=(line[6:].strip())
        #ss = ''.join(ss)
        #conf = ''.join(conf)
    # turn do dict
    ss_dict = {}
    counter = 1
    for i in ss:
        ss_dict[counter] = i
        counter+=1

    return ss_dict

def load_xl_data( xl_file ):
    file  = open(xl_file)
    from_site = 0
    to_site = 0
    score = 0
    gt_data = []
    xls = []
    is_decoy = 0   
    for line in file: 
        strline = str(line).strip().split(',')
        from_site = int(strline[7])-28
        to_site = int(strline[8])-28
        score = float(strline[9])
        is_decoy = strline[10]
        if from_site > 0 and to_site > 0 and abs(from_site-to_site) >= 12 and is_decoy == 'FALSE':
            xls.append(((from_site, to_site), score/30.0))
            gt_data.append((from_site, 'CA', to_site, 'CA', score))
    file.close()
    InputOutput.InputOutput.write_contact_file(gt_data, 'gt', upper_distance = 20)
   # gt_data.append(from_site, 'CA', to_site, 'CA', 1.0)
    return xls

def return_sorted_tuple(tuple):
    """Return the sorted tuple, such as the lower number is always first"""
    list_tup = list(tuple)
    list_tup.sort()
    tuple = (list_tup[0], list_tup[1])
    return tuple
#def toy_trans()

def helix_shift(tuple1, tuple2):

    anchor1 = tuple1[0]
    anchor2 = tuple1[1]

    if tuple2[0] == anchor1 - 4:
        if tuple2[1] == anchor2 - 4:
            return True
    if tuple2[0] == anchor1 + 4:
        if tuple2[1] == anchor2 + 4:
            return True
    if tuple2[0] == anchor1 - 4 or tuple2[1] == anchor2 - 4:
        return True
    if tuple2[0] == anchor1 + 4 or tuple2[1] == anchor2 + 4:
        return True
    return False

def share_neighbors(node_1, node_2, graph, loop_len = 1):
    loop = False
    for e in graph.neighbors(node_1[0]):
        for f in graph.neighbors(node_2[0]):
            if e == f:
            #if is_neighbourhood(graph.node[e]['xl'], graph.node[f]['xl'], delta=2):

                loop = True
        #print loop
    return loop

def is_neighbourhood(tuple_1, tuple_2, delta = 1, double = True):
    is_nei = False
    t_1 = return_sorted_tuple(tuple_1)
    t_2 = return_sorted_tuple(tuple_2)
    t_1_lower = t_1[0]
    t_2_lower =	t_2[0]
    t_1_upper =	t_1[1]
    t_2_upper =	t_2[1]
    
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

def gauss(x,a=1.0,b=1.0,c=1.0):
    return a * numpy.exp( -1.0 * ( (x-b)**2/2*c**2) )

def do_page_rank (xl_graph,pers, orig_scores,input_alpha):

    ranked_nodes = nx.pagerank(xl_graph,max_iter=10000, alpha=input_alpha, tol=1e-04,personalization=pers)#,weight=None)#, weight = 'weight')
    #ranked_nodes = nx.degree_centrality(xl_graph)
    #print ranked_nodes
    for_sorting = [ (score, node) for node, score in ranked_nodes.iteritems() if node <= options.length*999]
    for_comp = []
    for_sorting.sort()
    for_sorting.reverse()
    xl_ranked = []
    for score, n in for_sorting:
        #print node
        #print xl_graph.node[n]       
        res_lower = xl_graph.node[n]['xl'][0]
        res_upper = xl_graph.node[n]['xl'][1]
        #print res_upper
        xl_ranked.append((res_lower,'CA', res_upper,'CA', score))
        for_comp.append(((res_lower,res_upper),score))
    #print xl_ranked
    #return for_sorting
    #new_stuff = linear_combination(orig_scores, for_comp, 0.5)


    InputOutput.InputOutput.write_contact_file(xl_ranked, "%s/%sRRPAR_%s_%s"%(options.out_folder,options.pdb_id,input_alpha,options.top), upper_distance = 8)

def add_loops( xl_graph ):
    #pdb.set_trace()
    import itertools
    all_connections = []
    for n in xl_graph.nodes(data=True):
        connecting_nodes = []
        for o in xl_graph.nodes(data=True):
             if o[0] > n[0]:
                 print o, n
                 if is_neighbourhood(n[1]['xl'], o[1]['xl'], delta=2, double=False):
                     if len(xl_graph.neighbors(o[0])) > 0:
                         connecting_nodes.append(o)
                          
        if len(connecting_nodes) >= 1:
            for o in connecting_nodes:
                all_connections.append((n,o))
    for n,o in all_connections:
        xl_graph.add_edge(n[0],o[0])

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
                        xl_graph.add_edge(node_map[slice[0]], node_map[slice[1]],weight=1.0)#, weight=gauss(float(len(c)),b=6.0,c=0.5))
                else:
                    slice.append(link_tuples[0])
                    #link_tuples.append(return_sorted_tuple((slice[0],slice[1])))
                    if node_map.has_key ( slice[0] ) and  node_map.has_key ( slice[1] ):
                        xl_graph.add_edge(node_map[slice[0]], node_map[slice[1]],weight=1.0)#, weight=gauss(float(len(c)),b=6.0,c=0.5))
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

def vec_to_dict(vector,pos1,pos2):
    return_dict = {}
    for v in vector:
        return_dict[v[pos1]] = v[pos2]
    return return_dict

def linear_combination(original_vector, new_vector, alpha):
    orig_dict = vec_to_dict(original_vector,0,1)
    new_dict = vec_to_dict(new_vector,0,1)
    print orig_dict
    print new_dict
    new_scores = [ (alpha*values + (1.0-alpha)*new_dict[keys], keys) for keys, values in orig_dict.iteritems() ]
    new_scores.sort()
    new_scores.reverse()
    output_scores = [(keys[0], 'CA',keys[1],'CA', score) for score, keys in new_scores]
    return output_scores

def add_sec_struct_pseudo_nodes( ss_struct_dict, xl_graph,pers ):

    ss_string = [ ss_struct_dict.ss_dict[i] for i in  xrange(1,int(options.length+1))]

    helix = []
    sheet = []
    sec_structs = []
    #print ss_string
    #print ss_string.index('H')
    for i in xrange(0,len(ss_string)):
        if ss_string[i] == 'H':
            helix.append(i)
        if len(helix) > 1 and ss_string[i] != 'H':
            sec_structs.append(helix)
            helix = []
        #if ss_string[i] == 'E':
        #    sheet.append(i)
        #if len(sheet) > 1 and ss_string[i] != 'E':
        #    sec_structs.append(sheet)
        #    sheet = []
   #print sec_structs
    return sec_structs

    nodes_so_far = xl_graph.number_of_nodes()
    offset = 1
    for i in xrange(0,len(sec_structs)):
        for j in xrange(i+1,len(sec_structs)):
           # print i,j

            xl_graph.add_node(nodes_so_far+offset, res_lower=sec_structs[i], res_upper=sec_structs[j])
            pers[nodes_so_far+offset] = 1.0
            #print i,j
            offset+=1
    for i in xrange(nodes_so_far+1, nodes_so_far+offset):
        print xl_graph.node[i]
        for j in xrange(1,nodes_so_far):
            if xl_graph.node[j]['xl'][0] in xl_graph.node[i]['res_lower'] and xl_graph.node[j]['xl'][1] in xl_graph.node[i]['res_upper']:
                xl_graph.add_edge(i,j,weight = 1.0)


def add_loops_node_graph(xl_graph):
    ng = nx.Graph()
    counter = 1
    #for i in xrange(1,options.length+1):
    #    ng.add_node(i)
    #for i in xrange(1,options.length+1):
    #    if ng.has_node(i+1):
    #        ng.add_edge(i,i+1)
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
            


def toy_graph():
    y = nx.Graph()
    y.add_node(1)
    y.add_node(2)
    y.add_node(3)
    #y.add_node(4) 
    #y.add_node(5)
    y.add_edge(1,2)
    y.add_edge(2,3)
    #y.add_edge(1,5)
    return y

def get_relative_sec_struct_pos(ss_dict, i):
    anchor = ss_dict[i]
    pos = 1
    #print pos
    for j in xrange(1,20):
        if ss_dict.has_key(i-j):
            if anchor == ss_dict[i-j]:
                pos +=1
            else:
                return pos
        else:
            return pos
    return pos

def build_xl_graph( xl_data, length, shift_dict,sec_struct,sol ):
   
    g = nx.DiGraph()
    index = 1
    pers = {}
    for i, score in xl_data[:length]:
        g.add_node(index, xl=i, weight = score)
        pers[index] = score
        index += 1    
    #add_loops_node_graph( g )
    for n in g.nodes(data=True):
        sec_lower = sec_struct[n[1]['xl'][0]]
        sec_upper = sec_struct[n[1]['xl'][1]]
        """
        if sol[n[1]['xl'][0]] < 0.3:
            sol_lower = "B"
        else:
            sol_lower = "A"
        if sol[n[1]['xl'][1]] < 0.3:
            sol_upper = "B"
        else:
            sol_upper = "A"
        """
        #if (abs( n[1]['xl'][0]- n[1]['xl'][1]) >=12 and abs( n[1]['xl'][0]- n[1]['xl'][1]) <=24 ):
        #    seq_sep = (12,24)
        #else:
        #    seq_sep = (24,9999)
        #pos1 = get_relative_sec_struct_pos(sec_struct, n[1]['xl'][0])
        #pos2 = get_relative_sec_struct_pos(sec_struct, n[1]['xl'][1])

        ##if pos1 > 10:
           # pos1 = 10
        #if pos2 > 10:
        #    pos2 = 10
        for o in g.nodes(data=True):
            if o[0] != n[0]:


                    shift_tuple = (n[1]['xl'][0] - o[1]['xl'][0], n[1]['xl'][1] - o[1]['xl'][1])
                    sec_struct_shift_dict = shift_dict[(sec_lower,sec_upper)]

                    if sec_struct_shift_dict.has_key(shift_tuple):
                        g.add_edge(n[0],o[0], weight=sec_struct_shift_dict[shift_tuple])


                    #if helix_shift(n[1]['xl'], o[1]['xl']):
                    #    g.add_edge(n[0],o[0])

                 #if is_neighbourhood(n[1]['xl'], o[1]['xl'], delta=1):
                 #   pass
                    #g.add_edge(n[0],o[0],weight=0.2)#, weight =  numpy.min([n[1]['weight'], o[1]['weight']])  )
                     #g.add_edge(n[0],o[0], weight =   gauss(2.0))
                    # g.add_edge(n[0],o[0], weight = gauss(1.0))
                 #elif is_neighbourhood(n[1]['xl'], o[1]['xl'], delta=1):
                 #    g.add_edge(n[0],o[0], weight =   gauss(1.0))
                 #if is_neighbourhood(n[1]['xl'], o[1]['xl'], delta=2):
                 #    g.add_edge(n[0],o[0], weight =   gauss(2.0))

                    #g.add_edge(n[0],o[0], weight = gauss(1.0))
                 #elif is_neighbourhood(n[1]['xl'], o[1]['xl'], delta=1):
                 #    g.add_edge(n[0],o[0], weight =  gauss(1.0))
                 #elif is_neighbourhood(n[1]['xl'], o[1]['xl'], delta=3):
                 #   g.add_edge(n[0],o[0], weight = gauss(3.0))
                    #g.add_edge(n[0],o[0], weight =  numpy.min([n[1]['weight'], o[1]['weight']])  )
                    #print  ( n[1]['weight'] + o[1]['weight']) / 2.0
                #    print 'True'
                #else:
                #    print 'False' 
    """
    for i in xrange(0,3):
        to_add = []
        for n in g.nodes(data=True):
            for o in g.nodes(data=True):
                if o[0] > n[0]:
                    if share_neighbors( o, n,g ):
                        to_add.append((o[0],n[0]))
                    #g.add_edge(n[0],o[0])
        for i,j in to_add:
            g.add_edge(i,j, weight=1.0)
    """
    #add_loops( g )
    #for i in xrange(0,3)
    #add_loops_node_graph( g )

    return g, pers
def clean_sec_structs(sec_struct):

    for i in xrange(2,len(sec_struct)-1):
        if sec_struct[i-1] == sec_struct[i+1]:
            if sec_struct[i-1] == "H" or sec_struct[i-1] == "E":
                if sec_struct[i] == 'C':
                    sec_struct[i] = sec_struct[i-1]

def parse_scores(scores_file):
    ss_dict = {}
    residues = []
    for line in open(scores_file):
        if line.startswith('#') or not line.strip():
            continue
        else:
            line = line.split()
            rank = int(line[0])
            residue = line[1]
            coil, helix, strand = map(float, line[3:6])
            ss_dict[ rank] = ( coil,helix,strand )
    return ss_dict



def main():

   """Generic main function. Executes main functionality of program
   """
   #print is_neighbourhood((3,10), (4,11), delta = 1)
   #sys.exit(main())
   #xl_data = load_xl_data( options.example )
   #sec_struct = ResidueFeatureSecStruct.ResidueFeatureSecStruct(options.pdb_file)
   #buried_features = ResidueFeatureRelSasa.ResidueFeatureRelSasa(options.pdb_file)
   bur_dict = {}
   #for i in xrange(1,options.length+1):
   #    buried_features.calculate_feature(i,options.pdb_file)
   #    bur_dict[i] = buried_features.get_feature()
   #print bur_dict
   #return 0

   sec_struct = parse_psipred(options.psipred_file)
   clean_sec_structs(sec_struct)
   #parse_scores(scores_file)
   #print sec_struct
   #

   #print sec_struct
   #print
   #return 0
   #print sec_struct.ss_dict
   #return 0
   #print sec_struct.ss_dict
   #add_sec_struct_pseudo_nodes(sec_struct)
   #print
   #print i
   #return 0
   shift_dict = cPickle.load(open( "../probabilities/shifts.p", "rb" ))
   xl_data = InputOutput.InputOutput.load_restraints_pr(options.example,seq_sep_min=12)
   #print xl_data
   xl_graph,pers = build_xl_graph(xl_data,int(options.length*options.top), shift_dict, sec_struct,bur_dict)
   #add_sec_struct_pseudo_nodes(sec_struct, xl_graph,pers )

   for i in xl_graph.nodes(data=True):
       print i
   for e in xl_graph.edges(data=True):
       print e
   #return 0
   print xl_graph.number_of_edges()
   #xl_graph,pers = build_xl_graph(xl_data)
   do_page_rank(xl_graph,pers,xl_data[:int(options.length*options.top)],options.alpha)


if __name__ == '__main__':
    sys.exit(main())
else:
    print("Loaded as a module!")  

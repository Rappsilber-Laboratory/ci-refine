"""Author: Michael Schneider
   Very important notes
""" 
import os 
import sys 
import random
sys.path.append("../src")
sys.path.append("/scratch/schneider/libs/lib/python2.7/site-packages")
sys.path.append("/scratch/schneider/projects/contact_prediction_git/src/contact_git_code/contact_prediction/features/")
sys.path.append("/scratch/schneider/projects/contact_prediction_git/src/contact_git_code/contact_prediction/structure/")
import networkx as nx
import InputOutput
import numpy
from optparse import OptionParser
import pdb
import ResidueFeatureSecStruct
import ResidueFeatureRelSasa
import StructureContainer
import cPickle
import matplotlib.pyplot as plt
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
            conf += (line[6:].strip())
        elif line.startswith('Pred:'):
             ss += (line[6:].strip())
        #ss = ''.join(ss)
        #conf = ''.join(conf)
    # turn do dict
    ss_dict = {}
    counter = 1
    for i in ss:
        ss_dict[counter] = i
        counter += 1

    return ss_dict

def shift_matrix():
    matrix = []
    for i in xrange(-8,9):
        row = []
        for j in xrange(-8,9):
            if i == 0 and j == 0:
                pass
            else:
                row.append((i,j))
                matrix.append((i,j))
    return matrix


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
    #all_scores = {}
    #for i in xl_graph.nodes():
    #    all_scores[i] = 0.0
    #alphas = [0.85]
    #print pers
    ranked_nodes = nx.pagerank(xl_graph,max_iter=1000, alpha=input_alpha, tol=1e-04,personalization=pers)
    #for a in alphas:

       # ranked_nodes = nx.pagerank(xl_graph,max_iter=1000, alpha=a, tol=1e-04)#,personalization=pers)#,weight=None)#, weight = 'weight')
    #    print "Computing PageRank for Alpha:", a
    #    for node, score in ranked_nodes.iteritems():#

    #        all_scores[node] = all_scores[node] + score
    #for node, score in pers.iteritems():
    #    all_scores[node] = all_scores[node]*input_alpha + ( (1.0 -input_alpha)* score)
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

def write_edge_scores( graph, true_contacts ):
    true_dict = vec_to_dict(true_contacts,0,1)
    class_neg  = []
    class_pos  = []
    class_neg_pos = []
    all_class = []
    draw_graph(graph,true_dict)
    #print graph.nodes(data=True)
    #print graph.nodes()[96]
    for e in graph.edges(data=True):
        #print e
        #print graph.nodes(data=True)[e[0]-1][1]['xl']
        #print graph.nodes(data=True)[e[1]-1][1]['xl']
        all_class.append( (graph[e[0]][e[1]]['weight'], (e[0],e[1])))
        if true_dict.has_key(graph.nodes(data=True)[e[0]-1][1]['xl']) == False and true_dict.has_key(graph.nodes(data=True)[e[1]-1][1]['xl']) == False:
            class_neg.append( graph[e[0]][e[1]]['weight'] )
            #graph[e[0]][e[1]]['weight'] = 0.1
        elif true_dict.has_key(graph.nodes(data=True)[e[0]-1][1]['xl']) == True and true_dict.has_key(graph.nodes(data=True)[e[1]-1][1]['xl']) == True:
            class_pos.append( graph[e[0]][e[1]]['weight'] )

            #graph[e[0]][e[1]]['weight'] = 0.1
        else:
            #print true_dict.has_key(graph.nodes(data=True)[e[0]-1][1]['xl']), true_dict.has_key(graph.nodes(data=True)[e[1]-1][1]['xl'])
            class_neg_pos.append( graph[e[0]][e[1]]['weight'] )
    all_class.sort()
    """
    to_rm = []
    for n in graph.nodes():
        weight_list = []
        for nei in nx.neighbors(graph,n):

            weight_list.append((graph.edge[n][nei]['weight'],(n,nei)))
        weight_list.sort()
        weight_list.reverse()
        for i in weight_list[:10]:
            to_rm.append(i)

    for i,j in to_rm:
        if graph.has_edge(j[0],j[1]):
            graph.remove_edge(j[0],j[1])
    """
    #for e in graph.edges(data=True):
    #    if true_dict.has_key(graph.nodes(data=True)[e[0]-1][1]['xl']) == True and true_dict.has_key(graph.nodes(data=True)[e[1]-1][1]['xl']) == True:
    #        #class_pos.append( numpy.mean(class_neg)*2.0 )
    #        graph[e[0]][e[1]]['weight'] = numpy.mean(class_neg)*2.0
    #for i,j in all_class[:int(len(all_class)*0.25)]:
    #    graph.remove_edge(j[0],j[1])

    #all_class.reverse()
    #med = all_class[int(len(all_class)*0.5)][0]
    #to_rm = []
    ##for e in graph.edges(data=True):
    #    if graph[e[0]][e[1]]['weight'] <= med:
    #        to_rm.append()

    print "CLASS",  numpy.mean(class_pos), numpy.mean(class_neg), numpy.mean(class_neg_pos)
    #draw_graph(graph,true_dict)


def build_xl_graph( xl_data, length, shift_dict,sec_struct,sol, clust_aligns = None ):
    tmp_struct = StructureContainer.StructureContainer()
    tmp_struct.load_structure('xxxx', options.pdb_id[-1], options.pdb_file, seqsep =1)
    true_map = tmp_struct.get_contact_map().print_res_format()
    #xl_data = true_map
    g = nx.Graph()
    index = 1
    pers = {}
    for i, score in xl_data[:length]:
        g.add_node(index, xl=i, weight = score)
        pers[index] = score
        index += 1    
    #add_loops_node_graph( g )
    nodes_to_add = []
    """
    for n in g.nodes(data=True):
        test_vec = get_prediction_vector(xl_data[:length], n[1]['xl'][0], n[1]['xl'][1])

        #lowest_clust =  get_lowest_scoring_clust(test_vec, clust_aligns)
        sec_struct_shift_dict = get_averaged_dict(test_vec, clust_aligns)

        max_val = 0
        max_key = 0
        for keys, values in sec_struct_shift_dict.iteritems():
            if values > max_val:
                max_val = values
                max_key = keys
        nodes_to_add.append(( (n[1]['xl'][0] + max_key[0],n[1]['xl'][1] + max_key[1]), n[1]['weight']) )

    for i, score in nodes_to_add:
        g.add_node(index, xl=i, weight = score)
        pers[index] = score
        index += 1
    """
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
        #test_vec = get_prediction_vector(xl_data[:int(length*1.0)], n[1]['xl'][0], n[1]['xl'][1])
        #test_vec = get_prediction_vector(true_map, n[1]['xl'][0], n[1]['xl'][1])
        #lowest_clust =  get_lowest_scoring_clust(test_vec, clust_aligns)
        #print clust_aligns[(sec_lower, sec_upper)]
       # sec_struct_shift_dict = get_averaged_dict(test_vec, clust_aligns[(sec_lower, sec_upper)])
        sec_struct_shift_dict = shift_dict[(sec_lower,sec_upper)]
        #pseudo_shift = {}
        #pseudo_shift[0] = sec_struct_shift_dict
        #cPickle.dump(pseudo_shift, open( "../pseudo_shift.p", "wb" ),protocol=2 )
        #break
        if sec_struct_shift_dict != False:
            #sec_struct_shift_dict = shift_dict[(lowest_clust)]
            a = 'a'
            for o in g.nodes(data=True):
                if o[0] != n[0]:# and sec_lower == sec_struct[o[1]['xl'][0]] and sec_upper == sec_struct[o[1]['xl'][1]]:

                    shift_tuple = (n[1]['xl'][0] - o[1]['xl'][0], n[1]['xl'][1] - o[1]['xl'][1])

                    dist = numpy.sqrt(shift_tuple[0]**2+shift_tuple[1]**2)
                    if sec_struct_shift_dict.has_key(shift_tuple) and numpy.isnan(sec_struct_shift_dict[shift_tuple]) == False and sec_struct_shift_dict[shift_tuple] != 0.0:
                        if g.has_edge(n[0],o[0]):
                            #if g.edge[n[0]][o[0]]['weight'] < sec_struct_shift_dict[shift_tuple]:
                            old_weight = g.edge[n[0]][o[0]]['weight']
                            if old_weight > sec_struct_shift_dict[shift_tuple]:
                                g.add_edge(n[0],o[0], weight=sec_struct_shift_dict[shift_tuple] )
                        else:
                            g.add_edge(n[0],o[0], weight=sec_struct_shift_dict[shift_tuple])

    return g, pers


def which_clust(i, all_clust):
    for clust in all_clust:
        if i in clust:
            return clust

def clust_graph(graph):
    from sklearn import metrics
    from sklearn import cluster
    A = nx.adjacency_matrix(graph)
    max_score = 0
    max_labels = []
    max_num_clust = 0
    all_clust = []
    for i in xrange(2,8):
        spec = cluster.SpectralClustering(n_clusters=i)
        labels = spec.fit_predict(A)
        score = metrics.silhouette_score(A, labels, metric='euclidean')
        if score > max_score:
            max_score = score
            max_labels = labels
            max_num_clust = i
    #print max_labels
    for c in xrange(0,max_num_clust):
        clust = [i+1 for i,j in enumerate(max_labels) if j == c]

        all_clust.append(clust)

    #print all_clust

    print "CLUST", max_score, max_num_clust
    if max_num_clust == 2:
        return None
    else:
        return all_clust
    #print
    #print len(labels)
    #print labels
def draw_graph( graph, true_map ):
    true_nodes = []
    false_nodes = []
    #false_true_nodes = []
    for n in graph.nodes(data=True):
        if true_map.has_key(n[1]['xl']):
            true_nodes.append(n[0])
        else:
            false_nodes.append(n[0])
    pos=nx.spring_layout(graph)
    #nx.draw_networkx_nodes(graph,pos, nodelist=true_nodes, node_color='b', node_size=50, alpha=0.8)
    #nx.draw_networkx_nodes(graph,pos, nodelist=false_nodes, node_color='r', node_size=50, alpha=0.9)


    #clust_1 = [i+1 for i,j in enumerate(labels) if j == 0]
    #clust_2 = [i+1 for i,j in enumerate(labels) if j == 1]
    #clust_3 = [i+1 for i,j in enumerate(labels) if j == 2]
    #nx.draw_networkx_nodes(graph,pos, nodelist=clust_1, node_color='b', node_size=50, alpha=0.8)
    #nx.draw_networkx_nodes(graph,pos, nodelist=clust_2, node_color='r', node_size=50, alpha=0.9)
    #nx.draw_networkx_nodes(graph,pos, nodelist=clust_3, node_color='g', node_size=50, alpha=0.9)
    #all_clust = [clust_1,clust_2,clust_3]

    all_clust = clust_graph(graph)
    if all_clust != None:
        to_rm = []
    #for clust in all_clust:
        for e in graph.edges(data=True):
            print e
            clust = which_clust(e[0], all_clust )

            if e[1] in clust:
                pass
            else:
                to_rm.append((e[0],e[1]))
        for e in to_rm:
            graph.remove_edge(e[0],e[1])
        #same_clust = True
        #for clust in all_clust:
    #nx.draw_networkx_edges(graph,pos,width=0.2,alpha=0.5)
    #plt.show()
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



def get_prediction_vector( contact_list, i,j ):

    c_dict = vec_to_dict(contact_list,0,1)
    shift_mat = shift_matrix()

    pred_vec = []

    for i_shift, j_shift in shift_mat:
        if c_dict.has_key((i+i_shift, j+j_shift)):
            pred_vec.append(c_dict[(i+i_shift,j+j_shift)])
            #pred_vec.append(1.0)
        else:
            pred_vec.append(0.0)
    #sum_prob = numpy.sum(pred_vec)
    #for i in xrange(0,len(pred_vec)):
    #    pred_vec[i] = pred_vec[i]/sum_prob
    #print numpy.array(pred_vec)
    return numpy.array(pred_vec)

def get_clustered_aligns( shift_dict ):
    shift_mat = shift_matrix()
    vec_dict = {}
    for keys, values in shift_dict.iteritems():
        vec_dict[keys] = 1

        #print values
        pred_dict = {}
        for i in xrange(0, len(values)):
           # print values
            pred_vec = []
            #print i
            for i_shift, j_shift in shift_mat:
            #print i_shift
            #if i_shift != 0 and j_shift != 0:
                #print i_shift, j_shift, values[(i_shift,j_shift)]
              #  try:
             #       if values
                if values[i].has_key((i_shift,j_shift)):
                    pred_vec.append(values[i][(i_shift,j_shift)])
                else:
                    pred_vec.append(0.0)
               # except:
                #    pred_vec.append(0.0)
        #print keys, values
            pred_dict[i] = numpy.array(pred_vec)
        vec_dict[keys] = pred_dict
    return vec_dict

def get_lowest_scoring_clust(vec, clust_aligns):
    lowest_score = 0.0
    lowest_clust = 99
    for keys, values in clust_aligns.iteritems():
        if numpy.dot(vec,values) > lowest_score:
            lowest_clust = keys
            #print vec#, values
            #print values
            lowest_score = numpy.dot(vec,values)

        #print lowest_score

    if clust_aligns.has_key(lowest_clust):
        new_dict = {}
        shift_mat = shift_matrix()
        for shifts,val in zip(shift_mat, clust_aligns[lowest_clust]):
            new_dict[shifts] = val
        return new_dict
    else:
        return False

def normalize_per_position(clust_aligns):

    shift_mat = shift_matrix()
    shift_dict = {}

    dict_info = []
    for keys,values in clust_aligns.iteritems():
        for i in xrange(0,10):
            dict_info.append((keys,i))


    sum_per_stuff = [0]*len(clust_aligns[('H','H')][0])

    for keys in dict_info:

        for i in xrange(0,len(sum_per_stuff)):

            sum_per_stuff[i] = sum_per_stuff[i] + clust_aligns[keys[0]][keys[1]][i]

    for keys in dict_info:
        #print len(clust_aligns[('H','H')][0])
        new_val = [0]*len(clust_aligns[('H','H')][0])

        for i in xrange(0,len(sum_per_stuff)):

            new_val[i] = clust_aligns[keys[0]][keys[1]][i] / sum_per_stuff[i]

        clust_aligns[keys[0]][keys[1]] = new_val


def get_averaged_dict(vec, clust_aligns):
   # normalize_per_position(clust_aligns)
    shift_mat = shift_matrix()
    #all_vec = [0]*len(clust_aligns[('H','H')])
    all_vec = [0]*len(clust_aligns[0])
    new_dict = {}
    #pdb.set_trace()
    sum_score = 0.0
    scores = {}
    for keys, values in clust_aligns.iteritems():
        res = numpy.dot(vec,values)
        scores[keys] = res
        sum_score += res
        #print res, len(vec),len(values), vec
    if sum_score == 0:
        norm_scores = {}
        for keys, values in scores.iteritems():
            norm_scores[keys] = 0.1
    else:
        norm_scores = {}
        for keys, values in scores.iteritems():
            norm_scores[keys] = 0.2 #scores[keys] / sum_score

    print norm_scores

    for i, s in norm_scores.iteritems():

        for v_index in xrange(0,len(all_vec)):
            if s > 0.0:
                all_vec[v_index]= all_vec[v_index] + s * clust_aligns[i][v_index]

    n_sum = numpy.sum(all_vec)

    for shifts, val in zip(shift_mat, all_vec):
        new_dict[shifts] = val / n_sum
        #for shifts,values in zip(shift_mat,clust_aligns[i]):
           # all_vec
            #if new_dict.has_key((i_shift,j_shift)):
            #    new_dict[(i_shift,j_shift)] = new_dict(i_shift,j_shift) + s*

        #if numpy.dot(vec,values) > lowest_score:
            #lowest_clust = keys
            #print vec#, values
            #print values
            #lowest_score = numpy.dot(vec,values)
        #print lowest_score
    return new_dict


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
   #print shift_dict[('E','H')]
   #clust_aligns = get_clustered_aligns(shift_dict)
   #normalize_per_position(clust_aligns)
   #normalize_per_position(clust_aligns)
   xl_data = InputOutput.InputOutput.load_restraints_pr(options.example,seq_sep_min=12)
   #print xl_data
   #print xl_data[0][0]
   #test_vec = get_prediction_vector(xl_data[:int(options.length*options.top)], xl_data[0][0][0], xl_data[0][0][1])
   #print test_vec, clust_aligns[0]
   #print get_lowest_scoring_clust(test_vec, clust_aligns)

   #sys.exit()
   xl_graph,pers = build_xl_graph(xl_data,int(options.length*options.top), shift_dict, sec_struct,bur_dict)
   #add_sec_struct_pseudo_nodes(sec_struct, xl_graph,pers )

   #for i in xl_graph.nodes(data=True):
   #    print i
   #for e in xl_graph.edges(data=True):
   #    print e
   #return 0
   print xl_graph.number_of_edges()
   #xl_graph,pers = build_xl_graph(xl_data)
   do_page_rank(xl_graph,pers,xl_data[:int(options.length*options.top)],options.alpha)


if __name__ == '__main__':
    sys.exit(main())
else:
    print("Loaded as a module!")  

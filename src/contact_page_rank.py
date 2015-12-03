"""Author: Michael Schneider
""" 
import os 
import sys 
import random
import networkx as nx
import InputOutput
import numpy
import argparse
import cPickle
import matplotlib.pyplot as plt
import pdb

from features.ResidueFeatureSecStruct import ResidueFeatureSecStruct
from features.ResidueFeatureRelSasa  import ResidueFeatureRelSasa
from structure.StructureContainer import StructureContainer

options = {}

def main():
    parse_arguments()
    sec_struct = parse_psipred(options.psipred_file)
    shift_dict = cPickle.load(open("probabilities/shifts_sigma_0.05.txt", "rb"))
    xl_data = InputOutput.InputOutput.load_restraints_pr(options.contact_file, seq_sep_min=12)
    xl_graph, pers = build_ce_graph(xl_data, int(options.length*options.top), shift_dict, sec_struct)
    do_page_rank(xl_graph, pers, xl_data[:int(options.length*options.top)], options.alpha)


def parse_arguments():
    """Specify and parse command line inputs
    """
    global options
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", dest="contact_file", help="predicted contacts file in CASP format", required=True)
    parser.add_argument("-l", type=int, dest="length", help="number of residues in the protein", required=True)
    parser.add_argument("-p", dest="pdb_id", help="pdb id + chain id (5 letters)", required=True)
    parser.add_argument("-f", dest="pdb_file", help="pdb file. used for reference, not calculation", required=True)
    parser.add_argument("-s", dest="psipred_file", help="sequence and secondary structure file in psipred format", required=True)
    parser.add_argument("-t", type=float, dest="top", help="fraction of top probable contacts to use. 0 < x < 1", required=True)
    parser.add_argument("-a", type=float, dest="alpha", help="dampening parameter alpha", required=True)
    parser.add_argument("-o", dest="out_folder", help="output folder", default="../results/2015-12-04/")
    options = parser.parse_args()


def parse_psipred(psipred_file):
    ss = ''
    conf = ''
    for line in open(psipred_file):
        if line.startswith('Conf:'):
            conf += (line[6:].strip())
        elif line.startswith('Pred:'):
            ss += (line[6:].strip())

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
    t_2_lower =    t_2[0]
    t_1_upper =    t_1[1]
    t_2_upper =    t_2[1]
    
    t_1_low_nei = [ t_1_lower + i for i in xrange(-1*delta, 1*delta+1)]    
    t_2_low_nei = [ t_2_lower + i for i in xrange(-1*delta, 1*delta+1)]
    t_1_up_nei =  [ t_1_upper + i for i in xrange(-1*delta, 1*delta+1)]
    t_2_up_nei    = [ t_2_upper + i for i in xrange(-1*delta, 1*delta+1)]
    m_1 = 0
    for i in t_1_low_nei:
        for j in t_2_low_nei:
            if i == j:
                m_1 = 1
    m_2    = 0
    for    i in t_1_up_nei:
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

def gauss(x, a=1.0, b=1.0, c=1.0):
    return a * numpy.exp(-1.0 * ( (x-b)**2/2*c**2))


def do_page_rank (xl_graph, pers, orig_scores, input_alpha):
    #tmp_struct = StructureContainer()
    #tmp_struct.load_structure('xxxx', options.pdb_id[-1], options.pdb_file, seqsep =1)
    ranked_nodes = nx.pagerank(xl_graph,max_iter=1000, alpha=input_alpha, tol=1e-04,personalization=pers)#,weight=None)#, weight = 'weight')
    for_sorting = [ (score , node) for node, score in ranked_nodes.iteritems() if node <= options.length*999]
    for_comp = []
    for_sorting.sort()
    for_sorting.reverse()
    xl_ranked = []
    for score, n in for_sorting:
        res_lower = xl_graph.node[n]['xl'][0]
        res_upper = xl_graph.node[n]['xl'][1]
        xl_ranked.append((res_lower,'CB', res_upper,'CB', score))
        for_comp.append(((res_lower,res_upper),score))
    InputOutput.InputOutput.write_contact_file(xl_ranked, "%s%s_RRPAR_%s_%s"%(options.out_folder, options.pdb_id, input_alpha, options.top), upper_distance = 8)


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

            for i in xrange(0,len(c)):
                slice = c[i:i+2]

                if len(slice) == 2:
                    link_tuples.append(return_sorted_tuple((slice[0],slice[1])))

                else:
                    slice.append(c[0])
                    link_tuples.append(return_sorted_tuple((slice[0],slice[1])))

            for i in xrange(0,len(link_tuples)):
                slice = link_tuples[i:i+2]

                if len(slice) == 2:


                    if node_map.has_key ( slice[0] ) and  node_map.has_key ( slice[1] ):
                        xl_graph.add_edge(node_map[slice[0]], node_map[slice[1]],weight=1.0)#, weight=gauss(float(len(c)),b=6.0,c=0.5))
                else:
                    slice.append(link_tuples[0])
                    #link_tuples.append(return_sorted_tuple((slice[0],slice[1])))
                    if node_map.has_key ( slice[0] ) and  node_map.has_key ( slice[1] ):
                        xl_graph.add_edge(node_map[slice[0]], node_map[slice[1]],weight=1.0)#, weight=gauss(float(len(c)),b=6.0,c=0.5))
                print slice



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


def add_loops_node_graph(xl_graph):
    ng = nx.Graph()

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

    cycles = nx.cycle_basis(ng)

    add_cycles_to_graph( cycles, xl_graph, cycle_len=4 )


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
    for j in xrange(1 ,20):
        if ss_dict.has_key(i-j):
            if anchor == ss_dict[i-j]:
                pos += 1
            else:
                return pos
        else:
            return pos
    return pos


def is_same_sec_struct(tuple1,tuple2,ss_dict):

    ss_i_1 = get_sec_struct_limits(ss_dict,tuple1[0])
    ss_i_2 = get_sec_struct_limits(ss_dict,tuple1[1])

    if tuple2[0] in ss_i_1 or tuple2[1] in ss_i_1:
        if tuple2[1] in ss_i_2 or tuple2[0] in ss_i_2:
            return True

    return False


def write_edge_scores( graph, true_contacts, pers ):
    true_dict = vec_to_dict(true_contacts,0,1)
    class_neg  = []
    class_pos  = []
    class_neg_pos = []
    all_class = []
    #draw_graph(graph,true_dict, pers)
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


    print "CLASS",  numpy.mean(class_pos), numpy.mean(class_neg), numpy.mean(class_neg_pos)


def get_sec_struct_limits(ss_dict,i):
    anchor = ss_dict[i]
    pos = i
    lower_pos = i
    upper_pos = i
    #pdb.set_trace()
    for j in xrange(1,20):
        if ss_dict.has_key(i-j):
            if anchor == ss_dict[i-j]:
                lower_pos-=1
            else:
                break
        else:
            break
    for j in xrange(1,20):
        if ss_dict.has_key(i+j):
            if anchor == ss_dict[i+j]:
                upper_pos+=1
            else:
                break
        else:
            return [ i for i in xrange(lower_pos, upper_pos+1) ]
    return [ i for i in xrange(lower_pos, upper_pos+1) ]


def gauss_filter_probs(xl_data, length):
    x = numpy.array(numpy.zeros((options.length, options.length), numpy.float))
    for i, score in xl_data[:length]:
        x[i[0]-1][i[1]-1] = score
    import scipy.ndimage.filters as filters
    x = filters.gaussian_filter(x, 0.5)
    new_data = []
    for row in xrange(0,x.shape[0]):
        for col in  xrange(0,x.shape[1]):
            if col > row:
                new_data.append((x[row][col], (row+1,col+1)) )
    new_data.sort()
    new_data.reverse()
    return new_data[:length]


def build_ce_graph( xl_data, length, shift_dict, sec_struct):
    #tmp_struct = StructureContainer()
    #tmp_struct.load_structure('xxxx', options.pdb_id[-1], options.pdb_file, seqsep =1)
    g = nx.Graph()

    index = 1
    pers = {}
    for score, i in xl_data[:length]:
        g.add_node(index, xl=i, weight = score)
        pers[index] = score
        index += 1    

    for n in g.nodes(data=True):
        sec_lower = sec_struct[n[1]['xl'][0]]
        sec_upper = sec_struct[n[1]['xl'][1]]

        sec_struct_shift_dict = shift_dict[(sec_lower,sec_upper)]

        if sec_struct_shift_dict != False:
            for o in g.nodes(data=True):
                if o[0] != n[0]:# and sec_lower == sec_struct[o[1]['xl'][0]] and sec_upper == sec_struct[o[1]['xl'][1]]:
                    shift_tuple = (n[1]['xl'][0] - o[1]['xl'][0], n[1]['xl'][1] - o[1]['xl'][1])

                    if sec_struct_shift_dict.has_key(shift_tuple) and numpy.isnan(sec_struct_shift_dict[shift_tuple]) == False and sec_struct_shift_dict[shift_tuple] != 0.0:
                        if g.has_edge(n[0],o[0]):

                            old_weight = g.edge[n[0]][o[0]]['weight']
                            if old_weight > sec_struct_shift_dict[shift_tuple]:
                                g.add_edge(n[0],o[0], weight=sec_struct_shift_dict[shift_tuple] )
                        else:
                            g.add_edge(n[0],o[0], weight=sec_struct_shift_dict[shift_tuple])



    print len(g.edges())

    return g, pers

def average_weight(graph):
    average_weight = []
    for i in graph.edges(data=True):
        average_weight.append(i[2]['weight'])
    return numpy.mean(average_weight)

def remove_weight_percentile(graph):
    average_weight = []
    for i in graph.edges(data=True):
        average_weight.append(i[2]['weight'])

    perc =  numpy.percentile(average_weight,10)
    for i in graph.edges(data=True):
        if i[2]['weight'] < perc:
            graph.remove_edge(i[0],i[1])
    return graph

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
        spec = cluster.KMeans(n_clusters=i)
        labels = spec.fit_predict(A)
        score = metrics.silhouette_score(A, labels, metric='euclidean')
        if score > max_score:
            max_score = score
            max_labels = labels
            max_num_clust = i
    for c in xrange(0,max_num_clust):
        clust = [i+1 for i,j in enumerate(max_labels) if j == c]

        all_clust.append(clust)

    print "CLUST", max_score, max_num_clust
    if max_num_clust == 2:
        return None
    else:
        return all_clust


def draw_graph( graph, true_map, pers, clust=None ):
    true_nodes = []
    false_nodes = []

    true_pers = []
    false_pers = []
    for n in graph.nodes(data=True):
        if true_map.has_key(n[1]['xl']):
            true_nodes.append(n[0])
            true_pers.append( int(pers[n[0]]*10000))
        else:
            false_nodes.append(n[0])
            false_pers.append( int(pers[n[0]]*10000))
    pos=nx.spring_layout(graph)


    nx.draw_networkx_nodes(graph,pos, nodelist=true_nodes, node_color='b', node_size=true_pers, alpha=0.8)
    nx.draw_networkx_nodes(graph,pos, nodelist=false_nodes, node_color='r', node_size=false_pers, alpha=0.9)
    nx.draw_networkx_edges(graph,pos,width=0.2,alpha=0.5)
    plt.show()

def clean_sec_structs(sec_struct):

    for i in xrange(2,len(sec_struct)-1):
        if sec_struct[i-1] == sec_struct[i+1]:
            if sec_struct[i-1] == "H" or sec_struct[i-1] == "E":
                if sec_struct[i] == 'C':
                    sec_struct[i] = sec_struct[i-1]

def parse_scores(scores_file):
    ss_dict = {}
    for line in open(scores_file):
        if line.startswith('#') or not line.strip():
            continue
        else:
            line = line.split()
            rank = int(line[0])
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




if __name__ == '__main__':
    sys.exit(main())
else:
    print("Loaded as a module!")  

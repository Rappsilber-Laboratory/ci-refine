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
    parser.add_option("-o", type="int", dest="offset", help="offset")

    options, args = parser.parse_args()
    return options, args

options, args  = add_options( parser )

"""Add edges for all the "double" cross-links AND loops"""
def pagerank(G, alpha=0.85, personalization=None,
             max_iter=100, tol=1.0e-6, nstart=None, weight='weight',
             dangling=None):
    """Return the PageRank of the nodes in the graph.

    PageRank computes a ranking of the nodes in the graph G based on
    the structure of the incoming links. It was originally designed as
    an algorithm to rank web pages.

    Parameters
    -----------
    G : graph
      A NetworkX graph.  Undirected graphs will be converted to a directed
      graph with two directed edges for each undirected edge.

    alpha : float, optional
      Damping parameter for PageRank, default=0.85.

    personalization: dict, optional
      The "personalization vector" consisting of a dictionary with a
      key for every graph node and nonzero personalization value for each node.
      By default, a uniform distribution is used.

    max_iter : integer, optional
      Maximum number of iterations in power method eigenvalue solver.

    tol : float, optional
      Error tolerance used to check convergence in power method solver.

    nstart : dictionary, optional
      Starting value of PageRank iteration for each node.

    weight : key, optional
      Edge data key to use as weight.  If None weights are set to 1.

    dangling: dict, optional
      The outedges to be assigned to any "dangling" nodes, i.e., nodes without
      any outedges. The dict key is the node the outedge points to and the dict
      value is the weight of that outedge. By default, dangling nodes are given
      outedges according to the personalization vector (uniform if not
      specified). This must be selected to result in an irreducible transition
      matrix (see notes under google_matrix). It may be common to have the
      dangling dict to be the same as the personalization dict.

    Returns
    -------
    pagerank : dictionary
       Dictionary of nodes with PageRank as value

    Examples
    --------
    >>> G = nx.DiGraph(nx.path_graph(4))
    >>> pr = nx.pagerank(G, alpha=0.9)

    Notes
    -----
    The eigenvector calculation is done by the power iteration method
    and has no guarantee of convergence.  The iteration will stop
    after max_iter iterations or an error tolerance of
    number_of_nodes(G)*tol has been reached.

    The PageRank algorithm was designed for directed graphs but this
    algorithm does not check if the input graph is directed and will
    execute on undirected graphs by converting each edge in the
    directed graph to two edges.

    See Also
    --------
    pagerank_numpy, pagerank_scipy, google_matrix

    References
    ----------
    .. [1] A. Langville and C. Meyer,
       "A survey of eigenvector methods of web information retrieval."
       http://citeseer.ist.psu.edu/713792.html
    .. [2] Page, Lawrence; Brin, Sergey; Motwani, Rajeev and Winograd, Terry,
       The PageRank citation ranking: Bringing order to the Web. 1999
       http://dbpubs.stanford.edu:8090/pub/showDoc.Fulltext?lang=en&doc=1999-66&format=pdf
    """
    if len(G) == 0:
        return {}

    if not G.is_directed():
        D = G.to_directed()
    else:
        D = G

    # Create a copy in (right) stochastic form
    W = nx.stochastic_graph(D, weight=weight)
    N = W.number_of_nodes()

    # Choose fixed starting vector if not given
    if nstart is None:
        x = dict.fromkeys(W, 1.0 / N)
    else:
        # Normalized nstart vector
        s = float(sum(nstart.values()))
        x = dict((k, v / s) for k, v in nstart.items())

    if personalization is None:
        # Assign uniform personalization vector if not given
        p = dict.fromkeys(W, 1.0 / N)
    else:
        missing = set(G) - set(personalization)
        if missing:
            raise NetworkXError('Personalization dictionary '
                                'must have a value for every node. '
                                'Missing nodes %s' % missing)
        s = float(sum(personalization.values()))
        p = dict((k, v / s) for k, v in personalization.items())

    if dangling is None:
        # Use personalization vector if dangling vector not specified
        dangling_weights = p
    else:
        missing = set(G) - set(dangling)
        if missing:
            raise NetworkXError('Dangling node dictionary '
                                'must have a value for every node. '
                                'Missing nodes %s' % missing)
        s = float(sum(dangling.values()))
        dangling_weights = dict((k, v/s) for k, v in dangling.items())
    dangling_nodes = [n for n in W if W.out_degree(n, weight=weight) == 0.0]
    # power iteration: make up to max_iter iterations
    for _ in range(max_iter):
        xlast = x
        x = dict.fromkeys(xlast.keys(), 0)
        danglesum = alpha * sum(xlast[n] for n in dangling_nodes)
        for n in x:
            # this matrix multiply looks odd because it is
            # doing a left multiply x^T=xlast^T*W
            for nbr in W[n]:
                x[nbr] += alpha * xlast[n] * W[n][nbr][weight]
            x[n] += danglesum * dangling_weights[n] + (1.0 - alpha) * p[n]
        # check convergence, l1 norm
        err = sum([abs(x[n] - xlast[n]) for n in x])
        if err < N*tol:
            return x, err
    raise NetworkXError('pagerank: power iteration failed to converge '
                        'in %d iterations.' % max_iter)

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
    #print loop
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

   
def do_page_rank (xl_graph,pers):

    ranked_nodes, error = pagerank(xl_graph,max_iter=10000, alpha=0.85, tol=1e-06,personalization=pers)#, weight = 'weight')
    print (0.85 / (1-0.85))* error
    for_sorting = [ (score, node) for node, score in ranked_nodes.iteritems()]
    for_sorting.sort()
    for_sorting.reverse()
    xl_ranked = []
    for score, n in for_sorting:
        print score, n
        res_lower = xl_graph.node[n]['xl'][0]
        res_upper = xl_graph.node[n]['xl'][1]
        xl_ranked.append((res_lower,'CA', res_upper,'CA', score))
    InputOutput.InputOutput.write_contact_file(xl_ranked, 'foo', upper_distance = 20)

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
                #xl_graph.add_edge(n[0],o[0])
            #if is_neighbourhood(n[1]['xl'], o[1]['xl'], delta=2, double=False):
            #    connecting_nodes.append(o)
    #for n in xl_graph.nodes(data=True):
    #    connecting_nodes =[]
#	   for o in xl_graph.nodes(data=True):
#	      if is_neighbourhood(n[1]['xl'], o[1]['xl'], delta=2, double=False):#
#	         connecting_nodes.append(o)	
#	   if len(connecting_nodes) >= 2:
#	      for o in connecting_nodes:
#	         g.add_edge(n[0],o[0])

def toy_graph():
    y = nx.Graph()
    y.add_node(1)
    y.add_node(2)
    y.add_node(3)
    #y.add_node(4)
    y.add_edge(1,2)
    y.add_edge(2,3)
    #y.add_edge(3,4)
    return y

def build_xl_graph( xl_data ):
   
    g = nx.Graph()
    index = 1
    pers = {}
    # Add cross-links as nodes
    for i, score in xl_data:
        g.add_node(index, xl=i)#, weight = score)
        pers[index] = score
        index += 1
    # Add edges between if cross-links are within delta distance
    for n in g.nodes(data=True):
        for o in g.nodes(data=True):
            if o[0] > n[0]:
                if is_neighbourhood(n[1]['xl'], o[1]['xl'], delta=3):
                    g.add_edge(n[0],o[0])#, weight =  numpy.max([n[1]['weight'], o[1]['weight']])  )

    #add_loops( g )
    for i in xrange(0,2):
        to_add = []
        for n in g.nodes(data=True):
            for o in g.nodes(data=True):
                if o[0] > n[0]:
                    if has_loop(o, n, g):
                        to_add.append((o[0],n[0]))
                        #g.add_edge(n[0],o[0])
        for i,j in to_add:
            g.add_edge(i,j)
    return g, pers

def main():        
   """Generic main function. Executes main functionality of program
   """
   tg = toy_graph()
   #has_loop(2,4,tg)
   #sys.exit()
   #rint tg.edges()
   #pagerank(tg)
   #sys.exit()
   xl_data, gt_data = InputOutput.InputOutput.load_xl_data( options.example,options.offset )
   InputOutput.InputOutput.write_contact_file(gt_data, 'gt', upper_distance = 20)

   #sys.exit()
   #break
   xl_graph,pers = build_xl_graph(xl_data)
   for n in xl_graph.nodes(data=True):
       print xl_graph.degree(n[0]), pers[n[0]], n[0]
   #for n in xl_graph.edges(data=True):
   #    print n
   #print  xl_graph.number_of_edges()
   do_page_rank(xl_graph,pers)

if __name__ == '__main__':
    sys.exit(main())
else:
    print("Loaded as a module!")




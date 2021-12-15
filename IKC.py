import argparse
import networkit as nk
import csv
import heapq as hq

def main(args):

    edge_list = args.edgeList
    out_dir = args.outDir
    k = args.kvalue

    node_dict, edge_list = process_edgelist(edge_list, out_dir)

    edge_list_reader = nk.graphio.EdgeListReader('\t', 0, continuous=False, directed=True)
    graph1 = edge_list_reader.read(edge_list)

    graph, node_id_dict = format_graph(graph1)

    clusters = iterative_k_core_decomposition_MCS_ES(graph, k, node_id_dict)
    print ("IKC printing to", out_dir)
    print_clusters(clusters, out_dir, node_dict)


def print_clusters(clusters, out_dir, node_dict):
    '''
    This writes a csv containing lines with the:
    node Id, cluster nbr, and value of k for which cluster nbr was generated
    INPUT
    -----
    clusters : a list of clusters represented as lists of nodes in each cluster
    outDir : the file path and name of the ouput
    '''
    # the index indicates the order for when the cluster number was generated
    index = 0
    k = 0
    with open("{}".format(out_dir), "w") as output:
        for cluster_info in clusters:
            (cluster, k, modularity_score) = cluster_info

            # print a separate line for each node in each cluster
            for node in cluster:
                output.write('{}\t{}\n'.format(node_dict[node], index))
            
            index += 1


def iterative_k_core_decomposition_MCS_ES(graph, k, inverted_orig_node_ids):
    '''
    INPUT
    -----
    graph                  : the full networkit graph
    k                      : the minimum allowed value for k for valid clusters
    inverted_orig_node_ids : the dictionary mapping the compacted node IDs to the original node IDs
    OUTPUT
    ------
    final_clusters : the clustering output, a list of lists with clustered nodes
    '''
    orig_graph = nk.graphtools.subgraphFromNodes(graph, graph.iterNodes())
    L = orig_graph.numberOfEdges()
    singletons = []
    final_clusters = []

    nbr_failed_modularity = 0
    nbr_failed_k_valid = 0

    # continue finding clusters for different values of k until
    # a. there are no nodes left in the garph or
    # b. the maximum value of k is lower than the minumum allowed k for valid clusters
    while graph.numberOfNodes() > 0:

        # run kc on the graph to get the smallest kcore
        subgraph, max_k, kcore = kc(graph)
        if subgraph == None:
            print ("no available subgraph")
            break

        # if b. above is true, add all singletons and nodes left in the graph as individual clusters
        # and break
        if max_k < k:
            for node in graph.iterNodes():
               modularity = (-1)*(orig_graph.degree(inverted_orig_node_ids[node])/(2*L))**2
               final_clusters.append(([inverted_orig_node_ids[node]],0,modularity))
            for node in singletons:
               final_clusters.append(([node],0,0))
            break

        # compute the components
        cc = nk.components.WeaklyConnectedComponents(subgraph)
        cc.run()
        components = cc.getComponents()

        nodes_to_remove = set()

        # check components to make sure they are k-valid and m-valid 
        # then if so, add them to a cluster or break them up to make k-valid
        # finally add them to the final_clusters and remove those nodes from the graph
        for component in components:
            clusters = []
            valid_modularities = []

            # ensure the component is k_valid and modular
            # if not remode the nodes from the graph and add to the singletons
            if k_valid(component, subgraph, k):
                modularity = modular(component, orig_graph, inverted_orig_node_ids)
                if modularity > -1:
                    sub_components = [(component, modularity)]
                else:
                    print('failed modularity')
                    nbr_failed_modularity += 1
                    sub_components = []
                    nodes_to_remove.update(component)
                    singletons.extend(orig_id_component(component, inverted_orig_node_ids))
            else:
                print('failed k-valid')
                nbr_failed_k_valid += 1
                sub_components = []
                nodes_to_remove.update(component)
                singletons.extend(orig_id_component(component, inverted_orig_node_ids))

            # retreive the original node ID for the clustering output and add each new component as its own cluster
            for sub_component, modularity in sub_components:
                nodes_to_remove.update(sub_component)
                cluster = []
                for node in sub_component:
                    cluster.append(inverted_orig_node_ids[node])
                print ('adding cluster length', len(cluster))
                clusters.append((cluster, max_k, modularity))

            final_clusters.extend(clusters)

        # just prints information about the number of components to standard output
        component_sizes = cc.getComponentSizes()
        nbr_large_components = 0
        large_components = dict()
        for nbr,size in component_sizes.items():
            if size > 100:
                nbr_large_components += 1
                large_components[nbr] = size
        print ('nbr components:', len(components),
               ',  nbr components with more than 100 nodes:', nbr_large_components)

        # remove nodes marked for removal
        # (either already clustered or when a large cluster is not validly broken up)
        for node in nodes_to_remove:
            graph.removeNode(node)
            #print ("removing node: ", inverted_orig_node_ids[node])

        # compact the subgraph with continuous node IDs (needed for partitioning in kc)
        # and update inverted_orig_node_ids with the new IDs
        node_id_dict = nk.graphtools.getContinuousNodeIds(graph)
        inverted_node_id_dict = dict(map(reversed, node_id_dict.items()))
        temp_dict = dict()
        for new_id, old_id in inverted_node_id_dict.items():
            orig_id = inverted_orig_node_ids[old_id]
            temp_dict[new_id] = orig_id
        inverted_orig_node_ids = temp_dict
        graph = nk.graphtools.getCompactedGraph(graph, node_id_dict)
        print ('nodes left in graph: ', graph.numberOfNodes())

    print ("nbr of clusters which were rejected since they were not k-valid : ", nbr_failed_k_valid)
    print ("nbr of clusters which were rejected since they were not modular : ", nbr_failed_modularity)

    return final_clusters


def kc(graph, k=None):
    '''
    INPUT
    -----
    graph : networkit graph
    k     : the minimum node degree for all nodes in the subgraph
    OUTPUT
    ------
    subgraph : the subgraph containing only nodes of degree k or higher
    max_k    : the the largest value of k for which there are still nodes remaining in the subgraph
    kc       : the core decomposition networkit object
    '''

    # get the value of k associated with each node in the graph and store in a partition
    kc = nk.centrality.CoreDecomposition(graph, storeNodeOrder=True)
    kc.run()
    partition = kc.getPartition()

    kcore_members = []
    max_k = kc.maxCoreNumber()

    # default to the maximum k value
    if k == None:
       k = max_k

    save_k = k

    # if the value of k is lower than the minimum allowed exit
    if max_k < k:
        return None, max_k, kc

    # populate the kcore members with all nodes with k values between k and max_k
    while k <= max_k:
        kcore_members.extend(partition.getMembers(k))
        k += 1

    # return the subgraph with nodes from kcore members
    print ("k value", save_k,'nbr core members', len(kcore_members))
    return nk.graphtools.subgraphFromNodes(graph, kcore_members), max_k, kc


def k_valid(component, subgraph, k):
    #subgraph = nk.graphtools.subgraphFromNodes(graph, component)
    k_valid = True
    component_nodes = set(component)
    for node in subgraph.iterNodes():
        if node in component_nodes:
            if (subgraph.degreeIn(node) + subgraph.degreeOut(node)) < k:
                #print ("node", node)
                #print('fails k condition', subgraph.degree(node))
                k_valid = False
                break
    return k_valid


def modular(component, orig_graph, inverted_orig_node_ids):
    cluster = nk.graphtools.subgraphFromNodes(orig_graph, orig_id_component(component, inverted_orig_node_ids))

    l = orig_graph.numberOfEdges()
    ls = cluster.numberOfEdges()
    ds = 0

    for node in cluster.iterNodes():
       ds += orig_graph.degreeIn(node)
       ds += orig_graph.degreeOut(node)

    return (ls/l - (ds/(2*l))**2)


def orig_id_component(component, inverted_orig_node_ids):
    return [inverted_orig_node_ids[node] for node in component]


def format_graph(graph1):
    origNodeIdDict = nk.graphtools.getContinuousNodeIds(graph1)
    invertedOrigNodeIdDict = dict(map(reversed, origNodeIdDict.items()))
    graph = nk.graphtools.getCompactedGraph(graph1, origNodeIdDict)
    print("number of edges before", graph.numberOfEdges(), graph.isWeighted())

    edge_cnt = dict()
    for node1, node2, weight in graph.iterEdgesWeights():
        n = weight * 1
        #print (weight)
        edge_cnt[(node1, node2)] = n
    added = 0
    for (node1, node2), cnt in edge_cnt.items():
        graph.removeEdge(node1, node2)
        for x in range(int(round(cnt))):
            added += 1
            graph.addEdge(node1, node2, 1)

    graph.removeSelfLoops()

    print (graph.numberOfNodes(), graph.numberOfEdges(), added)

    return graph, invertedOrigNodeIdDict


def process_edgelist(edge_list, outpath):

    f = open(edge_list, "r")
    newEdgeListPath = "{}.newEdgeList".format(outpath)
    nf = open(newEdgeListPath, "w")
    index = 0
    node_dict = dict()
    for line in f:
        node1, node2, weight = line.split('\t')
        if node1 in node_dict:
            new_node1 = node_dict[node1]
        else:
            new_node1 = index
            node_dict[node1] = new_node1
            index += 1

        if node2 in node_dict:
            new_node2 = node_dict[node2]
        else:
            new_node2 = index
            node_dict[node2] = new_node2
            index += 1

        nf.write("{}\t{}\t{}\n".format(new_node1, new_node2, weight))
    nf.close()
    f.close()

    return dict(map(reversed, node_dict.items())), newEdgeListPath


def parseArgs():
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--edgeList", type=str,
                        help="Path to file containing edge lists",
                        required=True, default=None)

    parser.add_argument("-o", "--outDir", type=str,
                        help="Path to file containing output",
                        required=True, default=None)

    parser.add_argument("-k", "--kvalue", type=int,
                        help="non-negative integer value of the minimum required adjacent nodes for each node",
                        required=False, default=1)

    parser.add_argument("-v", "--version", action="version", version="1.0.0",
                        help="show the version number and exit")

    return parser.parse_args()

if __name__ == "__main__":
    main(parseArgs())

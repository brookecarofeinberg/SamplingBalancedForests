import networkx as nx
import json
from networkx.readwrite import json_graph
import matplotlib.pyplot as plt
import numpy as np
import random
import networkx as nx
from collections import deque

class WilsonsAlgorithm:
    def __init__(self, graph):
        self.graph = graph  #store networkX graph 
        self.adj_list = {v: list(graph.neighbors(v)) for v in graph.nodes() if list(graph.neighbors(v))}  #for each graph node return the neighbors and store it as a dictionary where node is the key and each value corresponds to a list of its neighbors 
        self.InTree = {node: False for node in self.adj_list}  #initalizes dictionary to keep track of whether nodes have been included in the spanning tree -- FALSE as initial value menaing node has yet to be included in the tree 
        self.Next = {node: None for node in self.adj_list}  #initalizes dictionary to track path each node takes in a random walk 
    def random_edge(self, v): #randomly selects an edge for a given node v 
        return random.choice(self.adj_list[v]) #retrieves the list of neighbors for node v in the adjacency list and randomly picks one -- simulating step ina  random walk 
    def RandomTreeWithRoot(self, r): #generate random spanning tree with r as the root node 
        self.InTree = {node: False for node in self.adj_list} #initalize all nodes as not in the spanning tree 
        self.InTree[r] = True  #mark root node as in the tree 
        self.Next[r] = -1  #root node has no parent aka no incoming edges 
        for node in self.adj_list: #loops through each nodein the graph, attempting to find a path to the root using a random walk 
            path = []  
            while not self.InTree[node]: #loop continues as long as node is not already in the spanning tree 
                path.append(node)  #record each step of the random walk 
                self.Next[node] = self.random_edge(node) #choose a random neighbor & record it as the next step in the walk (connects node to its next random neighbor in the path)
                node = self.Next[node]  #update next node and the new node so the loop continues till we reach a node already in the spanning tree 
            for vertex in path:  #loop has broken bc we've reached a node already part of the spanning tree 
                self.InTree[vertex] = True  #adds each vertex to the spanning tree from the root outward 
    def sample(self, seed=None):
        if seed is not None: #can include seed later on for reproducible results 
            random.seed(seed)
        root = random.choice(list(self.adj_list.keys()))  #picks a random node to be a root 
        self.RandomTreeWithRoot(root) #builds ranodm spannign tree with the ranodm root chosen 

        spanning_tree = []  #extract edges of spannign tree for plotting 
        for node in self.adj_list: 
            if self.Next[node] != -1 and self.Next[node] is not None: #find nodes connected to another node in the tree 
                spanning_tree.append((node, self.Next[node])) #addes connection or edge as a tuple of two nodes to the spanning tree list 
        return spanning_tree, root #returns all the edges in the form of (node1, node2) and the root node 
    def draw_tree(self, spanning_tree, root, node_locations): #plot the spanning tree on the graph
        plt.figure(figsize=(8, 8))
        nx.draw(self.graph, pos=node_locations, node_size=10, edge_color='lightgray') #draw original graph 
        nx.draw_networkx_edges(self.graph, pos=node_locations, edgelist=spanning_tree, edge_color='red', width=2) #draw spanning tree 
        nx.draw_networkx_nodes(self.graph, pos=node_locations, nodelist=[root], node_color='yellow', node_size=50)  #draw root node 
        plt.show()

class GraphPartitioner:
    def __init__(self, graph, k):
        self.graph = graph
        self.k = k  #choose k for confidence interval based on the Gamma Bernoulli approximation scheme 
    def create_weighted_partitions_and_draw_spanning_tree(self):
        while True:
            wilson = WilsonsAlgorithm(self.graph) #generate random spannign tree using wilson's algorithm 
            spanning_tree, root = wilson.sample() 
            tree_graph = nx.Graph() 
            tree_graph.add_edges_from(spanning_tree)  #take edges form spannign tree and add them to the new graph
            total_nodes = len(tree_graph) 
            ideal_size = total_nodes // self.k
            extra_nodes = total_nodes % self.k
            target_sizes = [ideal_size + 1 if i < extra_nodes else ideal_size for i in range(self.k)] #accomedate odd nodes 
            partitions = [[] for _ in range(self.k)] #initalize BFS queues 
            partition_sizes = [0] * self.k
            bfs_queues = [deque([root])] + [deque() for _ in range(self.k - 1)] 
            visited = set([root])
            partitions[0].append(root)
            partition_sizes[0] += 1
            for current_partition in range(self.k): #grow partitions with BFS
                while partition_sizes[current_partition] < target_sizes[current_partition] and bfs_queues[current_partition]:
                    node = bfs_queues[current_partition].popleft()
                    for neighbor in tree_graph.neighbors(node):
                        if neighbor not in visited:
                            partitions[current_partition].append(neighbor)
                            partition_sizes[current_partition] += 1
                            visited.add(neighbor)
                            bfs_queues[current_partition].append(neighbor)
                            if partition_sizes[current_partition] >= target_sizes[current_partition]:
                                break
                if partition_sizes[current_partition] < target_sizes[current_partition]: # check if partition has reached the right size--else continue the loop for a new tree
                    return False         
                unvisited_nodes = [n for n in tree_graph.nodes if n not in visited] # initalize next partition with an unvisited node 
                if unvisited_nodes and current_partition + 1 < self.k:
                    seed_node = unvisited_nodes[0]
                    partitions[current_partition + 1].append(seed_node)
                    partition_sizes[current_partition + 1] += 1
                    visited.add(seed_node)
                    bfs_queues[current_partition + 1].append(seed_node)
            if all(partition_sizes[i] == target_sizes[i] for i in range(self.k)): # check if all partitions are the right size
                return True

#plot of spannig tree for visualization 
node_locations = {v: (float(cnty.nodes()[v]["INTPTLON20"]) - 360 if float(cnty.nodes()[v]["INTPTLON20"]) > 0 else float(cnty.nodes()[v]["INTPTLON20"]),
        float(cnty.nodes()[v]["INTPTLAT20"]))for v in cnty.nodes()}
wilson = WilsonsAlgorithm(cnty)
spanning_tree, root = wilson.sample()
wilson.draw_tree(spanning_tree, root, node_locations)

def run_partition_simulations(graph, p, k): #Splitability analysis based on gamma bernoulli approx scheme 
    successful_splits = 0
    trials = 0
    while successful_splits < k:
        trials += 1
        partitioner = GraphPartitioner(graph, p)
        if partitioner.create_weighted_partitions_and_draw_spanning_tree():
            successful_splits += 1
    success_rate= (k/trials)*100
    print(f"It took {trials} trials to achieve {k} successful splits for {p} partitions.")
    print(f"average success rate is {success_rate}")
import networkx as nx
import numpy as np
from machine_state import MachineState
from utils import *
from route import *
from schedule import *
from machine import Trap, Segment, Junction
import random
import itertools as it

class RebalanceTraps:
    def __init__(self, machine, system_state):
        self.machine = machine
        self.ss = system_state
   
    def set_trap_check_order(self, demand, trap_free_space):
        """Checks candidate traps.
        A candidate trap has >1 free spaces and
        closest to traps with negative demand
        """
        m = self.machine
        # > 1 free space traps
        candidate_traps = [k for k in self.ss.trap_ions if trap_free_space[k] > 1]
        print(candidate_traps)

        # find closest
        neg_demand_traps = [k for k in self.ss.trap_ions if demand.get(m.traps[k], 0) == -1]
        print(neg_demand_traps)

        distances = {}
        # for dtrap in candidate_traps:
        #     dist = 0
        #     for strap in neg_demand_traps:
        #         dist += abs(dtrap - strap)

        #     distances[dtrap] = dist

        for src_trap in neg_demand_traps:
            d = {}
            for dest_trap in candidate_traps:
                dist = abs(src_trap - dest_trap)
                d[dest_trap] = dist

            print(d)
            dest = sorted(d, key=d.get)[0]
            print(dest)

            distances[src_trap] = dest

        prods = list(it.product(neg_demand_traps, candidate_traps))
        dists = {tup: abs(tup[0]-tup[1]) for tup in prods}

        sorted_src_dest = sorted(dists, key=dists.get)

        ftr = FreeTrapRoute(self.machine, self.ss)
        for source, destination in sorted_src_dest:
            status, route = ftr.find_route(source_trap=source, dest_trap=destination)
            if not status:
                shortest_rebal_path = (source, destination, route)
                break
            

        # min_dist = min(distances)

        # candidate_traps = [k for k, elem in enumerate(distances) if min_dist == elem]

        # candidate_traps = sorted(distances, key=distances.get)

        return distances, shortest_rebal_path

    
    def clear_all_blocks(self):
        m = self.machine
        graph = nx.DiGraph(m.graph)
        demand = {}
        weight = {}
        capacity = {}
        ss = self.ss
        trap_free_space = {}
        print(f"trap capacity {m.traps[0].capacity}")
        for k in self.ss.trap_ions:
            trap_free_space[k] = m.traps[k].capacity - len(ss.trap_ions[k])
        req_free_space = 0
        for k in self.ss.trap_ions:
            #If a trap is blocked, remove one ion from it
            if trap_free_space[k] == 0:
                req_free_space += 1
                demand[m.traps[k]] = -1
        
        # set iteration order
        distances, rebal_path = self.set_trap_check_order(demand, trap_free_space)
        print(f"distances {distances}")
        print(f"Shortest rebal path T{rebal_path[0]} -> T{rebal_path[1]}")

        return None, None, rebal_path

        """Old"""
        # for k in self.ss.trap_ions:
        # # for k in dest_traps:
        #     #If ions need to be moved, and this trap has 2 or more spaces, accepts ions
        #     if req_free_space != 0 and trap_free_space[k] > 1:
        #         offer = min(trap_free_space[k]-1, req_free_space)
        #         req_free_space -= offer
        #         demand[m.traps[k]] = offer #This trap accepts ions

        print(f"demand: {demand}")
        
        """New"""
        for src_traps, dest_trap in distances.items():
            # dest_trap = distances[src_traps]
            print(f"src: {src_traps} | dest: {dest_trap}")

            if m.traps[dest_trap] in demand:
                demand[m.traps[dest_trap]] += 1
            else:
                demand[m.traps[dest_trap]] = 1

        print("Demands")
        trap_ids_with_neg_demand = []
        for item in demand:
            print("T"+str(item.id), demand[item])
            if demand[item] < 0:
                trap_ids_with_neg_demand.append(item.id)
        nx.set_node_attributes(graph, demand, 'demand')
        for u, v in graph.edges:
            # print(f"u {u.id} v {v.id}")
            weight[(u,v)] = 1
            capacity[(u,v)] = 100
        nx.set_edge_attributes(graph, weight, 'weight')
        nx.set_edge_attributes(graph, capacity, 'capacity')

        flowCost, flowDict = nx.network_simplex(graph)
        print(f"flow cost {flowCost}")
        # print(f"Flow dict {flowDict}")
        # for i in flowDict:
        #     for j in flowDict[i]:
        #         print(i, j, flowDict[i][j])

        return flowDict, trap_ids_with_neg_demand, rebal_path

    def xclear_all_blocks(self):
        m = self.machine
        graph = nx.DiGraph(m.graph)
        demand = {}
        weight = {}
        capacity = {}
        ss = self.ss
        trap_free_space = {}
        for k in self.ss.trap_ions:
            trap_free_space[k] = m.traps[k].capacity - len(ss.trap_ions[k])
        req_free_space = 0
        for k in self.ss.trap_ions:
            #If a trap is blocked, remove one ion from it
            if trap_free_space[k] == 0:
                req_free_space += 1
                demand[m.traps[k]] = -1
        for k in self.ss.trap_ions:
            #If ions need to be moved, and this ion has 2 or more spaces, accepts ions
            if req_free_space != 0 and trap_free_space[k] > 1:
                offer = min(trap_free_space[k]-1, req_free_space)
                req_free_space -= offer
                demand[m.traps[k]] = offer #This trap accepts ions
        print("Demands")
        trap_ids_with_neg_demand = []
        for item in demand:
            print("T"+str(item.id), demand[item])
            if demand[item] < 0:
                trap_ids_with_neg_demand.append(item.id)
        nx.set_node_attributes(graph, demand, 'demand')
        num_edges = 0
        for u, v in graph.edges:
            weight[(u,v)] = 1
            capacity[(u,v)] = 100
            num_edges += 1 
        nx.set_edge_attributes(graph, weight, 'weight')
        nx.set_edge_attributes(graph, capacity, 'capacity')
        
        for node in graph.nodes:
            print(f"node type {type(node)} | id {node.id} | {graph.nodes[node]}")
        
        for edge in graph.edges:
            print(f"{type(edge[0])} {edge[0].id} {type(edge[1])} {edge[1].id} | {graph.edges[edge]}")
        flowCost, flowDict = nx.network_simplex(graph)

        for src, dest in flowDict.items():
            print(src.show(), {key.show(): val for key, val in dest.items()})

        return flowDict, trap_ids_with_neg_demand
    
    def clear_route(self, trap_list, route):
        #Set up node demands
        #Negative demand means node wants to send flow 
        #heuristic: For each blocked node on a path from i1 to i2, create negative demand
        #for free nodes outside this path if available, create positive demand
        assert 0 


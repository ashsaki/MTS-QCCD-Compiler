import networkx as nx
import numpy as np
from machine_state import MachineState
from utils import *
from route import *
from schedule import *
from machine import Trap, Segment, Junction
import random

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

        # find closest
        neg_demand_traps = [k for k in self.ss.trap_ions if demand.get(m.traps[k], 0) == -1]

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

            dest = sorted(d, key=d.get)[0]

            distances[src_trap] = dest



        # min_dist = min(distances)

        # candidate_traps = [k for k, elem in enumerate(distances) if min_dist == elem]

        # candidate_traps = sorted(distances, key=distances.get)

        return distances

    
    def clear_all_blocks(self):
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
        
        # set iteration order
        distances = self.set_trap_check_order(demand, trap_free_space)
        # print(f"candidate traps {distances}")

        """Old"""
        # for k in self.ss.trap_ions:
        # for k in dest_traps:
        #     #If ions need to be moved, and this trap has 2 or more spaces, accepts ions
        #     if req_free_space != 0 and trap_free_space[k] > 1:
        #         offer = min(trap_free_space[k]-1, req_free_space)
        #         req_free_space -= offer
        #         demand[m.traps[k]] = offer #This trap accepts ions

        """New"""
        for src_traps in distances:
            dest_trap = distances[src_traps]

            if m.traps[dest_trap] in demand:
                demand[m.traps[dest_trap]] += 1
            else:
                demand[m.traps[dest_trap]] = 1

        # print("Demands")
        # for item in demand:
        #    print("T"+str(item.id), demand[item])
        nx.set_node_attributes(graph, demand, 'demand')
        for u, v in graph.edges:
            weight[(u,v)] = 1
            capacity[(u,v)] = 100
        nx.set_edge_attributes(graph, weight, 'weight')
        nx.set_edge_attributes(graph, capacity, 'capacity')

        flowCost, flowDict = nx.network_simplex(graph) 
        return flowDict

    def clear_route(self, trap_list, route):
        #Set up node demands
        #Negative demand means node wants to send flow 
        #heuristic: For each blocked node on a path from i1 to i2, create negative demand
        #for free nodes outside this path if available, create positive demand
        assert 0 


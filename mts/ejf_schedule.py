'''
Schedules gates in topologically sorted order
Policy is called EJF Earliest Job First, because its similar to the famous job scheduling policy

done - handle case where both traps are full
done - handle junction traffic flow
todo - skip - handle preemption??
todo - skip - handle traffic based routing
'''
import copy
import random
import networkx as nx
import numpy as np
import math
from machine_state import MachineState
from utils import *
from route import *
from schedule import *
from machine import Trap, Segment, Junction
from rebalance import *
import matplotlib.pyplot as plt

class EJFSchedule:
    #Inputs are
    #1. gate dependency graph - IR
    #2. gate_info = what are the qubits used by a two-qubit gate?
    #3. M = machine object
    #4. init_map = initial qubit mapping
    def __init__(
        self,
        ir,
        gate_info,
        M,
        init_map,
        serial_trap_ops,
        serial_comm,
        global_serial_lock,
        **kwargs
    ):
        self.ir = ir
        self.gate_info = gate_info
        self.machine = M
        self.init_map = init_map
        self.gate_order = kwargs['gate_order']
        self.lookahead = True if kwargs['lookahead'] == 1 else False
        self.comm_cap = kwargs['comm_cap']
        self.moves = 0
        self.active_gate = None
        self.completed_gate = []
        self.layers = kwargs['layers']
        self.gate_to_layer = self.prep_gate_to_layer() # dict[gate id, layer]
        self.sx = True
        self.neighbors = {
            0: [1],
            1: [0, 2],
            2: [1, 3],
            3: [2, 4],
            4: [3, 5],
            5: [4]
        }
        self.free_up = True

        #Setup scheduler
        self.machine.add_comm_capacity(self.comm_cap) #Add space for 2 extra ions in each trap
        self.SerialTrapOps = serial_trap_ops # ---> serializes operations on a single trap zone
        self.SerialCommunication = serial_comm # ---> serializes all split/merge/move ops
        self.GlobalSerialLock = global_serial_lock #---> serialize gates and comm

        #SerialTrapOps enforces that all operations in a trap are serialized
        #i.e., no parallel gates in a single ion chain/region
        self.schedule = Schedule(M)
        self.router = BasicRoute(M)
        self.gate_finish_times = {}

        #Some scheduling statistics
        #Count the number of times we had to clear some traps because of traffic blocks
        self.count_rebalance = 0
        self.split_swap_counter = 0

        #Create the sys_stage object which is used to track system state
        #from the perspective of the scheduler
        trap_ions = {}
        seg_ions = {}
        for i in M.traps:
            if init_map[i.id]:
                trap_ions[i.id] = init_map[i.id][:]
            else:
                trap_ions[i.id] = []
        for i in M.segments:
            seg_ions[i.id] = []
        self.sys_state = MachineState(0, trap_ions, seg_ions)
        #self.sys_state.print_state()

    def prep_gate_to_layer(self):
        d = {}
        for layer_id, gates in enumerate(self.layers):
            for gate_id in gates:
                d[gate_id] = layer_id

        return d

    
    #Find the earliest time at which a gate can be scheduled
    #Earliest time = max(dependent gate times)
    def gate_ready_time(self, gate):
        ready_time = 0
        for in_edge in self.ir.in_edges(gate):
            #Each in edge is of the form (in_gate_id, this_gate_id)
            in_gate = in_edge[0]
            if in_gate in self.gate_finish_times:
                ready_time = max(ready_time, self.gate_finish_times[in_gate])
            else:
                print("Error: Finish time of depenedent gate not found", in_edge)
                assert 0
        return ready_time

    #Find the time at which a particular qubit/ion is ready for another operation
    def ion_ready_info(self, ion_id):
        s = self.schedule
        this_ion_ops = s.filter_by_ion(s.events, ion_id)
        this_ion_last_op_time = 0
        this_ion_trap = None

        #If there is some operation that has happened for this ion:
        if len(this_ion_ops):
            #The last operation on an ion is either a gate or a merge in a trap
            assert (this_ion_ops[-1][1] == Schedule.Gate) or (this_ion_ops[-1][1] == Schedule.Merge)
            #Pick up the time and location of last operation
            this_ion_last_op_time = this_ion_ops[-1][3]
            this_ion_trap = this_ion_ops[-1][4]['trap']
        else:
            #Find which trap originally held this ion
            #It shouldn't have changed because no ops have happened for this ion
            did_not_find = True
            for trap_id in self.init_map.keys():
                if ion_id in self.init_map[trap_id]:
                    this_ion_trap = trap_id
                    did_not_find = False
                    break
                    
            if did_not_find:
                print("Did not find:", ion_id)
            assert (did_not_find == False)

        #Double checking ion location from sys_state object
        if this_ion_trap != self.sys_state.find_trap_id_by_ion(ion_id):
            print(ion_id, this_ion_trap, self.sys_state.find_trap_id_by_ion(ion_id))
            self.sys_state.print_state()
            assert 0
        return this_ion_last_op_time, this_ion_trap

    #Add a split operation to the current schedule
    def add_split_op(self, clk, src_trap, dest_seg, ion):
        m = self.machine
        if self.SerialTrapOps == 1:
            last_event_time_on_trap = self.schedule.last_event_time_on_trap(src_trap.id)
            split_start = max(clk, last_event_time_on_trap)
        else:
            split_start = clk

        if self.SerialCommunication == 1:
            last_comm_time = self.schedule.last_comm_event_time()
            split_start = max(split_start, last_comm_time)

        if self.GlobalSerialLock == 1:
            last_event_time_in_system = self.schedule.get_last_event_ts()
            split_start = max(split_start, last_event_time_in_system)
        split_duration, split_swap_count, split_swap_hops, i1, i2, ion_swap_hops, self.sys_state = m.split_time(self.sys_state, src_trap.id, dest_seg.id, ion)
        self.split_swap_counter += split_swap_count
        split_end = split_start + split_duration
        self.schedule.add_split_or_merge(split_start, split_end, [ion], src_trap.id, dest_seg.id, Schedule.Split, split_swap_count, split_swap_hops, i1, i2, ion_swap_hops)
        return split_end

    #Add a merge operation to the current schedule
    def add_merge_op(self, clk, dest_trap, src_seg, ion):
        m = self.machine
        if self.SerialTrapOps == 1:
            last_event_time_on_trap = self.schedule.last_event_time_on_trap(dest_trap.id)
            merge_start = max(clk, last_event_time_on_trap)
        else:
            merge_start = clk

        if self.SerialCommunication == 1:
            last_comm_time = self.schedule.last_comm_event_time()
            merge_start = max(merge_start, last_comm_time)

        if self.GlobalSerialLock == 1:
            last_event_time_in_system = self.schedule.get_last_event_ts()
            merge_start = max(merge_start, last_event_time_in_system)
        merge_end = merge_start + m.merge_time(dest_trap.id)
        self.schedule.add_split_or_merge(merge_start, merge_end, [ion], dest_trap.id, src_seg.id, Schedule.Merge, 0, 0, 0, 0, 0)
        return merge_end

    #Add a move operation to the current schedule
    #This is one segment to segment move i.e., move ion from src_seg to dest_seg
    def add_move_op(self, clk, src_seg, dest_seg, junct, ion):
        m = self.machine
        move_start = clk
        if self.GlobalSerialLock == 1:
            last_event_time_in_system = self.schedule.get_last_event_ts()
            move_start = max(move_start, last_event_time_in_system)

        if self.SerialCommunication == 1:
            last_comm_time = self.schedule.last_comm_event_time()
            move_start = max(move_start, last_comm_time)

        move_end = move_start + m.move_time(src_seg.id, dest_seg.id) + m.junction_cross_time(junct)
        move_start, move_end = self.schedule.junction_traffic_crossing(src_seg, dest_seg, junct, move_start, move_end)
        self.schedule.add_move(move_start, move_end, [ion], src_seg.id, dest_seg.id)
        return move_end

    #Add a gate operation to the current schedule
    def add_gate_op(self, clk, trap_id, gate, ion1, ion2):
        fire_time = clk
        if self.SerialTrapOps == 1:
            last_event_time_on_trap = self.schedule.last_event_time_on_trap(trap_id)
            fire_time = max(clk, last_event_time_on_trap)

        if self.GlobalSerialLock == 1:
            last_event_time_in_system = self.schedule.get_last_event_ts()
            fire_time = max(fire_time, last_event_time_in_system)
        gate_duration = self.machine.gate_time(self.sys_state, trap_id, ion1, ion2)
        self.schedule.add_gate(fire_time, fire_time + gate_duration, [ion1, ion2], trap_id)
        self.gate_finish_times[gate] = fire_time + gate_duration
        return fire_time + gate_duration

    # Count future ops of each ion in current traps and decided the shuttle direction
    def count_future_ops(
        self,
        ion1,
        ion2,
        ion1_trap,
        ion2_trap,
        excess_cap1,
        excess_cap2,
        fire_time=None,
        rem_gates=None,
        mng_deadlock=True
    ):
        # if rem_gates is None:
        rem_gates = self.remaining_gates
        fin_time_new = 0


        ion1_trap1 = 0
        ion1_trap2 = 0

        ion2_trap1 = 0
        ion2_trap2 = 0

        w = 1
        check_gates = len(rem_gates) // w
        
        ntraps = len(self.machine.traps)
        ops_in_traps = np.zeros(ntraps, dtype=int)
        n_prev = 0
        
        
        threshold = 6 # window size for checking future gates
        print(f"gate [{ion1}, {ion2}] traps {ion1_trap}, {ion2_trap} th {threshold}")
        for n, g in enumerate(rem_gates[:check_gates]):
            idx = 0 # not used --> `idx` should be layer index. Need `def layer_idx(g)` also need `front_layer_idx`
            gate_ions = copy.deepcopy(self.gate_info[g])
            if ion1 in gate_ions:
                if (n-n_prev) > threshold:
                    print(f"breaking from ion1 n {n} n_prev {n_prev}")
                    break
                n_prev = n
                
                gate_ions.remove(ion1)
                other_ion = gate_ions[0]
                
                _, other_ion_trap = self.ion_ready_info(other_ion)
                print(f"{n} gate ions [{ion1}, {other_ion}] traps {ion1_trap} {other_ion_trap}")
                ops_in_traps[ion1_trap] += 1
                ops_in_traps[other_ion_trap] += 1

                if other_ion_trap == ion1_trap:
                    ion1_trap1 += 1
                    
                    idx += 1
                elif other_ion_trap == ion2_trap:
                    ion1_trap2 += 1
            elif ion2 in gate_ions:
                if (n-n_prev) > threshold:
                    print(f"breaking from ion2 n {n} n_prev {n_prev}")
                    break
                n_prev = n
                
                gate_ions.remove(ion2)
                
                other_ion = gate_ions[0]
                _, other_ion_trap = self.ion_ready_info(other_ion)
                print(f"{n} gate ions [{ion2}, {other_ion}] traps {ion2_trap} {other_ion_trap}")
                ops_in_traps[ion2_trap] += 1
                ops_in_traps[other_ion_trap] += 1

                if other_ion_trap == ion2_trap:
                    ion2_trap2 += 1
                    
                    idx += 1
                elif other_ion_trap == ion1_trap:
                    ion2_trap1 += 1

        ion2_move_score = ion1_trap1 + ion2_trap1
        ion1_move_score = ion2_trap2 + ion1_trap2
        
        print(f"ion1_trap1 {ion1_trap1} ion2_trap2 {ion2_trap2}")
        print(f"ion1 move score {ion1_move_score} ion2 move score {ion2_move_score}")

        
        if ion2_move_score > ion1_move_score:
            src = ion2_trap
            dest = ion1_trap
        elif ion2_move_score == ion1_move_score and not (excess_cap1 == excess_cap2):
            if excess_cap1 > excess_cap2:
                src, dest = ion2_trap, ion1_trap
            else:
                src, dest = ion1_trap, ion2_trap
        elif ion2_move_score == ion1_move_score and (excess_cap1 == excess_cap2):
            src, dest = ion2_trap, ion1_trap
        else:
            src = ion1_trap
            dest = ion2_trap


        """Logic for dead-lock resolution
        Check if dest trap cap is 0.
        if yes, check if scheduling another gate in the same  (or prev layers)
        can resolve the deadlock.
        """
        dest_excess_cap = (
            self.machine.traps[dest].capacity - 
            len(self.sys_state.trap_ions[dest])
        )
        if dest_excess_cap == 0 and mng_deadlock:
            print(f"managing deadlock")
            success = self.manage_deadlock(old_src=src, old_dest=dest, fire_time=fire_time)
            if success:
                print("Reordered | returning from count_future_ops deadlock")
                return None, None, False, fin_time_new
            
            biased_coin = True
            if not success and self.free_up and biased_coin:
                nbrs = self.neighbors[dest]
                
                nbr_excess_caps = {
                    nbr: self.excess_cap_calc(nbr)
                    for nbr in nbrs if self.excess_cap_calc(nbr) > 1
                }
                
                sorted_nbrs = sorted(nbr_excess_caps, key=nbr_excess_caps.get)
                if sorted_nbrs:
                    # print(f"sorted nbrs {sorted_nbrs} dict {nbr_excess_caps}")--
                    candidate_trap = sorted_nbrs[0]
                    candidate_ion = self.move_ion_for_rebalancing(
                        src_trap=dest,
                        dest_trap=candidate_trap
                    )
                    # print(f"move to trap {candidate_trap} move ion {candidate_ion}")--
                    fin_time_new = self.fire_shuttle(dest, candidate_trap, candidate_ion, fire_time)
                    return None, None, False, fin_time_new

        return src, dest, True, fin_time_new

    def excess_cap_calc(self, trap):
        excess_cap = self.machine.traps[trap].capacity - len(self.sys_state.trap_ions[trap])

        return excess_cap

    def manage_deadlock(self, old_src, old_dest, fire_time):
        """Deadlock will be resolved if any gate in the layer has
            filled `dest_trap` as `src_trap`.
        1. Check if a gate in the same layer requires a shuttle.
        2. if yes, find `src` and `dest` using the shuttle_direction function.
        3. if new_src == old_dest, then scheduling this gate earlier will resolve deadlock.
        4. if now gate in the `active_layer` resolves deadlock, try gate from immediate next
            layer?
            a. Maybe there are unscheduled gates from previous layers. Try those as well.
        """
        old_traps = [old_src, old_dest]
        active_layer = self.gate_to_layer[self.active_gate]
        

        gates_in_active_layer = []
        w = 0
        for k in range(active_layer+w):
            gates_in_active_layer += self.layers[k]

        for gate in gates_in_active_layer:
            if gate == self.active_gate:
                continue
            if not gate in self.remaining_gates:
                continue

            ion1 = self.gate_info[gate][0]
            ion2 = self.gate_info[gate][1]
            #Find time at which ions are ready
            _, ion1_trap = self.ion_ready_info(ion1)
            _, ion2_trap = self.ion_ready_info(ion2)

            if ion1_trap == ion2_trap:
                continue

            if not (ion1_trap == old_dest or ion2_trap == old_dest):
                continue

            excess_cap1 = (
                self.machine.traps[ion1_trap].capacity - 
                len(self.sys_state.trap_ions[ion1_trap])
            )
            excess_cap2 = (
                self.machine.traps[ion2_trap].capacity - 
                len(self.sys_state.trap_ions[ion2_trap])
            )  

            new_src, new_dest, _, _ = self.count_future_ops(ion1, ion2, ion1_trap, ion2_trap,
                                                 excess_cap1,
                                                 excess_cap2,
                                                 fire_time=fire_time,
                                                 mng_deadlock=False)

            new_dest_cap = self.machine.traps[new_dest].capacity - len(self.sys_state.trap_ions[new_dest])

            if new_dest_cap == 0:
                continue
            
            if new_src == old_dest:
                valid = self.check_valid_reorder(gate)
                if not valid:
                    continue
                self.reorder_remaining_gates(gate)
                return True

        return False
            
    def check_valid_reorder(self, gate):
        active_gate_idx = self.remaining_gates.index(self.active_gate)
        gate_index = self.remaining_gates.index(gate)
        gates_between_active_and_candidate = self.remaining_gates[active_gate_idx+1:gate_index]
        gate_ions = self.gate_info[gate]
        
        for interim_gate in gates_between_active_and_candidate:
            ion1, ion2 = self.gate_info[interim_gate]

            if ion1 in gate_ions or ion2 in gate_ions:
                return False

        return True

    
    def ops_in_trap(self, ion, trap):
        """Takes an ion and a trap.
        Returns how many future-ops that ion has in that trap
        """
        ops = 0
        for g in self.remaining_gates:
            ion1, ion2 = self.gate_info[g]

            if ion == ion1:
                _, ion2_trap = self.ion_ready_info(ion2)
                if trap == ion2_trap:
                    ops += 1
            elif ion == ion2:
                _, ion1_trap = self.ion_ready_info(ion1)
                if trap == ion1_trap:
                    ops += 1

        return ops


    
    def reorder_remaining_gates(self, gate):
        self.remaining_gates.remove(gate)
        insert_idx = self.remaining_gates.index(self.active_gate)
        self.remaining_gates.insert(insert_idx, gate)
        print(f"re-ordered remaining gates {self.remaining_gates}")

    #Heuristic to determine direction of shuttling for two traps
    def shuttling_direction(self, ion1, ion2, ion1_trap, ion2_trap, **kwargs):
        #Other Possible policies: lookahead, traffic/path based
        fire_time = kwargs['fire_time']
        m = self.machine
        ss = self.sys_state
        excess_cap1 = m.traps[ion1_trap].capacity - len(ss.trap_ions[ion1_trap])
        excess_cap2 = m.traps[ion2_trap].capacity - len(ss.trap_ions[ion2_trap])
        #both excess capacities can be 0 if the traps are full
        frac_empty = float(excess_cap1)/m.traps[ion1_trap].capacity
        fin_time = 0
        if self.lookahead:
            """New"""
            source_trap, dest_trap, success, fin_time = self.count_future_ops(
                ion1,
                ion2,
                ion1_trap,
                ion2_trap,
                excess_cap1,
                excess_cap2,
                fire_time=fire_time
            )
            if not success:
                return None, None, False, fin_time

            if source_trap == ion1_trap:
                src_cap = excess_cap1
                dest_cap = excess_cap2
            else:
                src_cap = excess_cap2
                dest_cap = excess_cap1
            
            if dest_cap == 0:
                # print(f"## dest trap {dest_trap}")--
                dest_trap, source_trap = source_trap, dest_trap
        else:
            """Old"""
            #Whichever trap has more excess capacity choose that as the destination
            if excess_cap1 > excess_cap2:
                dest_trap = ion1_trap
                source_trap = ion2_trap
            else:
                dest_trap = ion2_trap
                source_trap = ion1_trap

        if excess_cap1 <= 0 and excess_cap2 <= 0:
            print(ion1_trap, ion2_trap)
            print(ss.trap_ions)
            print("Both traps full", ion1_trap, m.traps[ion1_trap].capacity,  ss.trap_ions[ion1_trap])
            assert 0
        return source_trap, dest_trap, success, fin_time

    #Fire an end-to-end shuttle operation from src_trap to dest_trap
    def fire_shuttle(self, src_trap, dest_trap, ion, gate_fire_time, route=[]):
        s = self.schedule
        m = self.machine
        #If route is not specified in the function args, find a route using
        #the router object passed to the scheduler
        if len(route):
            rpath = route
        else:
            rpath = self.router.find_route(src_trap, dest_trap)

        #Find the time that it will take to do this entire shuttle
        #This is required to find a feasible time for scheduling this shuttle
        t_est = 0
        for i in range(len(rpath)-1):
            src = rpath[i]
            dest = rpath[i+1]
            if type(src) == Trap and type(dest) == Junction:
                my_seg = m.graph[src][dest]['seg']
                t_est += m.mparams.split_merge_time
                #split_time(self.sys_state, src.id, my_seg.id, ion)
            elif type(src) == Junction and type(dest) == Junction:
                t_est += m.move_time(src.id, dest.id)
            elif type(src) == Junction and type(dest) == Trap:
                t_est += m.merge_time(dest.id)

        #This is the traffic-unaware/conservative version where we wait for the full path to be available
        clk = self.schedule.identify_start_time(rpath, gate_fire_time, t_est)

        #Add the shuttling operations to the schedule based on the identified start time
        clk = self._add_shuttle_ops(rpath, ion, clk)

        #self.sys_state.trap_ions[src_trap].remove(ion)
        #self.sys_state.trap_ions[dest_trap].append(ion)
        return clk

    #Helper function to implement a shuttle
    def _add_shuttle_ops(self, spath, ion, clk):
        #Decompose into trap-trap paths
        #For each trap to trap path call a split-move*-merge sequence
        trap_pos = []
        for i in range(len(spath)):
            if type(spath[i]) == Trap:
                trap_pos.append(i)
        for i in range(len(trap_pos)-1):
            idx0 = trap_pos[i]
            idx1 = trap_pos[i+1]+1
            clk = self._add_partial_shuttle_ops(spath[idx0:idx1], ion, clk)
            self.sys_state.trap_ions[spath[trap_pos[i]].id].remove(ion)
            last_junct = spath[trap_pos[i+1]-1]
            dest_trap = spath[trap_pos[i+1]]
            last_seg = self.machine.graph[last_junct][dest_trap]['seg']
            orient = dest_trap.orientation[last_seg.id]
            if orient == 'R':
                self.sys_state.trap_ions[spath[trap_pos[i+1]].id].append(ion)
            else:
                self.sys_state.trap_ions[spath[trap_pos[i+1]].id].insert(0, ion)
        return clk

    #Helper function to implement a shuttle
    def _add_partial_shuttle_ops(self, spath, ion, clk):
        assert len([item for item in spath if type(item) == Trap]) == 2
        seg_list = []
        for i in range(len(spath)-1):
            u = spath[i]
            v = spath[i+1]
            seg_list.append(self.machine.graph[u][v]['seg'])
        clk = self.add_split_op(clk, spath[0], seg_list[0], ion)
        for i in range(len(seg_list)-1):
            u = seg_list[i]
            v = seg_list[i+1]
            junct = spath[1+i]
            clk = self.add_move_op(clk, u, v, junct, ion)
        clk = self.add_merge_op(clk, spath[-1], seg_list[-1], ion)
        return clk

    #Main scheduling function for a gate
    def schedule_gate(self, gate, specified_time=0):
        s = self.schedule
        #Find time at which the gate can be fired
        ready = self.gate_ready_time(gate)
        ion1 = self.gate_info[gate][0]
        ion2 = self.gate_info[gate][1]
        #Find time at which ions are ready
        ion1_time, ion1_trap = self.ion_ready_info(ion1)
        ion2_time, ion2_trap = self.ion_ready_info(ion2)
        fire_time = max(ready, ion1_time, ion2_time)
        fire_time = max(fire_time, specified_time)

        #print("Gate", gate, "I1", ion1, "I2", ion2, "IT1", ion1_trap, "IT2", ion2_trap, "Ready", ready, "FireTime:", fire_time)
        sc = True
        if ion1_trap == ion2_trap:
            #Ions are co-located in a trap, no shuttling required
            self.add_gate_op(fire_time, ion1_trap, gate, ion1, ion2)
        else:
            #Check if there is at least one path to shuttle from src to dest trap
            #rebalances the machine if needed i.e., clear traffic blocks
            rebal_flag, new_fin_time = self.rebalance_traps(focus_traps=[ion1_trap, ion2_trap], fire_time=fire_time)
            
            if rebal_flag:
                return False, new_fin_time
            
            if not rebal_flag:
                source_trap, dest_trap, sc, fin_time = self.shuttling_direction(ion1, ion2, ion1_trap, ion2_trap, fire_time=fire_time)

                if not sc:
                    return False, fin_time
                
                if source_trap == ion1_trap:
                    moving_ion = ion1
                else:
                    moving_ion = ion2
                n_ions = []
                n_traps = len(self.machine.traps)
                for trap in range(n_traps):
                    ni = len(self.sys_state.trap_ions[trap])
                    n_ions.append(ni)
                
                print(f"MOV {[moving_ion]} T{source_trap} ({len(self.sys_state.trap_ions[source_trap])}) --> T{dest_trap} ({len(self.sys_state.trap_ions[dest_trap])})")
                self.moves += abs(source_trap - dest_trap)
                clk = self.fire_shuttle(source_trap, dest_trap, moving_ion, fire_time)
                self.add_gate_op(clk, dest_trap, gate, ion1, ion2)
            # else:
            #     #This is for the rebalancing case, trap_ids compute till this point may be stale
            #     print("calling schedule gate from rebalancing")
            #     self.schedule_gate(gate, specified_time=new_fin_time)

        return True, 0

    #Checks and rebalances the machine if necessary using MCMF
    def rebalance_traps(self, focus_traps, fire_time):
        m = self.machine
        ss = self.sys_state
        t1 = focus_traps[0]
        t2 = focus_traps[1]
        excess_cap1 = m.traps[t1].capacity - len(ss.trap_ions[t1])
        excess_cap2 = m.traps[t2].capacity - len(ss.trap_ions[t2])
        need_rebalance = False

        ftr = FreeTrapRoute(m, ss)
        status12, route12 = ftr.find_route(t1, t2)
        status21, route21 = ftr.find_route(t2, t1)

        #If both traps are full
        if excess_cap1 == 0 and excess_cap2 == 0:
            need_rebalance = True
        else:
           #If no route exists either way
            if status12 == 1 and status21 == 1:
                need_rebalance = True
        if need_rebalance:
            #print("Rebalance procedure", "clk=", fire_time)
            finish_time = self.do_rebalance_traps(fire_time)
            #print("Rebalance procedure", "clk=", finish_time)
        if need_rebalance:
            return 1, finish_time
        else:
            return 0, fire_time

    def ops_by_trap_dict_prep(self, ion_list, trap):
        ops = {}
        for ion in ion_list:
            counts = self.ops_in_trap(ion, trap)
        
            ops[ion] = counts

        return ops
    
    def move_ion_for_rebalancing(self, src_trap, dest_trap):
        ions_src = self.sys_state.trap_ions[src_trap]

        """counts in dest"""
        ops_dest = {}
        for ion in ions_src:
            counts_in_dest = 0
            for gate in self.remaining_gates:
                gate_ions = copy.deepcopy(self.gate_info[gate])
                if ion in gate_ions:
                    gate_ions.remove(ion)
                    other_ion = gate_ions[0]
                    _, other_trap = self.ion_ready_info(other_ion)

                    if other_trap == dest_trap:
                        counts_in_dest += 1
        
            ops_dest[ion] = counts_in_dest
        

        """counts in src"""
        ops_src = {}
        for ion in ions_src:
            counts_in_src = 0
            for gate in self.remaining_gates:
                gate_ions = copy.deepcopy(self.gate_info[gate])
                if ion in gate_ions:
                    gate_ions.remove(ion)
                    other_ion = gate_ions[0]
                    _, other_trap = self.ion_ready_info(other_ion)

                    if other_trap == src_trap:
                        counts_in_src += 1
        
            ops_src[ion] = counts_in_src

        ops = {}
        # w_dest = random.choice([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
        w_dest = 0.49
        w_src = 1.0 - w_dest
        for ion in ions_src:
            try:
                if ops_dest[ion] == ops_src[ion]:
                    w_dest = 0.49
                else:
                    w_dest = 0.50
                w_src = 1 - w_dest
                
                score = w_dest * ops_dest[ion] - w_src * ops_src[ion]
            except ZeroDivisionError: # TODO: remove. not needed anymore.
                score = abs(ops_dest[ion]-ops_src[ion])

            ops[ion] = score

        
        sorted_ions = sorted(ops, key=ops.get)
        move_ion = sorted_ions[-1]
        sorted_ops = [round(ops[ion], 2) for ion in sorted_ions]

        return move_ion

    def ion_w_min_ops(self, ions, trap):
        ops_list = []
        for ion in ions:
            ops = self.ops_in_trap(ion, trap)
            ops_list.append((ion, ops))

        ops_list_sorted = sorted(ops_list, key=lambda x:x[1])

        return ops_list_sorted[0][0]



    def do_rebalance_traps(self, fire_time):
        self.count_rebalance += 1
        rebal = RebalanceTraps(self.machine, self.sys_state)
        flow_dict = rebal.clear_all_blocks()
        shuttle_graph = nx.DiGraph()
        used_flow = {}
        clk = fire_time
        for i in flow_dict:
            for j in flow_dict[i]:
                if flow_dict[i][j] != 0:
                    shuttle_graph.add_edge(i, j, weight=flow_dict[i][j])
                    used_flow[(i, j)] = 0
        fin_time = fire_time
        
        for node in shuttle_graph.nodes():
            if (shuttle_graph.in_degree(node) == 0) and type(node) == Trap:
                g_updated = False
                updated_graph = shuttle_graph.copy()
                for edge in used_flow:
                    if used_flow[edge] == updated_graph[edge[0]][edge[1]]['weight']:
                        updated_graph.remove_edge(edge[0], edge[1])
                T = nx.dfs_tree(updated_graph, source=node)
                for tnode in T:
                    if T.out_degree(tnode) == 0:
                        shuttle_route = nx.shortest_path(T, node, tnode)
                        break
                for i in range(len(shuttle_route)-1):
                    e0 = shuttle_route[i]
                    e1 = shuttle_route[i+1]
                    if (e0, e1) in used_flow:
                        used_flow[(e0, e1)] += 1
                    elif (e1, e0) in used_flow:
                        used_flow[(e1, e0)] += 1
                
                moving_ion = self.move_ion_for_rebalancing(node.id, tnode.id)
                
                ion_time, _ = self.ion_ready_info(moving_ion)
                fire_time = max(fire_time, ion_time)
                
                fin_time_new = self.fire_shuttle(node.id, tnode.id, moving_ion, fire_time, route=shuttle_route)
                fin_time = max(fin_time, fin_time_new)
        return fin_time

    def run(self):
        self.gates = list(nx.topological_sort(self.ir))
        
        last_gate = []
        if self.gate_order:
            # print(f"gate order {self.gate_order}")
            self.gates = []
            for item in self.gate_order:
                if type(item) == list:
                    self.gates += item
                    last_gate.append(item[-1])
                else:
                    self.gates = self.gate_order
                    break
            # self.gates = self.gate_order
        
        cnt = 0
        self.remaining_gates = copy.deepcopy(self.gates)
        print(f"gate order {self.gates}")
        # self.sys_state.print_state()
        # for g in self.gates:
        specified_time = 0
        while self.remaining_gates:
            g = self.remaining_gates[0]
            gions = self.gate_info[g]
            
            self.active_gate = g
            _, ion1_trap = self.ion_ready_info(gions[0])
            _, ion2_trap = self.ion_ready_info(gions[1])
            
            sx, new_fin_time = self.schedule_gate(g, specified_time)
            
            if new_fin_time != 0:
                specified_time = new_fin_time
            else:
                specified_time = 0

            if sx == False:
                sx = True
                continue
            
            self.remaining_gates.remove(g)
            active_layer = self.gate_to_layer[g]
            self.layers[active_layer].remove(g)
            
            cnt += 1
            
            if g in last_gate:
                self.moves = 0
        # self.schedule.print_events()
        # self.sys_state.print_state()

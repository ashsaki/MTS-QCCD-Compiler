import copy
import random

TRAP_CAP = 15
COMM_CAP = 2

class CostCalculator:
    def __init__(self, init_qubit_layout):
        self.init_qubit_layout = init_qubit_layout
        self.sys_state = copy.deepcopy(init_qubit_layout)
        self.qubit_location_dict = self.prep_qubit_location_dict()
        self.excess_caps = self.prep_excess_caps()

    def prep_excess_caps(self):
        excess_caps = {}

        for trap in self.init_qubit_layout:
            excess_caps[trap] = TRAP_CAP - len(self.init_qubit_layout[trap]) + COMM_CAP

        return excess_caps

    def prep_qubit_location_dict(self):
        location = {}

        for trap, ions in self.init_qubit_layout.items():
            for ion in ions:
                location[ion] = trap

        return location

    
class LayerSplitter:
    def __init__(self, layer, cx_gate_map, qubit_location_dict):
        self.layer = layer
        self.gate_info = cx_gate_map
        self.location = qubit_location_dict

    def split_layer(self):
        no_shuttle_gates = []
        shuttle_gates = []
        for gate in self.layer:
            ion1, ion2 = self.gate_info[gate]

            if self.location[ion1] == self.location[ion2]:
                no_shuttle_gates.append(gate)
            else:
                shuttle_gates.append(gate)

        return no_shuttle_gates, shuttle_gates

def calc_n_moves(ion1, ion2, qubit_location_dict):
    trap1 = qubit_location_dict[ion1]
    trap2 = qubit_location_dict[ion2]

    return abs(trap1 - trap2)

def h_cost(gates, cx_gate_map, qubit_location_dict):
    cost = 0
    flat_list = [item for sublist in gates for item in sublist]
    for g in flat_list:
        ion1, ion2 = cx_gate_map[g]

        dist = calc_n_moves(ion1, ion2, qubit_location_dict)

        cost += dist

    return cost

def updater(move_ion, src_trap, dest_trap, location, ecc, idx, gate_order, cx_gate_map, min_vals):
    gx = calc_n_moves(ion1, ion2, location)
    assert gx > 0
    location[move_ion] = dest_trap
    ecc[dest_trap] -= 1
    ecc[src_trap] += 1
    # print(f"ion1 {ion1} ion2 {ion2} trap1 {trap1} trap2 {trap2} gate {gate} gx {gx}")
    hx = h_cost(gate_order[idx+1:], cx_gate_map, location)
    total_cost = gx + hx

    if total_cost < min_cost:
        min_cost = total_cost
        min_cost_gate = gate
        min_cost_location = location
        min_cost_excs_cap = ecc
        min_cost_gx = gx
    elif total_cost == min_cost:
        update = random.choice([True, False])

        if update:
            min_cost = total_cost
            min_cost_gate = gate
            min_cost_location = location
            min_cost_excs_cap = ecc
            min_cost_gx = gx

    return min_cost, min_cost_gate, min_cost_location, min_cost_excs_cap, min_cost_gx, location, ecc


def gate_reorder(gate_order, cx_gate_map, init_qubit_layout):
    gate_reorder = []
    total_shuttles = 0
    cc = CostCalculator(init_qubit_layout)
    qubit_location_dict = cc.qubit_location_dict
    excess_caps = cc.excess_caps
    
    for idx, layer in enumerate(gate_order):
        new_layer = []
        # print(f"Layer {layer}")
        ls = LayerSplitter(layer, cx_gate_map, qubit_location_dict)
        no_shuttle_gates, shuttle_gates = ls.split_layer()
        # print(f"No Shuttle gates {no_shuttle_gates} Shuttle gates {shuttle_gates}")
        new_layer += no_shuttle_gates
        
        while shuttle_gates:
            min_cost = float("inf")
            min_cost_gate = None
            min_cost_location = None
            min_cost_excs_cap = None
            min_cost_gx = None
            w = 5
            
            for gate in shuttle_gates:
                # print(f"Checking gate {gate}")
                ion1, ion2 = cx_gate_map[gate]
                trap1 = qubit_location_dict[ion1]
                trap2 = qubit_location_dict[ion2]

                # location = copy.deepcopy(qubit_location_dict)
                # ecc = copy.deepcopy(excess_caps)
                # # print(f"Extra Cap {excess_caps} trap1 {trap1} trap2 {trap2}")
                
                # if ecc[trap1] == 1 and ecc[trap2] == 1:
                #     assert False, f"Both traps have extra cap of 1"
                # elif ecc[trap1] == 1:
                #     # moved ion1 trap1 --> trap2
                # elif ecc[trap2] == 1:

                # # ion1 moved trap1 --> trap2
                # if ecc[trap2] > 1:
                #     min_cost, min_cost_gate, min_cost_location, min_cost_excs_cap, min_cost_gx, location, ecc = updater(
                #         ion1, trap1, trap2, location, ecc, idx, gate_order, cx_gate_map, min_cost
                #     )
                    
                
                location = copy.deepcopy(qubit_location_dict)
                ecc = copy.deepcopy(excess_caps)
                if ecc[trap1] >= ecc[trap2]:
                    gx = calc_n_moves(ion1, ion2, location)
                    assert gx > 0
                    location[ion2] = trap1
                    ecc[trap1] -= 1
                    ecc[trap2] += 1
                    # print(f"ion1 {ion1} ion2 {ion2} trap1 {trap1} trap2 {trap2} gate {gate} gx {gx}")
                    hx = h_cost(gate_order[idx+1:], cx_gate_map, location)
                    total_cost = gx + hx

                    if total_cost < min_cost:
                        min_cost = total_cost
                        min_cost_gate = gate
                        min_cost_location = location
                        min_cost_excs_cap = ecc
                        min_cost_gx = gx
                    elif total_cost == min_cost:
                        update = random.choice([True, False])

                        if update:
                            min_cost = total_cost
                            min_cost_gate = gate
                            min_cost_location = location
                            min_cost_excs_cap = ecc
                            min_cost_gx = gx
                else:
                    gx = calc_n_moves(ion1, ion2, location)
                    assert gx > 0
                    location[ion1] = trap2
                    ecc[trap2] -= 1
                    ecc[trap1] += 1
                    # print(f"ion1 {ion1} ion2 {ion2} trap1 {trap1} trap2 {trap2} gate {gate} gx {gx}")
                    hx = h_cost(gate_order[idx+1:], cx_gate_map, location)
                    total_cost = gx + hx

                    if total_cost < min_cost:
                        min_cost = total_cost
                        min_cost_gate = gate
                        min_cost_location = location
                        min_cost_excs_cap = ecc
                        min_cost_gx = gx
                    elif total_cost == min_cost:
                        update = random.choice([True, False])

                        if update:
                            min_cost = total_cost
                            min_cost_gate = gate
                            min_cost_location = location
                            min_cost_excs_cap = ecc
                            min_cost_gx = gx


            
            # print(f"Gate with min cost {min_cost_gate} gx {min_cost_gx} cap {min_cost_excs_cap}")
            new_layer.append(min_cost_gate)
            shuttle_gates.remove(min_cost_gate)
            qubit_location_dict = min_cost_location
            excess_caps = min_cost_excs_cap
            total_shuttles += min_cost_gx

        
        gate_reorder.append(new_layer)

    return gate_reorder, total_shuttles

                    

                    


import sys
import random
from parse import InputParse
from mappers import *
from machine import Machine, MachineParams, Trap, Segment
from ejf_schedule import Schedule, EJFSchedule
from analyzer import *
from test_machines import *
import numpy as np
from qiskit import QuantumCircuit
from qiskit.converters import circuit_to_dag, dag_to_circuit
from qiskit.visualization import dag_drawer
from cost_calculator import gate_reorder

np.random.seed(12345)

#Command line args
#Machine attributes
openqasm_file_name = sys.argv[1]
machine_type = sys.argv[2]
num_ions_per_region = int(sys.argv[3])
mapper_choice = sys.argv[4]
reorder_choice = sys.argv[5]
serial_trap_ops = int(sys.argv[6])
serial_comm = int(sys.argv[7])
serial_all = int(sys.argv[8])
gate_type = sys.argv[9]
swap_type = sys.argv[10]
lookahead = int(sys.argv[11])

##########################################################
mpar_model1 = MachineParams()
mpar_model1.alpha = 0.003680029
mpar_model1.beta = 39.996319971
mpar_model1.split_merge_time = 80
mpar_model1.shuttle_time = 5
mpar_model1.junction2_cross_time = 5
mpar_model1.junction3_cross_time = 100
mpar_model1.junction4_cross_time = 120
mpar_model1.gate_type = gate_type
mpar_model1.swap_type = swap_type
mpar_model1.ion_swap_time = 42
machine_model = "MPar1"

'''
mpar_model2 = MachineParams()
mpar_model2.alpha = 0.003680029
mpar_model2.beta = 39.996319971
mpar_model2.split_merge_time = 80
mpar_model2.shuttle_time = 5
mpar_model2.junction2_cross_time = 5
mpar_model2.junction3_cross_time = 100
mpar_model2.junction4_cross_time = 120
mpar_model2.alpha
machine_model = "MPar2"
'''

print("Simulation")
print("Program:", openqasm_file_name)
print("Machine:", machine_type)
print("Model:", machine_model)
print("Ions:", num_ions_per_region)
print("Mapper:", mapper_choice)
print("Reorder:", reorder_choice)
print("SerialTrap:", serial_trap_ops)
print("SerialComm:", serial_comm)
print("SerialAll:", serial_all)
print("Gatetype:", gate_type)
print("Swaptype:", swap_type)

#Create a test machine
if machine_type == "G2x3":
    m = test_trap_2x3(num_ions_per_region, mpar_model1)
elif machine_type == "L6":
    m = make_linear_machine(6, num_ions_per_region, mpar_model1)
elif machine_type == "H6":
    m = make_single_hexagon_machine(num_ions_per_region, mpar_model1)
else:
    assert 0

#Parse the input program DAG
ip = InputParse()
ip.parse_ir(openqasm_file_name)
ip.visualize_graph("visualize_graph_2.gexf") # dumps parser graph into file


#Map the program onto the machine regions
#For every program qubit, this gives a region id
if mapper_choice == "LPFS":
    qm = QubitMapLPFS(ip,m)
elif mapper_choice == "Agg":
    qm = QubitMapAgg(ip, m)
elif mapper_choice == "Random":
    qm = QubitMapRandom(ip, m)
elif mapper_choice == "PO":
    qm = QubitMapPO(ip, m)
elif mapper_choice == "Greedy":
    qm = QubitMapGreedy(ip, m)
else:
    assert 0
mapping = qm.compute_mapping()

#Reorder qubits within a region to increse the use of high fidelity operations
if mapper_choice == "Greedy":
    init_qubit_layout = mapping
else:
    qo = QubitOrdering(ip, m, mapping)
    if reorder_choice == "Naive":
        init_qubit_layout = qo.reorder_naive()
    elif reorder_choice == "Fidelity":
        init_qubit_layout = qo.reorder_fidelity()
    else:
        assert 0

print(init_qubit_layout)



# DAG-based gate order
cx_gate_map_copy = copy.deepcopy(ip.cx_gate_map)
def prep_gate_order():
    order = []
    with open("temp.tmp", "r") as f_rd:
        for line in f_rd:
            if line.startswith("cx"):
                base = ''.join(line.split()).split(',')
                qbit1 = int(base[0].split('[')[1].split(']')[0])
                qbit2 = int(base[1].split('[')[1].split(']')[0])

                global_gates_ids = sorted(cx_gate_map_copy.keys())
                for global_gate_id in global_gates_ids:
                    if cx_gate_map_copy[global_gate_id] == [qbit1, qbit2]:
                        order.append(global_gate_id)
                        del cx_gate_map_copy[global_gate_id]
                        break

    random.shuffle(order)
    return order

dag_based_ordering = True
gate_order = []

if dag_based_ordering:
    qc = QuantumCircuit.from_qasm_file(openqasm_file_name)
    dag = circuit_to_dag(qc)
    layers = dag.layers()

    for idx, layer in enumerate(layers):
        graph = layer['graph']
        circ = dag_to_circuit(graph)
        circ.qasm(filename='temp.tmp')
        circ.qasm(filename=f'./tmp/temp{idx}.tmp')
        order = prep_gate_order()
        # print(f"Order {len(order)} {order}")
        # gate_order += order
        gate_order.append(order)

# The above function is only computing which gate belongs to which layer
# the DAG-based gate ordering is not used in the main code. Instead,
# original topological sort-based ordering was used.
# that is why after getting the list of gates in the same layer
# `gate_order` is forcefully set to an empty list.
layers = copy.deepcopy(gate_order)
gate_order = []

#Schedule gates in the prorgam in topological sorted order
#EJF = earliest job first, here it refers to earliest gate first
#This step performs the shuttling
ejfs = EJFSchedule(
    ip.gate_graph,
    ip.cx_gate_map,
    m,
    init_qubit_layout,
    serial_trap_ops,
    serial_comm,
    serial_all,
    gate_order=gate_order,
    lookahead=lookahead,
    comm_cap=2,
    layers=layers
)
ejfs.run()

#Analyze the output schedule and print statistics
analyzer = Analyzer(ejfs.schedule, m, init_qubit_layout)
analyzer.move_check()
print("SplitSWAP:", ejfs.split_swap_counter)
#analyzer.print_events()
print("----------------")

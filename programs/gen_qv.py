import qiskit
import qiskit.ignis.verification.quantum_volume as qv
from qiskit.compiler import transpile

nqubits = 64
qubit_lists=[list(range(nqubits))]
ntrials = 1

_, qv_circs_nomeas = qv.qv_circuits(qubit_lists, ntrials)

qv_circs_nomeas[0] = qiskit.compiler.transpile(qv_circs_nomeas[0], basis_gates=['id','rz','sx', 'x', 'cx'])
qv_circs_nomeas[0][0].qasm(filename='qv64.qasm')

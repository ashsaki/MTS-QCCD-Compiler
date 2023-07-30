from qiskit import Aer, QuantumCircuit
from qiskit.circuit.library import RealAmplitudes, ZZFeatureMap
from qiskit_machine_learning.neural_networks import TwoLayerQNN, CircuitQNN
from qiskit_machine_learning.algorithms.classifiers import NeuralNetworkClassifier, VQC
from qiskit_machine_learning.algorithms.regressors import NeuralNetworkRegressor, VQR
from qiskit.compiler import transpile
#
from qiskit import QuantumRegister, QuantumCircuit, BasicAer
from qiskit.circuit.library import TwoLocal, UniformDistribution


# qnn
num_inputs = 64

feature_map = ZZFeatureMap(num_inputs)

# construct ansatz
ansatz = RealAmplitudes(num_inputs, reps=1)

# construct quantum circuit
qc = QuantumCircuit(num_inputs)
qc.append(feature_map, range(num_inputs))
qc.append(ansatz, range(num_inputs))

qct = transpile(qc, basis_gates=['id', 'rz', 'sx', 'x', 'cx'])
print(qct.count_ops()['cx'])

# qgan
ansatz = TwoLocal(64, 'x', 'cx', reps=1)
ansatz = transpile(ansatz, basis_gates=['id', 'rz', 'sx', 'x', 'cx'])
ansatz.qasm(filename='twolocal.qasm')

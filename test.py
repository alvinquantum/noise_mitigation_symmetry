# Needed for functions
import numpy as np
import time
from copy import deepcopy
import mymodule

# Import Qiskit classes
import qiskit
import qiskit.quantum_info as qi
from qiskit import QuantumRegister, QuantumCircuit, ClassicalRegister, Aer
from qiskit.providers.aer import noise
from qiskit.compiler import assemble

# Tomography functions
from qiskit.ignis.verification.tomography import state_tomography_circuits, StateTomographyFitter
from qiskit.ignis.verification.tomography import process_tomography_circuits, ProcessTomographyFitter
from qiskit.ignis.verification.tomography import gateset_tomography_circuits, GatesetTomographyFitter
import qiskit.ignis.mitigation.measurement as mc

# Auxiliary methods
from qiskit.quantum_info import Choi, Kraus
from qiskit.extensions import HGate, XGate

# Create the actual circuit
q2 = QuantumRegister(3)
ghz = QuantumCircuit(q2)
ghz.h(q2[0])
ghz.cx(q2[0], q2[1])
ghz.cx(q2[1], q2[2])
ghz.h(q2[2])
print(ghz)

qst_ghz = state_tomography_circuits(ghz, [q2[0],q2[1]])
print(qst_ghz[0])

#Make a copy without the ancilla register
qst_ghz_no_anc = deepcopy(qst_ghz)
ca = ClassicalRegister(1)
for qst_ghz_circ in qst_ghz:
    qst_ghz_circ.add_register(ca)
    qst_ghz_circ.measure(q2[2],ca[0])

#Run in Aer
job = qiskit.execute(qst_ghz, Aer.get_backend('qasm_simulator'), shots=10000)
raw_results = job.result()

new_result = mymodule.post_select_on_ancilla(raw_results, "1", 2)

tomo_bell = StateTomographyFitter(new_result, qst_ghz_no_anc)
# Perform the tomography fit
# which outputs a density matrix
rho_fit_bell = tomo_bell.fit(method='lstsq')


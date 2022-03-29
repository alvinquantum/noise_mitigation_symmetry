from circtester import CircuitRunner
import circgenerator
from qiskit import QuantumCircuit
import checksfinder_new
import utilities

def test_circrunner_getsso():
    dist1={"1111": .9, "0000": 0.1}
    dist2={"1111": .9, "0000": 0.1}
    assert CircuitRunner.get_sso(dist1, dist2) == 1.0

def test_count_gate():
    gate_count=5
    circ = QuantumCircuit(gate_count)
    for idx in range(gate_count):
        circ.x(idx)
        circ.y(idx)
        circ.z(idx)
        circ.h(idx)
        circ.cx(0, 1)
        circ.swap(0, 1)
    assert circgenerator.count_gate(circ, "x") == gate_count
    assert circgenerator.count_gate(circ, "y") == gate_count
    assert circgenerator.count_gate(circ, "z") == gate_count
    assert circgenerator.count_gate(circ, "h") == gate_count
    assert circgenerator.count_gate(circ, "cx") == gate_count
    assert circgenerator.count_gate(circ, "swap") == gate_count
    assert circgenerator.count_gate(circ, "not_a_gate") == 0

def test_get_single_weight_pauli_table_at_qubit_idx():
    qubit_idx=0
    number_of_qubits=5
    single_weights=checksfinder_new.get_single_weight_pauli_table_at_qubit_idx(qubit_idx, number_of_qubits)
    assert single_weights[0]=="XIIII"
    assert single_weights[1]=="YIIII"
    assert single_weights[2]=="ZIIII"
    qubit_idx=4
    number_of_qubits=5
    single_weights=checksfinder_new.get_single_weight_pauli_table_at_qubit_idx(qubit_idx, number_of_qubits)
    assert single_weights[0]=="IIIIX"
    assert single_weights[1]=="IIIIY"
    assert single_weights[2]=="IIIIZ"

def test_generate_checks_file_name_from_raw_qasm():
    file_name="some_name_raw_.qasm"
    assert checksfinder_new.generate_checks_file_name_from_raw_qasm(file_name) == "some_name_.qasm"

def test_generate_obj_file_name_from_checks_file_name():
    file_name="some_name_.qasm"
    assert checksfinder_new.generate_obj_file_name_from_checks_file_name(file_name) == "some_name_.obj"
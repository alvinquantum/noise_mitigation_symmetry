from circtester import CircuitRunner
import circgenerator
from qiskit import QuantumCircuit
import utilities
import checksfinder_new as checksfinder

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
    single_weights=checksfinder.get_single_weight_pauli_table_at_qubit_idx(qubit_idx, number_of_qubits)
    assert single_weights[0]=="XIIII"
    assert single_weights[1]=="ZIIII"
    qubit_idx=4
    number_of_qubits=5
    single_weights=checksfinder.get_single_weight_pauli_table_at_qubit_idx(qubit_idx, number_of_qubits)
    assert single_weights[0]=="IIIIX"
    assert single_weights[1]=="IIIIZ"

def test_get_min_checks_linear():
    circ=QuantumCircuit(2)
    # Test pushing through X.
    circ.x(0)
    p2s=["IX", "IY", "IZ"]
    finder= checksfinder.ChecksFinder(2, circ)
    _, _, pauli_str_p1s, pauli_str_p2s=finder.get_min_checks_linear(p2s)
    assert pauli_str_p1s==["+1IX", "-1IY", "-1IZ"] and pauli_str_p2s==["+1IX", "+1IY", "+1IZ"]
    # Test pushing through Y.
    circ=QuantumCircuit(2)
    circ.y(0)
    finder= checksfinder.ChecksFinder(2, circ)
    _, _, pauli_str_p1s, pauli_str_p2s=finder.get_min_checks_linear(p2s)
    assert pauli_str_p1s==["-1IX", "+1IY", "-1IZ"] and pauli_str_p2s==["+1IX", "+1IY", "+1IZ"]
    # Test pushing through Z.
    circ=QuantumCircuit(2)
    circ.z(0)
    finder= checksfinder.ChecksFinder(2, circ)
    _, _, pauli_str_p1s, pauli_str_p2s=finder.get_min_checks_linear(p2s)
    assert pauli_str_p1s==["-1IX", "-1IY", "+1IZ"] and pauli_str_p2s==["+1IX", "+1IY", "+1IZ"]
    # Test pushing through H.
    circ=QuantumCircuit(2)
    circ.h(0)
    finder= checksfinder.ChecksFinder(2, circ)
    _, _, pauli_str_p1s, pauli_str_p2s=finder.get_min_checks_linear(p2s)
    assert pauli_str_p1s==["+1IZ", "-1IY", "+1IX"] and pauli_str_p2s==["+1IX", "+1IY", "+1IZ"]
    # Test pushing through S.
    circ=QuantumCircuit(2)
    circ.s(0)
    finder= checksfinder.ChecksFinder(2, circ)
    _, _, pauli_str_p1s, pauli_str_p2s=finder.get_min_checks_linear(p2s)
    assert pauli_str_p1s==["-1IY", "+1IX", "+1IZ"] and pauli_str_p2s==["+1IX", "+1IY", "+1IZ"]
    # Test pushing through Sdg.
    circ=QuantumCircuit(2)
    circ.sdg(0)
    finder= checksfinder.ChecksFinder(2, circ)
    _, _, pauli_str_p1s, pauli_str_p2s=finder.get_min_checks_linear(p2s)
    assert pauli_str_p1s==["+1IY", "-1IX", "+1IZ"] and pauli_str_p2s==["+1IX", "+1IY", "+1IZ"]
    # Test pushing through CX. #Note that we are using little endian because Qiskit uses it.
    circ=QuantumCircuit(2)
    circ.cx(0,1)
    finder= checksfinder.ChecksFinder(2, circ)
    p2s=["IX", "IY", "IZ", "XI"] 
    _, _, pauli_str_p1s, pauli_str_p2s=finder.get_min_checks_linear(p2s)
    assert pauli_str_p1s==["+1XX", "+1XY", "+1IZ", "+1XI"] and pauli_str_p2s==["+1IX", "+1IY", "+1IZ", "+1XI"]
    p2s=["XX", "XY", "XZ"]
    _, _, pauli_str_p1s, pauli_str_p2s=finder.get_min_checks_linear(p2s)
    assert pauli_str_p1s==["+1IX", "+1IY", "+1XZ"] and pauli_str_p2s==["+1XX", "+1XY", "+1XZ"]
    p2s=["YI", "YX", "YY", "YZ"]
    _, _, pauli_str_p1s, pauli_str_p2s=finder.get_min_checks_linear(p2s)
    assert pauli_str_p1s==["+1YZ", "+1ZY", "-1ZX", "+1YI"] and pauli_str_p2s==["+1YI", "+1YX", "+1YY", "+1YZ"]
    p2s=["ZI", "ZX", "ZY", "ZZ"]
    _, _, pauli_str_p1s, pauli_str_p2s=finder.get_min_checks_linear(p2s)
    assert pauli_str_p1s==["+1ZZ", "-1YY", "+1YX", "+1ZI"] and pauli_str_p2s==["+1ZI", "+1ZX", "+1ZY", "+1ZZ"]

# def test_generate_checks_file_name_from_raw_qasm():
#     file_name="some_name_raw_.qasm"
#     assert checksfinder_new.generate_checks_file_name_from_raw_qasm(file_name) == "some_name_.qasm"

# def test_generate_obj_file_name_from_checks_file_name():
#     file_name="some_name_.qasm"
#     assert checksfinder_new.generate_obj_file_name_from_checks_file_name(file_name) == "some_name_.obj"

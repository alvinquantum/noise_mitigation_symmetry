from circtester import CircuitRunner
import circgenerator
from qiskit import QuantumCircuit

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
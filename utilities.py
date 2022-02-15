from qiskit import QuantumCircuit

def add_controlU(circ, pauli_str, number_of_qubits, quantum_register, ancilla_register):
    '''Testing cricuits: Adds a controlled Pauli to circ.'''
    print(pauli_str)
    c_minus_x = QuantumCircuit.from_qasm_str("""
    OPENQASM 2.0;
    include "qelib1.inc";
    qreg q[2];
    u3(1.0592554e-13,3.288277,2.8116442) q[0];
    u3(pi,5.8991791,2.7575864) q[1];
    cx q[0],q[1];
    u3(2*pi,5.4995268,4.1085153) q[0];
    u3(pi,0.38400625,3.5255989) q[1];
    """.strip())

    c_minus_j_x = QuantumCircuit.from_qasm_str("""
    OPENQASM 2.0;
    include "qelib1.inc";
    qreg q[2];
    u3(6.9759301e-13,4.3371506,0.076820611) q[0];
    u3(6.9089241e-09,5.477472,3.947306) q[1];
    cx q[0],q[1];
    u3(2*pi,0.60772643,2.832284) q[0];
    u3(5.5644668e-09,5.4477338,3.9770441) q[1];
    """.strip())

    c_plus_j_x =QuantumCircuit.from_qasm_str("""
    OPENQASM 2.0;
    include "qelib1.inc";
    qreg q[2];
    u3(pi,4.6477413,6.258122) q[0];
    u3(1.7599424,3*pi/2,3*pi/2) q[1];
    cx q[0],q[1];
    u3(pi,3.1666559,3.2062404) q[0];
    u3(4.9015351,3*pi/2,3*pi/2) q[1];
    """.strip())

    #Check the length of the pauli string. It maybe empty when there were no matches found.
    ancilla_qpos=ancilla_register[0]
    if len(pauli_str)>number_of_qubits:
        phase=pauli_str[:2]
        pauli_str=pauli_str[2:]
        qubit_pos=number_of_qubits-1
        phase_added=False
        for _, char in enumerate(pauli_str):
            # First case is most common
            if phase_added or phase=="+1":
                if char=="X":
                    # print(phase_added, phase, char)
                    circ.cx(ancilla_qpos, quantum_register[qubit_pos])
                elif char=="Y":
                    # print(phase_added, phase, char)
                    circ.sdg(quantum_register[qubit_pos])
                    circ.cx(ancilla_qpos, quantum_register[qubit_pos])
                    circ.s(quantum_register[qubit_pos])
                elif char=="Z":
                    # print(phase_added, phase, char)
                    circ.h(quantum_register[qubit_pos])
                    circ.cx(ancilla_qpos, quantum_register[qubit_pos])
                    circ.h(quantum_register[qubit_pos])
            # -1 phase
            else:
                if char=="X":
                    # print(phase_added, phase, char)
                    circ.compose(c_minus_x, qubits=[ancilla_qpos, quantum_register[qubit_pos]], inplace=True)
                elif char=="Y":
                    # print(phase_added, phase, char)
                    circ.sdg(quantum_register[qubit_pos])
                    circ.compose(c_minus_x, qubits=[ancilla_qpos, quantum_register[qubit_pos]], inplace=True)
                    circ.s(quantum_register[qubit_pos])
                elif char=="Z":
                    # print(phase_added, phase, char)
                    circ.h(quantum_register[qubit_pos])
                    circ.compose(c_minus_x, qubits=[ancilla_qpos, quantum_register[qubit_pos]], inplace=True)
                    circ.h(quantum_register[qubit_pos])
                # Note the phase added needs this check because of identity terms in the pauli strings.
                # With no check, if there is an identity it will change phase_added to true and we won't get the
                # added necessary phase on of the pauli terms.
                if char!="I":
                    phase_added=True
            qubit_pos-=1
    return circ
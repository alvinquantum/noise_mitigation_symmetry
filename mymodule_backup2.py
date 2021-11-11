# This code has been modified from Qiskit

# This code is part of Qiskit.
#
# (C) Copyright IBM 2017.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.


from copy import deepcopy
import numpy as np

from qiskit.circuit import QuantumRegister, ClassicalRegister, QuantumCircuit
from qiskit.circuit import Reset
from qiskit.circuit.library.standard_gates import (IGate, U1Gate, U2Gate, U3Gate, XGate,
                                                   YGate, ZGate, HGate, SGate, SdgGate, TGate,
                                                   TdgGate, RXGate, RYGate, RZGate, CXGate,
                                                   CYGate, CZGate, CHGate, CRZGate, CU1Gate,
                                                   CU3Gate, SwapGate, RZZGate,
                                                   CCXGate, CSwapGate)
from qiskit.circuit.exceptions import CircuitError

def split_circ():
    pass

def my_random_circuit(num_qubits, depth, seed=None):
    """ Generates a random circuit with num_qubits and depth (each wire has the specified depth). 
    The circuit uses h, s, t, cx gates and no measurements.
    Derived from https://qiskit.org/documentation/_modules/qiskit/circuit/random/utils.html#random_circuit.

    Returns:
        QuantumCircuit: constructed circuit
    """

    # This code is part of Qiskit.
    #
    # (C) Copyright IBM 2017.
    #
    # This code is licensed under the Apache License, Version 2.0. You may
    # obtain a copy of this license in the LICENSE.txt file in the root directory
    # of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
    #
    # Any modifications or derivative works of this code must retain this
    # copyright notice, and modified files need to carry a notice indicating
    # that they have been altered from the originals.

    # if max_operands < 1 or max_operands > 3:
    #     raise CircuitError("max_operands must be between 1 and 3")

    max_operands=2 
    # measure=False
    # conditional=False 
    # reset=False 

    one_q_ops = [HGate, SGate, TGate]
    # one_param = [U1Gate, RXGate, RYGate, RZGate, RZZGate, CU1Gate, CRZGate]
    # two_param = [U2Gate]
    # three_param = [U3Gate, CU3Gate]
    two_q_ops = [CXGate]
    # three_q_ops = [CCXGate, CSwapGate]

    qr = QuantumRegister(num_qubits, 'q')
    qc = QuantumCircuit(num_qubits)

    # if measure or conditional:
    #     cr = ClassicalRegister(num_qubits, 'c')
    #     qc.add_register(cr)

    # if reset:
    #     one_q_ops += [Reset]

    if seed is None:
        seed = np.random.randint(0, np.iinfo(np.int32).max)
    rng = np.random.default_rng(seed)

    # apply arbitrary random operations at every depth
    for _ in range(depth):
        # choose either 1, 2, or 3 qubits for the operation
        remaining_qubits = list(range(num_qubits))
        while remaining_qubits:
            max_possible_operands = min(len(remaining_qubits), max_operands)
            num_operands = rng.choice(range(max_possible_operands)) + 1
            rng.shuffle(remaining_qubits) # We shuffle and then apply the operations to the necessary number of qubits
            #  at the beginning.
            operands = remaining_qubits[:num_operands]
            remaining_qubits = [q for q in remaining_qubits if q not in operands]
            if num_operands == 1:
                operation = rng.choice(one_q_ops)
            else:
                operation = rng.choice(two_q_ops)
            # elif num_operands == 3:
            #     operation = rng.choice(three_q_ops)
            # if operation in one_param:
            #     num_angles = 1
            # elif operation in two_param:
            #     num_angles = 2
            # elif operation in three_param:
            #     num_angles = 3
            # else:
            # num_angles = 0
            # angles = [rng.uniform(0, 2 * np.pi) for x in range(num_angles)]
            # register_operands = [qr[i] for i in operands]
            # op = operation(*angles)

            # # with some low probability, condition on classical bit values
            # if conditional and rng.choice(range(10)) == 0:
            #     value = rng.integers(0, np.power(2, num_qubits))
            #     op.condition = (cr, value)

            register_operands = [qr[i] for i in operands]
            op = operation()

            qc.append(op, register_operands)

    # if measure:
    #     qc.measure(qr, cr)

    return qc

def post_select_on_ancilla(res, ancilla_value, new_nqubits):
    """
    strip the results where ancilla was not equal to `ancilla_value`
    This is some voodoo copied from 
    https://qiskit.org/documentation/tutorials/noise/8_tomography.html#2-Qubit-Conditional-State-Tomography  
    """
    assert(isinstance(ancilla_value, str))
    res_new = deepcopy(res)
    for resultidx, _ in enumerate(res.results):
        old_counts = res.get_counts(resultidx)
        new_counts = {}
        # print(res_new.results[resultidx].header.clbit_labels)
        res_new.results[resultidx].header.creg_sizes = [res_new.results[resultidx].header.creg_sizes[0]]
        res_new.results[resultidx].header.clbit_labels = res_new.results[resultidx].header.clbit_labels[0:-1]
        # print(res_new.results[resultidx].header.clbit_labels)
        res_new.results[resultidx].header.memory_slots = new_nqubits 
        for reg_key in old_counts:
            # print(reg_key)
            reg_bits = reg_key.split(' ')
            assert(len(reg_bits) == 2)
            assert(len(reg_bits[0]) == 1)
            if reg_bits[0]==ancilla_value:
                new_counts[reg_bits[1]]=old_counts[reg_key]
            res_new.results[resultidx].data.counts = new_counts
    return res_new

if __name__ == "__main__":
    pass
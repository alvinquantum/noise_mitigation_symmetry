#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
# The number of threads is capped on Bebop.
# os.environ['OPENBLAS_NUM_THREADS'] = '1'

from qiskit.quantum_info.operators.symplectic.pauli_utils import pauli_basis
import circgenerator
import sys

if __name__ == "__main__":
    NUMBER_OF_QUBITS=int(sys.argv[1])
    CNOT_COUNT=int(sys.argv[2])
    NUMBER_OF_CIRCUITS=int(sys.argv[3])
    # circgenerator.generate_circuits_with_checks_to_files(NUMBER_OF_QUBITS, CNOT_COUNT, NUMBER_OF_CIRCUITS)
    circgenerator.write_circs_no_checks_to_qasm_file(circgenerator.generate_random_circuits(NUMBER_OF_QUBITS, CNOT_COUNT, NUMBER_OF_CIRCUITS), NUMBER_OF_QUBITS, CNOT_COUNT)
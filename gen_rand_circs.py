#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
from tkinter import N
# The number of threads is capped on Bebop.
# os.environ['OPENBLAS_NUM_THREADS'] = '1'

from qiskit.quantum_info.operators.symplectic.pauli_utils import pauli_basis
from qiskit import QuantumCircuit
import circgenerator
import checksfinder_new as checksfinder
import circtester
import sys

if __name__ == "__main__":
    NUMBER_OF_QUBITS=int(sys.argv[1])
    CNOT_COUNT=int(sys.argv[2])
    NUMBER_OF_CIRCUITS=int(sys.argv[3])
    CODE_DIR=os.path.abspath(os.path.dirname(__file__))
    # subdir="qubits_10"
    # output_folder="raw_circuits"
    # input_folder="qubits_10_results"
    # BASE_INPUT_PATH=os.path.join(CODE_DIR,subdir, input_folder)
    # BASE_OUTPUT_PATH=os.path.join(CODE_DIR, subdir, output_folder)
    # # circgenerator.generate_circuits_with_checks_to_files(NUMBER_OF_QUBITS, CNOT_COUNT, NUMBER_OF_CIRCUITS)
    # circ=circgenerator.generate_a_random_circuit(NUMBER_OF_QUBITS, CNOT_COUNT, NUMBER_OF_CIRCUITS)
    # circ_with_checks=checksfinder.find_and_add_checks_to_circ(circ, CNOT_COUNT)
    # circs=circgenerator.generate_random_circuits(NUMBER_OF_QUBITS, CNOT_COUNT, NUMBER_OF_CIRCUITS)
    # circgenerator.write_circs_no_checks_to_qasm_file(circs, NUMBER_OF_QUBITS, CNOT_COUNT)
    # circs_with_checks=checksfinder.find_and_append_checks_to_all_circuits(circs, CNOT_COUNT)
    # # checksfinder.save_circuits_qasm(circs_with_checks, NUMBER_OF_QUBITS, CNOT_COUNT)
    # # circtester.simulate_all_circuits_concise(NUMBER_OF_QUBITS, CNOT_COUNT, 0, NUMBER_OF_CIRCUITS-1)
    # file_manipulator=circtester.FilesManipulator(BASE_INPUT_PATH, NUMBER_OF_QUBITS, CNOT_COUNT)
    # qasm_files=file_manipulator.get_all_qasm_files_with_checks(BASE_INPUT_PATH)
    # # print(len(qasm_files))
    # # print(qasm_files)
    # circtester.FilesManipulator.get_main_circuit_from_qasm_to_file(BASE_INPUT_PATH, BASE_OUTPUT_PATH, qasm_files, NUMBER_OF_QUBITS)
    circgenerator.generate_random_circuits_to_files(NUMBER_OF_QUBITS, CNOT_COUNT, NUMBER_OF_CIRCUITS)

    # # subdir="data"
    # subdir="qubits_10"
    # raw_subdir="raw_circuits"
    # initial_states_subdir="qubits_10_results"
    # split_results_subdir="split_circs_result"
    # BASE_PATH=os.path.join(CODE_DIR,subdir)
    # BASE_PATH_RAW_QASM=os.path.join(BASE_PATH, raw_subdir)
    # BASE_PATH_INITIAL_STATES=os.path.join(BASE_PATH, initial_states_subdir)
    # BASE_PATH_SPLIT_RESULTS=os.path.join(BASE_PATH, split_results_subdir)
    # # circs=circgenerator.generate_random_circuits(NUMBER_OF_QUBITS, CNOT_COUNT, NUMBER_OF_CIRCUITS)
    # # circs_with_checks=checksfinder.find_and_append_checks_to_all_circuits(circs, CNOT_COUNT)
    # # circgenerator.write_circs_no_checks_to_qasm_file(circs, NUMBER_OF_QUBITS, CNOT_COUNT)
    # raw_file_names=circtester.FilesManipulator.get_all_qasm_files_by_cnot(BASE_PATH_RAW_QASM, CNOT_COUNT)
    # initial_state_file_names=circgenerator.get_initial_state_file_names(raw_file_names)
    # assert len(raw_file_names)==len(initial_state_file_names)
    # # print(initial_state_file_names)
    # raw_circs_split=circgenerator.CircuitSplitter.split_circs_by_qubits(BASE_PATH_RAW_QASM, raw_file_names[0], 5, 5)
    # # initial_states=circtester.get_initial_states_from_qasm_files(BASE_PATH_INITIAL_STATES, initial_state_file_names)
    # # full_combined_circs=circtester.combine_circuits_for_simulate()
    # # circtester.simulate_all_given_circuits_concise(NUMBER_OF_QUBITS, CNOT_COUNT, full_combined_circs, initial_state_file_names, BASE_PATH_SPLIT_RESULTS)

    # # checksfinder.save_circuits_qasm(circs_with_checks, NUMBER_OF_QUBITS, CNOT_COUNT)
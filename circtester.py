from copy import deepcopy
from logging import error
from tkinter import N
from cirq.sim import density_matrix_simulator, simulator
import numpy as np
import psutil
from qiskit import transpile, QuantumCircuit
from qiskit.circuit import QuantumRegister
from qiskit.circuit.classicalregister import ClassicalRegister
from qiskit.visualization import circuit_drawer
from qiskit.opflow import X,Y,Z
import os
import pickle
import cirq
import json
from cirq.contrib.qasm_import import circuit_from_qasm
from os import listdir
from os.path import isfile
import scipy
from math import sqrt
from pandas import DataFrame as df
import time
from multiprocessing import Pool

class NoiselessCircuits:
    '''Testing Circuits'''
    def __init__(self, qiskit_circ_with_checks, cirq_circ_with_checks, cirq_circ_no_checks, qubits_label):
        self.qiskit_circ_with_checks = qiskit_circ_with_checks
        self.cirq_circ_with_checks = cirq_circ_with_checks
        self.cirq_circ_no_checks = cirq_circ_no_checks
        #This is useful for getting a handle on specific qubits for cirq versions of the circuits.
        self.qubits_label=qubits_label

class NoiselessMeasurementCircuits:
    def __init__(self, cirq_circ_measurements_with_checks, cirq_circ_measurements_no_checks, qubits_label):
        self.cirq_circ_measurements_with_checks=cirq_circ_measurements_with_checks
        self.cirq_circ_measurements_no_checks= cirq_circ_measurements_no_checks
        self.qubits_label=qubits_label

class CircuitMaker:
    '''Contains methods for generating the full circuits.'''
    def __init__(self, qasm_file, number_of_qubits):
        '''circ_pieces: iterable containing p1, main circuit, p2.'''
        circ_pieces=self.split_qasmfile_by_barrier(qasm_file)
        # assert len(circ_pieces)==3, f"{qasm_file} does not have the propper qasm format. There should be barriers between p1, main circuit, and p2."
        self.circ_pieces = circ_pieces
        self.number_of_qubits = number_of_qubits

    @staticmethod
    def qasm_file_to_cirq_circ(base_path, file_name):
        '''Outputs a cirq.circuit. Cirq can't import qasm with barriers.'''
        circ_pieces=CircuitMaker.split_qasmfile_by_barrier(os.path.join(base_path, file_name))
        new_circ=circ_pieces[0]
        if len(circ_pieces)>1:
            for piece in circ_pieces[1:]:
                new_circ.compose(piece, inplace=True)
        return circuit_from_qasm(new_circ.qasm())

    @staticmethod
    def qasm_string_to_cirq_circ(qasm_string):
        '''Outputs a cirq.circuit. Cirq can't import qasm with barriers.'''
        circ_pieces=CircuitMaker.split_qasm_str_by_barrier(qasm_string)
        new_circ=circ_pieces[0]
        if len(circ_pieces)>1:
            for piece in circ_pieces[1:]:
                new_circ.compose(piece, inplace=True)
        return circuit_from_qasm(new_circ.qasm())

    @staticmethod
    def split_qasmfile_by_barrier(qasm_file_path):
        '''Testing circuits: Split circuits by barrier.'''
        with open(qasm_file_path, "r") as file:
            qasm = file.read()
            prelude = []
            circuits = [[]]
            for line in qasm.splitlines():
                if any([line.startswith(t) for t in ['OPENQASM', 'include', 'qreg', 'creg']]):
                    prelude.append(line)
                elif line.startswith('barrier'):
                    circuits.append([])
                else:
                    circuits[-1].append(line)
            circuits_with_prelude = [prelude+circuit for circuit in circuits]
            return list(map(lambda x: QuantumCircuit.from_qasm_str("\n".join(x)), circuits_with_prelude))

    @staticmethod
    def split_qasm_str_by_barrier(qasm):
        '''Testing circuits: Split circuits by barrier.'''
        prelude = []
        circuits = [[]]
        for line in qasm.splitlines():
            if any([line.startswith(t) for t in ['OPENQASM', 'include', 'qreg', 'creg']]):
                prelude.append(line)
            elif line.startswith('barrier'):
                circuits.append([])
            else:
                circuits[-1].append(line)
        circuits_with_prelude = [prelude+circuit for circuit in circuits]
        return list(map(lambda x: QuantumCircuit.from_qasm_str("\n".join(x)), circuits_with_prelude))

    @staticmethod
    def add_rand_input_state(number_of_qubits, quantum_register, circ_with_checks, circ_no_checks):
        '''Testing circuits: Create a random state. Need to send both circs at the same time so they have the same random initial state.'''
        #Insert random state generator
        random_params = np.random.uniform(size=(number_of_qubits, 3))
        for i,qreg in enumerate(quantum_register):
            for j, pauli in enumerate([X,Y,Z]):
                rand_tuple=(random_params[i][j])
                circ_with_checks.compose((pauli * rand_tuple).exp_i().to_circuit(), [qreg], inplace=True)
                circ_no_checks.compose((pauli * rand_tuple).exp_i().to_circuit(), [qreg], inplace=True)

    def make_noiseless_circs(self):
        '''Testing circuits: Creates circs no checks and with checks. The circs have the same random initial state.
        Return type: list'''
        # The size of the circuit is number_of_qubit+1 since we have an ancila. We label the qubits
        # so that we can access the ancilla qubit later, i.e., the ancilla label is "q_{num}".format(number_of_qubits)
        number_of_qubits=self.number_of_qubits 
        qubits_label="q"
        quantum_register=QuantumRegister(number_of_qubits+1, qubits_label)
        ancilla_qreg=quantum_register[-1]
        qiskit_circ_with_checks=QuantumCircuit(quantum_register)
        qiskit_circ_no_checks=QuantumCircuit(quantum_register)
        self.add_rand_input_state(number_of_qubits, quantum_register[:number_of_qubits:], qiskit_circ_with_checks, qiskit_circ_no_checks)
        # For saving the initial qiskit circuit with added initial state. We do this because cirq cannot handle barriers. Thus,
        # the cirq printing of circuits is not so good.
        qiskit_circ_with_checks_store=deepcopy(qiskit_circ_no_checks)
        # # Initial hadamard.
        # qiskit_circ_with_checks.h(ancilla_qreg)
        # qiskit_circ_with_checks_store.barrier()
        # qiskit_circ_with_checks_store.h(ancilla_qreg)
        # Copy the pieces in the circuit.
        for elem in self.circ_pieces:
            qiskit_circ_with_checks.compose(elem, inplace=True)
            qiskit_circ_with_checks_store.barrier()
            qiskit_circ_with_checks_store.compose(elem, inplace=True)
        # # The final hadamard.
        # qiskit_circ_with_checks.h(ancilla_qreg)
        # qiskit_circ_with_checks_store.barrier()
        # qiskit_circ_with_checks_store.h(ancilla_qreg)
        # The no checks only uses the main compute circuit.
        qiskit_circ_no_checks.compose(self.circ_pieces[2], inplace=True)
        # print(qiskit_circ_no_checks)
        # Add the measurement for the qiskit circuit that we will print out.
        classical_register=ClassicalRegister(1, "c")
        qiskit_circ_with_checks_store.add_register(classical_register)
        qiskit_circ_with_checks_store.barrier()
        qiskit_circ_with_checks_store.measure(quantum_register[-1], classical_register[0])

        # We should transpile to a basis.
        basis_gates=['u1', 'u2', 'u3', 'cx']
        qiskit_circ_with_checks=transpile(qiskit_circ_with_checks, basis_gates=basis_gates, optimization_level=0)
        qiskit_circ_no_checks=transpile(qiskit_circ_no_checks, basis_gates=basis_gates, optimization_level=0)
        cirq_circ_with_checks =circuit_from_qasm(qiskit_circ_with_checks.qasm())
        cirq_circ_measurements_with_checks=deepcopy(cirq_circ_with_checks)
        self.add_measurements(cirq_circ_measurements_with_checks)
        cirq_circ_no_checks=circuit_from_qasm(qiskit_circ_no_checks.qasm())
        cirq_circ_measurements_no_checks=deepcopy(cirq_circ_no_checks)
        self.add_measurements(cirq_circ_measurements_no_checks)

        ancilla_qubit=cirq.NamedQubit(f"{qubits_label}_{number_of_qubits}")
        # Creates a channel that applies the zero projector. We use this to get the measurement zero outcome of the
        # density matrix. Since the resulting trial density matrix is unormalized we can get the percentages of outcomes that
        # we keep. In the protocol, we keep the zero measurement outcome results for the ancilla.
        projector0_channel=cirq.KrausChannel(
            kraus_ops=(np.array([[1,0],[0,0]]),),
            validate=False
        )
        cirq_circ_with_checks.append([projector0_channel.on(ancilla_qubit)]) 

        return (
            NoiselessCircuits(qiskit_circ_with_checks_store, cirq_circ_with_checks, cirq_circ_no_checks, qubits_label),
            NoiselessMeasurementCircuits(cirq_circ_measurements_with_checks, cirq_circ_measurements_no_checks, qubits_label))

    @staticmethod
    def add_measurements(circ):
        '''Adds measurements on all the qubits. Does this in place because 
        cirq.Circuit.append does append inplace.'''
        all_qubits=circ.all_qubits()
        ops=[]
        for qubit in all_qubits:
            ops.append(cirq.measure(qubit))
        moment=cirq.Moment(ops)
        # Does this inplace
        circ.append(moment, strategy=cirq.InsertStrategy.NEW)        

class CircuitSimulator:
    '''Testing circuits: For running the simulations.'''
    def __init__(self, noiseless_circs, number_of_qubits, keep_qubits):
        self.qubits_label=noiseless_circs.qubits_label
        self.cirq_circ_with_checks=noiseless_circs.cirq_circ_with_checks
        self.cirq_circ_no_checks=noiseless_circs.cirq_circ_no_checks
        self.number_of_qubits=number_of_qubits
        self.rho_correct=self.get_result_rho(self.cirq_circ_no_checks, number_of_qubits, keep_qubits)
        self.rho_checks=self.get_result_rho(self.cirq_circ_with_checks, number_of_qubits+1, keep_qubits)
        self.keep_qubits=keep_qubits
        self.sanity_check_fidelity=self.get_fidelity(self.rho_checks, self.rho_correct)

    def simulate_test(self, error_info):
        '''Uses cirq.simulate'''
        cirq_circ_with_checks=self.cirq_circ_with_checks
        cirq_circ_no_checks=self.cirq_circ_no_checks
        rho_correct=self.rho_correct
        number_of_qubits=self.number_of_qubits
        keep_qubits=self.keep_qubits
        sanity_check_fidelity=self.sanity_check_fidelity
        error_idx=error_info[0]
        single_qubit_error=error_info[1]

        print("getting result rho...")
        noisy_cirq_circ_with_checks=self.add_noise(cirq_circ_with_checks, single_qubit_error)
        noisy_rho_with_checks=self.get_result_rho(noisy_cirq_circ_with_checks, number_of_qubits+1, keep_qubits)
        print("taking trace...")
        ancilla_zero_outcome_probability=np.real(np.trace(noisy_rho_with_checks))
        fidelity_noisy_rho_with_check=self.get_fidelity(noisy_rho_with_checks* 1/ancilla_zero_outcome_probability, rho_correct)
        
        noisy_cirq_circ_no_checks=self.add_noise(cirq_circ_no_checks, single_qubit_error)
        noisy_rho_no_checks=self.get_result_rho(noisy_cirq_circ_no_checks, number_of_qubits, keep_qubits)
        fidelity_noisy_rho_no_check=self.get_fidelity(noisy_rho_no_checks, rho_correct)
        print(f"single qubit error rate: {single_qubit_error}")
        print(f"ancilla probability of 0 outcome: {ancilla_zero_outcome_probability}")
        print(f"fidelity no check: {fidelity_noisy_rho_no_check}")
        print(f"fidelity with check: {fidelity_noisy_rho_with_check}")
        print()

        return {"percent_results_before_postselect": 1, "percent_results_after_postselect": ancilla_zero_outcome_probability, "error_idx": error_idx, 
            "one_qubit_err": single_qubit_error, "two_qubit_err": 10*single_qubit_error, 
            "state_fidelity_no_checks_with_errors": fidelity_noisy_rho_no_check, 
            "state_fidelity_with_checks_with_errors": fidelity_noisy_rho_with_check, 
            "state_fidelity_with_checks_no_errors": sanity_check_fidelity}

    def add_noise(self, circ, single_qubit_error):
        '''Testing circuits: Uses Google Cirq. Adds noise to circ.'''
        two_qubit_error=10*single_qubit_error
        twoqubit_noise_model=cirq.DepolarizingChannel(p=two_qubit_error, n_qubits=2)
        singlequbit_noise_model=cirq.DepolarizingChannel(p=single_qubit_error)
        # all_qubits=circ.all_qubits()
        moments=[]
        noisy_cirq=cirq.Circuit()
        # Iterate through the moments. For two qubit gates we add a two qubit depolarization at a two qubit error rate.
        # print(all_qubits)
        for idx, moment in enumerate(circ):
            # Skip the initial state. Which consists of the first 3 moments.
            if idx <3:
                moments+=[moment]
                # We're dealing with the check circuit. Add a depolarizing noise after the hadamard gate.
                # W have to add the noise here because the hadamard goes into the first moment.
                if len(circ.all_qubits())==self.number_of_qubits+1 and idx == 2:
                    ancila=cirq.NamedQubit(f"{self.qubits_label}_{self.number_of_qubits}")
                    error_op=singlequbit_noise_model.on_each(ancila)
                    moments+=[cirq.ops.Moment(error_op)]
            else:
                error_ops=[]
                # TODO: Should we do this? 
                # For the moment find all the qubits with no operations and noise to them. 
                # for qubit in all_qubits:
                #     if qubit not in moment.qubits:
                #         # print(qubit)
                #         error_ops += singlequbit_noise_model.on_each(qubit)

                # Go through the operations in the moment. For 2 qubit operations use the
                # two qubit error. For everything else, i.e. single qubit gates, use single
                # qubit error. 
                for operation in moment.operations:
                    # print(operations)
                    # print(type(operation.gate))
                    # print(type(moment.qubits))
                    # print("all moment qubits: ", moment.qubits)
                    if isinstance(operation.gate, cirq.ops.common_gates.CXPowGate) or isinstance(operation.gate, cirq.ops.SwapPowGate):
                        # print("2 qubit ", operation)
                        # print(twoqubit_noise_model.on_each(operation.qubits))
                        error_ops += twoqubit_noise_model.on_each(operation.qubits)
                        # error_ops += singlequbit_noise_model.on_each(operation.qubits)

                    else:
                        # print("1 qubit ", operation)
                        # print(singlequbit_noise_model.on_each(operation.qubits))
                        error_ops += singlequbit_noise_model.on_each(operation.qubits)
                # print()
                # print(error_ops)
                moments+=[moment, cirq.ops.Moment(error_ops)]
            # print(moments)
        noisy_cirq+=cirq.Circuit(moments)
        # noisy_cirq=circ.with_noise(cirq.depolarize(p=single_qubit_error))
        return noisy_cirq

    def get_result_rho(self, circ, number_of_qubits, keep_qubits):
        '''Testing circuits: Uses Google Cirq. 
        number_of_qubits: total number of qubits in circ. 
        keep_qubits: is a list. 
        Returns resulting rho from simulation of circ.'''
        simulator=cirq.DensityMatrixSimulator()
        print("simulating....")
        zero_state = np.array([[1, 0], [0, 0]], dtype=np.complex64)
        initial_state = zero_state
        for _ in range(number_of_qubits-1):
            initial_state = np.kron(initial_state, zero_state)
        trial_result = simulator.simulate(circ, initial_state=initial_state)
        rho=trial_result.final_density_matrix
        final_size=len(keep_qubits)
        if number_of_qubits!=final_size:
            # Have to expand the indices for each qubit, e.g., for 2 qubits a_ij|i><j|-->a_ijkl|ij><kl|. Then 
            # to get the first state sum over indices j and l which are indices 1 and 3. For the second state
            # sum over i and k which are indices 0 and 2. In the cirq package, we specify which indices to keep.
            # But we only need to specify the first index not both. Thus, we can think of keep_indices as a list of
            # qubits to keep.
            print("taking partial trace...")
            # print(rho)
            rho_reduced=cirq.partial_trace(np.reshape(rho, [2,2]*(number_of_qubits)), keep_indices=keep_qubits)
            # We have to reshape back to a square matrix. We reasign to rho so we can just return rho.
            rho=np.reshape(rho_reduced, (2**final_size, 2**final_size))
        return rho

    def simulate_all_tests_parallel(self, pool, single_qubit_error_space):
        '''Uses cirq.simulate
        single_qubit_error_space: iterable containing the error numbers.'''
        #In some cases pool.imap_unordered needs to be wrapped in list in order to return properly. 
        #see: https://stackoverflow.com/questions/5481104/multiprocessing-pool-imap-broken
        return list(pool.imap(self.simulate_test, enumerate(single_qubit_error_space), chunksize=1))

    def simulate_all_tests(self, single_qubit_error_space):
        '''Uses cirq.simulate. Non parallel tests. 
        single_qubit_error_space: iterable containing the error numbers.'''
        results=[]
        for error_info in enumerate(single_qubit_error_space):
                results.append(self.simulate_test(error_info))
        return results

    @staticmethod
    def get_fidelity(rho1, rho2):
        '''Returns the fidelity through the qutip package.'''
        # Sums up the singular value of sqrt(rho1)*sqrt(rho2). The square of
        # this sum is equal to the fidelity. This calcuation is better than the builtin
        # functions in cirq and qutip. Credit to https://github.com/qutip/qutip/issues/925#issuecomment-542318121
        rho1_sqrt=scipy.linalg.sqrtm(rho1)
        rho2_sqrt=scipy.linalg.sqrtm(rho2)
        return scipy.linalg.svdvals(rho1_sqrt @ rho2_sqrt).sum()**2

class CircuitRunner:
    '''Testing circuits: For running the simulations.'''
    def __init__(self, noiseless_circs_measurements, number_of_qubits, number_of_runs):
        self.qubits_label=noiseless_circs_measurements.qubits_label
        self.cirq_circ_measurements_with_checks=noiseless_circs_measurements.cirq_circ_measurements_with_checks
        self.cirq_circ_measurements_no_checks=noiseless_circs_measurements.cirq_circ_measurements_no_checks
        self.number_of_qubits=number_of_qubits
        self.number_of_runs=number_of_runs
        # Get distribution returns the final counts.
        _, self.distribution_correct=self.get_distribution(self.cirq_circ_measurements_no_checks, self.number_of_runs)
        _, self.distribution_checks=self.get_distribution(self.cirq_circ_measurements_with_checks, self.number_of_runs)
        self.sanity_check_sso=self.get_sso(self.distribution_correct, self.distribution_checks)

    def run_test(self, error_info):
        '''Uses the probability distribution of outcomes.'''
        cirq_circ_measurements_with_checks=self.cirq_circ_measurements_with_checks
        cirq_circ_measurements_no_checks=self.cirq_circ_measurements_no_checks
        # number_of_qubits=self.number_of_qubits
        error_idx=error_info[0]
        single_qubit_error=error_info[1]
        distribution_correct=self.distribution_correct
        sanity_check_sso=self.sanity_check_sso
        number_of_runs=self.number_of_runs

        print("adding noise...")
        # print(cirq_circ_measurements_with_checks)
        noisy_cirq_circ_with_checks=self.add_noise(cirq_circ_measurements_with_checks, single_qubit_error)
        print("getting distribution...")
        # print(noisy_cirq_circ_with_checks)
        counts_checks, noisy_distribution_with_checks=self.get_distribution(noisy_cirq_circ_with_checks, number_of_runs)
        sso_noisy_distribution_with_checks=self.get_sso(noisy_distribution_with_checks, distribution_correct)
        print("taking trace...")

        # print(noisy_cirq_circ_with_checks)
        noisy_cirq_circ_no_checks=self.add_noise(cirq_circ_measurements_no_checks, single_qubit_error)
        _, noisy_distribution_no_checks=self.get_distribution(noisy_cirq_circ_no_checks, number_of_runs)
        sso_noisy_distribution_no_checks=self.get_sso(noisy_distribution_no_checks, distribution_correct)
        percent_keep=counts_checks/number_of_runs
        print(f"single qubit error rate: {single_qubit_error}")
        print(f"ancilla probability of 0 outcome: {percent_keep}")
        print(f"SSO no check: {sso_noisy_distribution_no_checks}")
        print(f"SSO with check: {sso_noisy_distribution_with_checks}")
        print()

        return {"percent_results_before_postselect": 1, "percent_results_after_postselect": percent_keep, "error_idx": error_idx, 
            "one_qubit_err": single_qubit_error, "two_qubit_err": 10*single_qubit_error, 
            "sso_no_checks_with_errors": sso_noisy_distribution_no_checks, 
            "sso_with_checks_with_errors": sso_noisy_distribution_with_checks, 
            "sany_check_sso": sanity_check_sso,
            "correct_distribution": self.distribution_correct}
    
    def add_noise(self, circ, single_qubit_error):
        '''Testing circuits: Uses Google Cirq. Adds noise to circ.'''
        two_qubit_error=10*single_qubit_error
        twoqubit_noise_model=cirq.DepolarizingChannel(p=two_qubit_error, n_qubits=2)
        singlequbit_noise_model=cirq.DepolarizingChannel(p=single_qubit_error)
        # all_qubits=circ.all_qubits()
        moments=[]
        noisy_cirq=cirq.Circuit()
        # Iterate through the moments. For two qubit gates we add a two qubit depolarization at a two qubit error rate.
        # print(all_qubits)
        for idx, moment in enumerate(circ):
            if idx==len(circ)-1:
                #The circuits have measurments. Ignore adding noise gates after measurements.
                moments+=[moment]
            # Skip the initial state. Which consists of the first 3 moments.
            elif idx <3:
                moments+=[moment]
                # We're dealing with the check circuit. Add a depolarizing noise after the hadamard gate.
                # We have to add the noise here because the hadamard goes into the first moment.
                if len(circ.all_qubits())==self.number_of_qubits+1 and idx == 2:
                    ancila=cirq.NamedQubit(f"{self.qubits_label}_{self.number_of_qubits}")
                    error_op=singlequbit_noise_model.on_each(ancila)
                    moments+=[cirq.ops.Moment(error_op)]
            else:
                error_ops=[]
                # TODO: Should we do this? 
                # For the moment find all the qubits with no operations and noise to them. 
                # for qubit in all_qubits:
                #     if qubit not in moment.qubits:
                #         # print(qubit)
                #         error_ops += singlequbit_noise_model.on_each(qubit)

                # Go through the operations in the moment. For 2 qubit operations use the
                # two qubit error. For everything else, i.e. single qubit gates, use single
                # qubit error. 
                for operation in moment.operations:
                    # print(operations)
                    # print(type(operation.gate))
                    # print(type(moment.qubits))
                    # print("all moment qubits: ", moment.qubits)
                    if (
                        isinstance(operation.gate, cirq.ops.common_gates.CXPowGate) 
                        or isinstance(operation.gate, cirq.ops.SwapPowGate)):
                        # print("2 qubit ", operation)
                        # print(twoqubit_noise_model.on_each(operation.qubits))
                        error_ops += twoqubit_noise_model.on_each(operation.qubits)
                        # error_ops += singlequbit_noise_model.on_each(operation.qubits)

                    else:
                        # print("1 qubit ", operation)
                        # print(singlequbit_noise_model.on_each(operation.qubits))
                        error_ops += singlequbit_noise_model.on_each(operation.qubits)
                # print()
                # print(error_ops)
                moments+=[moment, cirq.ops.Moment(error_ops)]
            # print(moments)
        noisy_cirq+=cirq.Circuit(moments)
        # noisy_cirq=circ.with_noise(cirq.depolarize(p=single_qubit_error))
        return noisy_cirq

    def get_distribution(self, circ, number_of_runs):
        '''Uses cirq.run'''
        t0=time.time()
        # print(circ)
        simulator=cirq.Simulator()
        # Returns a pnadas DataFrame
        distribution=simulator.sample(circ, repetitions=number_of_runs)
        print("finished running...")
        # print(f"running time: {time.time()-t0}")
        # print(f"original distribution: {distribution}")
        distribution.sort_index(axis="columns", ascending=True, inplace=True)
        # print("sorted", distribution)
        circ_num_qubits=len(circ.all_qubits())
        # Sort the columns
        if circ_num_qubits!=self.number_of_qubits:
            ancilla=f"{self.qubits_label}_{self.number_of_qubits}"
            distribution=distribution[distribution[ancilla] ==0]
            distribution=distribution.drop(ancilla, axis="columns")
            # print("dropped", distribution)
        # Concatenate the columns
        output_col="outcomes"
        concatenated=df(distribution.astype(str).agg("".join, axis="columns"),columns=[output_col])
        # print(f"combined: {concatenated}")
        histogram=df(concatenated[output_col].value_counts(), columns=[output_col])
        # print("counts", histogram)
        counts_final=histogram[output_col].sum()
        histogram=(histogram[output_col]/counts_final).to_dict()
        # print("histogram", histogram)
        print(f"counts final: {counts_final}")
        return counts_final, histogram

    def run_all_tests_parallel(self, pool, single_qubit_error_space):
        '''Uses probablity distribution of outcomes.'''
        #In some cases pool.imap_unordered needs to be wrapped in list in order to return properly. 
        #see: https://stackoverflow.com/questions/5481104/multiprocessing-pool-imap-broken
        return list(pool.imap(self.run_test, enumerate(single_qubit_error_space), chunksize=1))

    def run_all_tests(self, single_qubit_error_space):
        '''Uses probablity distribution of outcomes.'''
        results=[]
        for error_info in enumerate(single_qubit_error_space):
            results.append(self.run_test(error_info))
        return results

    @staticmethod
    def get_sso(dist1, dist2):
        '''Returns the square of the statistical overlap. 
        dist1 and dist2 are probability distributions.
        dist1: list
        dits2: list'''
        sum=0
        common_keys=dist1.keys() & dist2.keys()
        # print(f"dist1: {dist1}")
        # print(f"dist2: {dist2}")
        for key in common_keys:
            # print(f"key: {key}")
            # print(f"dist1[key]: {dist1[key]}")
            # print(f"dist2[key]: {dist2[key]}")
            sum+=sqrt(dist1[key]*dist2[key])
        return sum**2


class FilesManipulator:
    '''Class for dealing with files.'''
    def __init__(self, base_path, number_of_qubits, cnot_count):
        '''
        base_path: file path to the output folder.
        '''
        self.base_path=base_path
        self.number_of_qubits=number_of_qubits
        self.cnot_count=cnot_count

    @staticmethod
    def get_all_qasm_files_by_cnot(base_path_raw_qasm, cnot_count):
        all_files=[f for f in listdir(base_path_raw_qasm) if isfile(os.path.join(base_path_raw_qasm, f))]
        found_qasm_files=[]
        for file_name in all_files:
            if file_name.endswith(".qasm") and f"CNOTS_{cnot_count}_" in file_name:
                # split_name_nums=[int(elem) for elem in file_name.split("_") if elem.isdigit()]
                # if split_name_nums[1]==cnot_count:
                found_qasm_files.append(file_name)
        return found_qasm_files

    def result_exists(self, file_name, result_num=0):
        '''Testing circuits: Checks if the result file exists.'''
        #Strip file extension. This assumes that the file ends with "_.someextention".
        name_split=file_name.split("_")
        return os.path.isfile(os.path.join(self.base_path, f"{'_'.join(name_split[:-1:])}_result_{result_num}_.txt"))

    def result_sso_exists(self, file_name, result_num=0):
        '''Testing circuits: Checks if the result file exists.'''
        #Strip file extension. This assumes that the file ends with "_.someextention".
        name_split=file_name.split("_")
        name_split_nums=[int(elem) for elem in name_split if elem.isdigit()]
        # print(f"qubits_{name_split_nums[0]}_CNOTS_{name_split_nums[1]}_circuit_{name_split_nums[2]}_resultsso_{result_num}_.txt")
        return os.path.isfile(os.path.join(
            self.base_path, 
            f"qubits_{name_split_nums[0]}_CNOTS_{name_split_nums[1]}_circuit_{name_split_nums[2]}_resultsso_{result_num}_.txt"))

    def store_fidelity_results_verbose(self, circ_porp_file_name, noiseless_circs, results):
        '''Testing circuits: Stores all the results.'''
        base_path=self.base_path
        with open(os.path.join(base_path, circ_porp_file_name), "rb") as circ_file:
            circ_info=pickle.load(circ_file)

        #File naming stuff. 
        #Strip the extension.
        split_circ_file_name=circ_porp_file_name.split("_")
        # Pick the elements of the list that are numbers and then choose the last one.
        circ_num=[int(elem) for elem in split_circ_file_name if elem.isdigit()][-1]
        file_name_no_extension="_".join(split_circ_file_name[:-1])

        # print(output_file_name_obj)
        # print(os.path.join(base_path,output_file_name_obj))
        temp_file_number=0
        output_file_name_obj=f"{file_name_no_extension}_result_{temp_file_number}_.obj"
        while os.path.isfile(os.path.join(base_path,output_file_name_obj)):
            temp_file_number+=1
            output_file_name_obj=f"{file_name_no_extension}_result_{temp_file_number}_.obj"
        output_file_name_txt=f"{file_name_no_extension}_result_{temp_file_number}_.txt"
        # circ_full_no_measure=noiseless_circs.qiskit_circ_with_checks
        # qiskit_circ_full_no_measure=QuantumCircuit.from_qasm_str(cirq.qasm(circ_full_no_measure))

        # Dump all the results into a pickle
        with open(os.path.join(base_path, output_file_name_obj), "wb") as circ_file:
            pickle.dump({
                    "cx": circ_info["cx"], "rz": circ_info["rz"], "qubits": circ_info["qubits"], 
                    "circuit_num" : circ_num, "found_matches: ": circ_info["found_matches"],
                    "max_pauli_weight": circ_info["max_pauli_weight"], "max_pauli_str_p1": circ_info["max_pauli_str_p1"],
                    "max_pauli_str_p2": circ_info["max_pauli_str_p2"] ,"results": results}, circ_file)

        #Print text results to file
        output_file_txt_full_path=os.path.join(base_path, output_file_name_txt)
        qiskit_circ_with_checks= noiseless_circs.qiskit_circ_with_checks
        output_file_name_qasm=f"{file_name_no_extension}_result_{temp_file_number}_.qasm"
        qiskit_circ_with_checks.qasm(filename=os.path.join(base_path, output_file_name_qasm))
        circuit_drawer(qiskit_circ_with_checks, filename=output_file_txt_full_path)
        with open(output_file_txt_full_path, "a") as output_file_txt:
            output_file_txt.write("\n")
            output_file_txt.write(json.dumps(qiskit_circ_with_checks.count_ops()))
            for result in results:
                output_file_txt.write("\n")
                output_file_txt.write(f"Error idx: {result['error_idx']}\n")
                output_file_txt.write(f"One qubit error: {result['one_qubit_err']}\n")
                # print("One_qubit_err:", one_qubit_err)
                output_file_txt.write(f"Two qubit error: {result['two_qubit_err']}\n")
                # print("Two qubit error:", two_qubit_err)
                output_file_txt.write(f"ancilla probability of 0 outcome: {result['percent_results_after_postselect']}\n")
                output_file_txt.write(f"State fidelity no checks and with errors: {result['state_fidelity_no_checks_with_errors']}\n")
                output_file_txt.write(f"State fidelity with checks and with errors: {result['state_fidelity_with_checks_with_errors']}\n")
                output_file_txt.write(f"Sanity check fidelity with checks and no errors: {result['state_fidelity_with_checks_no_errors']}\n")

    def store_fidelity_results_concise(self, noiseless_circs, results, file_name_no_extension, extension=""):
        '''Testing circuits: Stores all the results concisely.'''
        base_path=self.base_path
        temp_file_number=0
        output_file_name_obj=f"{file_name_no_extension}_result{extension}_{temp_file_number}_.obj"
        while os.path.isfile(os.path.join(base_path,output_file_name_obj)):
            temp_file_number+=1
            output_file_name_obj=f"{file_name_no_extension}_result{extension}_{temp_file_number}_.obj"
        output_file_name_txt=f"{file_name_no_extension}_result{extension}_{temp_file_number}_.txt"
        # circ_full_no_measure=noiseless_circs.qiskit_circ_with_checks
        # qiskit_circ_full_no_measure=QuantumCircuit.from_qasm_str(cirq.qasm(circ_full_no_measure))

        # Dump all the results into a pickle
        with open(os.path.join(base_path, output_file_name_obj), "wb") as circ_file:
            pickle.dump({"results": results}, circ_file)

        #Print text results to file
        output_file_txt_full_path=os.path.join(base_path, output_file_name_txt)
        qiskit_circ_with_checks= noiseless_circs.qiskit_circ_with_checks
        output_file_name_qasm=f"{file_name_no_extension}_result{extension}_{temp_file_number}_.qasm"
        qiskit_circ_with_checks.qasm(filename=os.path.join(base_path, output_file_name_qasm))
        circuit_drawer(qiskit_circ_with_checks, filename=output_file_txt_full_path)
        with open(output_file_txt_full_path, "a") as output_file_txt:
            output_file_txt.write("\n")
            output_file_txt.write(json.dumps(qiskit_circ_with_checks.count_ops()))
            for result in results:
                output_file_txt.write("\n")
                output_file_txt.write(f"Error idx: {result['error_idx']}\n")
                output_file_txt.write(f"One qubit error: {result['one_qubit_err']}\n")
                # print("One_qubit_err:", one_qubit_err)
                output_file_txt.write(f"Two qubit error: {result['two_qubit_err']}\n")
                # print("Two qubit error:", two_qubit_err)
                output_file_txt.write(f"ancilla probability of 0 outcome: {result['percent_results_after_postselect']}\n")
                output_file_txt.write(f"State fidelity no checks and with errors: {result['state_fidelity_no_checks_with_errors']}\n")
                output_file_txt.write(f"State fidelity with checks and with errors: {result['state_fidelity_with_checks_with_errors']}\n")
                output_file_txt.write(f"Sanity check fidelity with checks and no errors: {result['state_fidelity_with_checks_no_errors']}\n")

    def store_sso_results(self, circ_porp_file_name, noiseless_circs, results):
        '''Testing circuits: Stores all the results.'''
        base_path=self.base_path
        with open(os.path.join(base_path, circ_porp_file_name), "rb") as circ_file:
            circ_info=pickle.load(circ_file)

        #File naming stuff. 
        #Strip the extension.
        split_circ_file_name=circ_porp_file_name.split("_")
        # Pick the elements of the list that are numbers and then choose the last one.
        circ_num=[int(elem) for elem in split_circ_file_name if elem.isdigit()][-1]
        file_name_no_extension="_".join(split_circ_file_name[:-1])

        # print(output_file_name_obj)
        # print(os.path.join(base_path,output_file_name_obj))
        temp_file_number=0
        output_file_name_obj=f"{file_name_no_extension}_resultsso_{temp_file_number}_.obj"
        while os.path.isfile(os.path.join(base_path,output_file_name_obj)):
            temp_file_number+=1
            output_file_name_obj=f"{file_name_no_extension}_resultsso_{temp_file_number}_.obj"
        output_file_name_txt=f"{file_name_no_extension}_resultsso_{temp_file_number}_.txt"
        # circ_full_no_measure=noiseless_circs.qiskit_circ_with_checks
        # qiskit_circ_full_no_measure=QuantumCircuit.from_qasm_str(cirq.qasm(circ_full_no_measure))

        # Dump all the results into a pickle
        with open(os.path.join(base_path, output_file_name_obj), "wb") as circ_file:
            pickle.dump({
                    "cx": circ_info["cx"], "rz": circ_info["rz"], "qubits": circ_info["qubits"], 
                    "circuit_num" : circ_num, "found_matches: ": circ_info["found_matches"],
                    "max_pauli_weight": circ_info["max_pauli_weight"], "max_pauli_str_p1": circ_info["max_pauli_str_p1"],
                    "max_pauli_str_p2": circ_info["max_pauli_str_p2"] ,"results": results}, circ_file)

        #Print text results to file
        output_file_txt_full_path=os.path.join(base_path, output_file_name_txt)
        qiskit_circ_with_checks= noiseless_circs.qiskit_circ_with_checks
        output_file_name_qasm=f"{file_name_no_extension}_resultsso_{temp_file_number}_.qasm"
        qiskit_circ_with_checks.qasm(filename=os.path.join(base_path, output_file_name_qasm))
        circuit_drawer(qiskit_circ_with_checks, filename=output_file_txt_full_path)
        with open(output_file_txt_full_path, "a") as output_file_txt:
            output_file_txt.write("\n")
            output_file_txt.write(json.dumps(qiskit_circ_with_checks.count_ops()))
            for result in results:
                output_file_txt.write("\n")
                output_file_txt.write(f"Error idx: {result['error_idx']}\n")
                output_file_txt.write(f"One qubit error: {result['one_qubit_err']}\n")
                # print("One_qubit_err:", one_qubit_err)
                output_file_txt.write(f"Two qubit error: {result['two_qubit_err']}\n")
                # print("Two qubit error:", two_qubit_err)
                output_file_txt.write(f"ancilla probability of 0 outcome: {result['percent_results_after_postselect']}\n")
                output_file_txt.write(f"State sso no checks and with errors: {result['sso_no_checks_with_errors']}\n")
                output_file_txt.write(f"State sso with checks and with errors: {result['sso_with_checks_with_errors']}\n")
                output_file_txt.write(f"Sanity check sso with checks and no errors: {result['sany_check_sso']}\n") 

    def get_files(self, start_circ_number, end_circ_number):
        '''Testing circuits: Get the desired files for testing. Grabs qasm files
        that have checks.'''
        base_path=self.base_path
        # Gets the files that match the string. Files include the path
        all_files=[f for f in listdir(base_path) if isfile(os.path.join(base_path, f))]
        rand_circ_files=[]
        circ_properties_files=[]
        for file in all_files:
            name_split=file.split("_")
            name_split_nums=[int(num) for num in name_split if num.isdigit()]
            # print(name_split)
            # print(name_split_nums)
            if (
                ".qasm" in name_split and "result" not in name_split
                and "raw" not in name_split 
                and name_split_nums[1]==self.cnot_count 
                and name_split_nums[0]==self.number_of_qubits 
                and start_circ_number<=name_split_nums[2]<=end_circ_number):
                rand_circ_files.append(file)
                circ_properties_files.append(f"{'_'.join(name_split[:-1])}_.obj")

        return rand_circ_files, circ_properties_files

    @staticmethod
    def store_basic_sso_result(base_path, circ_file_name, results):
        '''Testing circuits: Stores all the results.'''
        #File naming stuff. 
        #Strip the extension.
        split_circ_file_name=circ_file_name.split("_")
        # Pick the elements of the list that are numbers and then choose the last one.
        file_name_no_extension="_".join(split_circ_file_name[:-1])

        temp_file_number=0
        output_file_name_obj=f"{file_name_no_extension}_resultsso_{temp_file_number}_.obj"
        while os.path.isfile(os.path.join(base_path,output_file_name_obj)):
            temp_file_number+=1
            output_file_name_obj=f"{file_name_no_extension}_resultsso_{temp_file_number}_.obj"
        output_file_name_txt=f"{file_name_no_extension}_resultsso_{temp_file_number}_.txt"

        # Dump all the results into a pickle
        with open(os.path.join(base_path, output_file_name_obj), "wb") as circ_file:
            pickle.dump({"results": results}, circ_file)

        #Print text results to file
        output_file_txt_full_path=os.path.join(base_path, output_file_name_txt)
        with open(output_file_txt_full_path, "w") as output_file_txt:
            for result in results:
                output_file_txt.write("\n")
                output_file_txt.write(f"Error idx: {result['error_idx']}\n")
                output_file_txt.write(f"One qubit error: {result['one_qubit_err']}\n")
                # print("One_qubit_err:", one_qubit_err)
                output_file_txt.write(f"Two qubit error: {result['two_qubit_err']}\n")
                # print("Two qubit error:", two_qubit_err)
                output_file_txt.write(f"ancilla probability of 0 outcome: {result['percent_results_after_postselect']}\n")
                output_file_txt.write(f"State sso no checks and with errors: {result['sso_no_checks_with_errors']}\n")
                output_file_txt.write(f"State sso with checks and with errors: {result['sso_with_checks_with_errors']}\n")
                output_file_txt.write(f"Sanity check sso with checks and no errors: {result['sany_check_sso']}\n") 

    @staticmethod
    def remove_checks_in_complete_qasm(base_path, input_file_name, output_file_name):
        '''Initially the circuit has init|h|p1|main|p2|h. This returns
        init|main.'''
        circ_pieces=CircuitMaker.split_qasmfile_by_barrier(os.path.join(base_path, input_file_name))
        assert len(circ_pieces)>=6, f"number of pieces in {input_file_name} is wrong."
        initstate=circ_pieces[0]
        main_circ=circ_pieces[3]
        circ_new=deepcopy(initstate)
        circ_new.barrier()
        circ_new.compose(main_circ, inplace=True)
        circ_new.qasm(filename=os.path.join(base_path, output_file_name))

    @staticmethod
    def get_all_qasm_files(base_path):
        '''Gets all qasm files in base_path.'''
        all_files=[f for f in listdir(base_path) if isfile(os.path.join(base_path, f))]
        return [file_name for file_name in all_files if file_name.endswith(".qasm")]

    @staticmethod
    def get_all_qasm_files_with_checks(base_path):
        '''Grabs all qasm files with checks and no initial state from base_path.'''
        all_files=[f for f in listdir(base_path) if isfile(os.path.join(base_path, f))]
        circ_files=[]
        for file in all_files:
            name_split=file.split("_")
            name_split_nums=[int(num) for num in name_split if num.isdigit()]
            # print(name_split)
            # print(name_split_nums)
            if (
                ".qasm" in name_split and "result" not in name_split
                and "raw" not in name_split):
                circ_files.append(file)
        return circ_files

    @staticmethod
    def get_main_circuit_from_qasm_to_file(input_path, output_path, qasm_files, number_of_qubits):
        '''Starts with circuits: p1|main|p2
        Gets the main circuit from the qasm files.'''
        for file_name in qasm_files:
            # print(QuantumCircuit.from_qasm_file(os.path.join(input_path, file_name)))
            circ_pieces=CircuitMaker.split_qasmfile_by_barrier(os.path.join(input_path, file_name))
            assert len(circ_pieces)==3, f"number of pieces in {file_name} is wrong."
            main_circ=circ_pieces[1]
            output_file_name=f"{file_name[:-6]}_raw_.qasm"
            print(output_file_name)
            main_circ.qasm(filename=os.path.join(output_path, output_file_name))


    @staticmethod
    def get_result_qasm_and_original_prop(base_path, search_file):
        '''Testing circuits: Get the desired files for testing.'''
        # Gets the files that match the string. Files include the path
        all_files=[f for f in listdir(base_path) if isfile(os.path.join(base_path, f))]
        for current_file in all_files:
            name_split=search_file.split("_")
            name_split_nums=[int(elem) for elem in name_split if elem.isdigit()]
            name_split.remove("result")
            if current_file==search_file:
                return current_file, f"qubits_{name_split_nums[0]}_CNOTS_{name_split_nums[1]}_circuit_{name_split_nums[2]}_.obj"

    @staticmethod
    def remove_measurements_in_qasm_file(base_path, file_name):
        '''Strips measurements from qasm file.'''
        absolute_file_name=os.path.join(base_path, file_name)
        circ=QuantumCircuit.from_qasm_file(absolute_file_name)
        circ.remove_final_measurements(inplace=True)
        circ.qasm(filename=absolute_file_name)

    @staticmethod
    def add_all_measurements_in_qasm_file(base_path, file_name):
        '''Strips measurements from qasm file.'''
        absolute_file_name=os.path.join(base_path, file_name)
        circ=QuantumCircuit.from_qasm_file(absolute_file_name)
        circ.measure_all()
        circ.qasm(filename=absolute_file_name)

    @staticmethod
    def transpile_qasm_file(base_path, file_name, basis_gates):
        absolute_file_name=os.path.join(base_path, file_name)
        qiskit_circ=QuantumCircuit.from_qasm_file(absolute_file_name)
        qiskit_circ=transpile(qiskit_circ, basis_gates=basis_gates, optimization_level=0)
        qiskit_circ.qasm(filename=absolute_file_name)

def simulate_all_circuits_verbose(number_of_qubits, cnot_count, start_circ_number, end_circ_number, subdir="data"):
    '''Runs simulations on all the circuits in the given subdir (default is data).'''
    time0=time.time()
    print("running...")
    #Program parameters.
    # NUMBER_OF_QUBITS=int(sys.argv[1])
    # CNOT_COUNT=int(sys.argv[2])
    # START_CIRC_NUMBER=int(sys.argv[3])
    # END_CIRC_NUMBER=int(sys.argv[4])
    #Determines if we run parallel or not.
    PARALLEL=True
    # Gets the file path of the script
    CODE_DIR=os.path.abspath(os.path.dirname(__file__))
    BASE_PATH=os.path.join(CODE_DIR,subdir)
    #Error space
    SINGLE_QUBIT_ERROR_SPACE=np.logspace(-5, -2, num=21) #goes from 10^-5 to 10^-2
    files_manipulator=FilesManipulator(BASE_PATH, number_of_qubits, cnot_count)
    # Gets the files that match the string. Files include the path.
    # The returned files correspond accordingly, e.g., circ_file[0] and circ_properties_files[0] refer
    # to the same circuit.
    circ_files, circ_properties_files=files_manipulator.get_files(start_circ_number, end_circ_number)
    if PARALLEL:
        pool=Pool(psutil.cpu_count(logical=False))

    #Get the qasm and pickle info
    for file_idx, file_name in enumerate(circ_files):
        #Stopwatch
        time1=time.time()
        #If the file exists we already did this so just skip. Later on we can remove this for other initial states.
        if files_manipulator.result_exists(file_name):
            continue

        circuit_maker=CircuitMaker(os.path.join(BASE_PATH, file_name), number_of_qubits)
        noiseless_circs, _=circuit_maker.make_noiseless_circs()

        keep_qubits=list(range(number_of_qubits))
        circ_tester=CircuitSimulator(noiseless_circs, number_of_qubits, keep_qubits)
        assert circ_tester.sanity_check_fidelity>0.98, f"Sanity check fidelity {circ_tester.sanity_check_fidelity} failed for circuit {file_name}"
        print(f"sanity check fidelity: {circ_tester.sanity_check_fidelity}")

        if PARALLEL:
            results=circ_tester.simulate_all_tests_parallel(pool, SINGLE_QUBIT_ERROR_SPACE)
        else:
            results=circ_tester.simulate_all_tests(SINGLE_QUBIT_ERROR_SPACE)
        files_manipulator.store_fidelity_results_verbose(circ_properties_files[file_idx], noiseless_circs, results)

        print(f"file execution time {time.time()-time1}")
    if PARALLEL:
        pool.close()
        pool.join()
    print(f"total execution time {time.time()-time0}")
    print("Finished.")

def simulate_all_circuits_concise(number_of_qubits, cnot_count, start_circ_number, end_circ_number, base_path_to_circs):
    '''Runs simulations on all the circuits in the given subdir (default is data).'''
    time0=time.time()
    print("running...")
    #Program parameters.
    # NUMBER_OF_QUBITS=int(sys.argv[1])
    # CNOT_COUNT=int(sys.argv[2])
    # START_CIRC_NUMBER=int(sys.argv[3])
    # END_CIRC_NUMBER=int(sys.argv[4])
    #Determines if we run parallel or not.
    PARALLEL=True
    # Gets the file path of the script
    CODE_DIR=os.path.abspath(os.path.dirname(__file__))
    # BASE_PATH=os.path.join(CODE_DIR,subdir)
    #Error space
    SINGLE_QUBIT_ERROR_SPACE=np.logspace(-5, -2, num=21) #goes from 10^-5 to 10^-2
    files_manipulator=FilesManipulator(base_path_to_circs, number_of_qubits, cnot_count)
    # Gets the files that match the string. Files include the path.
    # The returned files correspond accordingly, e.g., circ_file[0] and circ_properties_files[0] refer
    # to the same circuit.
    circ_files, _=files_manipulator.get_files(start_circ_number, end_circ_number)
    if PARALLEL:
        pool=Pool(psutil.cpu_count(logical=False))

    #Get the qasm and pickle info
    for file_name in circ_files:
        #Stopwatch
        time1=time.time()
        #If the file exists we already did this so just skip. Later on we can remove this for other initial states.
        # if files_manipulator.result_exists(file_name):
        #     continue

        circuit_maker=CircuitMaker(os.path.join(base_path_to_circs, file_name), number_of_qubits)
        noiseless_circs, _=circuit_maker.make_noiseless_circs()

        keep_qubits=list(range(number_of_qubits))
        circ_tester=CircuitSimulator(noiseless_circs, number_of_qubits, keep_qubits)
        assert circ_tester.sanity_check_fidelity>0.98, f"Sanity check fidelity {circ_tester.sanity_check_fidelity} failed for circuit {file_name}"
        print(f"sanity check fidelity: {circ_tester.sanity_check_fidelity}")

        if PARALLEL:
            results=circ_tester.simulate_all_tests_parallel(pool, SINGLE_QUBIT_ERROR_SPACE)
        else:
            results=circ_tester.simulate_all_tests(SINGLE_QUBIT_ERROR_SPACE)
        file_name_no_extension=file_name[:-6] #remove _.qasm
        files_manipulator.store_fidelity_results_concise(noiseless_circs, results, file_name_no_extension,)

        print(f"file execution time {time.time()-time1}")
    if PARALLEL:
        pool.close()
        pool.join()
    print(f"total execution time {time.time()-time0}")
    print("Finished.")

def simulate_all_given_circuits_concise(number_of_qubits, cnot_count, circs, circ_file_names, path_to_store_results):
    '''Runs simulations on all the circuits in the given subdir (default is data).'''
    time0=time.time()
    print("running...")
    #Program parameters.
    # NUMBER_OF_QUBITS=int(sys.argv[1])
    # CNOT_COUNT=int(sys.argv[2])
    # START_CIRC_NUMBER=int(sys.argv[3])
    # END_CIRC_NUMBER=int(sys.argv[4])
    #Determines if we run parallel or not.
    PARALLEL=True
    # Gets the file path of the script
    CODE_DIR=os.path.abspath(os.path.dirname(__file__))
    # BASE_PATH=os.path.join(CODE_DIR,subdir)
    #Error space
    SINGLE_QUBIT_ERROR_SPACE=np.logspace(-5, -2, num=21) #goes from 10^-5 to 10^-2
    files_manipulator=FilesManipulator(path_to_store_results, number_of_qubits, cnot_count)
    # Gets the files that match the string. Files include the path.
    # The returned files correspond accordingly, e.g., circ_file[0] and circ_properties_files[0] refer
    # to the same circuit.
    circ_files, _=files_manipulator.get_files(start_circ_number, end_circ_number)
    if PARALLEL:
        pool=Pool(psutil.cpu_count(logical=False))

    #Get the qasm and pickle info
    for file_name in circ_files:
        #Stopwatch
        time1=time.time()
        #If the file exists we already did this so just skip. Later on we can remove this for other initial states.
        # if files_manipulator.result_exists(file_name):
        #     continue

        circuit_maker=CircuitMaker(os.path.join(path_to_store_results, file_name), number_of_qubits)
        noiseless_circs, _=circuit_maker.make_noiseless_circs()

        keep_qubits=list(range(number_of_qubits))
        circ_tester=CircuitSimulator(noiseless_circs, number_of_qubits, keep_qubits)
        assert circ_tester.sanity_check_fidelity>0.98, f"Sanity check fidelity {circ_tester.sanity_check_fidelity} failed for circuit {file_name}"
        print(f"sanity check fidelity: {circ_tester.sanity_check_fidelity}")

        if PARALLEL:
            results=circ_tester.simulate_all_tests_parallel(pool, SINGLE_QUBIT_ERROR_SPACE)
        else:
            results=circ_tester.simulate_all_tests(SINGLE_QUBIT_ERROR_SPACE)
        file_name_no_extension=file_name[:-6] #remove _.qasm
        files_manipulator.store_fidelity_results_concise(noiseless_circs, results, file_name_no_extension,)

        print(f"file execution time {time.time()-time1}")
    if PARALLEL:
        pool.close()
        pool.join()
    print(f"total execution time {time.time()-time0}")
    print("Finished.")
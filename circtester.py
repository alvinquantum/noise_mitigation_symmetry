from copy import deepcopy
import numpy as np
from qiskit import transpile, QuantumCircuit
from qiskit.circuit import QuantumRegister
from qiskit.visualization import circuit_drawer
from qiskit.opflow import X,Y,Z
import os
import pickle
import cirq
import json
from cirq.contrib.qasm_import import circuit_from_qasm
from os import listdir
from os.path import isfile
import utilities

class CircuitTester:
    '''Testing circuits: For running the simulations.'''
    def __init__(self, cirq_circ_with_checks, cirq_circ_no_checks, number_of_qubits, rho_correct, sanity_check_fidelity, keep_qubits):
        self.cirq_circ_with_checks=cirq_circ_with_checks
        self.cirq_circ_no_checks=cirq_circ_no_checks
        self.number_of_qubits=number_of_qubits
        self.rho_correct=rho_correct
        self.keep_qubits=keep_qubits
        self.sanity_check_fidelity=sanity_check_fidelity

    def run_test(self, error_info):
        '''Testing circuits: runs the test.'''
        cirq_circ_with_checks=self.cirq_circ_with_checks
        cirq_circ_no_checks=self.cirq_circ_no_checks
        rho_correct=self.rho_correct
        number_of_qubits=self.number_of_qubits
        keep_qubits=self.keep_qubits
        sanity_check_fidelity=self.sanity_check_fidelity
        error_idx=error_info[0]
        single_qubit_error=error_info[1]

        '''Testing circuits: get testing results'''
        print("adding noise...")
        noisy_cirq_circ_with_checks=add_noise(cirq_circ_with_checks, single_qubit_error) 
        # print(noisy_cirq_circ_with_checks)   
        print("getting result rho...")    
        noisy_rho_with_checks=get_result_rho(noisy_cirq_circ_with_checks, number_of_qubits+1, keep_qubits)
        print("taking trace...")
        ancilla_zero_outcome_probability=np.trace(noisy_rho_with_checks)
        fidelity_noisy_rho_with_check=cirq.fidelity(np.around(noisy_rho_with_checks* 1/ancilla_zero_outcome_probability, 5), 
            rho_correct, qid_shape=(2**(number_of_qubits),), validate=False)
        
        noisy_cirq_circ_no_checks=add_noise(cirq_circ_no_checks, single_qubit_error)
        # print(noisy_cirq_circ_with_checks)
        noisy_rho_no_checks=get_result_rho(noisy_cirq_circ_no_checks, number_of_qubits, keep_qubits)
        fidelity_noisy_rho_no_check=cirq.fidelity(noisy_rho_no_checks, rho_correct, qid_shape=(2**(number_of_qubits),), validate=False)
        print(f"single qubit error rate: {single_qubit_error}")
        print(f"ancilla 0 prob outcome: {ancilla_zero_outcome_probability}")
        print(f"fidelity no check: {fidelity_noisy_rho_no_check}")
        print(f"fidelity with check: {fidelity_noisy_rho_with_check}")
        print()

        return {"num_results_before_postselect": 1, "num_results_after_postselect": ancilla_zero_outcome_probability, "error_idx": error_idx, 
            "one_qubit_err": single_qubit_error, "two_qubit_err": 10*single_qubit_error, 
            "state_fidelity_no_checks_with_errors": fidelity_noisy_rho_no_check, 
            "state_fidelity_with_checks_with_errors": fidelity_noisy_rho_with_check, 
            "state_fidelity_with_checks_no_errors": sanity_check_fidelity}

def result_exists(base_path, file_name, result_num):
    '''Testing circuits: Checks if the result file exists.'''
    name_split=file_name.split("_")
    result_name=f"{'_'.join(name_split[:-1])}_result_{result_num}_.txt"
    # print(file_name)
    # print(result_name)
    output_file_txt_path=os.path.join(base_path, result_name)
    # print(output_file_txt_path)
    return os.path.isfile(output_file_txt_path)

def store_results(circ, circ_full_with_ancilla, base_path, circ_file_name, number_of_qubits, results, cnot_count, rz_count, circuit_num,
    count, p2_weight, pauli_str_p1, pauli_str_p2):
    '''Testing circuits: Stores all the results.'''
    #File naming stuff. 
    #Strip the extension.
    split_circ_file_name=circ_file_name.split("_")
    file_name_no_extension="_".join(split_circ_file_name[:-1])

    # print(output_file_name_obj)
    # print(os.path.join(base_path,output_file_name_obj))
    temp_file_number=0
    output_file_name_obj=f"{file_name_no_extension}_result_{temp_file_number}_.obj"
    while os.path.isfile(os.path.join(base_path,output_file_name_obj)):
        temp_file_number+=1
        output_file_name_obj=f"{file_name_no_extension}_result_{temp_file_number}_.obj"
    output_file_name_txt=f"{file_name_no_extension}_result_{temp_file_number}_.txt"
    output_file_name_qasm=f"{file_name_no_extension}_result_{temp_file_number}_.qasm"
    circ_full_no_measure=circ_full_with_ancilla[:-1:]
    qiskit_circ_full_no_measure=QuantumCircuit.from_qasm_str(cirq.qasm(circ_full_no_measure))
    qiskit_circ_full_no_measure.qasm(filename=os.path.join(base_path, output_file_name_qasm))

    # Dump all the results into a pickle
    with open(os.path.join(base_path, output_file_name_obj), "wb") as circ_file:
        pickle.dump(
            {"cx": cnot_count, "rz": rz_count, "qubits": number_of_qubits, "circuit_num" : circuit_num, "found_matches: ": count,
            "max_pauli_weight": p2_weight, "max_pauli_str_p1": pauli_str_p1, "max_pauli_str_p2": pauli_str_p2 ,"results": results}, circ_file)

    #Print text results to file
    output_file_txt_path=os.path.join(base_path, output_file_name_txt)
    circuit_drawer(qiskit_circ_full_no_measure, filename=output_file_txt_path)
    with open(output_file_txt_path, "a") as output_file_txt:
        output_file_txt.write("\n")
        qiskit_circ=QuantumCircuit.from_qasm_str(cirq.qasm(circ))
        output_file_txt.write(json.dumps(qiskit_circ.count_ops()))
        for result in results:
            output_file_txt.write("\n")
            output_file_txt.write(f"Error idx: {result['error_idx']}\n")
            output_file_txt.write(f"One qubit error: {result['one_qubit_err']}\n")
            # print("One_qubit_err:", one_qubit_err)
            output_file_txt.write(f"Two qubit error: {result['two_qubit_err']}\n")
            # print("Two qubit error:", two_qubit_err)
            output_file_txt.write(f"State fidelity no checks and with errors: {result['state_fidelity_no_checks_with_errors']}\n")
            output_file_txt.write(f"State fidelity with checks and with errors: {result['state_fidelity_with_checks_with_errors']}\n")
            output_file_txt.write(f"Sanity check fidelity with checks and no errors: {result['state_fidelity_with_checks_no_errors']}\n")

def get_files(base_path, number_of_qubits, cnot_count, start_circ_number, end_circ_number):
    '''Testing circuits: Get the desired files for testing.'''

    # Gets the files that match the string. Files include the path
    all_files=[f for f in listdir(base_path) if isfile(os.path.join(base_path, f))]
    rand_circ_files=[]
    circ_properties_files=[]
    for file in all_files:
        name_split=file.split("_")
        name_split_nums=[int(num) for num in name_split if num.isdigit()]
        # print(name_split)
        # print(name_split_nums)
        if (".qasm" in name_split and "result" not in name_split and name_split_nums[1]==cnot_count and name_split_nums[0]==number_of_qubits 
            and start_circ_number<=name_split_nums[2]<=end_circ_number):
            rand_circ_files.append(file)
            circ_properties_files.append(f"{'_'.join(name_split[:-1])}_.obj")

    return rand_circ_files, circ_properties_files

def add_rand_input_state(number_of_qubits, quantum_register, circ_with_checks, circ_no_checks):
    '''Testing circuits: Create a random state. Need to send both circs at the same time so they have the same random initial state.'''
    #Insert random state generator
    random_params = np.random.uniform(size=(number_of_qubits, 3))
    for i,qreg in enumerate(quantum_register):
        for j, pauli in enumerate([X,Y,Z]):
            rand_tuple=(random_params[i][j])
            circ_with_checks.compose((pauli * rand_tuple).exp_i().to_circuit(), [qreg], inplace=True)
            circ_no_checks.compose((pauli * rand_tuple).exp_i().to_circuit(), [qreg], inplace=True)

def get_checks_properties(base_path, file_name):
    '''Testing circuits: Gets the check properties. Return type: CheckProperties'''
    with open(os.path.join(base_path, file_name), "rb") as circ_file:
        circ_info=pickle.load(circ_file)
            # circ=circ_info["circ"]
    return utilities.ChecksProperties(
        circ_info["found_matches"], circ_info["max_pauli_weight"], 
        circ_info["max_pauli_str_p1"], circ_info["max_pauli_str_p2"])

class NoiselessCircuits:
    '''Testing Circuits'''
    def __init__(self, qiskit_circ_with_checks, cirq_circ_with_checks, cirq_circ_no_checks):
        self.qiskit_circ_with_checks = qiskit_circ_with_checks
        self.cirq_circ_with_checks = cirq_circ_with_checks
        self.cirq_circ_no_checks = cirq_circ_no_checks

class NoisyCircuits:
    '''Testing Circuits'''
    def __init__(self, ):
        self.cirq_circ_with_checks_noisy
        self.cirq_circ_no_checks_noisy

def create_circs(number_of_qubits, circ_pieces):
    '''Testing circuits: Creates circs no checks and with checks. The circs have the same random initial state.
    Return type: list'''
    # The size of the circuit is number_of_qubit+1 since we have an ancila. We label the qubits
    # so that we can access the ancilla qubit later, i.e., the ancilla label is "q_{num}".format(number_of_qubits) 
    qubits_label="q"
    quantum_register=QuantumRegister(number_of_qubits+1, qubits_label)
    qiskit_circ_with_checks=QuantumCircuit(quantum_register)
    #Create the circuit with no checks
    circ_no_checks=QuantumCircuit(quantum_register)
    add_rand_input_state(number_of_qubits, quantum_register[:number_of_qubits:], qiskit_circ_with_checks, circ_no_checks)
    # Save the initial state circuit.
    circ_initial_state=deepcopy(circ_no_checks)
    # print("random initial state", rand_initial_state)
    # circ_with_checks.barrier()
    qiskit_circ_with_checks.h(number_of_qubits)
    for elem in circ_pieces:
        # circ_with_checks.barrier()
        qiskit_circ_with_checks.compose(elem, inplace=True)
    # circ_with_checks.barrier()
    circ_no_checks.compose(circ_pieces[1], inplace=True)
    qiskit_circ_with_checks.h(number_of_qubits)

    #Test
    basis_gates=['u1', 'u2', 'u3', 'cx']
    qiskit_circ_with_checks=transpile(qiskit_circ_with_checks, basis_gates=basis_gates, optimization_level=0)
    circ_no_checks=transpile(circ_no_checks, basis_gates=basis_gates, optimization_level=0)
    circ_initial_state=transpile(circ_initial_state, basis_gates=basis_gates, optimization_level=0)
    cirq_circ_with_checks =circuit_from_qasm(qiskit_circ_with_checks.qasm())
    cirq_circ_no_checks=circuit_from_qasm(circ_no_checks.qasm())
    cirq_circ_initial_state=circuit_from_qasm(circ_initial_state.qasm())
    print(qiskit_circ_with_checks)
    print(circ_no_checks)

    ancilla_qubit=cirq.NamedQubit(f"{qubits_label}_{number_of_qubits}")
    # Creates a channel that applies the zero projector. We use this to get the measurement zero outcome of the
    # density matrix. Since the resulting trial density matrix is unormalized we can get the percentages of outcomes that
    # we discard.
    projector0_channel=cirq.KrausChannel(
        kraus_ops=(np.array([[1,0],[0,0]]),),
        validate=False
    )
    cirq_circ_with_checks.append([projector0_channel.on(ancilla_qubit)])

    return NoiselessCircuits(qiskit_circ_with_checks, cirq_circ_with_checks, cirq_circ_no_checks)

def split_circuit_by_barrier(qasm_file_path):
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

def get_result_rho(circ, number_of_qubits, keep_qubits):
    '''Testing circuits: Uses Google Cirq. 
    number_of_qubits: total number of qubits in circ. 
    keep_qubits: is a list. 
    Returns resulting rho from simulation of circ.'''
    simulator=cirq.DensityMatrixSimulator()
    print("simulating....")
    trial_result=simulator.simulate(circ)
    rho= np.around(trial_result.final_density_matrix, 5)
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
        # print("reduced....")
        # print(rho)
        # print(rho.shape)
    return rho

def add_noise(circ, single_qubit_error):
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
        # Skip the initial state.
        if idx <3:
            moments+=[moment]
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
                #     # print(twoqubit_noise_model.on_each(operation.qubits))
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

if __name__ == "__main__":
    pass
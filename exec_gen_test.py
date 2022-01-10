#!/usr/bin/env python3
import os

# For automating the process on a local machine. Don't use for HPC.
if __name__ == "__main__":
    num_circs=50
    qubits=5
    start_circ=0
    end_circ=num_circs-1
    cnots=[1,5,10,15,20,25,30,35,40]
    dir=os.path.abspath(os.path.dirname(__file__))
    # 5 qubit process
    # for cnot_count in cnots:
    #     os.system(f"python {os.path.join(dir, 'gen_rand_circs.py')} {qubits} {cnot_count} {num_circs}")
    # for cnot_count in cnots:
    #     os.system(f"python {os.path.join(dir, 'test_circuits_cirq.py')} {qubits} {cnot_count} {start_circ} {end_circ}")
    # os.system(f"python {os.path.join(dir, 'create_graphs.py')}")

    # 10 qubit. Do a step number of circiuts for each cnot count at a time.
    step=5
    for start_circ in list(range(0, num_circs, step)):
        for cnot_count in cnots:
            os.system(f"python {os.path.join(dir, 'test_circuits_cirq.py')} {qubits} {cnot_count} {start_circ} {start_circ+step-1}")
    os.system(f"python {os.path.join(dir, 'create_graphs.py')}")
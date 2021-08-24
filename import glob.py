import glob
import sys
import uuid
from pathlib import Path
import qiskit
from qsearch import gates, unitaries, Project, leap_compiler, post_processing, parallelizers,multistart_solvers, gatesets, assemblers, utils
from itertools import combinations
import numpy as np

if __name__ == '__main__':
    # @Ashish: update output file name as needed
    outpath = "example.qasm"
    if Path(outpath).exists():
        print(f"{outpath} already present, exiting")
        sys.exit(0)
    # Simple example: compile CNOT
    U1 = np.array([[1,0,0,0],[0,1,0,0],[0,0,0,-1j],[0,0,-1j,0]], dtype=complex)
    num_qubits = int(np.log2(U1.shape[0]))
    print(f"Compiling unitary on {num_qubits} qubits:")
    print(U1)
    circ_name= "Perm12"
    proj_name = 'Permutation_SV' + '-' + str(uuid.uuid4().hex)
    with Project(proj_name) as project:
        project.add_compilation(circ_name, U1)
        project["compiler_class"] = leap_compiler.LeapCompiler
        project["min_depth"] = 3
        project["gateset"] = gatesets.QubitCNOTAdjacencyList([[i,j] for i, j in combinations(range(num_qubits), r=2)])
        project["verbosity"] = 2
        project.run()
        print('Done running')
        project.post_process(post_processing.LEAPReoptimizing_PostProcessor(), solver=multistart_solvers.MultiStart_Solver(24),  parallelizer=parallelizers.ProcessPoolParallelizer, weight_limit=5)
        print('Done post processing')
        qasm = project.assemble(circ_name, assembler=assemblers.ASSEMBLER_IBMOPENQASM)
    qc = qiskit.QuantumCircuit.from_qasm_str(qasm)
    print(qc.draw())
    print(f"Generated circuit with {qc.num_nonlocal_gates()} CNOTs")
    qc.qasm(filename=outpath)
#!/usr/bin/env python3

import os, sys, glob, pickle, operator
import matplotlib.pyplot as plt
import numpy as np
from os import listdir
from os.path import isfile

from qiskit.quantum_info.synthesis import one_qubit_decompose

def get_lines(file_path):
    '''Return the lines of file but in reverse order.'''
    output=[]
    with open(file_path, "r", errors="ignore") as file:
        lines=file.readlines()
        for line in reversed(lines):
            if isinstance(line, str):
                # print(line)
                # print(type(line))
                output.append(line)
    return output

def get_line_str_stripped(file_lines, str):
    '''Get count associated with str'''
    for line in file_lines:
        if line.startswith(str):
            return line[len(str):]
    return "0"

def create_cnot_plot():
    cnot_counts=[]
    max_weight=[]
    for file in os.listdir():
        if file.endswith(".txt") and not file.endswith("result0.txt"):
            print(file)
            file_lines_reverse=get_lines(file)
            for item in file_lines_reverse:
                if item.startswith("CNOT count: "):
                    cnot_counts.append(item.split(" ")[2].strip("\n"))
            # for item in file_lines_reverse:
                elif item.startswith("Max Weight: "):
                    max_weight.append(item.split(" ")[2].strip("\n"))
    # print(len(cnot_counts))
    # print(len(max_weight))
            # print(file_lines_reverse)
            # print(next(obj for obj in file_lines_reverse if obj.startswith("CNOT count: ")))
    #         cnot_counts.append(next(obj for obj in file_lines_reverse if obj.startswith("CNOT count: ")).split(" ")[1])
    #         max_weight.append(next(obj for obj in file_lines_reverse if obj.startswith("Max weight: ")).split(" ")[1])
    cnot_counts=list(map(int, cnot_counts))
    max_weight=list(map(int, max_weight))
    # print(cnot_counts)
    # print(max_weight)
    # X axis is cnot_counts and y is max_weight of Pauli P2.
    plt.plot(cnot_counts,max_weight)
    plt.xlabel("CNOT count")
    plt.ylabel("Pauli Weight")
    plt.savefig(PATH+"weight_vs_cnot_count.png")
    plt.show()

def get_fidelities(erroridx_files, check_string):
    '''Get the number after the string. Returns list'''
    #same depth and error so average the gain
    fidelities=[]
    for file in erroridx_files:
        lines=get_lines(file)
        for item in lines:
            if item.startswith(check_string):
                fidelity_no_checks_temp=item.split(": ")[1].strip("\n")
                fidelities.append(float(fidelity_no_checks_temp))
            #     # if depth==25:
            #     #     fidelity_no_checks_temp=item.split(": ")[1].strip("\n")
            #     #     fidelity_no_checks.append(float(fidelity_no_checks_temp))
            #     #     # print("file: "+ file+ "error: "+ str(ONE_QUBIT_ERROR_SPACE[erroridx])+" State fidelity no checks with errors: "+ fidelity_no_checks_temp)
            #     #     # print(float(fidelity_no_checks_temp))
    return fidelities

def get_file_info(fn):
    '''returns a list of qubit_no, depth_no, circuit_no, [erroridx, result_no]'''
    temp_file=[elem for elem in fn.split("_") if elem.isdigit()]
    return temp_file

def create_fidelity_plot_depth():
    global NUMBER_OF_QUBITS, SUBDIR, PATH, ONE_QUBIT_ERROR_SPACE, NUMBER_OF_ERROR_POINTS
    # Gets the file path of the script.
    code_dir=sys.path[0]
    # Get all the result txt files.
    rand_circ_files=glob.glob("*.txt")
    for file_name in rand_circ_files:
        file_info=get_file_info(file_name)
        print(file_info)
    # print(rand_circ_files)
    # depths=(1,5,10,15,20,25)
    # #For each depth we plot fidelities vs single qubit error.
    # for depth in depths:
    #     plot_avg_fidelity_gain=[]
    #     plot_avg_fidelity_checks=[]
    #     plot_avg_fidelity_no_checks=[]
    #     erroridx_list=[]
    #     error_points=[]
    #     depth_str_check="depth"+str(depth)+"_"
    #     # Separate the current depth files.
    #     depth_files=[elem for elem in rand_circ_files if depth_str_check in elem]    
    #     for erroridx in range (NUMBER_OF_ERROR_POINTS):
    #         erroridx_str_check="erroridx"+str(erroridx)+"_"
    #         # Separate the current erroridx files.
    #         erroridx_files=[elem for elem in depth_files if erroridx_str_check in elem]
    #         # print(len(erroridx_files))
    #         if len(erroridx_files)!=0:
    #             # The error point has data.
    #             erroridx_list.append(erroridx)
    #             # get average gain
    #             fidelity_checks_temp=get_fidelities(erroridx_files, "State fidelity with checks and errors: ")
    #             fidelity_no_checks_temp=get_fidelities(erroridx_files, "State fidelity no checks with errors: ")

    #             assert len(fidelity_checks_temp)==len(fidelity_no_checks_temp), "Number of fidelities do not match."
    #             plot_avg_fidelity_checks.append(sum(fidelity_checks_temp)/len(fidelity_checks_temp))
    #             plot_avg_fidelity_no_checks.append(sum(fidelity_no_checks_temp)/len(fidelity_no_checks_temp))
    #             plot_avg_fidelity_gain.append(plot_avg_fidelity_checks[-1]-plot_avg_fidelity_no_checks[-1])

    #     #Convert the erroridx to error values.
    #     for elem in erroridx_list:
    #         error_points.append(ONE_QUBIT_ERROR_SPACE[elem])
    #     # print(error_points)        

    #     plt.plot(error_points, plot_avg_fidelity_gain)
    #     plt.title("Qubits: "+str(NUMBER_OF_QUBITS)+ " Depth: " +str(depth))
    #     plt.xscale("log")
    #     plt.xlabel("1 Qubit Error")
    #     plt.ylabel("Average Fidelity Gain")
    #     plt.savefig(PATH+"Qubits"+str(NUMBER_OF_QUBITS)+ "_Depth" +str(depth)+".png")
    #     plt.show()


def create_fidelity_plot_cnots():
    global CNOT_COUNT, NUMBER_OF_QUBITS, SUBDIR, BASE_PATH, ONE_QUBIT_ERROR_SPACE, NUMBER_OF_ERROR_POINTS, ONE_QUBIT_ERROR_SPACE
    # Gets the file path of the script.
    code_dir=sys.path[0]
    all_files=[f for f in listdir(BASE_PATH) if isfile(os.path.join(BASE_PATH, f))]
    # print(all_files)
    circ_result_property_files=[]
    for file in all_files:
        name_split=file.split("_")
        # Get the digits in the names
        name_split_nums=[int(num) for num in name_split if num.isdigit()]
        # print(name_split)
        # print(name_split_nums)
        if ".obj" in name_split and "result" in name_split and name_split_nums[1]==CNOT_COUNT: #and start_circ_number<=name_split_nums[2]<=end_circ_number:
            circ_result_property_files.append(file)
    
    # Load the results
    all_circs_results_sorted=[]
    for file_name in circ_result_property_files:
        #We're looking at one circuit
        # Load data from each pickle and get the results object for the circuit. 
        with open(os.path.join(BASE_PATH, file_name), "rb") as circ_file:
                circ_results=pickle.load(circ_file)["results"]
        
        # for info in circ_infos:
        #     circ_results.append(info["results"])
        
        # # Get the results in order. 
        circ_results_sorted=sorted(circ_results, key=lambda d: d['error_idx'])
        assert len(circ_results_sorted)==len(circ_results), "error when sorting outcomes."
        all_circs_results_sorted.append(circ_results_sorted)

    assert len(all_circs_results_sorted)==len(circ_result_property_files), "lost some files from sorting."

    results=[]
    for circ in all_circs_results_sorted:
        circ_outcomes_temp=[]
        # print(len(circ))
        for idx, rnd in enumerate(circ):
            # print(idx)
            circ_outcomes_temp.append(rnd["state_fidelity_with_checks_with_errors"]-rnd["state_fidelity_no_checks_with_errors"])
        results.append(circ_outcomes_temp)
    
    results=[sum(x)/len(x) for x in zip(*results)]

    plt.plot(ONE_QUBIT_ERROR_SPACE, results)
    plt.title("Qubits: "+str(NUMBER_OF_QUBITS)+ " CNOTS: " + str(CNOT_COUNT))
    plt.xlabel("1 Qubit Error")
    plt.ylabel("Average Fidelity Gain")
    plt.xscale("log")
    plt.savefig(BASE_PATH+"Qubits"+str(NUMBER_OF_QUBITS)+ "_CNOTS_" + str(CNOT_COUNT)+".png")
    plt.show()

                
    
    # rand_circ_files=[elem for elem in glob.glob(code_dir+SUBDIR+"qubits"+str(NUMBER_OF_QUBITS)+"*result*.txt")]
    # cnots1to5=[] #idx 0
    # cnots6to10=[] #idx 1
    # cnots11to15=[] #idx 2
    # cnots16to20=[] #idx 3
    # cnots21to25=[] #idx 4
    # cnots26to30=[] #idx 5
    # cnots31to35=[] #idx 6
    # cnots36to40=[] #idx 7
    # cnots41to45=[] #idx 8
    # grouped_files=[cnots1to5,cnots6to10,cnots11to15,cnots16to20,cnots21to25,cnots26to30,cnots31to35,cnots36to40,cnots41to45]
    # for file in rand_circ_files:
    #     lines=get_lines(file)
    #     gate_count=get_line_str_stripped(lines , "{").split(",")
    #     cnot_count=get_cnot_count(gate_count)
    #     if 1<=cnot_count<=5: 
    #         cnots1to5.append(file)
    #     elif 6<=cnot_count<=10:
    #         cnots6to10.append(file)
    #     elif 11<=cnot_count<=15:
    #         cnots11to15.append(file)
    #     elif 16<=cnot_count<=20:
    #         cnots16to20.append(file)
    #     elif 21<=cnot_count<=25:
    #         cnots21to25.append(file)
    #     elif 26<=cnot_count<=30:
    #         cnots26to30.append(file)
    #     elif 31<=cnot_count<=35:
    #         cnots31to35.append(file)
    #     elif 36<=cnot_count<=40:
    #         cnots36to40.append(file)
    #     elif 41<=cnot_count<=45:
    #         cnots41to45.append(file)
    #     elif cnot_count>45:
    #         print("cnot count of "+ str(cnot_count) + " wasn't included.")
    # plot_from_grouped_cnots(grouped_files)

def plot_from_grouped_cnots(grouped_files):
    for group_num, temp_list in enumerate(grouped_files):
        plot_avg_fidelity_gain=[]
        plot_avg_fidelity_checks=[]
        plot_avg_fidelity_no_checks=[]
        erroridx_list=[]
        error_points=[]
        num_circs=0
        for erroridx in range (NUMBER_OF_ERROR_POINTS):
            erroridx_str_check="erroridx"+str(erroridx)+"_"
            # Separate the current erroridx files.
            erroridx_files=[elem for elem in temp_list if erroridx_str_check in elem]
            
            if len(erroridx_files)!=0:
                num_circs=len(erroridx_files)
                # The error point has data.
                erroridx_list.append(erroridx)
                # get average gain
                fidelity_checks_temp=get_fidelities(erroridx_files, "State fidelity with checks and errors: ")
                fidelity_no_checks_temp=get_fidelities(erroridx_files, "State fidelity no checks with errors: ")
                if group_num == 1:
                    print(len(erroridx_files))
                    print("error: ", ONE_QUBIT_ERROR_SPACE[erroridx])
                    print(len(fidelity_checks_temp))
                    print(len(fidelity_no_checks_temp))

                assert len(fidelity_checks_temp)==len(fidelity_no_checks_temp), "Number of fidelities do not match."
                plot_avg_fidelity_checks.append(sum(fidelity_checks_temp)/len(fidelity_checks_temp))
                plot_avg_fidelity_no_checks.append(sum(fidelity_no_checks_temp)/len(fidelity_no_checks_temp))
                plot_avg_fidelity_gain.append(plot_avg_fidelity_checks[-1]-plot_avg_fidelity_no_checks[-1])
        for elem in erroridx_list:
            error_points.append(ONE_QUBIT_ERROR_SPACE[elem])
        
        cnot_count=0
        if group_num==0:
            cnot_count="1-5"
        elif group_num==1:
            cnot_count="6-10"
            # print(sorted(fidelity_checks_temp))
            # print(sorted(fidelity_no_checks_temp))
        elif group_num==2:
            cnot_count="11-15"
        elif group_num==3:
            cnot_count="16-20"
        elif group_num==4:
            cnot_count="21-25"
        elif group_num==5:
            cnot_count="26-30"
        elif group_num==6:
            cnot_count="31-35"
        elif group_num==7:
            cnot_count="36-40"
        elif group_num==8:
            cnot_count="41-45"
        
        plt.plot(error_points, plot_avg_fidelity_gain)
        plt.title("Qubits: "+str(NUMBER_OF_QUBITS)+ " CNOTS: " +cnot_count + " Number of Circuits: " + str(num_circs))
        plt.xlabel("1 Qubit Error")
        plt.ylabel("Average Fidelity Gain")
        plt.xscale("log")
        plt.savefig(PATH+"Qubits"+str(NUMBER_OF_QUBITS)+ "_CNOTS" +cnot_count+".png")
        plt.show()



def get_cnot_count(gate_count):
    for item in gate_count:
        if item.startswith("\"cx\": "):
            return int(item[6:])
        else:
            return 0


if __name__ == "__main__":
    # File path stuff
    CODE_DIR=sys.path[0]
    SUBDIR="/data/"
    BASE_PATH=CODE_DIR+SUBDIR
    os.chdir(BASE_PATH)
    NUMBER_OF_QUBITS=5
    CNOT_COUNT=20
    NUMBER_OF_ERROR_POINTS=21
    ONE_QUBIT_ERROR_SPACE=np.logspace(-5, -2, num=NUMBER_OF_ERROR_POINTS)
    # create_cnot_plot()
    # create_fidelity_plot_depth()
    create_fidelity_plot_cnots()
#!/usr/bin/env python3

import os, sys, glob
import matplotlib.pyplot as plt
import operator
import numpy as np

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

# def get_cnot_counts(file_lines):
#     '''Get CNOT counts'''
#     for line in file_lines:
#         if line.startswith("CNOT count"):
#             print(line)
#             return

def get_count(file_lines, str):
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

def create_fidelity_plot():
    global NUMBER_OF_QUBITS, SUBDIR, PATH, ONE_QUBIT_ERROR_SPACE
    # Gets the file path of the script
    code_dir=sys.path[0]
    rand_circ_files=[elem for elem in glob.glob(code_dir+SUBDIR+"qubits"+str(NUMBER_OF_QUBITS)+"*result0.txt")]
    depths=(1,5,10,15,20,25)
    for depth in depths:
        plot_fidelity_gain=[]
        plot_fidelity_checks=[]
        plot_fidelity_no_checks=[]
        erroridx_list=[]
        depth_str_check="depth"+str(depth)
        depth_files=[elem for elem in rand_circ_files if depth_str_check in elem]    
        # fidelity_no_checks=[]
        # fidelity_with_checks=[]
        for erroridx in range (21):
            erroridx_str_check="erroridx"+str(erroridx)
            idx_files=[elem for elem in depth_files if erroridx_str_check in elem]
            if len(idx_files)!=0:
                # print(idx_files)
                #same depth and error so average the gain
                fidelity_gain=[]
                fidelity_checks=[]
                fidelity_no_checks=[]
                average_fidelity_gain=0
                for file in idx_files:
                    lines=get_lines(file)
                    for item in lines:
                        fidelity_checks_temp=0
                        fidelity_no_checks_temp=0
                        if item.startswith("State fidelity with checks and errors: "):
                            fidelity_checks_temp=item.split(": ")[1].strip("\n")
                            fidelity_checks.append(float(fidelity_checks_temp))
                            # if depth==25:
                            #     fidelity_checks_temp=item.split(": ")[1].strip("\n")
                            #     fidelity_checks.append(float(fidelity_checks_temp))
                            #     # print("file: "+ file+ " error: "+ str(ONE_QUBIT_ERROR_SPACE[erroridx])+" State fidelity with checks and errors: "+ fidelity_checks_temp)
                            #     # print(float(fidelity_checks_temp))
                        elif item.startswith("State fidelity no checks with errors: "):
                            fidelity_no_checks_temp=item.split(": ")[1].strip("\n")
                            fidelity_no_checks.append(float(fidelity_no_checks_temp))
                            # if depth==25:
                            #     fidelity_no_checks_temp=item.split(": ")[1].strip("\n")
                            #     fidelity_no_checks.append(float(fidelity_no_checks_temp))
                            #     # print("file: "+ file+ "error: "+ str(ONE_QUBIT_ERROR_SPACE[erroridx])+" State fidelity no checks with errors: "+ fidelity_no_checks_temp)
                            #     # print(float(fidelity_no_checks_temp))
                zip_object=zip(fidelity_checks, fidelity_no_checks)
                for item1, item2 in zip_object:
                    fidelity_gain.append(item1-item2)    
                average_fidelity_gain=sum(fidelity_gain)/len(fidelity_gain)
                print(average_fidelity_gain)
                average_fidelity_checks=sum(fidelity_checks)/len(fidelity_checks)
                average_fidelity_no_checks=sum(fidelity_no_checks)/len(fidelity_no_checks)
                plot_fidelity_gain.append(average_fidelity_gain)
                plot_fidelity_no_checks.append(average_fidelity_no_checks)
                plot_fidelity_checks.append(average_fidelity_checks)
                erroridx_list.append(erroridx)
                # print("Depth:" + str(depth)+ "; Error: " + str(ONE_QUBIT_ERROR_SPACE[erroridx]))
                # print(plot_fidelity_gain)
                # print()
                # print(plot_fidelity_checks)
                # print()
                # print(plot_fidelity_no_checks)
        #Convert the erroridx to error values.
        error_points=[]
        for elem in erroridx_list:
            error_points.append(ONE_QUBIT_ERROR_SPACE[elem])
        

        plt.plot(error_points, plot_fidelity_gain)
        plt.title("Qubits: "+str(NUMBER_OF_QUBITS)+ " Depth: " +str(depth))
        plt.xlabel("1 Qubit Error")
        plt.ylabel("Average Fidelity Gain")
        plt.savefig(PATH+"Qubits"+str(NUMBER_OF_QUBITS)+ "_Depth" +str(depth)+".png")
        plt.show()

if __name__ == "__main__":
    # File path stuff
    CODE_DIR=sys.path[0]
    SUBDIR="/data/"
    PATH=CODE_DIR+SUBDIR
    os.chdir(PATH)
    NUMBER_OF_QUBITS=5
    ONE_QUBIT_ERROR_SPACE=np.logspace(-5, -2, num=21)
    # create_cnot_plot()
    create_fidelity_plot()
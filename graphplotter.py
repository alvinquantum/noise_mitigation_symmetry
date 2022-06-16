#!/usr/bin/env python3

import os, sys, pickle
import matplotlib.pyplot as plt
import numpy as np
from os import listdir
from os.path import isfile
import seaborn as sns
import matplotlib.transforms as mtransforms
from matplotlib.pyplot import figure

def get_all_sorted_circ_results(base_path, cnot_count, number_of_qubits, result_str="result"):
    '''Gets the sorted (by error idx, i.e., smallest single qubit error to largest) test results for 
    the given cnot count and number of qubits.'''
    all_files=[f for f in listdir(base_path) if isfile(os.path.join(base_path, f))]
    # print(all_files)
    target_files=[]
    for file in all_files:
        name_split=file.split("_")
        # Get the digits in the names for checking
        name_split_nums=[int(num) for num in name_split if num.isdigit()]
        if (
            ".obj" in name_split and result_str in name_split 
            and name_split_nums[1]==cnot_count 
            and name_split_nums[0]==number_of_qubits): #and start_circ_number<=name_split_nums[2]<=end_circ_number:
            target_files.append(file)

    all_circ_results_sorted=[]
    if target_files:
        # Load the results
        for file_name in target_files:
            #We're looking at one circuit
            # Load data from each pickle and get the results object for the circuit. 
            with open(os.path.join(base_path, file_name), "rb") as circ_file:
                    pickle_info=pickle.load(circ_file)
                    circ_result=pickle_info["results"]
                    # assert pickle_info["qubits"]==number_of_qubits, f"{file_name} does not have the required qubit count."
                    # assert pickle_info["cx"]==cnot_count, f"{file_name} does not have the required CNOT count."
            
            # Get the results in order. 
            circ_result_sorted=sorted(circ_result, key=lambda d: d['error_idx'])
            assert len(circ_result_sorted)==len(circ_result), "error when sorting outcomes."
            all_circ_results_sorted.append(circ_result_sorted)

        assert len(all_circ_results_sorted)==len(target_files), "lost some files from sorting."
    return all_circ_results_sorted

def get_all_sorted_circ_results_for_layer(base_path, cnot_count, number_of_qubits, layers_count, result_str="result"):
    '''Gets the sorted (by error idx, i.e., smallest single qubit error to largest) test results for 
    the given cnot count and number of qubits.'''
    all_files=[f for f in listdir(base_path) if isfile(os.path.join(base_path, f))]
    # print(all_files)
    target_files=[]
    for file in all_files:
        name_split=file.split("_")
        # Get the digits in the names for checking
        name_split_nums=[int(num) for num in name_split if num.isdigit()]
        if (f"layers_{layers_count}_" in file and
            ".obj" in name_split and result_str in name_split 
            and name_split_nums[1]==cnot_count 
            and name_split_nums[0]==number_of_qubits): #and start_circ_number<=name_split_nums[2]<=end_circ_number:
            target_files.append(file)

    all_circ_results_sorted=[]
    if target_files:
        # print(target_files)
        # Load the results
        for file_name in target_files:
            #We're looking at one circuit
            # Load data from each pickle and get the results object for the circuit. 
            with open(os.path.join(base_path, file_name), "rb") as circ_file:
                pickle_info=pickle.load(circ_file)
                circ_result=pickle_info["results"]
                # assert pickle_info["qubits"]==number_of_qubits, f"{file_name} does not have the required qubit count."
                # assert pickle_info["cx"]==cnot_count, f"{file_name} does not have the required CNOT count."
            
            # Get the results in order. 
            circ_result_sorted=sorted(circ_result, key=lambda d: d['error_idx'])
            assert len(circ_result_sorted)==len(circ_result), "error when sorting outcomes."
            all_circ_results_sorted.append(circ_result_sorted)

        assert len(all_circ_results_sorted)==len(target_files), "lost some files from sorting."
    return all_circ_results_sorted

def create_fidelity_plot_cnot(base_path, number_of_qubits, cnot_count, one_qubit_error_space):
    '''Creates plots based on CNOT counts.'''
    # Gets the file path of the script.
    all_sorted_circ_results=get_all_sorted_circ_results(base_path, cnot_count, number_of_qubits)
    if all_sorted_circ_results:    
        average_fidelity_gains=[]
        fidelity_gains=[]
        for circ_result in all_sorted_circ_results:
            circ_outcomes_temp=[]

            for result in circ_result:
                circ_outcomes_temp.append(result["state_fidelity_with_checks_with_errors"]-result["state_fidelity_no_checks_with_errors"])
            fidelity_gains.append(circ_outcomes_temp)
        
        # print("raw fidelity", fidelity_gains)
        average_fidelity_gains=[sum(x)/len(x) for x in zip(*fidelity_gains)]
        # print(average_fidelity_gains)
        std_deviations=[np.std(x) for x in zip(*fidelity_gains)]
        # print(std_deviations)

        plt.clf()
        plt.errorbar(one_qubit_error_space, average_fidelity_gains, yerr=std_deviations, capsize=3.0)
        plt.title("Qubits: "+str(number_of_qubits)+ " CNOTS: " + str(cnot_count))
        plt.xlabel("Single Qubit Error")
        plt.ylabel("Average Fidelity Gain")
        plt.xscale("log")
        # plt.savefig(os.path.join(base_path, f"qubits_{number_of_qubits}_CNOTS_{cnot_count}_fidelity_.png"))
        return (cnot_count, average_fidelity_gains, std_deviations)

def create_fidelity_plot_cnot_layers(base_path, number_of_qubits, cnot_count, one_qubit_error_space, layers_count, plot_result):
    '''Creates plots based on CNOT counts and layers.'''
    # Gets the file path of the script.
    all_sorted_circ_results=get_all_sorted_circ_results_for_layer(base_path, cnot_count, number_of_qubits, layers_count)
    if all_sorted_circ_results:
        average_fidelity_gains=[]
        fidelity_gains=[]
        # Looking at one circuit.
        for circ_result in all_sorted_circ_results:
            circ_outcomes_temp=[]
            # Looking at one error idx for this specific circuit.
            for result in circ_result:
                circ_outcomes_temp.append(result["state_fidelity_with_checks_with_errors"]-result["state_fidelity_no_checks_with_errors"])
            fidelity_gains.append(circ_outcomes_temp)
        
        # print("raw fidelity", fidelity_gains)
        average_fidelity_gains=[sum(x)/len(x) for x in zip(*fidelity_gains)]
        # print("average gains ", average_fidelity_gains)
        # print(len(average_fidelity_gains))
        std_deviations=[np.std(x) for x in zip(*fidelity_gains)]
        # print(len(std_deviations))
        # print(len(one_qubit_error_space))
        # print(std_deviations)

        if plot_result:
            plt.clf()
            plt.errorbar(one_qubit_error_space, average_fidelity_gains, yerr=std_deviations, capsize=3.0)
            plt.title(f"Compute Qubits: {str(number_of_qubits)}, CNOTS: {cnot_count}, Layers: {layers_count}")
            plt.xlabel("Single Qubit Error")
            plt.ylabel("Average Fidelity Gain")
            plt.xscale("log")
            plt.savefig(os.path.join(base_path, f"qubits_{number_of_qubits}_CNOTS_{cnot_count}_fidelity_layers_{layers_count}_.png"))
        return (cnot_count, average_fidelity_gains, std_deviations)

def create_all_fidelity_plot_cnots_layers_rz(base_path, number_of_qubits, cnot_counts, one_qubit_error_space, palette, layers_count, rz_count, ax):
    # plt.clf()
    all_results=[]
    for cnot_count in cnot_counts:
        result=create_fidelity_plot_cnot_layers(base_path, number_of_qubits, cnot_count, one_qubit_error_space, layers_count, plot_result= False)
        if result:
            all_results.append(result)


    for idx, elem in enumerate(all_results):
        # print(f"fidelity gain for cnots {cnot_counts[idx]}: {elem[1]}")
        ax.errorbar(one_qubit_error_space, elem[1], yerr=elem[2], capsize=3.0, color=palette[idx*5], label="CNOT Count: "+str(elem[0]))
    ax.set_title(f"Compute Qubits: {str(number_of_qubits)}, Layers: {layers_count}, Rz: {rz_count}", fontsize=12)
    # ax.set_xlabel("Single Qubit Error")
    ax.set_ylabel("Average Fidelity Gain", fontsize=12)
    ax.set_xscale("log")
    # ax.legend()    
    # ax.savefig(os.path.join(base_path, f"qubits_{number_of_qubits}_allCNOTS_fidelity_gain_layers_{layers_count}_.png"))

def create_2qubit_plot_example_fidelity(base_path, number_of_qubits, cnot_count, one_qubit_error_space, layers, ax):
    '''Creates plots based on CNOT counts.'''
    NLINES = len(layers)+1 #plus one for the zero checks.
    palette = sns.color_palette("viridis",NLINES)
    for layer_idx, layer in enumerate(layers):
        all_sorted_circ_results=get_all_sorted_circ_results_for_layer(base_path, cnot_count, number_of_qubits, layer)

        if all_sorted_circ_results:
            # Looking at one circuit.
            for circ_result in all_sorted_circ_results:

                # Looking at one error value for that circuit.
                temp_fidelity_no_checks=[]
                temp_fidelity_with_checks=[]
                for result in circ_result:
                    temp_fidelity_no_checks.append(result["state_fidelity_no_checks_with_errors"])
                    temp_fidelity_with_checks.append(result["state_fidelity_with_checks_with_errors"])

                if layer_idx==0:
                    ax.errorbar(one_qubit_error_space, temp_fidelity_no_checks, capsize=3.0, color=palette[0], label="Number of layers: 0")
                # Each temp_fidelity_with_checks is a complete circuit result. So plot each separately.
                ax.errorbar(one_qubit_error_space, temp_fidelity_with_checks, capsize=3.0, color=palette[layer_idx+1], label=f"Number of layers: {layer}")
    ax.set_title(f"Compute Qubits: {number_of_qubits}, CNOTS: {cnot_count}", fontsize=12)
    ax.set_ylabel("Fidelity", fontsize=12)
    ax.set_xscale("log")

def create_count_percentage_plot(base_path, number_of_qubits, cnot_count, one_qubit_error_space):
    '''Plots the average zero outcome for ancilla vs the single qubit error rate.'''
    '''Creates plots based on CNOT counts.'''
    # Gets the file path of the script.
    all_sorted_circ_results=get_all_sorted_circ_results(base_path, cnot_count, number_of_qubits)

    if all_sorted_circ_results:
        count_diffs=[]
        for circ_result in all_sorted_circ_results:
            circ_outcomes_temp=[]
            for result in circ_result:
                # In our simulations we have the full density matrix so the percent results before post selecting is actually 1.
                # We do it this way so the graphing can accomidate physical runs where these two values are actually counts.
                circ_outcomes_temp.append(result["percent_results_after_postselect"]/result["percent_results_before_postselect"])
            count_diffs.append(circ_outcomes_temp)
        
        average_counts_diffs=[sum(x)/len(x) for x in zip(*count_diffs)]
        # print(average_counts_diffs)
        std_deviations=[np.std(x) for x in zip(*count_diffs)]
        # print(std_deviations)

        plt.clf()
        plt.errorbar(one_qubit_error_space, average_counts_diffs, yerr=std_deviations, capsize=3.0)
        plt.title(f"Compute Qubits: {number_of_qubits}, CNOTS: {cnot_count}")
        plt.xlabel("Single Qubit Error")
        plt.ylabel("Average Probability of Zero Measurement on Ancilla")
        plt.xscale("log")
        plt.savefig(os.path.join(base_path, f"qubits_{number_of_qubits}_CNOTS_{cnot_count}_counts_.png"))
        return (cnot_count, average_counts_diffs, std_deviations)

def create_count_percentage_plot_layers(base_path, number_of_qubits, cnot_count, one_qubit_error_space, layers_count, plot_result):
    '''Plots the average zero outcome for ancilla vs the single qubit error rate.'''
    '''Creates plots based on CNOT counts.'''
    # Gets the file path of the script.
    all_sorted_circ_results=get_all_sorted_circ_results_for_layer(base_path, cnot_count, number_of_qubits, layers_count = layers_count)

    if all_sorted_circ_results:
        count_diffs=[]
        for circ_result in all_sorted_circ_results:
            circ_outcomes_temp=[]
            for result in circ_result:
                # In our simulations we have the full density matrix so the percent results before post selecting is actually 1.
                # We do it this way so the graphing can accomidate physical runs where these two values are actually counts.
                circ_outcomes_temp.append(result["percent_results_after_postselect"]/result["percent_results_before_postselect"])
            count_diffs.append(circ_outcomes_temp)
        
        average_counts_diffs=[sum(x)/len(x) for x in zip(*count_diffs)]
        # print(average_counts_diffs)
        std_deviations=[np.std(x) for x in zip(*count_diffs)]
        # print(std_deviations)

        if plot_result:
            plt.clf()
            plt.errorbar(one_qubit_error_space, average_counts_diffs, yerr=std_deviations, capsize=3.0)
            plt.title(f"Compute Qubits: {number_of_qubits}, CNOTS: {cnot_count}, Layers: {layers_count}")
            plt.xlabel("Single Qubit Error")
            plt.ylabel("Average Probability of Zero Measurement on Ancilla")
            plt.xscale("log")
            plt.savefig(os.path.join(base_path, f"qubits_{number_of_qubits}_CNOTS_{cnot_count}_counts_layers_{layers_count}_.png"))
        return (cnot_count, average_counts_diffs, std_deviations)

def create_2qubit_plot_example_count_percentage(base_path, number_of_qubits, cnot_count, one_qubit_error_space, layers, ax):
    '''Plots the average zero outcome for ancilla vs the single qubit error rate.'''
    '''Creates plots based on CNOT counts.'''
    NLINES = len(layers)+1 #plus one for the zero checks.
    palette = sns.color_palette("viridis", NLINES)
    for layer_idx, layer in enumerate(layers):
        all_sorted_circ_results=get_all_sorted_circ_results_for_layer(base_path, cnot_count, number_of_qubits, layer)
        if all_sorted_circ_results:
            # each element is a circuit.
            for circ_result in all_sorted_circ_results:
                counts_after_temp=[]
                counts_before_temp=[]
                # each element here is one error idx for the circuit.
                for result in circ_result:
                    # In our simulations we have the full density matrix so the percent results before post selecting is actually 1.
                    # We do it this way so the graphing can accomidate physical runs where these two values are actually counts.
                    counts_after_temp.append(result["percent_results_after_postselect"])
                    counts_before_temp.append(result["percent_results_before_postselect"])

                if layer_idx==0:
                    ax.errorbar(one_qubit_error_space, counts_before_temp, color=palette[0], capsize=3.0, label="Number of layers: 0")
                ax.errorbar(one_qubit_error_space, counts_after_temp,
                    capsize=3.0, color=palette[layer_idx+1], label=f"Number of layers: {layer}")
        ax.set_xlabel("Single Qubit Error", fontsize=12)
        ax.set_ylabel("Probability of Zero\n Measurement on All Ancillas", fontsize=12)
        ax.set_xscale("log")
        ax.legend()


def create_all_count_percentage_plots_layers(base_path, number_of_qubits, cnot_counts, one_qubit_error_space, palette, layers_count, ax):
    '''Wrapper for create_count_percentage_plot(). Makes individual plots with all the lines and a single plot'''
    # plt.clf()
    all_results=[]
    for cnot_count in cnot_counts:
        result=create_count_percentage_plot_layers(base_path, number_of_qubits, cnot_count, one_qubit_error_space, layers_count=layers_count, plot_result=False)
        if result:
            all_results.append(result)

    for idx, elem in enumerate(all_results):
        ax.errorbar(one_qubit_error_space, elem[1], yerr=elem[2], capsize=3.0, color=palette[idx*5], label="CNOT Count: "+str(elem[0]))
    # ax.set_title(f"Compute Qubits: {str(number_of_qubits)}, Layers: {layers_count}")
    ax.set_xlabel("Single Qubit Error", fontsize=12)
    ax.set_ylabel("Average Probability of Zero\n Measurement on All Ancillas", fontsize=12)
    ax.set_xscale("log")
    ax.legend()    
    # ax.savefig(os.path.join(base_path, f"qubits_{number_of_qubits}_allCNOTS_counts_layers_{layers_count}_.png"))

def create_all_count_percentage_plots(base_path, number_of_qubits, cnot_counts, one_qubit_error_space, palette, NLINES):
    '''Wrapper for create_count_percentage_plot(). Makes individual plots with all the lines and a single plot'''
    plt.clf()
    all_results=[]
    for cnot_count in cnot_counts:
        result=create_count_percentage_plot(base_path, number_of_qubits, cnot_count, one_qubit_error_space)
        if result:
            all_results.append(result)

    for idx, elem in enumerate(all_results):
        plt.errorbar(one_qubit_error_space, elem[1], yerr=elem[2], capsize=3.0, color=palette[idx*5], label="CNOT Count: "+str(elem[0]))
    plt.title(f"Compute Qubits: {number_of_qubits}, Layers: 1")
    plt.xlabel("Single Qubit Error")
    plt.ylabel("Average Probability of Zero Measurement on Ancilla")
    plt.xscale("log")
    plt.legend()    
    plt.savefig(os.path.join(base_path, f"qubits_{number_of_qubits}_allCNOTS_counts_.png"))

def create_split_layers_fidelity_plot_cnots_subplots(base_path, number_of_qubits, cnot_count, one_qubit_error_space, layers, ax, rz_count):
    '''Creates plots based on CNOT counts.'''
    # Gets the file path of the script.
    NLINES = len(layers)+1 #plus one for the zero checks.
    # NLINES = 45
    palette = sns.color_palette("viridis",NLINES)
    for layer_idx, layer in enumerate(layers):
        all_sorted_circ_results=get_all_sorted_circ_results_for_layer(base_path, cnot_count, number_of_qubits, layer)
        # print(all_sorted_circ_results)
        if all_sorted_circ_results:
            #Regular plots
            # average_fidelity_gains=[]
            # fidelity_gains=[]

            # Each element contains results for one error value
            fidelity_no_checks=[]
            fidelity_with_checks=[]
            # Looking at one circuit.
            for circ_result in all_sorted_circ_results:
                # circ_outcomes_temp=[]

                # Looking at one error value for that circuit.
                temp_fidelity_no_checks=[]
                temp_fidelity_with_checks=[]
                for result in circ_result:
                    temp_fidelity_no_checks.append(result["state_fidelity_no_checks_with_errors"])
                    temp_fidelity_with_checks.append(result["state_fidelity_with_checks_with_errors"])
                #     circ_outcomes_temp.append(result["state_fidelity_with_checks_with_errors"]-result["state_fidelity_no_checks_with_errors"])
                # fidelity_gains.append(circ_outcomes_temp)
                fidelity_no_checks.append(temp_fidelity_no_checks)
                fidelity_with_checks.append(temp_fidelity_with_checks)
            average_fidelity_with_checks=[sum(x)/len(x) for x in zip(*fidelity_with_checks)]
            # if cnot_count==40:
                # print(f"layer: {layer}, fidelity with check: {average_fidelity_with_checks}")
            # print("raw fidelity", fidelity_gains)
            # average_fidelity_with_checks=[sum(x)/len(x) for x in zip(*fidelity_with_checks)]
            # print(f"average fidelity with checks {average_fidelity_with_checks}")
            std_deviations_with_checks=[np.std(x) for x in zip(*fidelity_with_checks)]
            # print(f"std deviation with checks {std_deviations_with_checks}")

            average_fidelity_no_checks=[sum(x)/len(x) for x in zip(*fidelity_no_checks)]
            # if cnot_count==40:
                # print(f"layer: {layer}, fidelity with no check: {average_fidelity_no_checks}")

            # print(f"average fidelity no checks {average_fidelity_no_checks}")
            std_deviations_no_checks=[np.std(x) for x in zip(*fidelity_no_checks)]
            # print(f"std deviation no checks {std_deviations_no_checks}")
            if layer_idx==0:
                ax.errorbar(one_qubit_error_space, average_fidelity_no_checks, std_deviations_no_checks, capsize=3.0, color=palette[0], label="Number of layers: 0")
            # # Each element of fidelity_with_checks is a complete circuit result. So plot each separately.
            # for circ_idx, circ_result in enumerate(fidelity_with_checks):
            #     layer_count=circ_idx+1
            ax.errorbar(one_qubit_error_space, average_fidelity_with_checks, std_deviations_with_checks, capsize=3.0, color=palette[layer_idx+1], label=f"Number of layers: {layer}")
    ax.set_title(f"Compute Qubits: {number_of_qubits}, CNOTS: {cnot_count}, Rz: {rz_count}", fontsize=12)
    ax.set_ylabel("Average Fidelity", fontsize=12)
    ax.set_xscale("log")
    # ax.legend()
    # return (cnot_count, average_fidelity_gains, std_deviations)

def create_fidelity_vs_cnots_plots_filtered_zero_prob(base_path, number_of_qubits, cnot_counts, one_qubit_err, one_qubit_err_idx, layers, ax):
    # plt.clf()
    nlines = len(layers)+1 #plus one for the zero checks.
    palette = sns.color_palette("viridis", nlines)
    # Get the results.
    for layer_idx, layer in enumerate(layers):
        average_fidelity_forlayers=[]
        stddev_of_avg_fidelities_forlayer=[]
        #No mitigation stuff.
        average_fidelity_fornolayers=[]
        stddev_of_avg_fidelities_fornolayers=[]
        for cnot_idx, cnot_count in enumerate(cnot_counts):
            # looking at current layer and cnot count.
            # Each element in contains one circuit.
            all_sorted_circ_results=get_all_sorted_circ_results_for_layer(base_path, cnot_count, number_of_qubits, layer)
            #Each element contains the result of one circuit.
            temp_fidelities=[result[one_qubit_err_idx]["state_fidelity_with_checks_with_errors"] for result in all_sorted_circ_results
                        if result[0]["percent_results_after_postselect"]>0.001] 
            #Get the zero layers stuff.
            if layer_idx==0:
                 #Each element contains the result of one circuit.
                temp_fidelities_nolayers=[result[one_qubit_err_idx]["state_fidelity_no_checks_with_errors"] for result in all_sorted_circ_results
                        if result[0]["percent_results_after_postselect"]>0.001] 
                # Average over the results and append.
                if temp_fidelities_nolayers:
                    average_fidelity_fornolayers.append(sum(temp_fidelities_nolayers)/len(temp_fidelities_nolayers))
                    stddev_of_avg_fidelities_fornolayers.append(np.std(temp_fidelities_nolayers))


            #Average over the results and append.
            if temp_fidelities:
                average_fidelity_forlayers.append(sum(temp_fidelities)/len(temp_fidelities))
                stddev_of_avg_fidelities_forlayer.append(np.std(temp_fidelities))

        if layer_idx==0 and average_fidelity_fornolayers:
            ax.errorbar(cnot_counts, average_fidelity_fornolayers, stddev_of_avg_fidelities_fornolayers, capsize=3.0, color=palette[layer_idx], label=f"Layers: {0}")
        if average_fidelity_forlayers:
            # layers_plot=list(range(1, len(average_fidelity_gains_forlayer)+1, 1))
            ax.errorbar(cnot_counts, average_fidelity_forlayers, stddev_of_avg_fidelities_forlayer, capsize=3.0, color=palette[layer_idx+1], label=f"Layers: {layer}")
        # if layer_idx==0:#len(layers)-1:
            # print(f"layers {layer_idx}: has avg fidelities of: {average_fidelity_fornolayers}")
            # # Each element of fidelity_with_checks is a complete circuit result. So plot each separately.
            # for circ_idx, circ_result in enumerate(fidelity_with_checks):
            #     layer_count=circ_idx+1
            # ax.errorbar(one_qubit_error_space, average_fidelity_with_checks, std_deviations_with_checks, capsize=3.0, color=palette[layer_idx+1], label=f"Number of layers: {layer}")

    ax.set_title(f"Compute Qubits: {number_of_qubits}, Single Qubit Error: {round(one_qubit_err, 5)}", fontsize=12)
    ax.set_ylabel("Average Fidelity", fontsize=12)
    # ax.set_xscale("log")
    # ax.legend()

def create_counts_vs_cnots_plots_filtered_zero_prob(base_path, number_of_qubits, cnot_counts, one_qubit_err, one_qubit_err_idx, layers, ax):
     # plt.clf()
    nlines = len(layers)+1 #plus one for the zero checks.
    palette = sns.color_palette("viridis", nlines)
    # Get the results.
    for layer_idx, layer in enumerate(layers):
        average_counts_forlayer=[]
        stddev_of_avg_counts_forlayer=[]
        #No mitigation stuff.
        average_counts_fornolayers=[]
        stddev_of_avg_counts_fornolayers=[]
        for cnot_idx, cnot_count in enumerate(cnot_counts):
            # looking at current layer and cnot count.
            # Each element in contains one circuit.
            all_sorted_circ_results=get_all_sorted_circ_results_for_layer(base_path, cnot_count, number_of_qubits, layer)
            #Each element contains the result of one circuit.
            temp_fidelities=[result[one_qubit_err_idx]["percent_results_after_postselect"] for result in all_sorted_circ_results
                    if result[0]["percent_results_after_postselect"]>0.001] 
            
            #Get the zero layers stuff.
            if layer_idx==0:
                 #Each element contains the result of one circuit.
                temp_counts_nolayers=[result[one_qubit_err_idx]["percent_results_before_postselect"] for result in all_sorted_circ_results
                        if result[0]["percent_results_after_postselect"]>0.001]
                # Average over the results and append.
                if temp_counts_nolayers:
                    average_counts_fornolayers.append(sum(temp_counts_nolayers)/len(temp_counts_nolayers))
                    stddev_of_avg_counts_fornolayers.append(np.std(temp_counts_nolayers))

            if temp_fidelities:
                average_counts_forlayer.append(sum(temp_fidelities)/len(temp_fidelities))
                stddev_of_avg_counts_forlayer.append(np.std(temp_fidelities))

        if layer_idx==0:
            ax.errorbar(cnot_counts, average_counts_fornolayers, stddev_of_avg_counts_fornolayers, capsize=3.0, color=palette[layer_idx], label=f"Layers: {0}")
        if average_counts_forlayer:
            # layers_plot=list(range(1, len(average_counts_forlayer)+1, 1))
            ax.errorbar(cnot_counts, average_counts_forlayer, stddev_of_avg_counts_forlayer, capsize=3.0, color=palette[layer_idx+1], label=f"Layers: {layer}")

        # if layer_idx==len(layers)-1:
            # print(f"layers {layer_idx+1}: has avg post selection prob of: {average_counts_forlayer}")
            # # Each element of fidelity_with_checks is a complete circuit result. So plot each separately.
            # for circ_idx, circ_result in enumerate(fidelity_with_checks):
            #     layer_count=circ_idx+1
            # ax.errorbar(one_qubit_error_space, average_fidelity_with_checks, std_deviations_with_checks, capsize=3.0, color=palette[layer_idx+1], label=f"Number of layers: {layer}")

    # ax.set_title(f"Compute Qubits: {number_of_qubits}, Error: {round(one_qubit_err, 5)}", fontsize=12)
    ax.set_ylabel("Average Probability of Zero\n Measurement on All Ancillas", fontsize=12)
    ax.set_xlabel("Number of CNOTs", fontsize=12)
    # ax.set_xscale("log")
    ax.legend()

def create_fidelity_vs_cnots_plots(base_path, number_of_qubits, cnot_counts, one_qubit_err, one_qubit_err_idx, layers, ax):
    # plt.clf()
    nlines = len(layers)+1 #plus one for the zero checks.
    palette = sns.color_palette("viridis", nlines)
    # Get the results.
    for layer_idx, layer in enumerate(layers):
        average_fidelity_forlayers=[]
        stddev_of_avg_fidelities_forlayer=[]
        #No mitigation stuff.
        average_fidelity_fornolayers=[]
        stddev_of_avg_fidelities_fornolayers=[]
        for cnot_idx, cnot_count in enumerate(cnot_counts):
            # looking at current layer and cnot count.
            # Each element in contains one circuit.
            all_sorted_circ_results=get_all_sorted_circ_results_for_layer(base_path, cnot_count, number_of_qubits, layer)
            #Each element contains the result of one circuit.
            # temp_fidelities=[result[one_qubit_err_idx]["state_fidelity_with_checks_with_errors"] for result in all_sorted_circ_results] 
            temp_fidelities=[result[0]["state_fidelity_with_checks_with_errors"] for result in all_sorted_circ_results] 

            #Get the zero layers stuff.
            if layer_idx==0:
                 #Each element contains the result of one circuit.
                temp_fidelities_nolayers=[result[one_qubit_err_idx]["state_fidelity_no_checks_with_errors"] for result in all_sorted_circ_results] 
                temp_fidelities_nolayers=[result[one_qubit_err_idx]["state_fidelity_no_checks_with_errors"] for result in all_sorted_circ_results] 
                # Average over the results and append.
                if temp_fidelities_nolayers:
                    average_fidelity_fornolayers.append(sum(temp_fidelities_nolayers)/len(temp_fidelities_nolayers))
                    stddev_of_avg_fidelities_fornolayers.append(np.std(temp_fidelities_nolayers))


            #Average over the results and append.
            if temp_fidelities:
                average_fidelity_forlayers.append(sum(temp_fidelities)/len(temp_fidelities))
                stddev_of_avg_fidelities_forlayer.append(np.std(temp_fidelities))

        if layer_idx==0 and average_fidelity_fornolayers:
            ax.errorbar(cnot_counts, average_fidelity_fornolayers, stddev_of_avg_fidelities_fornolayers, capsize=3.0, color=palette[layer_idx], label=f"Layers: {0}")
        if average_fidelity_forlayers:
            # layers_plot=list(range(1, len(average_fidelity_gains_forlayer)+1, 1))
            ax.errorbar(cnot_counts, average_fidelity_forlayers, stddev_of_avg_fidelities_forlayer, capsize=3.0, color=palette[layer_idx+1], label=f"Layers: {layer}")
        # if layer_idx==0:#len(layers)-1:
            # print(f"layers {layer_idx}: has avg fidelities of: {average_fidelity_fornolayers}")
            # # Each element of fidelity_with_checks is a complete circuit result. So plot each separately.
            # for circ_idx, circ_result in enumerate(fidelity_with_checks):
            #     layer_count=circ_idx+1
            # ax.errorbar(one_qubit_error_space, average_fidelity_with_checks, std_deviations_with_checks, capsize=3.0, color=palette[layer_idx+1], label=f"Number of layers: {layer}")

    ax.set_title(f"Compute Qubits: {number_of_qubits}, Single Qubit Error: {round(one_qubit_err, 5)}", fontsize=12)
    ax.set_ylabel("Average Fidelity", fontsize=12)
    # ax.set_xscale("log")
    # ax.legend()

def create_counts_vs_cnots_plots(base_path, number_of_qubits, cnot_counts, one_qubit_err, one_qubit_err_idx, layers, ax):
     # plt.clf()
    nlines = len(layers)+1 #plus one for the zero checks.
    palette = sns.color_palette("viridis", nlines)
    # Get the results.
    for layer_idx, layer in enumerate(layers):
        average_counts_forlayer=[]
        stddev_of_avg_counts_forlayer=[]
        #No mitigation stuff.
        average_counts_fornolayers=[]
        stddev_of_avg_counts_fornolayers=[]
        for cnot_idx, cnot_count in enumerate(cnot_counts):
            # looking at current layer and cnot count.
            # Each element in contains one circuit.
            all_sorted_circ_results=get_all_sorted_circ_results_for_layer(base_path, cnot_count, number_of_qubits, layer)
            #Each element contains the result of one circuit.
            temp_fidelities=[result[one_qubit_err_idx]["percent_results_after_postselect"] for result in all_sorted_circ_results] 
            
            #Get the zero layers stuff.
            if layer_idx==0:
                 #Each element contains the result of one circuit.
                temp_counts_nolayers=[result[one_qubit_err_idx]["percent_results_before_postselect"] for result in all_sorted_circ_results]
                # Average over the results and append.
                if temp_counts_nolayers:
                    average_counts_fornolayers.append(sum(temp_counts_nolayers)/len(temp_counts_nolayers))
                    stddev_of_avg_counts_fornolayers.append(np.std(temp_counts_nolayers))

            if temp_fidelities:
                average_counts_forlayer.append(sum(temp_fidelities)/len(temp_fidelities))
                stddev_of_avg_counts_forlayer.append(np.std(temp_fidelities))

        if layer_idx==0:
            ax.errorbar(cnot_counts, average_counts_fornolayers, stddev_of_avg_counts_fornolayers, capsize=3.0, color=palette[layer_idx], label=f"Layers: {0}")
        if average_counts_forlayer:
            # layers_plot=list(range(1, len(average_counts_forlayer)+1, 1))
            ax.errorbar(cnot_counts, average_counts_forlayer, stddev_of_avg_counts_forlayer, capsize=3.0, color=palette[layer_idx+1], label=f"Layers: {layer}")

        # if layer_idx==len(layers)-1:
            # print(f"layers {layer_idx+1}: has avg post selection prob of: {average_counts_forlayer}")
            # # Each element of fidelity_with_checks is a complete circuit result. So plot each separately.
            # for circ_idx, circ_result in enumerate(fidelity_with_checks):
            #     layer_count=circ_idx+1
            # ax.errorbar(one_qubit_error_space, average_fidelity_with_checks, std_deviations_with_checks, capsize=3.0, color=palette[layer_idx+1], label=f"Number of layers: {layer}")

    # ax.set_title(f"Compute Qubits: {number_of_qubits}, Error: {round(one_qubit_err, 5)}", fontsize=12)
    ax.set_ylabel("Average Probability of Zero\n Measurement on All Ancillas", fontsize=12)
    ax.set_xlabel("Number of CNOTs", fontsize=12)
    # ax.set_xscale("log")
    ax.legend()

def create_allcnot_fidelity_vs_layers_plots(base_path, number_of_qubits, cnot_counts, one_qubit_err, one_qubit_err_idx, layers, ax, rz_count):
    # plt.clf()
    nlines = len(cnot_counts)+1 #plus one for the zero checks.
    palette = sns.color_palette("viridis", nlines)
    # Get the results.
    for cnot_idx, cnot_count in enumerate(cnot_counts):
        average_fidelity_gains_forcnot=[]
        stddev_of_avg_fidelities_forcnot=[]
        for layer_idx, layer in enumerate(layers):
            # looking at current layer and cnot count.
            # Each element in contains one circuit.
            all_sorted_circ_results=get_all_sorted_circ_results_for_layer(base_path, cnot_count, number_of_qubits, layer)
            temp_fidelities=[result[one_qubit_err_idx]["state_fidelity_with_checks_with_errors"]-
                result[one_qubit_err_idx]["state_fidelity_no_checks_with_errors"] for result in all_sorted_circ_results] #Each element contains the result of one circuit.
            if layer_idx==0 and all_sorted_circ_results:
                no_layers_fidelity_gain=[result[one_qubit_err_idx]["state_fidelity_no_checks_with_errors"]-
                result[one_qubit_err_idx]["state_fidelity_no_checks_with_errors"] for result in all_sorted_circ_results]
                average_fidelity_gains_forcnot.insert(0,sum(no_layers_fidelity_gain)/len(no_layers_fidelity_gain))
                stddev_of_avg_fidelities_forcnot.insert(0, np.std(no_layers_fidelity_gain))
            # for result in all_sorted_circ_results:
            #     #We only append the result for the specific error value.
            #     temp_fidelities.append(result[one_qubit_err_idx]["state_fidelity_with_checks_with_errors"])
                # print(result[one_qubit_err_idx])
            #Average over the results and append.
            
            if temp_fidelities:
                average_fidelity_gains_forcnot.append(sum(temp_fidelities)/len(temp_fidelities))
                stddev_of_avg_fidelities_forcnot.append(np.std(temp_fidelities))
            # print(len(average_fidelity_gains_forcnot))
            # print(average_fidelity_gains_forcnot)
            # print(stddev_of_avg_fidelities_forcnot)
        if average_fidelity_gains_forcnot:
            layers_plot=list(range(0, len(average_fidelity_gains_forcnot), 1))
            # print(layers_plot)
            ax.errorbar(layers_plot, average_fidelity_gains_forcnot, stddev_of_avg_fidelities_forcnot, capsize=3.0, color=palette[cnot_idx], label=f"CNOT Count: {cnot_count}")
            # # Each element of fidelity_with_checks is a complete circuit result. So plot each separately.
            # for circ_idx, circ_result in enumerate(fidelity_with_checks):
            #     layer_count=circ_idx+1
            # ax.errorbar(one_qubit_error_space, average_fidelity_with_checks, std_deviations_with_checks, capsize=3.0, color=palette[layer_idx+1], label=f"Number of layers: {layer}")

    ax.set_title(f"Compute Qubits: {number_of_qubits}, Rz: {rz_count}, Error: {round(one_qubit_err, 5)}", fontsize=12)
    ax.set_ylabel("Average Fidelity Gain", fontsize=12)
    # ax.set_xscale("log")
    # ax.legend()

def create_allcnot_counts_vs_layers_plots(base_path, number_of_qubits, cnot_counts, one_qubit_err, one_qubit_err_idx, layers, ax):
     # plt.clf()
    nlines = len(cnot_counts)+1 #plus one for the zero checks.
    palette = sns.color_palette("viridis", nlines)
    # Get the results.
    for cnot_idx, cnot_count in enumerate(cnot_counts):
        average_counts_forcnot=[]
        stddev_of_avg_counts_forcnot=[]
        for layer_idx, layer in enumerate(layers):
            # looking at current layer and cnot count.
            # Each element in contains one circuit.
            all_sorted_circ_results=get_all_sorted_circ_results_for_layer(base_path, cnot_count, number_of_qubits, layer)
            temp_fidelities=[result[one_qubit_err_idx]["percent_results_after_postselect"] for result in all_sorted_circ_results] #Each element contains the result of one circuit.
            # for result in all_sorted_circ_results:
            #     #We only append the result for the specific error value.
            #     temp_fidelities.append(result[one_qubit_err_idx]["state_fidelity_with_checks_with_errors"])
                # print(result[one_qubit_err_idx])
            #Average over the results and append.
            if layer_idx==0 and all_sorted_circ_results:
                no_layers_postselect=[result[one_qubit_err_idx]["percent_results_before_postselect"] for result in all_sorted_circ_results]
                average_counts_forcnot.insert(0,sum(no_layers_postselect)/len(no_layers_postselect))
                stddev_of_avg_counts_forcnot.insert(0, np.std(no_layers_postselect))

            if temp_fidelities:
                average_counts_forcnot.append(sum(temp_fidelities)/len(temp_fidelities))
                stddev_of_avg_counts_forcnot.append(np.std(temp_fidelities))

        # if layer_idx==0:
        if average_counts_forcnot:
            layers_plot=list(range(0, len(average_counts_forcnot), 1))
            ax.errorbar(layers_plot, average_counts_forcnot, stddev_of_avg_counts_forcnot, capsize=3.0, color=palette[cnot_idx], label=f"CNOT Count: {cnot_count}")
            # # Each element of fidelity_with_checks is a complete circuit result. So plot each separately.
            # for circ_idx, circ_result in enumerate(fidelity_with_checks):
            #     layer_count=circ_idx+1
            # ax.errorbar(one_qubit_error_space, average_fidelity_with_checks, std_deviations_with_checks, capsize=3.0, color=palette[layer_idx+1], label=f"Number of layers: {layer}")

    # ax.set_title(f"Compute Qubits: {number_of_qubits}, Error: {round(one_qubit_err, 5)}", fontsize=12)
    ax.set_ylabel("Average Probability of Zero\n Measurement on All Ancillas", fontsize=12)
    ax.set_xlabel("Number of Layers", fontsize=12)
    # ax.set_xscale("log")
    ax.legend()


def create_split_layers_count_percentage_subplots(base_path, number_of_qubits, cnot_count, one_qubit_error_space, layers, ax):
    '''Plots the average zero outcome for ancilla vs the single qubit error rate.'''
    '''Creates plots based on CNOT counts.'''
    # Gets the file path of the script.
    NLINES = len(layers)+1 #plus one for the zero checks.
    palette = sns.color_palette("viridis", NLINES)
    for layer_idx, layer in enumerate(layers):
        all_sorted_circ_results=get_all_sorted_circ_results_for_layer(base_path, cnot_count, number_of_qubits, layer)
        # print(all_sorted_circ_results)
        if all_sorted_circ_results:
            #Regular plots
            # NLINES = 45
            # each element is a circuit.
            counts_after=[]
            counts_before=[]
            for circ_result in all_sorted_circ_results:
                counts_after_temp=[]
                counts_before_temp=[]
                # each element here is one error idx for the circuit.
                for result in circ_result:
                    # In our simulations we have the full density matrix so the percent results before post selecting is actually 1.
                    # We do it this way so the graphing can accomidate physical runs where these two values are actually counts.
                    counts_after_temp.append(result["percent_results_after_postselect"])
                    counts_before_temp.append(result["percent_results_before_postselect"])
                counts_after.append(counts_after_temp)
                counts_before.append(counts_before_temp)
            
            average_counts_after=[sum(x)/len(x) for x in zip(*counts_after)]
            # print(average_counts_diffs)
            std_deviations_after=[np.std(x) for x in zip(*counts_after)]
            # print(std_deviations)

            if layer_idx==0:
                ax.errorbar(one_qubit_error_space, counts_before[0], color=palette[0], capsize=3.0, label="Number of layers: 0")
            # for layer_idx, circ_count_after in enumerate(counts_after):
            #     layer_count=layer_idx+1
            ax.errorbar(one_qubit_error_space, average_counts_after, std_deviations_after,
                capsize=3.0, color=palette[layer_idx+1], label=f"Number of layers: {layer}")
        ax.set_xlabel("Single Qubit Error", fontsize=12)
        ax.set_ylabel("Average Probability of Zero\n Measurement on All Ancillas", fontsize=12)
        ax.set_xscale("log")
        ax.legend()

def get_5qubits_fid_and_counts_vs_error_low_weightchks(rz_count, layers_all):
    '''Creates fidelity vs error rate and post select prob vs error rate
    on a single plot. Low weight checks.'''
    cnot_counts=[1,5,10,15,20,25,30,35,40]
    NUMBER_OF_ERROR_POINTS=21
    ONE_QUBIT_ERROR_SPACE=np.logspace(-5, -2, num=NUMBER_OF_ERROR_POINTS)
    NLINES = 45
    palette = sns.color_palette("viridis",NLINES)
    CODE_DIR=sys.path[0]
    NUMBER_OF_QUBITS=5
    SUBDIR=os.path.join(f"qubits_{NUMBER_OF_QUBITS}_rz_{rz_count}", "low_weight_checks", "results")
    BASE_PATH=os.path.join(CODE_DIR,SUBDIR)
    for layers in layers_all:
    # if layers == 4:
        # fig, axs = plt.subplots(2,1, figsize=(6.5,9.75), sharex=True)
        fig, axs = plt.subplots(2,1, figsize=(7,10.5), sharex=True)

        create_all_fidelity_plot_cnots_layers_rz(BASE_PATH, NUMBER_OF_QUBITS, cnot_counts, ONE_QUBIT_ERROR_SPACE, palette, layers, rz_count, axs[0])
        create_all_count_percentage_plots_layers(BASE_PATH, NUMBER_OF_QUBITS, cnot_counts, ONE_QUBIT_ERROR_SPACE, palette, layers, axs[1])
        plt.subplots_adjust(hspace=0.05) 
        # create_all_count_percentage_plots(BASE_PATH, NUMBER_OF_QUBITS, cnot_counts, ONE_QUBIT_ERROR_SPACE, palette, NLINES)
        str_labels=["(a)", "(b)"]
        for idx, ax in enumerate(axs):
            str_label=str_labels[idx]
            # label physical distance to the left and up:
            # trans = mtransforms.ScaledTranslation(-55/72, -10/72, fig.dpi_scale_trans)
            trans = mtransforms.ScaledTranslation(-65/72, -10/72, fig.dpi_scale_trans)

            ax.text(0.0, 1.0, str_label, transform=ax.transAxes + trans,
                    fontsize=16, weight="bold", va='center')
            ax.autoscale(tight=False)
        fig.savefig(os.path.join(BASE_PATH, f"qubits_{NUMBER_OF_QUBITS}_allCNOTS_fidelity_and_counts_layers_{layers}_rz_{rz_count}_.png"))

def create_all_maxweight_fidelity_vs_error(base_path, number_of_qubits, cnot_counts, one_qubit_error_space, palette, layers_count, rz_count):
    plt.cla()
    figure(figsize=(6.4, 4.8))
    all_results=[]
    for cnot_count in cnot_counts:
        result=create_fidelity_plot_cnot(base_path, number_of_qubits, cnot_count, one_qubit_error_space)
        if result:
            all_results.append(result)

    for idx, elem in enumerate(all_results):
        # print(f"fidelity gain for cnots {cnot_counts[idx]}: {elem[1]}")
        plt.errorbar(one_qubit_error_space, elem[1], yerr=elem[2], capsize=3.0, color=palette[idx*5], label="CNOT Count: "+str(elem[0]))
    plt.title(f"Qubits: {str(number_of_qubits)}, Layers: {layers_count}, Rz: {rz_count}", fontsize=12)
    plt.xlabel("Single Qubit Error", fontsize=12)
    plt.ylabel("Average Fidelity Gain", fontsize=12)
    plt.xscale("log")
    plt.legend()    
    plt.savefig(os.path.join(base_path, f"qubits_{number_of_qubits}_allCNOTS_fidelity_rz_{rz_count}_maxweight_.png"))

def get_5qubits_fid_vs_error_max_weightchks(rz_count):
    '''Creates fidelity vs error rate and post select prob vs error rate
    on a single plot. High weight checks.'''
    cnot_counts=[1,5,10,15,20,25,30,35,40]
    NUMBER_OF_ERROR_POINTS=21
    ONE_QUBIT_ERROR_SPACE=np.logspace(-5, -2, num=NUMBER_OF_ERROR_POINTS)
    NLINES = 50
    palette = sns.color_palette("viridis",NLINES)
    layers_all=[1]
    CODE_DIR=sys.path[0]
    NUMBER_OF_QUBITS=5
    SUBDIR=os.path.join(f"qubits_{NUMBER_OF_QUBITS}_rz_{rz_count}", "high_weight_checks", "results")
    BASE_PATH=os.path.join(CODE_DIR,SUBDIR)
    for layers in layers_all:
        create_all_maxweight_fidelity_vs_error(BASE_PATH, NUMBER_OF_QUBITS, cnot_counts, ONE_QUBIT_ERROR_SPACE, palette, layers, rz_count)

def get_5qubits_fidelity_vs_layers(rz_count):
    '''5 qubits and 10 qubits experiments.
    Plots the average fidelity gain vs number of layers. 
    The single qubit error is fixed. All CNOTs are on each plot.'''
    cnot_counts=[1,5,10,15,20,25,30,35,40]
    NUMBER_OF_ERROR_POINTS=21
    ONE_QUBIT_ERROR_SPACE=np.logspace(-5, -2, num=NUMBER_OF_ERROR_POINTS)
    # NLINES = 45
    # palette = sns.color_palette("viridis",NLINES)
    layers_all=[1, 2, 3, 4, 5, 6]
    CODE_DIR=sys.path[0]
    NUMBER_OF_QUBITS=5
    SUBDIR=os.path.join(f"qubits_{NUMBER_OF_QUBITS}_rz_{rz_count}", "low_weight_checks", "results")
    BASE_PATH=os.path.join(CODE_DIR,SUBDIR)
    for error_idx, error in enumerate(ONE_QUBIT_ERROR_SPACE):
        if error_idx==16:
            # if layers == 4:
            fig, axs = plt.subplots(2,1, figsize=(6.5,9.75), sharex=True)
            create_allcnot_fidelity_vs_layers_plots(BASE_PATH, NUMBER_OF_QUBITS, cnot_counts, ONE_QUBIT_ERROR_SPACE[error_idx], error_idx, layers_all, axs[0], rz_count)
            create_allcnot_counts_vs_layers_plots(BASE_PATH, NUMBER_OF_QUBITS, cnot_counts, ONE_QUBIT_ERROR_SPACE[error_idx], error_idx, layers_all, axs[1])
            plt.subplots_adjust(hspace=0.05) 

            str_labels=["(a)", "(b)"]
            for idx, ax in enumerate(axs):
                str_label=str_labels[idx]
                # label physical distance to the left and up:
                trans = mtransforms.ScaledTranslation(-55/72, -10/72, fig.dpi_scale_trans)
                ax.text(0.0, 1.0, str_label, transform=ax.transAxes + trans,
                        fontsize=16, weight="bold", va='center')
                ax.autoscale(tight=False)
            fig.savefig(os.path.join(BASE_PATH, f"qubits_{NUMBER_OF_QUBITS}_allCNOTS_fidelity_and_counts_vs_layers_rz_{rz_count}_.png"))

def get_2qubits_unit_fidel_ex():
    '''Generates example plots for 2 qubits that give fidelity 1 for noise only effecting compute qubits.'''
    cnot_count=30#[1,5,10,15,20,25,30,35,40]
    NUMBER_OF_ERROR_POINTS=21
    ONE_QUBIT_ERROR_SPACE_2QUBIT_EXAMPLE=np.logspace(-5, -1, num=NUMBER_OF_ERROR_POINTS)
    layers_all=[1, 2, 3, 4, 5, 6]
    CODE_DIR=sys.path[0]
    NUMBER_OF_QUBITS=2
    SUBDIR_2QUBITS_EXAMPLE=os.path.join("qubits_2_rz_0", "noiseless_checks_ex", "low_weight_checks", "results")
    BASE_PATH_2QUBITS=os.path.join(CODE_DIR,SUBDIR_2QUBITS_EXAMPLE)
    fig, axs = plt.subplots(2,1, figsize=(6.5,9.75), sharex=True)
    create_2qubit_plot_example_fidelity(BASE_PATH_2QUBITS, NUMBER_OF_QUBITS, cnot_count, ONE_QUBIT_ERROR_SPACE_2QUBIT_EXAMPLE, layers_all, axs[0])
    create_2qubit_plot_example_count_percentage(BASE_PATH_2QUBITS, NUMBER_OF_QUBITS, cnot_count, ONE_QUBIT_ERROR_SPACE_2QUBIT_EXAMPLE, layers_all, axs[1])
    plt.subplots_adjust(hspace=0.05) 
    str_labels=["(a)", "(b)"]
    for idx, ax in enumerate(axs):
        str_label=str_labels[idx]
        # label physical distance to the left and up:
        trans = mtransforms.ScaledTranslation(-55/72, -10/72, fig.dpi_scale_trans)
        ax.text(0.0, 1.0, str_label, transform=ax.transAxes + trans,
                fontsize=16, weight="bold", va='center')
        ax.autoscale(tight=False)
    fig.savefig(os.path.join(BASE_PATH_2QUBITS, f"qubits_{NUMBER_OF_QUBITS}_example_fidelity_and_counts_layers_all_.png"))

def get_5qubits_raw_fidel_and_counts_vs_error(rz_count):
    '''5 qubits
    Plots the average fidelity and post selection probability vs single qubit error. 
    CNOT count is fixed for each plot. All layers are on each plot.'''
    cnot_counts=[40]
    layers_all=[1, 2, 3, 4, 5, 6]
    CODE_DIR=sys.path[0]
    NUMBER_OF_QUBITS=5
    SUBDIR=os.path.join(f"qubits_{NUMBER_OF_QUBITS}_rz_{rz_count}", "low_weight_checks", "results")
    BASE_PATH=os.path.join(CODE_DIR,SUBDIR)
    NUMBER_OF_ERROR_POINTS=21
    ONE_QUBIT_ERROR_SPACE=np.logspace(-5, -2, num=NUMBER_OF_ERROR_POINTS)
    for cnot_count in cnot_counts:
        plt.clf()
        fig, axs = plt.subplots(2,1, figsize=(6.5,9.75), sharex=True)
        create_split_layers_fidelity_plot_cnots_subplots(BASE_PATH, NUMBER_OF_QUBITS, cnot_count, ONE_QUBIT_ERROR_SPACE, layers_all, axs[0], rz_count)
        create_split_layers_count_percentage_subplots(BASE_PATH, NUMBER_OF_QUBITS, cnot_count, ONE_QUBIT_ERROR_SPACE, layers_all, axs[1])
        plt.subplots_adjust(hspace=0.05) 

        str_labels=["(a)", "(b)"]
        for idx, ax in enumerate(axs):
            str_label=str_labels[idx]
            # label physical distance to the left and up:
            trans = mtransforms.ScaledTranslation(-55/72, -10/72, fig.dpi_scale_trans)
            ax.text(0.0, 1.0, str_label, transform=ax.transAxes + trans,
                    fontsize=16, weight="bold", va='center')
            ax.autoscale(tight=False)
        fig.savefig(os.path.join(BASE_PATH, f"qubits_{NUMBER_OF_QUBITS}_CNOTS_{cnot_count}_fidelity_and_counts_layers_{layers_all[-1]}_rz_{rz_count}_.png"))

def get_2qubits_fidel_counts_vs_cnots():
    ''' Plots of fidelity vs CNOTs.
    Plots the average fidelity gain vs number of layers. 
    The single qubit error is fixed.'''
    cnot_counts=[1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]
    layers_all=[1, 2, 3, 4]
    CODE_DIR=sys.path[0]
    NUMBER_OF_ERROR_POINTS=21
    ONE_QUBIT_ERROR_SPACE=np.logspace(-5, -2, num=NUMBER_OF_ERROR_POINTS)
    NUMBER_OF_QUBITS=2
    SUBDIR_2QUBITS_EXAMPLE=os.path.join("qubits_2_rz_0", "low_weight_checks", "results")
    BASE_PATH=os.path.join(CODE_DIR,SUBDIR_2QUBITS_EXAMPLE)
    for error_idx, _ in enumerate(ONE_QUBIT_ERROR_SPACE):
        if error_idx==14:
            # if layers == 4:
            fig, axs = plt.subplots(2,1, figsize=(6.5,9.75), sharex=True)
            # create_fidelity_vs_cnots_plots(BASE_PATH, NUMBER_OF_QUBITS, cnot_counts, ONE_QUBIT_ERROR_SPACE[error_idx], error_idx, layers_all, axs[0])
            # create_counts_vs_cnots_plots(BASE_PATH, NUMBER_OF_QUBITS, cnot_counts, ONE_QUBIT_ERROR_SPACE[error_idx], error_idx, layers_all, axs[1])
            create_fidelity_vs_cnots_plots_filtered_zero_prob(BASE_PATH, NUMBER_OF_QUBITS, cnot_counts, ONE_QUBIT_ERROR_SPACE[error_idx], error_idx, layers_all, axs[0])
            create_counts_vs_cnots_plots_filtered_zero_prob(BASE_PATH, NUMBER_OF_QUBITS, cnot_counts, ONE_QUBIT_ERROR_SPACE[error_idx], error_idx, layers_all, axs[1])

            plt.subplots_adjust(hspace=0.05) 

            str_labels=["(a)", "(b)"]
            for idx, ax in enumerate(axs):
                str_label=str_labels[idx]
                # label physical distance to the left and up:
                trans = mtransforms.ScaledTranslation(-55/72, -10/72, fig.dpi_scale_trans)
                ax.text(0.0, 1.0, str_label, transform=ax.transAxes + trans,
                        fontsize=16, weight="bold", va='center')
                ax.autoscale(tight=False)
            fig.savefig(os.path.join(BASE_PATH, f"qubits_{NUMBER_OF_QUBITS}_fidelity_and_counts_vs_cnots_.png"))

def get_10qubits_fid_and_counts_vs_error_low_weightchks(rz_count, layers_all):
    '''Creates fidelity vs error rate and post select prob vs error rate
    on a single plot. Low weight checks.'''
    cnot_counts=[1,5,10,15,20,25,30,35,40, 80]
    NUMBER_OF_ERROR_POINTS=21
    ONE_QUBIT_ERROR_SPACE=np.logspace(-5, -2, num=NUMBER_OF_ERROR_POINTS)
    NLINES = 50
    palette = sns.color_palette("viridis",NLINES)
    CODE_DIR=sys.path[0]
    NUMBER_OF_QUBITS=10
    SUBDIR=os.path.join(f"qubits_{NUMBER_OF_QUBITS}_rz_{rz_count}", "low_weight_checks", "results")
    BASE_PATH=os.path.join(CODE_DIR,SUBDIR)
    for layers in layers_all:
    # if layers == 4:
        # fig, axs = plt.subplots(2,1, figsize=(6.5,9.75), sharex=True)
        fig, axs = plt.subplots(2,1, figsize=(7,10.5), sharex=True)

        create_all_fidelity_plot_cnots_layers_rz(BASE_PATH, NUMBER_OF_QUBITS, cnot_counts, ONE_QUBIT_ERROR_SPACE, palette, layers, rz_count, axs[0])
        create_all_count_percentage_plots_layers(BASE_PATH, NUMBER_OF_QUBITS, cnot_counts, ONE_QUBIT_ERROR_SPACE, palette, layers, axs[1])
        plt.subplots_adjust(hspace=0.05) 
        # create_all_count_percentage_plots(BASE_PATH, NUMBER_OF_QUBITS, cnot_counts, ONE_QUBIT_ERROR_SPACE, palette, NLINES)
        str_labels=["(a)", "(b)"]
        for idx, ax in enumerate(axs):
            str_label=str_labels[idx]
            # label physical distance to the left and up:
            # trans = mtransforms.ScaledTranslation(-55/72, -10/72, fig.dpi_scale_trans)
            trans = mtransforms.ScaledTranslation(-65/72, -10/72, fig.dpi_scale_trans)

            ax.text(0.0, 1.0, str_label, transform=ax.transAxes + trans,
                    fontsize=16, weight="bold", va='center')
            ax.autoscale(tight=False)
        fig.savefig(os.path.join(BASE_PATH, f"qubits_{NUMBER_OF_QUBITS}_allCNOTS_fidelity_and_counts_layers_{layers}_rz_{rz_count}_.png"))

def get_10qubits_fid_vs_error_max_weightchks(rz_count):
    '''Creates fidelity vs error rate and post select prob vs error rate
    on a single plot. High weight checks.'''
    cnot_counts=[1,5,10,15,20,25,30,35,40, 80]
    NUMBER_OF_ERROR_POINTS=21
    ONE_QUBIT_ERROR_SPACE=np.logspace(-5, -2, num=NUMBER_OF_ERROR_POINTS)
    NLINES = 50
    palette = sns.color_palette("viridis",NLINES)
    layers_all=[1]
    CODE_DIR=sys.path[0]
    NUMBER_OF_QUBITS=10
    SUBDIR=os.path.join(f"qubits_{NUMBER_OF_QUBITS}_rz_{rz_count}", "high_weight_checks", "results")
    BASE_PATH=os.path.join(CODE_DIR,SUBDIR)
    for layers in layers_all:
        create_all_maxweight_fidelity_vs_error(BASE_PATH, NUMBER_OF_QUBITS, cnot_counts, ONE_QUBIT_ERROR_SPACE, palette, layers, rz_count)
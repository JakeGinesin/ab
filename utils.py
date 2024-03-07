#!/usr/bin/env python3

from automata.nfa import NFA
from automata.buchi import Buchi
from automata.dfa import DFA
from typing import Set, Tuple, List, Union, Dict

import random
import os
import subprocess
import re

def generate_random_fa(
    num_states : int = 10, # Number of states
    num_acceptance_states : int = 1, # Number of acceptance states
    num_symbols : int = None, # Number of (transition) symbols
    max_transitions : int = None, # Max number of transitions generated
    custom_alphabet : Set[str] = {"a","b","c"}, # Can give a custom alphabet if desired
    no_epsilon : bool = True,
    fa : str = "Buchi") -> Union[NFA, Buchi, DFA]:  # type of finite automata to generate 
    """ Generates a random finite automata """

    assert fa == "Buchi" or fa == "NFA" or fa == "DFA", "finite automata inputted is not found"

    states = {'q' + str(i) for i in range(num_states)}
    assert 0 < num_acceptance_states < num_states, "Number of acceptance states must be between 1 and the number of states minus 1"
    initial_state = random.choice(tuple(states))

    # Select a random subset of Q as the acceptance states F, ensuring initial_state is not included
    remaining_states = list(states - {initial_state})
    acceptance_states = set(random.sample(remaining_states, num_acceptance_states))

    # Create a random alphabet Sigma if alphabet not provided
    alphabet = custom_alphabet if custom_alphabet is not None else {chr(i + 97) for i in range(num_symbols)}

    # Generate a list of transitions T, ensure max_transitions is at least the number of states
    if max_transitions is None:
        max_transitions = len(states) * len(alphabet)

    transitions = []
    for _ in range(max_transitions):
        start_state = random.choice(tuple(states))
        end_state = random.choice(tuple(states))
        if fa == "Buchi" or fa == "NFA": 
            if no_epsilon : input_symbol = random.choice(tuple(alphabet))
            else : input_symbol = random.choice(tuple(alphabet) + ('',))  # including epsilon transitions

        else : input_symbol = random.choice(tuple(alphabet))
        transitions.append((start_state, input_symbol, end_state))

    if fa == "Buchi" :  return Buchi(states, initial_state, acceptance_states, alphabet, transitions)
    elif fa == "NFA" : return NFA(states, initial_state, acceptance_states, alphabet, transitions)
    else : return DFA(states, initial_state, acceptance_states, alphabet, transitions)

def read_ba_to_buchi(filename : str) -> Buchi:
    """ reads a .ba file and converts it to a Buchi Automata """
    with open(filename, 'r') as file:
        lines = file.readlines()

    # Initialize parameters for the Buchi object
    initial_state = None
    states = set()
    acceptance_states = set()
    transitions = []
    alphabet = set()

    # Parse initial state
    initial_state = lines[0].strip()
    states.add(initial_state)

    for line in lines[1:]:
        line = line.strip()

        # Parse transition
        if "->" in line:
            label, transition = line.split(",")
            src, dest = transition.split("->")
            src, dest = src.strip(), dest.strip()

            states.add(src)
            states.add(dest)
            alphabet.add(label)

            transitions.append((src, label, dest))

        # Parse acceptance state
        elif line:
            acceptance_states.add(line)

    return Buchi(states, initial_state, acceptance_states, alphabet, transitions)

# we can always interpret a buchi automata as a NFA, so... there's no problem with this
def read_ba_to_NFA(filename : str) -> NFA:
    """ reads a .ba file and converts it to a NFA """
    with open(filename, 'r') as file:
        lines = file.readlines()

    initial_state = None
    states = set()
    acceptance_states = set()
    transitions = []
    alphabet = set()

    # Parse initial state
    initial_state = lines[0].strip()
    states.add(initial_state)

    for line in lines[1:]:
        line = line.strip()

        # Parse transition
        if "->" in line:
            label, transition = line.split(",")
            src, dest = transition.split("->")
            src, dest = src.strip(), dest.strip()

            states.add(src)
            states.add(dest)
            alphabet.add(label)

            transitions.append((src, label, dest))

        # Parse acceptance state
        elif line:
            acceptance_states.add(line)

    return NFA(states, initial_state, acceptance_states, alphabet, transitions)

def ltl_to_ba(LTL : str) -> Buchi:
    """ converts an LTL formula to a Buchi Automata """
    cmd = ['/bin/ltl2tgba', '-8B', str(LTL), '-d']
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()

    # convert bytes to string
    stdout = stdout.decode()
    stderr = stderr.decode()

    if proc.returncode != 0:
        print(f"The command '{' '.join(cmd)}' failed with error:\n{stdout}")

    try:
        result_pan = subprocess.run(cmd, check=True, text=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        print(f"The command '{' '.join(cmd_pan)}' failed with error:\n{e.output}")
        return None

    lines = result_pan.stdout.splitlines()
    init_passed = False
    buchi = {}
    buchi["states"] = set()
    buchi["acceptance_states"] = set()
    buchi["transitions"] = []
    buchi["alphabet"] = set()
    buchi["initial_state"] = None
    for line in lines:
        if not init_passed:
            res_init = re.search("^\s*I\s->\s(\d+).*$", line)
            if res_init: 
                init_passed = True
                buchi["initial_state"] = res_init.group(1)
                buchi["states"].add(res_init.group(1))
                continue
        else:
            res_transition = re.search("^\s*(\d+) -> (\d+) \[label=\"(.*)\"\]\s*$", line)
            if res_transition:
                a1 = res_transition.group(1)
                a2 = res_transition.group(2)
                t = res_transition.group(3)
                buchi["states"].add(a1)
                buchi["states"].add(a2)
                buchi["transitions"].append((a1, t, a2))
                buchi["alphabet"].add(t)
                # add to language
                continue

            res_accept_state = re.search("^\s*(\d+) \[label=\"(\d+)\", peripheries=2]\s*$", line)
            if res_accept_state:
                a1 = res_accept_state.group(1)
                buchi["states"].add(a1)
                buchi["acceptance_states"].add(a1)
                continue

            res_state = re.search("^\s*(\d+) \[label=\"(\d+)\"]\s*$", line)
            if res_state:
                a1 = res_state.group(1)
                buchi["states"].add(a1)
                continue

    return Buchi(buchi["states"], buchi["initial_state"], buchi["acceptance_states"], buchi["alphabet"], buchi["transitions"])

def construct_buchi_from_promela(path : str) -> Dict[str, Buchi]:
    """ construct multiple Buchi Automatas from a Promela file """

    cmd = ['spin', '-a', '-run', '-DNOREDUCE', path]
    cmd_pan = ['./pan', '-d']

    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    filename = os.path.basename(path)
    userdir = os.getcwd()

    # Convert bytes to string
    stdout = stdout.decode()
    stderr = stderr.decode()

    if proc.returncode != 0:
        print(f"The command '{' '.join(cmd)}' failed with error:\n{stdout}")

    try:
        result_pan = subprocess.run(cmd_pan, check=True, text=True, capture_output=True, cwd=userdir)
    except subprocess.CalledProcessError as e:
        print(f"The command '{' '.join(cmd_pan)}' failed with error:\n{e.output}")
        return None

    lines = result_pan.stdout.splitlines()

    current_process = {}
    current_proctype = None
    for line in lines:
        # if we hit a match, we stop adding states to the previous proctype
        res_proctype = re.search("^proctype (\S*)$", line)
        res_init = re.search("^init$", line)
        res_claim = re.search("^claim (\S*)$", line) 
        if res_proctype:
            proctype = res_proctype.group(1)
            current_process[proctype] = {}
            current_process[proctype]["states"] = set()
            # current_process[proctype]["states"].add('1')
            current_process[proctype]["acceptance_states"] = set()
            current_process[proctype]["transitions"] = []
            current_process[proctype]["initial_state"] = None
            current_process[proctype]["alphabet"] = set()
            current_proctype = proctype

        elif res_init:  
            proctype = "init"
            current_process[proctype] = {}
            current_process[proctype]["states"] = set()
            # current_process[proctype]["states"].add('1')
            current_process[proctype]["acceptance_states"] = set()
            current_process[proctype]["transitions"] = []
            current_process[proctype]["initial_state"] = None
            current_process[proctype]["alphabet"] = set()
            current_proctype = proctype
        
        elif res_claim:
            proctype = res_claim.group(1)
            current_process[proctype] = {}
            current_process[proctype]["states"] = set()
            # current_process[proctype]["states"].add('1')
            current_process[proctype]["acceptance_states"] = set()
            current_process[proctype]["transitions"] = []
            current_process[proctype]["initial_state"] = None
            current_process[proctype]["alphabet"] = set()
            current_proctype = proctype
        else:
            # we're either in a transition line, or we're at the end
            #data = re.search("^ *state +([0-9]+) +\-[(]tr +([0-9]*)[)]\-\> +state +([0-9]*) +\[id +[0-9]* tp +[0-9]*\] \[(\S\S\S\S\S)\] +\S+ => (\S+) *$", line)
            data = re.search("^\s*state\s+([0-9]+)\s+\-[(]tr\s+([0-9]*)[)]\-\>\s+state +([0-9]*)\s+\[id +[0-9]*\stp +[0-9]*\]\s\[(\S\S\S\S\S)\]\s+\S+\s=>\s(.+)$", line)
            if data:
                if current_process[current_proctype]["states"] == set():
                    current_process[proctype]["initial_state"] = data.group(1)
                else:
                    current_process[proctype]["acceptance_states"].add(data.group(1))
                    current_process[proctype]["acceptance_states"].add(data.group(3))

                current_process[current_proctype]["states"].add(data.group(1))
                current_process[current_proctype]["states"].add(data.group(3))
                current_process[current_proctype]["transitions"].append((data.group(1), data.group(5), data.group(3)) )
                current_process[current_proctype]["alphabet"].add(data.group(5))
            
    ret = {} 
    for (name, d) in current_process.items():
        print(name)
        ret[name] = Buchi(current_process[name]["states"], current_process[name]["initial_state"], current_process[name]["acceptance_states"],  current_process[name]["alphabet"], current_process[name]["transitions"])

    panloc = os.path.join(userdir, "pan") 

    if os.path.exists(panloc):
        os.remove(panloc)

    trailoc = os.path.join(userdir, filename + ".trail")

    if os.path.exists(trailoc):
        os.remove(trailoc)

    return ret

def construct_buchi_from_hoa_file(filename : str) -> Buchi:
    """ 
    construct a Buchi Automata from a HOA file 
    HOA files have lots of detail: https://adl.github.io/hoaf/examples.html
    """
    with open(filename, 'r') as file:
        content = file.read()
        lines = content.splitlines()

    initial_state = re.search(r'Start:\s*(\d+)', content).group(1)
    states = set()
    acceptance_states = set()
    transitions = []

    alphabet_parse = re.search(r'AP: (\d+)(( \"[\w\d\_]+\")*)', content).group(2)
    alphabet_list = [s.replace('"', '') for s in alphabet_parse[1:].split(" ")]
    alphabet = set(alphabet_list)

    parse = False
    current_state = ""
    for line in lines:
        if line == "--BODY--" : parse = True
        elif line == "--END--" : parse = False

        if parse:
            if "State" in line:
                current_state = re.search(r'State:\s*(\d+)', line).group(1)
                states.add(current_state)
                if re.search(r'{\d+}', line) : acceptance_states.add(current_state)
            elif "[" in line:
                parse = re.search(r'(\[[\!\d\&\s\|\]+]+)\s(\d+)', line)
                labels_unparsed = parse.group(1)
                endstate = parse.group(2)

                # parse the horrible label format
                for item in labels_unparsed.replace("[","").replace("]","").split(" | "):
                    labels = [s for s in item.split("&") if "!" not in s]
                    if labels == [] and len(item.split("&")) == len(alphabet):
                        # epsilon transition
                        transitions.append((str(current_state), '', str(endstate)))
                    elif len(labels) == 1:
                        transitions.append((str(current_state), str(alphabet_list[int(labels[0])]), str(endstate)))
                        # normal transition
                    elif labels == [] and len(item.split("&")) < len(alphabet):
                        continue
                        # nothing

    return Buchi(states, initial_state, acceptance_states, alphabet, transitions)

def construct_buchi_from_hoa_string(content : str) -> Buchi:
    """ 
    construct a Buchi Automata from a HOA string
    HOA files have lots of detail: https://adl.github.io/hoaf/examples.html
    """
    lines = content.splitlines()

    initial_state = re.search(r'Start:\s*(\d+)', content).group(1)
    states = set()
    acceptance_states = set()
    transitions = []

    alphabet_parse = re.search(r'AP: (\d+)(( \"[\w\d\_]+\")*)', content).group(2)
    alphabet_list = [s.replace('"', '') for s in alphabet_parse[1:].split(" ")]
    alphabet = set(alphabet_list)

    parse = False
    current_state = ""
    for line in lines:
        if line == "--BODY--" : parse = True
        elif line == "--END--" : parse = False

        if parse:
            if "State" in line:
                current_state = re.search(r'State:\s*(\d+)', line).group(1)
                states.add(current_state)
                if re.search(r'{\d+}', line) : acceptance_states.add(current_state)
            elif "[" in line:
                parse = re.search(r'(\[[\!\d\&\s\|\]+]+)\s(\d+)', line)
                labels_unparsed = parse.group(1)
                endstate = parse.group(2)

                # parse the horrible label format
                for item in labels_unparsed.replace("[","").replace("]","").split(" | "):
                    labels = [s for s in item.split("&") if "!" not in s]
                    if labels == [] and len(item.split("&")) == len(alphabet):
                        # epsilon transition
                        transitions.append((str(current_state), '', str(endstate)))
                    elif len(labels) == 1:
                        transitions.append((str(current_state), str(alphabet_list[int(labels[0])]), str(endstate)))
                        # normal transition
                    elif labels == [] and len(item.split("&")) < len(alphabet):
                        continue
                        # nothing

    return Buchi(states, initial_state, acceptance_states, alphabet, transitions)

#!/usr/bin/env python3
from __future__ import annotations

from automata.finite_automata import Finite_Automata
from typing import Set, Tuple, List, Dict
import heapq
from z3 import *
import os
import subprocess
import re

class Buchi(Finite_Automata):
    def __init__(
        self,
        states : set,
        initial_state : str, 
        acceptance_states : set,
        alphabet : set,
        transitions : List[Tuple[str, str, str]] = []):

        Finite_Automata.__init__(self, states, initial_state, acceptance_states, alphabet, transitions)

        assert "" not in self.alphabet, "epsilon transitions (defined just as empty transitions, "") not allowed"
        self.alphabet.append('')

        self.construct_transition_function_hashmap()
        self.construct_reachability_hashmap()

    def minimal_parikh_image(self, init_state : str = None, usable_states : set = set()) -> Dict[str, int]:
        """
        a parikh image is simply a vector representing the number of characters in any given accepted string.
        finding the *minimal* parikh image for a NFA/buchi automata is trivially done using SMT/ILP. 
        see https://arxiv.org/abs/1002.1464 for details

        this is also in relation to parikh automata 
        """
        if init_state == None : init_state = self.initial_state
        # note: we only consider states relevant to the language of the automata
        if usable_states == set() : usable_states = self.get_wellconnected_states()

        # if init state is an accepting state, our parikh image is empty
        if init_state in self.acceptance_states : return {char: 0 for char in self.alphabet} 

        incoming_transitions, outgoing_transitions = {}, {}
        for state in usable_states:
            incoming_transitions[state] = set()
            outgoing_transitions[state] = set()

        z3_transitions = []
        for t in self.transitions:
            if t[0] in usable_states and t[2] in usable_states : z3_transitions.append(t)

        # maintain a mapping between transition labels and actual transitions
        transition_mapping = {}
        transition_labels = ["t" + str(i) for i in range(len(z3_transitions))]

        count = 0
        for (sA, i, sB) in z3_transitions:
            transition_mapping[transition_labels[count]] = (sA, i, sB)
            outgoing_transitions[sA].add("t" + str(count))
            incoming_transitions[sB].add("t" + str(count))

            count = count + 1

        # plug flow constraints into z3

        solver = Optimize()
        tosolve = {}
        accepting_states_to_solve = set()
        for t in transition_labels : tosolve[t] = Int(str(t))
        for t in tosolve.values() : solver.add(t >= 0) # transitions must not have a negative value

        solver.minimize(Sum(list(tosolve.values())))

        for state in usable_states:
            if state in self.acceptance_states : accepting_states_to_solve.add(state)
            elif state == init_state:
                # initial state -> incoming + 1 = outgoing (i.e. outgoing is one higher than incoming)
                solver.add(Sum([tosolve[t] for t in incoming_transitions[state]]) + 1 == Sum([tosolve[t] for t in outgoing_transitions[state]]))
            else:
                # incoming = outgoing
                solver.add(Sum([tosolve[t] for t in incoming_transitions[state]]) == Sum([tosolve[t] for t in outgoing_transitions[state]]))

        # assert acceptance states have at least one flow difference between input and output
        solver.add(
            Or([Sum([tosolve[t] for t in incoming_transitions[state]]) + 0 == 
                Sum([tosolve[t] for t in outgoing_transitions[state]]) + 1 
                for state in accepting_states_to_solve])
        )

        # assert acceptance states have at least two incoming flow
        solver.add(Or([ Sum([tosolve[t] for t in incoming_transitions[state]]) == 2 for state in accepting_states_to_solve]))

        # assert acceptance states have at least one outgoing flow
        solver.add(Or([ Sum([tosolve[t] for t in outgoing_transitions[state]]) == 1 for state in accepting_states_to_solve]))

        if solver.check() == sat:
            m = solver.model()
            solns = {var: m[var] for var in tosolve.values()}
            
            char_total = {}
            for char in self.alphabet:
                char_total[char] = 0

            for transition, value in solns.items():
                sA, i, sB = transition_mapping[str(transition)]
                char_total[i] += int(str(value))

            return char_total

        else:
            return {}

    def parikh_image_inclusion(self, parikh_image = Dict[str, int], init_state : str = None, usable_states : set = set()) -> bool:
        """
        determines if a given parikh image is included in the language of the buchi automata. 
        this problem can nicely be encoded as an ILP problem, and is NP hard via: https://arxiv.org/abs/1002.1464
        """
        if init_state == None : init_state = self.initial_state
        # note: we only consider states relevant to the language of the automata
        if usable_states == set() : usable_states = self.get_wellconnected_states()
        for char in self.alphabet : assert char in parikh_image, "parikh image doens't have all the required characters"

        # if init state is an accepting state, our parikh image is empty
        if init_state in self.acceptance_states : return {char: 0 for char in self.alphabet} 

        incoming_transitions, outgoing_transitions = {}, {}
        for state in usable_states:
            incoming_transitions[state] = set()
            outgoing_transitions[state] = set()

        z3_transitions = []
        for t in self.transitions:
            if t[0] in usable_states and t[2] in usable_states : z3_transitions.append(t)

        # maintain a mapping between transition labels and actual transitions
        transition_mapping = {}
        transition_labels = ["t" + str(i) for i in range(len(z3_transitions))]

        # map characters to certain transitions
        character_mapping = {}
        for char in self.alphabet : character_mapping[char] = set()

        count = 0
        for (sA, i, sB) in z3_transitions:
            transition_mapping[transition_labels[count]] = (sA, i, sB)
            outgoing_transitions[sA].add("t" + str(count))
            incoming_transitions[sB].add("t" + str(count))
            character_mapping[i].add(transition_labels[count])

            count = count + 1

        # plug flow constraints into z3

        solver = Solver()
        tosolve = {}
        accepting_states_to_solve = set()
        for t in transition_labels : tosolve[t] = Int(str(t))
        for t in tosolve.values() : solver.add(t >= 0) # transitions must not have a negative value

        for char in self.alphabet:
            solver.add(
                Sum([tosolve[t] for t in character_mapping[char]]) == parikh_image[char]
            )

        for state in usable_states:
            if state in self.acceptance_states : accepting_states_to_solve.add(state)
            elif state == init_state:
                # initial state -> incoming + 1 = outgoing (i.e. outgoing is one higher than incoming)
                solver.add(Sum([tosolve[t] for t in incoming_transitions[state]]) + 1 == Sum([tosolve[t] for t in outgoing_transitions[state]]))
            else:
                # incoming = outgoing
                solver.add(Sum([tosolve[t] for t in incoming_transitions[state]]) == Sum([tosolve[t] for t in outgoing_transitions[state]]))

        # assert acceptance states have at least one flow difference between input and output
        solver.add(
            Or([Sum([tosolve[t] for t in incoming_transitions[state]]) + 0 == 
                Sum([tosolve[t] for t in outgoing_transitions[state]]) + 1 
                for state in accepting_states_to_solve])
        )

        # assert acceptance states have at least two incoming flow
        solver.add(Or([ Sum([tosolve[t] for t in incoming_transitions[state]]) == 2 for state in accepting_states_to_solve]))

        # assert acceptance states have at least one outgoing flow
        solver.add(Or([ Sum([tosolve[t] for t in outgoing_transitions[state]]) == 1 for state in accepting_states_to_solve]))

        if solver.check() == sat:
            m = solver.model()
            solns = {var: m[var] for var in tosolve.values()}
            
            char_total = {}
            for char in self.alphabet:
                char_total[char] = 0

            for transition, value in solns.items():
                sA, i, sB = transition_mapping[str(transition)]
                char_total[i] += int(str(value))

            return char_total

        else:
            return {}

    def to_ba_file(self, path : str = "out.ba") -> None:
        """
        turn this Buchi Automata into a .ba formatted automata file
        see: https://languageinclusion.org/doku.php?id=tools#the_ba_format
        """
        outba = ""
        if all(isinstance(str(item), (int)) for item in self.states):
            # if this is the case, then we don't have to enumerate states
            outba+= "[" + str(self.initial_state) + "]\n"
            for (sA,i,sB) in self.transitions:
                outba += str(i) + "," + "[" + str(sA) + "]->[" + str(sB) + "]\n"

            for state in self.acceptance_states:
                outba += "[" + str(state) + "]\n"
        else:
            states_mapping = {}
            list_states = list(self.states)
            for i in range(len(list_states)):
                states_mapping[list_states[i]] = i

            outba+= "[" + str(states_mapping[self.initial_state]) + "]\n"
            for (sA,i,sB) in self.transitions:
                outba += str(i) + "," + "[" + str(states_mapping[sA]) + "]->[" + str(states_mapping[sB]) + "]\n"

            for state in self.acceptance_states:
                outba += "[" + str(states_mapping[state]) + "]\n"

        with open(path, 'w') as file:
            file.write(outba)

    def complement(self) -> Buchi:
        """
        take the complement of this buchi automata. 
        uses ranker: https://github.com/vhavlena/ranker?tab=readme-ov-file
        """
        self.to_ba_file("/tmp/out.ba")
        cmd = ['../bin/ranker', "/tmp/out.ba"]

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

        content = result_pan.stdout

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

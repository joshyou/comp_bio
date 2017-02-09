import math
import random
import re
import sys


'''
This class implements the Viterbi algorithm for finding the best state path of a 
Hidden Markov Model. find_prob recursively computes the probability that 
the sequence up to the given index ended in a given state. This updates the 
table state_probs storing these results, which find_best_path can use to 
reconstruct the best path to the end of the sequence.
'''
class HMM:

    def __init__(self, sequence, states, state_matrix, output_matrix, begin_matrix):
        self.sequence = sequence
        self.states = states

        #contains transition probabilities between each pair of states
        self.state_matrix = state_matrix

        #contains initial probability that the sequence will begin in any state
        self.begin_matrix = begin_matrix

        #contains probabilities that each output will result from each state
        self.output_matrix = output_matrix

        #Dynamic programming table which, for each state, stores
        #the probability of the optimal state-path ending at that state
        #and each index, as well as the preceding state
        #Each entry in the dictionary is a list of pairs [probability, last state]
        self.state_probs = {}

        #this initialization makes it so the indices start at 1
        for state in self.states:
            self.state_probs[state] = [[0, "begin"]]


    def find_best_path(self):
        
        best_k = self.states[0]

        #finds the state with the most probable path to the end of the sequence 

        for k in self.states:
            if self.state_probs[k][-1][0] > self.state_probs[best_k][-1][0]:
                best_k = k

        best_path = [best_k]

        #uses backpointers to reconstruct the state path
        index = 1
        current_state = best_k
        while len(best_path) < len(self.sequence):
            current_state = self.state_probs[current_state][-1 * index][1]
            best_path.append(current_state)
            index += 1
        
        #reverses the state path
        reversed = []
        for i in range(1, len(best_path)+1):
            reversed.append(best_path[-i])

        return reversed


    #returns probability of the best path ending at that index
    def find_best_prob(self, index):
        if index > len(self.sequence):
            return None

        if index == 0:
            return 1

        best_prob = -float("inf")

        for state in self.states:
            state_prob = self.find_prob(state, index)
            if state_prob > best_prob:
                best_prob = state_prob

        return best_prob

    #Recursively computes and returns probability of best path ending in given
    #state and index. Indices start at 1, so index = 1 means looking at the first 
    #character in the sequence
    def find_prob(self, state, index):
        #if the probability of a given state and index have already been
        #computed, return the stored value
        if len(self.state_probs[state]) > index:
            return self.state_probs[state][index][0]

        best_prob = -float("inf")
        out = self.prob_output(self.sequence[index-1], state)
        best_k = state
        #loop through states and find the most likely preceding state
        for k in self.states:
            #change_prob is probability k will transition to state
            #for first index, look at begin_matrix since it is transitioning from 
            #the "begin" state to state
            if index == 1:
                change_prob = self.begin_matrix[state]
            else:
                change_prob = self.prob_state_change(k, state)

            #state_prob is probability of best path ending in k and index - 1
            if index < 2:
                state_prob = 0
            else:
                state_prob = self.find_prob(k, index - 1)

            #state_prob * change_prob is probability the HMM followed the most likely path
            #ending at state k at index - 1 and then transitioned to state
            #probabilities are logarithms so they are added instead of multiplied
            if change_prob + state_prob > best_prob:
                best_prob = change_prob + state_prob
                best_k = k

        best_prob += self.prob_output(self.sequence[index-1], state)
        #update dictionary
        self.state_probs[state].append([best_prob, best_k])
        return best_prob
            
    #probability that state 1 transitions to state 2
    def prob_state_change(self, state1, state2):
        return self.state_matrix[state1][state2]

    #probability that a state produces a given output
    def prob_output(self, output, state):
        return self.output_matrix[state][output]

    def get_state_probs(self):
        return self.state_probs

def ln(n):
    return math.log(n)

def find_CGI(sequence):


    CPG_states = ["A+", "C+", "G+", "T+", "A-", "C-", "G-", "T-"]

    CPG_begin_matrix = {"A-": ln(.2465), "C-": ln(.2465), "G-": ln(.2465), "T-": ln(.2465),
        "A+": ln(.0035), "C+": ln(.0035), "G+": ln(.0035), "T+": ln(.0035)}

    CPG_state_matrix = {"A+": {"A+": ln(0.176), "C+": ln(.268), "G+": ln(.417), "T+": ln(.117), "A-": ln(.0037), "C-": ln(.0056), "G-": ln(.0086), "T-": ln(.0025)},
    "C+": {"A+": ln(.167), "C+":ln(.36), "G+":ln(.268), "T+":ln(.184), "A-":ln(.00354), "C-":ln(.00747), "G-":ln(.00559), "T-":ln(.00387)},
    "G+": {"A+": ln(.157), "C+":ln(.332), "G+":ln(.367), "T+":ln(.112), "A-":ln(.0034), "C-":ln(.0069), "G-":ln(.0076), "T-":ln(.0026)},
    "T+": {"A+": ln(.077), "C+":ln(.348), "G+":ln(.376), "T+":ln(.178), "A-":ln(.0017), "C-":ln(.0072), "G-":ln(.0078), "T-":ln(.00376)},
    "A-": {"A+": ln(.00042), "C+":ln(.00033), "G+":ln(.000408), "T+":ln(.00033), "A-":ln(.299), "C-":ln(.2047), "G-":ln(.285), "T-":ln(.2097)},
    "C-": {"A+": ln(.000447), "C+":ln(.00042), "G+":ln(.0002), "T+":ln(.000427), "A-":ln(.321), "C-":ln(.2975), "G-":ln(.078), "T-":ln(.301)},
    "G-": {"A+": ln(.0003), "C+":ln(.00036), "G+":ln(.000417), "T+":ln(.000417), "A-":ln(.177), "C-":ln(.239), "G-":ln(.2915), "T-":ln(.2915)},
    "T-": {"A+": ln(.000372), "C+":ln(.00037), "G+":ln(.000423), "T+":ln(.00033), "A-":ln(.2476), "C-":ln(.2456), "G-":ln(.2975), "T-":ln(.2077)}
    }

    CPG_output_matrix = {"A+": {"A":0, "C":-float("inf"), "G":-float("inf"), "T":-float("inf")},
    "C+": {"A":-float("inf"), "C":0, "G":-float("inf"), "T":-float("inf")},
    "G+": {"A":-float("inf"), "C":-float("inf"), "G":0, "T":-float("inf")},
    "T+": {"A":-float("inf"), "C":-float("inf"), "G":-float("inf"), "T":0},
    "A-": {"A":0, "C":-float("inf"), "G":-float("inf"), "T":-float("inf")},
    "C-": {"A":-float("inf"), "C":0, "G":-float("inf"), "T":-float("inf")},
    "G-": {"A":-float("inf"), "C":-float("inf"), "G":0, "T":-float("inf")},
    "T-": {"A":-float("inf"), "C":-float("inf"), "G":-float("inf"), "T":0},
    }

    CPG_hmm = HMM(sequence, CPG_states, CPG_state_matrix, CPG_output_matrix, CPG_begin_matrix)
    CPG_hmm.find_best_prob(len(sequence))

    path = CPG_hmm.find_best_path()
    
    #once the path is found, the start and endpoints of CG islands are located
    #'+' states are in CGI islands, so search for sequences of + states
    index = 0
    islands = []
    
    for k in path:
        #beginning of CGI
        if "+" in k and (index == 0 or "-" in path[index-1]):
            new_island = [index]
            islands.append(new_island)
        elif "-" in k:
            if index > 0 and "+" in path[index-1]:
                islands[-1].append(index)

        if (index == len(path) - 1) and len(islands) > 0 and len(islands[-1]) < 2:
            islands[-1].append(index)
        index += 1

    for island in islands:
        print island

    if len(islands) == 0:
        print "no CGIs detected"
    return islands

def main():
    filename = sys.argv[1]
    f = open(filename, 'r')
    sequence = ""
    for line in f:
        #removes characters that aren't A, G, C, or T
        line = re.sub(r'[^ACGT]', "", line)
        line = line.replace("\n", "")
        sequence += line
    sys.setrecursionlimit(10000)
    find_CGI(sequence)

if __name__ == "__main__":
    main()
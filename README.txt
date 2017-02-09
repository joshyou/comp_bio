hmm.py implements the Viterbi algorithm, which is a dynamic programming algorithm that identifies the most 
probable series of hidden states that a Hidden Markov Model went through to produce a set of results. A 
Hidden Markov Model is a model with a set of hidden states, where each state has a probability distribution
for what output it will produce, as well as which state it will transition to. The state transition model was 
determined empirically from Han et. al. (2008).

Han, Leng, et al. "CpG island density and its correlations with genomic features in mammalian genomes." Genome Biol 9.5 (2008): R79.


Usage: 

>python hmm.py [filename]

hmm.py takes in a filename as a command line argument and returns the start and endpoints of any CGIs it identifies 
from the sequence. The program strips any characters that aren't ACGT so the file can be any format.
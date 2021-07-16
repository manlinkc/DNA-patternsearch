"""M345SC Homework 1, part 1
Manlin Chawla 01205586
"""
def ksearchdictionary(S,k,f,x):
    """
    Dictionary implementation that is referenced to in the docstring of ksearch.
    """
    # Intialize lists for efficiency, will append to these later
    L1,L2,L3=[],[],[]

    # Initialize empty dictionary
    d = dict()

    # Define the bases to use for x-point mutations
    bases = ['A','C','G','T']

    # Get length of S
    lenS=len(S)

    # Iterate across the string S adding each k-mer to the dictionary
    for i in range(lenS-k+1):

        # Condition to check whether k-mer is in the dictionary
        # 0th entry of dictionary value records frequency of k-mers
        # 1st, 2nd, 3rd ... entry of dictionary value records location of k-mer

        # If kmer is in the dictionary
        if S[i:k+i] in d:

            # +1 counter to record frequency of k-mer
            d[S[i:k+i]][0] += 1
            # Add location of k-mer to the value
            d[S[i:k+i]].append(i)

        # If k-mer is not in the dictionary
        else:
            # Set frequency to 1 and record location of k-mer
            d[S[i:k+i]] = [1,i]

    # L1 contains kmers that have occured at least f times
    L1 = [key for key in d if (d[key][0] >= f)]

    # L2 contains the location of the frequently occuring kmers
    L2 = [d[key][1:] for key in d if (d[key][0] >= f)]

    # Iterate of frequently occuring kmers
    for kmer in L1:
        # Change kmer string into a list
        #kmerlist = list(kmer)
        # Initialize sum to count number of times x-mutations occur for each kmer
        sum = 0
        # Iterate over the bases
        for base in bases:
            # If xth position of kmer is the same as base, continue to next base
            if (kmer[x]==base):
                continue
            # If xth position of kmer is not the same as base
            else:
                # Replace xth position of kmer with the new base letter
                #kmerlist[x]=base
                # Join list back to a string
                #newkmer=''.join(kmerlist)
                newkmer=kmer[:x]+base+kmer[x+1:]

                # If newkmer exists in dictionary, add frequency of newkmer to sum
                if newkmer in d:
                    sum = sum + d[newkmer][0]

        # Append frequency of x-mutation of the current kmer to L3
        L3.append(sum)

    return L1,L2,L3

if __name__=='__main__':
    #Sample input and function call. Use/modify if/as needed
    #infile = open('test_sequence.txt', 'r')
    #S=infile.read()
    #infile.close()
    S='CCCTATGTTTGTGTGCACGGACACGTGCTAAAGCCAAATATTTTGCTAGT'
    k=3
    x=1
    f=2

    L1,L2,L3=ksearch(S,k,f,x)

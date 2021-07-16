# Rolling Hash Function Implementation

def initialhash(kmer,k,d,q,power4):
    """
    This function computes the rolling hash value of a k-length string
    Input:
    kmer  : k-length string consisting of A,C,G and Ts
    k     : The size of k-mer
    d     : Dictionary containing genebases and the corresponding numbers assigned used for computing rolling hash
    q     : Large prime used for modulu operations
    power4: List containing the powers of 4 in descending order of magnitude, precomputed in the ksearch
    """
    # Initialize sum to compute hash of kmer
    hashvalue = 0

    # Iterate over the characters in kmer
    for i in range(k):
        # Add character*powerof4 modulo prime
        hashvalue = hashvalue + ((d[kmer[i]]*power4[i]) % q)

    return hashvalue

def iteratedhash(oldkmer,newkmer,d,bm,q,rhash):
    """
    This function computes the iterated rolling hash value of a k-length string
    based on the hash value of the previous k-length string.
    Input:
    oldkmer : Previous k-length string consisting of A,C,G and Ts
    newkmer : New k-length string consisting of A,C,G and Ts
    rhash   : Hash value of oldkmer that is used to compute hash value of newkmer
    d       : Dictionary containing genebases and the corresponding numbers assigned used for computing rolling hash
    q       : Large prime used for modulu operations
    bm      : (4**k)%q value precomputed in ksearch to use for new hash value
    """
    # Compute hash of kmer based on hash of previous kmer
    iterhash = (4*rhash - d[oldkmer[0]]*bm+d[newkmer[-1]]) % q

    return iterhash

def ksearchhash(S,k,f,x):
    """
    Search for frequently-occurring k-mers within a DNA sequence
    and find the number of point-x mutations for each frequently
    occurring k-mer.
    Input:
    S: A string consisting of A,C,G, and Ts
    k: The size of k-mer to search for (an integer)
    f: frequency parameter -- the search should identify k-mers which
    occur at least f times in S
    x: the location within each frequently-occurring k-mer where
    point-x mutations will differ.

    Output:
    L1: list containing the strings corresponding to the frequently-occurring k-mers
    L2: list containing the locations of the frequent k-mers
    L3: list containing the number of point-x mutations for each k-mer in L1.

    Discussion: Add analysis here
    """

    # Intialize lists for efficiency, will append to these later
    L1,L2,L3=[],[],[]

    # Get length of S
    lenS=len(S)

    # q is a large prime that is roughly the same magnitude as N (the length of S)
    #q=997
    q=1500450271

    # Set up dictionary for the bases
    d={'A':0,'C':1,'G':2,'T':3}

    # Define the bases to use for x-point mutations
    bases=['A','C','G','T']

    # Initialize an empty dictionary to store hash values for each distinct kmer
    dhash=dict()

    # Compute powers of 4 in decreasing magnitude to use for hash operations, avoids having to recompute in each loop
    # i.e 4**(k-1), 4**(k-2),....,4**2, 4, 1
    power4=[pow(4,k-(i+1)) for i in range(0,k)]

    # Compute value to use in iterated hash operation, avoids having to recompute in each loop
    # Feed this value into the iteratedhash function
    bm=(4**k) % q

    # Use function initial hash to compute the hash of the first kmer in S
    # Compute the initial hash to be used in iterations of rolling hash process
    rhash=initialhash(S[:k],k,d,q,power4)

    # Since this is the first hash computed it won't be in the empty hash dictionary
    # Can assign a value for this hashkey and add it to the hash dictionary
    # 0th entry of dictionary value records frequency of k-mers
    # 1st, 2nd, 3rd ... entry of dictionary value records location of k-mer
    # Initial hash => frequency = 1, location = 0
    dhash[rhash] = [1,0]

    # Iterate over the remaining kmers in S
    for i in range(1,lenS-k+1):

        # Use the iteratedhash function to compute hash based on previous hash
        rhash=iteratedhash(S[i-1:k+i-1],S[i:k+i],d,bm,q,rhash)

        # Condition to check whether k-mer is in the dictionary
        # If kmer is in the dictionary
        if rhash in dhash:
            #Check if this is a collision
            #If this is a collision
            if S[i:k+i]!=S[dhash[rhash][1]:k+dhash[rhash][1]]:
                print('Hash Collision')
            #If it isn't a collision
            else:# +1 counter to record frequency of k-mer
                dhash[rhash][0] += 1
                # Add location of k-mer to the value
                dhash[rhash].append(i)

        # If k-mer is not in the dictionary
        else:
            # Set frequency to 1 and record location of k-mer
            dhash[rhash] = [1,i]


    # Iterate over keys in the hash dictionary
    for key in dhash:
        # Initialize sum to count number of times x-mutations occur for each kmer
        sum=0

        # Check if the key has appeared atleast f times
        if (dhash[key][0]>=f):
            # L1 contains kmers that have occured at least f times
            L1.append(S[dhash[key][1]:k+dhash[key][1]])
            # L2 contains the location of the frequently occuring kmers
            L2.append(dhash[key][1:])

            # Pick out the frequently occuring kmer
            kmer=S[dhash[key][1]:k+dhash[key][1]]

            # Iterate over the bases
            for base in bases:
                # If xth position of kmer is the same as base, continue to next base
                if (kmer[x]==base):
                    continue
                # If xth position of kmer is not the same as base
                else:
                    # Compute new hash function of a kmer with xth position replaced with the new base letter
                    newkey = (key - d[kmer[x]]*power4[x] + d[base]*power4[x]) % q

                    # If newkmer exists in dictionary, add frequency of newkmer to sum
                    if newkey in dhash:
                        sum = sum+dhash[newkey][0]

            # Append frequency of x-mutation of the current kmer to L3
            L3.append(sum)

    return L1,L2,L3


if __name__=='__main__':
    #Sample input and function call. Use/modify if/as needed
    infile = open('test_sequence.txt', 'r')
    S=infile.read()
    infile.close()
    #S='CCCTATGTTTGTGTGCACGGACACGTGCTAAAGCCAAATATTTTGCTAGT'
    k=10
    x=2
    f=2

    #L1,L2,L3=ksearch(S,k,f,x):

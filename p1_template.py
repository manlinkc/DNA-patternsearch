"""M345SC Homework 1, part 1
Manlin Chawla 01205586
"""

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

    # Take modulus of final value
    hashvalue = hashvalue % q

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

def ksearch(S,k,f,x):
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

    Analysis:
    My algorithm is based on the Rabin Karp Algorithm. The algorithm works by
    taking the first kmer in S and rewriting it in base 4 digits. To do this, I
    first set up a dictionary containing the gene bases and the corresponding digit.
    Setting up this dictionary has O(1) complexity.

    For the rolling hash I generated the hash by evaluating a K-1th order polynomial
    that converts the kmer from base 4 to base 10. To do this, I wrote a function
    called initialhash that computes first hash which will be used in the iterations
    of the rolling hash process. Computing the initial hash has a time complexity of O(K).

    To store each hash value, I initialised a separate dictionary at the beginning
    of the algorithm. As this is the first kmer it won’t be in the empty dictionary
    so this hash is added to the dictionary. Initializing the dictionary and storing
    the hash has O(1) complexity.

    When dealing with long strings S there can be a problem with hash collisions.
    To deal with the hash collisions I used a list of lists to build the key to
    store with hash. Each kmer is stored using its hash as the index, it is stored
    with a corresponding list of the form [frequency, location1, location2, ….].
    Python dictionary applies a hash function to the keys to know where to sort it.
    Since the computed hash is an integer the inbuilt hash function maps the integer
    to the same integer. This means there aren't extra computations involved in the
    storing step.

    Next, the algorithm iterates all of the kmers in S consecutively and computes
    the rolling hash by using iteration. I have written a function to do this called
    iteratedhash. Iterating across all kmers has complexity O(N-K) and the iteratedhash
    function consists of 4 operations so has O(1) complexity.

    Now the algorithm checks for hash collision. To do this it iterates across the
    list of lists to see if kmer is already being stored at that index. For each list
    inside the list the kmer is identified by picking the  location and evaluating
    the kmer at that location to the current kmer of interest. If the kmers match
    then the new location is appended to its list. If not, a new list containing
    [1,location] is appended to that index. Assuming that the hash table is well
    designed there will be minimal collision so in the best case can assume there
    is O(1) complexity for these steps.

    Next, the algorithm iterates of whole dictionary to find frequently occurring
    kmers and appends frequency and location to L1 and L2. These steps have complexity
    O(N-K) in the worst case. For each frequently occurring kmer, the hash of the
    x-mutation is easy to compute as the hash of the kmer has already been computed.
    Then this hash is looked up in the dictionary and if it exists the frequency
    is appended to L3. For looking up kmers to add to L1,L2, L3 and assuming there
    are minimal hash collisions, these steps have O(1) complexity.

    When k is large, calculations become to big for python to handle so I have used
    the modulo operator to deal with this. As there might be patterns in the gene sequence,
    I used a large prime to keep hash collisions to a minimum (probability of collision 1/P).
    I choose the size of the prime to be large enough to work for a wide range of
    values of k but I have hardwired other primes into the code that can be used.

    I also coded up an implementation that uses dictionaries only, I have added the
    code for this in the file in a function called ksearchdictionary.  An advantage
    of this approach is that it is simpler to implement, it works well for small
    values of K and produces the correct answer. As Python dictionaries are implemented
    as hash tables, Python works in the background to avoid hash collisions. A
    disadvantage of this approach is that as k increases a lot of memory is used,
    there is asignificant slowdown and the kernel crashes (crashes for k=100 on my
    laptop).On the other hand, the rolling hash approach is harder to implement and
    you have to explicitly check for hash collisions and need to create lists within
    lists to deal with this effectively. However, the advantage is that this approach
    works for a broader range of values of k especially works for large values of k
    that the dictionary implementation could not handle.

    b)The main contributions to the complexity come from the following steps:
    O(K): Computing a hash by using the initialhash function (i.e computing an k-1th polynomial). This is only used once.
    O(N-K): iterating over S
    O(1): computing hash using rolling hash iteration function
    O(1): Storing Hash, assuming prime is suitable enough for there to be minimal hash collisions
    O(K): matching the strings
    O(N-K): worst case for iterating over hash dictionary

    Therefore leading order time complexity becomes: O(K(N-K))
    """

    # Initialize lists for efficiency, will append to these later
    L1,L2,L3=[],[],[]

    # Get length of S
    lenS=len(S)

    # q is a large prime that is roughly the same magnitude as N (the length of S)
    # Can pick a suitable prime from the list below to use depending on size of k:
    #q=997 (3 digits)
    #q=12,239 (5 digits)
    q=1500450271
    #q=12764787846358441471 (20 digits)
    #q=115756986668303657898962467957 (30 digits)
    #q=5992830235524142758386850633773258681119 (40 digits)

    # Set up dictionary for the bases
    d={'A':0,'C':1,'G':2,'T':3}

    # Define the bases to use for x-point mutations
    bases=['A','C','G','T']

    # Initialize an empty dictionary to store hash values for each distinct kmer
    dhash=dict()

    # Compute powers of 4 in decreasing magnitude to use for hash operations, avoids having to recompute in each loop
    # i.e 4**(k-1), 4**(k-2),....,4**2, 4, 1
    power4=[pow(4,k-(i+1),q) for i in range(0,k)]

    # Compute value to use in iterated hash operation, avoids having to recompute in each loop
    # Feed this value into the iteratedhash function
    bm=pow(4,k,q)

    # Use function initial hash to compute the hash of the first kmer in S
    # Compute the initial hash to be used in iterations of rolling hash process
    rhash=initialhash(S[:k],k,d,q,power4)

    # Since this is the first hash computed it won't be in the empty hash dictionary
    # Can assign a value for this hashkey and add it to the hash dictionary
    # 0th entry of dictionary value records frequency of k-mers
    # 1st, 2nd, 3rd ... entry of dictionary value records location of k-mer
    # Initial hash => frequency = 1, location = 0
    dhash[rhash] = [[1,0]]

    # Iterate over the remaining kmers in S
    for i in range(1,lenS-k+1):

        # Use the iteratedhash function to compute hash based on previous hash
        rhash=iteratedhash(S[i-1:k+i-1],S[i:k+i],d,bm,q,rhash)

        # Condition to check whether k-mer is in the dictionary
        # If kmer is in the dictionary
        if rhash in dhash:

            # Condition to check if there is a hash collision
            if S[i:k+i]==S[dhash[rhash][0][1]:k+dhash[rhash][0][1]]:

            # If there isn't a collision
            # +1 counter to record frequency of k-mer
                dhash[rhash][0][0] += 1
                # Add location of k-mer to the value
                dhash[rhash][0].append(i)

            # If there is a collision
            else:
                # Switch will change if k-mer is already stored at this index
                hasappeared='no'
                # Iterate over lists alrady stored at this index
                for j in range(len(dhash[rhash])):
                    # If the string appears then
                    if S[i:k+i]==S[dhash[rhash][j][1]:k+dhash[rhash][j][1]]:
                        # Add one to the frequency value in it's list
                        dhash[rhash][j][0] += 1
                        # Append the new location to it's list
                        dhash[rhash][j].append(i)
                        # Flip the switch
                        hasappeared='yes'
                        # Already found so won't appear again, can break loop
                        break

                # If the kmer doesn't already exist at this index
                if hasappeared=='no':
                    # Store a corresponding list at this index
                    dhash[rhash].append([1,i])

        # If k-mer is not in the dictionary
        else:
            # Set frequency to 1 and record location of k-mer
            dhash[rhash] = [[1,i]]

    # Iterate over keys in the hash dictionary
    for key in dhash:
        # Iterate over the list of lists
        for j in range(len(dhash[rhash])):

            # Check if the key has appeared atleast f times
            if (dhash[key][j][0]>=f):

                # L1 contains kmers that have occured at least f times
                L1.append(S[dhash[key][j][1]:k+dhash[key][j][1]])
                # L2 contains the location of the frequently occuring kmers
                L2.append(dhash[key][j][1:])

                # Pick out the frequently occuring kmer
                kmer=S[dhash[key][j][1]:k+dhash[key][j][1]]
                # Initialize sum to count number of times x-mutations occur for each kmer
                sum=0
                # Iterate over the bases
                for base in bases:
                    # If xth position of kmer is the same as base, continue to next base
                    if (kmer[x]==base):
                        continue
                    # If xth position of kmer is not the same as base
                    else:
                        newkmer=kmer[:x]+base+kmer[x+1:]
                        # Compute new hash function of a kmer with xth position replaced with the new base letter
                        newkey = (key - d[kmer[x]]*power4[x] + d[base]*power4[x]) % q

                        # If newkmer exists in dictionary, add frequency of newkmer to sum
                        if newkey in dhash:
                            # Iterate through list of lists
                            for j in range(len(dhash[newkey])):
                                # Find the list that matches to the x-mutated kmer
                                if newkmer==S[dhash[newkey][j][1]:k+dhash[newkey][j][1]]:
                                    # Add frequency to the sum
                                    sum = sum + dhash[newkey][j][0]
                # Append final frequency of x-mutation of the current kmer to L3
                L3.append(sum)

    return L1,L2,L3

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
    x=2
    f=2

    #L1,L2,L3=ksearch(S,k,f,x)

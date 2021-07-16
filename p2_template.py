"""M345SC Homework 1, part 2
Manlin Chawla 01205586
"""
import random
import numpy as np
import matplotlib.pyplot as plt
import time

def listgenerator(P,M,N):
    """Function that generates a list containing N partially sorted lists of size
    M. The first P elements of each list is sorted into ascending order, the remaining
    M-P elements are random and unsorted.
    Inputs:
    P: Length of sorted list
    M: Length of entire lists
    N: Number of lists

    a)	My approach to an algorithm for nsearch works as follows:

    Nsearch takes a list L as its input. To generate this list I wrote a function
    called listgenerator that generates N partially sorted lists of length M as its
    output. This list L is used as an input for nsearch.
    For the nsearch algorithm, I first initialized an empty list for Lout to which
    I later append sublists of the form [list, position] for each occurrence the
    target is found. I also computed the value of N which is the number of lists
    of length M contained within L. These steps have a time complexity O(1).I then
    iterate over all N lists, this step has a time complexity of O(N). As the first
    P elements of each list are sorted, I have set an if condition to check the
    first value and the last value of the sorted sublist. If the target is in the
    range of the lowest value and the highest value, then the algorithm carries
    out the binary search. If the binary search is carried out and the target appears
    in the sorted list then the algorithm continues to search around either side
    to find the range of indexes where the target appears. This needed as binary
    search only returns one index where the target is found in the sorted list.
    To find all occurrences I wrote a function called searchrange. Searchrange works
    by using the index returned by the binarysearch as a starting point. It then uses
    linear search to check element by element to find the first elements to the
    left and to the right of the index where the target doesn’t appear. There are
    some boundary conditions and I have used if statements to deal with that which
    add a complexity of O(1).If target is not found by initial binary search a search
    around is not needed. I have a included an if statement that first checks if the
    binary check has returned anything, if it has then the searchrange function is
    used. If it hasn’t then the algorithm proceeds to checking the unsorted portion
    of the list. This if statement step is has O(1) complexity. The algorithm then
    proceeds to searching through the unsorted portion of the list. For the unsorted
    portion of the list I used a linear search to iterate through each element of
    the list and check if it matches against the target. As the unsorted list has
    length M-P, the time complexity of the linear search is O(M-P).  If the value
    matches the target then the list number and position are appended to the Lout,
     which is an O(1) operation.

     b) For the nsearch algorithm the running time depends mainly on the following steps:

     O(1): Initializing Lout and finding the value of N
     O(N): Iterating over the N lists.
     O(log(P)): Within this for loop, binary search on sorted portion of the list
     O(X):  If target is located by binary search, need to use searchrange function.
            This involves doing a linear search on either side and the number of
            searches is the same as the number of repeat appearances. X=no of searches carried out.
    X*O(1): Appending the range of locations found by the searchrange function.
    O(M-P): Linear search over the unsorted portion of the list

    Overall the time complexity of nsearch is of the form Time Complexity = O(N(log(P)*x+ x + M-P)).
    In the best case the target will only appear once in the sorted list so the best
    case time complexity becomes O(N(log(P)+M-P)). In the worst case the whole sorted
    list consists of the target so the worst case time complexity becomes
    O(N(P*log(P)+P+M-P))=O(N(P*log(P)+M)).

    c)I have tried to make my algorithm efficient by taking the following steps:

    For the sorted portion of the list, as this list is sorted into ascending order,
    I used this fact to improve the efficiency of my algorithm. I first have an if
    statement condition that checks the first value and last value of the partially
    sorted list of length P. As the sublist is in ascending order, if the last value
    of the sublist is smaller than the target we are looking for then it follows
    that the target cannot appear in the sublist. Likewise, if the first value is
    larger than the target we are looking for then it follows that the target cannot
    appear in the target in the sublist. If this is the case, then implementing the
    binary search on the partially sorted list is useless so the algorithm skips
    to looking for the target in the unsorted portion of the list.

    If the binary search step is carried and the target is found then I have setup
    a function searchrange to search around. I have used linear search for this.
    For a list of length P searching around using linear search has a worst case
    time complexity of O(P). However, as the lists are lists are randomly generated
    it is likely that there will only be a few repeats of the target. Assuming
    there are a low number of repetitions then then the linear search reduces
    to only checking a few values and has close to O(1)/small time complexity.

    e)
    Figure 1:
    This is a plot of the runtime of nsearch as P (the length of the sorted portion
    of each list) is varied whilst keeping M and N fixed. The parameters I have used
    are M=10000, N=5 and target=442. This plot shows that as P is increased the
    runtime decrease and the relationship is roughly linear. There is more variation
    in the runtimes for smaller values of P and less variation for larger values of P.
    P corresponds to the length of the sorted section of the list. Keeping M fixed,
    as the input variable P is increased, a higher percentage of the entire list is
    sorted and lower percentage is unsorted. As P increases a larger part of the
    operations of the algorithm consists of the binary search and lower part of
    the operations consists of the linear search. As binary search has complexity
    O(log(P)), this explains that as P increases there is logarithmic growth in the
    runtime due to the growing size of the sorted list. However, the linear search
    has a complexity of O(M-P). As P increases this complexity M-P becomes smaller
    and smaller, asymptotically the complexity tends to O(1). The decrease in the
    complexity of the unsorted portion of the list is faster than the logarithmic
    increase in the complexity corresponding to the sorted portion of the list. This
    helps to explain the negative trend of runtime as P increases. The last point in
    this plot shows when P=M meaning the entire list is sorted and there is no linear
    search, this point has the lowest runtime.

    Figure 2:
    This is a plot of the runtime of nsearch as M (the length of each list within L)
    is varied whilst keeping P and N fixed. The parameters I have used are P=100, N=5
    and target=442. This plot shows that as M is increased the runtime increases
    the relationship is roughly linear. There is some variation in the runtimes
    for larger values of M and and less variation for smaller values of P. M corresponds
    to the length of the list within L. Keeping P fixed, as M is increased, a higher
    percentage of the entire list is unsorted and a lower percentage is sorted. As
    M increases a larger part of the operations of the algorithm consists of the
    linear search and lower part of the operations consists of the binary search.
    The algorithm spends more time completing the linear search operations and
    less time on the binary search operation. Linear search has a complexity of
    O(M-P) in this algorithm. As M is increased asymptotically this complexity
    goes to O(M) so there is a linear increase in the time taken on linear search
    as expected. Binary search has complexity O(log(P)) and as M increases this
    complexity contributes a smaller amount to total runtime relative to the contribution
    of the linear search. This helps to explain the positive trend of runtime as M increases.
    The first point in this plot shows when P=M meaning the entire list is sorted and there is no
    linear search, this point has the lowest runtime.

    Figure 3:
    This is a plot of the runtime of nsearch as N (the number of lists within L)
    is varied whilst keeping P and M fixed. The parameters I have used are P=100,
    M=1000 and target=442. This plot shows that as N is increased the runtime increases.
    The plot shows that the relationship between N and the runtime is roughly linear.
    There is some variation in the runtimes but in the most part the data fits a straight
    line. N, the number of lists is independently set with respect to M and N. The plot
    shows that a N increases as expected the runtime increases as there are more lists
    to search for and this increase is linear.

    Figure 4:
    This is a plot showing the runtime as N is increases for different ratios of P/M.
    To generate this plot I have kept M fixed and plotted five different lines corresponding
    to where the value P is 0%, 25%,50%,75% and 100% of the value M. The parameters I
    have used are M=1000 and target=442.The plot shows that when P is 100% the value
    of M this line has the lowest runtime as N increases. As the percentage P relative
    to M increases the runtime of each line also increases higher for all values of N.
    The plot shows that when the value of P is 0% the value of M then this line has
    the highest runtime as N increases.  There is some fluctuations in runtime within
    all lines but there is generally a linear positive relationship between N and runtime.
    When the value of P is 100% the value of M then this means that P=M and the entire
    list is sorted. The algorithm only uses binary search to find the target in this
    fully sorted list. This line acts as a baseline for comparison, it has the lowest
    runtime compared to other lines as there is no component of linear search involved
    so complexity is just O(log(M)). On the other end when the value of P is 0% the value
    of M then this means that the entire list is unsorted and the algorithm only uses
    linear search. This line has the highest run time as there is no component of binary
    search involved so time complexity is just O(M). For large problem sizes O(M) is
    generally larger than O(log(M)) so this helps to explain why decreasing the sorted
    portion of the list increases the runtime in this plot.

    Figure 5:
    This is a plot showing the runtime as M increases for different ratios of P/M.
    To generate this plot I have kept N fixed and plotted five different lines corresponding
    to where the value P is 0%, 25%,50%,75% and 100% of the value M. The parameters I have
    used are N=5 and target=442. This plot shows that as M increases, as P is just a
    fixed proportion of M this means P also increases. This explains the positive relationship
    between runtime and M increasing. Similar to the previous plot when P is a 100% of the
    value of M the entire list is sorted and only binary search is used. The complexity
    of this is lower that the complexity of linear search. The difference in complexity
    explains why when P is 0% the value of M then the entire list is unsorted and only
    linear search is used and hence this line has the highest runtime as M is increased.
    As the ratio increases all of the lines follow this trend with higher runtimes as the
    ratio of P/M decreases.

    """
    # Preallocate for efficienyc
    L = np.zeros(N).tolist()

    # Iterate over N to make each list separately
    for i in range(N):
        # Generate a list of length P containing random integers
        listP = [random.randint(0,999) for i in range(P)]
        # Sort Plist into ascending order
        listP = np.sort(listP).tolist()
        # Generate a list of length M-P containing random integers
        listM = [random.randint(0,999) for i in range(M-P)]
        # Concatenate to form a list of length P
        L[i] = listP+listM

    return L

def binarysearch(sublist, target):
    #Set an arbitrary negative index, if target is found index will change
    index=-1

    # Set start and end points
    start, end = 0, len(sublist)-1

    # Set up condition
    while start <= end:
        # Compute mid value
        mid = int(0.5*(start + end))

        # Check conditions
        if target == sublist[mid]:
            index=mid
            break
        elif target < sublist[mid]:
            end = mid - 1
        else:
            start = mid + 1

    return index

def searchrange(sublist,target,index):
    # Intialize indexes for searching around
    firstindex=index
    lastindex=index

    # Search around to the left
    while target==sublist[firstindex]:
        # Check boundary point, can't -1 as index is not in list range
        if firstindex==0:
            break
        # Iterate to the left
        else:
            firstindex=firstindex-1

    # Search around to the right
    while target==sublist[lastindex]:
        # Check boundary point, can't +1 as index is not in list range
        if lastindex==(len(sublist)-1):
            break
        else:
            lastindex=lastindex+1

    # Final check to compensate for boundary cases
    # Check right index
    if target!=sublist[lastindex]:
        lastindex=lastindex-1
    # Check left index
    if target!=sublist[firstindex]:
        firstindex=firstindex+1

    return firstindex,lastindex

def nsearch(L,P,target):
    """Input:
    L: list containing *N* sub-lists of length M. Each sub-list
    contains M numbers (floats or ints), and the first P elements
    of each sub-list can be assumed to have been sorted in
    ascending order (assume that P<M). L[i][:p] contains the sorted elements in
    the i+1th sub-list of L
    P: The first P elements in each sub-list of L are assumed
    to be sorted in ascending order
    target: the number to be searched for in L

    Output:
    Lout: A list consisting of Q 2-element sub-lists where Q is the number of
    times target occurs in L. Each sub-list should contain 1) the index of
    the sublist of L where target was found and 2) the index within the sublist
    where the target was found. So, Lout = [[0,5],[0,6],[1,3]] indicates
    that the target can be found at L[0][5],L[0][6],L[1][3]. If target
    is not found in L, simply return an empty list (as in the code below)
    """
    # Initalize output for efficiency
    Lout=[]

    # Find number of lists to iterate over
    N=len(L)

    # Iterate over each list
    for i in range(N):
        # Checking if target could possibly appear in P List
        if L[i][P-1] >= target >= L[i][0] :
            # Apply binary search, uses function I wrote which has code included above
            index=binarysearch(L[i][:P],target)
            # If target is found in list by binary search
            if index != -1:
                # Search around to the left and right
                firstindex,lastindex=searchrange(L[i][:P],target,index)
                # If target only appears once
                if firstindex==lastindex:
                    # Append to Lout
                    Lout.append([i,firstindex])
                # If target appears multiple times then
                else:
                    # For entire range of index target appears
                    for k in range(firstindex,lastindex+1):
                        # Append to Lout
                        Lout.append([i,k])

        # Checking unsorted M-P list
        for y,x in enumerate(L[i][P:]):
            # Use linear search to check element by element
            if target==x:
                # If target is found append to list
                Lout.append([i,y+P])

    return Lout

def nsearch_time(figurenum, display=False):
    """Analyze the running time of nsearch.
    Add input/output as needed, add a call to this function below to generate
    the figures you are submitting with your codes.

    Discussion: (add your discussion here)

    """
    # Plot of runtime as P (length of sorted portion of each list) is varied whilst keeping M and N fixed
    # Plot is a linear line going downwards
    if figurenum==1:
        # Set parameters M, N, number of points to be plotted and target to be searched for
        M=10000
        N=5
        numpoints=100
        target=442

        # Preallocate for efficiency, list to contain the runtimes
        runtimevalues=np.zeros(numpoints)

        # Generate a list containing values of P to iterate over
        Pvalues=np.linspace(10,10000,numpoints,dtype=int)

        for i,P in enumerate(Pvalues):
            # Generate L using the function listgenerator function with M,N as defined above and the value of P in current iteration
            L=listgenerator(P,M,N)

            t1=time.time()
            nsearch(L,P,target)
            t2=time.time()
            runtimevalues[i]=t2-t1

        # Plotting graph and format
        plt.hold(True)
        plt.plot(Pvalues,runtimevalues,'.',c='r')
        plt.xlabel('P')
        plt.ylabel('Runtime')
        plt.title('Manlin Chawla: nsearch_time(1) \n Runtime of nsearch for different values of P (keeping M, N fixed)')
        plt.legend()
        plt.hold(False)

        if display==True:
            plt.show()
    #-----------------------------------------------------------------------------------------------------------------------------
    # Plot of runtime as M (length of each list within L) is varied whilst keeping P and N fixed
    # Plot linear line going upwards
    if figurenum==2:
        # Set parameters P, N, number of points to be plotted and target to be searched for
        P=100
        N=5
        numpoints=100
        target=442

        # Preallocate for efficiency, list to contain the runtimes
        runtimevalues=np.zeros(numpoints)

        # Generate a list containing values of M to iterate over
        Mvalues=np.linspace(100,10000,numpoints,dtype=int)

        for i,M in enumerate(Mvalues):
            # Generate L using the function listgenerator function with P,N as defined above and the value of M in current iteration
            L=listgenerator(P,M,N)

            t1=time.time()
            nsearch(L,P,target)
            t2=time.time()
            runtimevalues[i]=t2-t1

        # Plotting graph and format
        plt.hold(True)
        plt.plot(Mvalues,runtimevalues,'.',c='b')
        plt.xlabel('M')
        plt.ylabel('Runtime')
        plt.title('Manlin Chawla: nsearch_time(2) \n Runtime of nsearch for different values of M (keeping P, N fixed)')
        plt.legend()
        plt.hold(False)

        if display==True:
            plt.show()
    #---------------------------------------------------------------------------------------------------------------------------------------
    # Plot of runtime as N (the number of lists within L) is varied whilst keeping P and M fixed
    # Plot description
    if figurenum==3:
        # Set parameters P, N, number of points to be plotted and target to be searched for
        P=100
        M=1000
        numpoints=100
        target=442

        # Preallocate for efficiency, list to contain the runtimes
        runtimevalues=np.zeros(numpoints)

        # Generate a list containing values of N to iterate over
        Nvalues=np.linspace(0,1000,numpoints,dtype=int)
        Nvalues[0]=1

        for i,N in enumerate(Nvalues):
            # Generate L using the function listgenerator function with P,M as defined above and the value of N in current iteration
            L=listgenerator(P,M,N)

            # Time function
            t1=time.time()
            nsearch(L,P,target)
            t2=time.time()
            runtimevalues[i]=t2-t1

        # Plotting graph and format
        plt.hold(True)
        plt.plot(Nvalues,runtimevalues,'.',c='k')
        plt.xlabel('N')
        plt.ylabel('Runtime')
        plt.title('Manlin Chawla: nsearch_time(3) \n Runtime of nsearch for different values of N (keeping P, M fixed)')
        plt.legend()
        plt.hold(False)

        if display==True:
            plt.show()
    #--------------------------------------------------------------------------------------------------------------------------------

    # Plot of runtime as N (the number of lists within L) is varied whilst varying the portion of list that is sorted,
    #ie. (different raios of P/M)
    # Plot description
    if figurenum==4:
        # Set parameters P, N, number of points to be plotted and target to be searched for
        M=1000
        numpoints=100
        target=442

        # Preallocate for efficiency, list to contain the runtimes
        Pvalues=np.linspace(0,M,5,dtype=int)
        runtimevalues=np.zeros(numpoints)

        # Generate a list containing values of N to iterate over
        Nvalues=np.linspace(0,1000,numpoints,dtype=int)
        Nvalues[0]=1

        plt.hold(True)

        for P in Pvalues:
            for i, N in enumerate(Nvalues):
                    L=listgenerator(P,M,N)
                    t1=time.time()
                    nsearch(L,P,target)
                    t2=time.time()
                    runtimevalues[i]=t2-t1

            plt.plot(Nvalues,runtimevalues,label='P='+str(P/10)+'% of M')

        plt.xlabel('N')
        plt.ylabel('Runtime')
        plt.title('Manlin Chawla: nsearch_time(4) \n Runtime of nsearch for different values of N and different percentages of sorted list P')
        plt.legend()
        plt.hold(False)

        if display==True:
            plt.show()

    #--------------------------------------------------------------------------------------------------------------------------------

    # Plot of runtime as M (the size of lists) is varied whilst varying the portion of list that is sorted,
    #ie. (different raios of P/M)
    # Plot description
    if figurenum==5:
        # Set parameters P, N, number of points to be plotted and target to be searched for
        N=5

        numpoints=100
        target=442

        # Preallocate for efficiency, list to contain the runtimes
        Mvalues=np.linspace(100,1000,numpoints,dtype=int)
        runtimevalues=np.zeros(numpoints)

        Pratio=np.linspace(0,1,5)
        plt.hold(True)

        for ratio in Pratio:
            for i,M in enumerate(Mvalues):
                P=int(M*ratio)
                L=listgenerator(P,M,N)
                t1=time.time()
                nsearch(L,P,target)
                t2=time.time()
                runtimevalues[i]=t2-t1

            plt.plot(Mvalues,runtimevalues,label='P='+str(P/10)+'% of M')

        plt.xlabel('M')
        plt.ylabel('Runtime')
        plt.title('Manlin Chawla: nsearch_time(5) \n Runtime of nsearch for different values of M and different percentages of sorted list P')
        plt.legend()
        plt.hold(False)

        if display==True:
            plt.show()


    return None #Modify as needed


if __name__ == '__main__':
    #add call(s) to nsearch here which generate the figure(s) you are submitting

    output_a = nsearch_time(1)
    plt.savefig('fig1FINAL.png', bbox_inches="tight")
    plt.clf()
    print('nsearch_time(1) plot saved')

    output_a = nsearch_time(2)
    plt.savefig('fig2FINAL.png', bbox_inches="tight")
    plt.clf()
    print('nsearch_time(2) plot saved')

    output_a = nsearch_time(3)
    plt.savefig('fig3FINAL.png', bbox_inches="tight")
    plt.clf()
    print('nsearch_time(3) plot saved')

    output_a = nsearch_time(4)
    plt.savefig('fig4FINAL.png', bbox_inches="tight")
    plt.clf()
    print('nsearch_time(4) plot saved')

    output_a = nsearch_time(5)
    plt.savefig('fig5FINAL.png', bbox_inches="tight")
    plt.clf()
    print('nsearch_time(5) plot saved')

import os.path
import time
from multiprocessing import Pool as pool
import numpy as np

increasing_boolean_functions = {}

k_max = 7
num_processors = 7

def duplicate_list(L):
    """
    This function takes as input a two-dimensional boolean array of size N x n,   n = 2^k
    representing a list of N boolean functions on {0, 1}^k.
    Returns the list obtained by concatenating all vectors of the shape (l1, l2)
    where l1, l2 are rows of L, with l1 <= l2 for the lexicographical order.
    """
    N, n = L.shape
    Lsorted = np.take(L, np.lexsort(np.flip(np.transpose(L), axis=0)), axis=0)
    # Lsorted is L sorted lexicographically so that l1 <= l2 if and only if the index of l1 is at most the index of l2
    Lsquare = np.repeat(Lsorted, repeats=N, axis=0).reshape((N, N, n))
    Lsquare_rev = np.transpose(Lsquare, (1, 0, 2))
    Lconcat = np.concatenate((Lsquare, Lsquare_rev), axis=2)
    # now Lconcat[i, j] is  (l1, l2) where l1 is of index i and l2 is of index j (in L sorted)
    # We extract the flattened upper-triangular part of this matrix
    Lextracted = Lconcat[np.triu_indices(N)]
    return Lextracted

def comparables(n):
    """
    Returns the list of pairs (x, y) of vectors of {0, 1}^n such that x != y and y <= y for the dominant order,
    and also x[0] = 0 and y[0] = 1.
    """
    if n == 1:
        return np.array([[[0],[0]],[[0],[1]],[[1],[1]]], dtype = np.bool)
    res = []
    previous_comparables = comparables(n-1)
    for l1, l2 in previous_comparables:
        res.append([np.append(l1, [np.False_]), np.append(l2, [np.False_])])
        res.append([np.append(l1, [np.False_]), np.append(l2, [np.True_])])
        res.append([np.append(l1, [np.True_]), np.append(l2, [np.True_])])
        if np.sum(l1) < np.sum(l2):
            res.append([np.append(l1, [np.True_]), np.append(l2, [np.False_])])
    res2 = [u for u in res if not np.all(u[0] == u[1]) and not u[0][0] == u[1][0]] 
    return np.array(res2, dtype = np.bool)

def vlist_to_int(vlist):
    """
    Takes as input vlist a list of vectors of {0, 1}^n.
    Returns the list of integers represented by these vectors in base 2,
    which for a single vector v is the sum of 2^{n-1-i} * v[i] for i in 0, ..., n-1.
    """
    N, n = vlist.shape
    l = np.array(vlist.transpose(), dtype = np.uint16)
    l = np.diag(np.array([2**(n-1-k) for k in range(n)], dtype = np.uint16)) @ l
    return np.sum(l, axis = 0)

def eval_funlist_on_vlist(flist, vlist):
    """
    Takes as input :
      flist, a list of boolean functions (a boolean array N x 2^n)
      vlist, a list of vectors (a boolean array M x n)
    Returns a N x M array of boolean, representing the values F(v) for F in flist, and v in vlist.
    """
    vlist_int = vlist_to_int(vlist)
    return np.take(flist, vlist_int, axis=1)
    

def list_increasing_lorenz_boolean_functions(N):
    """
    This function computes the list of all boolean functions which are increasing for the dominant order.
    It stores the result into increasing_boolean_functions[N].
    """
    global increasing_boolean_functions, compar

    mem_max = 4 * 1024**3

    # If it has already been computed, returns the cached list.
    if N in increasing_boolean_functions:
        return increasing_boolean_functions[N]
        
    if N==0:
        increasing_boolean_functions[N] = np.array([[0],[1]], dtype = np.bool)
        return increasing_boolean_functions[N]

    # We first compute the (large) list of all fonctions F obtained by concatenating monotone functions of dimension N-1.
    # This is because if F in increasing on {0, 1}^N, then certainly F(0, v) and F(1, v) are increasing on {0, 1}^(N-1).
    # This also means that we only need to compare values F(0, v) with F(1, w) for various v, w in order to assess
    # whether F is increasing (whence the condition x[0] = 0 and y[0] = 1 in comparables.
    L = duplicate_list(list_increasing_lorenz_boolean_functions(N-1))

    # We compute the list of all comparable vectors of dimension N.
    compar = comparables(N).transpose((1, 0, 2))
    v1list = compar[0]
    v2list = compar[1]

    i = 0
    # Next we sieve the boolean functions L to keep only those F which are increasing, meaning F[v1] <= F[v2]
    # for every pair (v1, v2) in v1list x v2list.
    # To keep memory usage under control, we compare in batches of size at most mem_max.
    ni = mem_max // v1list.shape[0]
    f_total = np.zeros((0, 2**N), dtype = np.bool)
    while i<L.shape[0]:
        Lwork = L[i:i+ni]
        # We compute the values of functions in Lwork in v1list and then v2list.
        fval1 = eval_funlist_on_vlist(Lwork, v1list)
        fval2 = eval_funlist_on_vlist(Lwork, v2list)
        # We compare them for each v.
        fcompare = fval1 <= fval2
        # Next we keep only those F for which fcompare[F, v] is true for every v.
        fsieved = np.all(fcompare, axis=1)
        # Finally we extract from Lwork those functions detected by fcrible.
        f_total = np.append(f_total, Lwork[fsieved], axis=0)
        i += ni
    
    if N not in increasing_boolean_functions:
        increasing_boolean_functions[N] = f_total
    return increasing_boolean_functions[N]
    
    

def write_Fk(k):
    """
    Stores the list of increasing boolean functions of dimension k.
    """
    with open(f"BooleanFunctions/data/Fk-k={k}", "w") as f:
        for l in increasing_boolean_functions[k]:
            f.write("".join([str(int(i)) for i in l]) + "\n")

def read_Fk(k):
    """
    Recovers the list of increasing boolean functions of dimension k.
    """
    with open(f"BooleanFunctions/data/Fk-k={k}", "r") as f:
        s = f.read()
        Ls = s.split("\n")
        t = [[int(i) for i in x] for x in Ls[:-1]]
        increasing_boolean_functions[k] = np.array(t, dtype = np.bool).reshape((len(t), 2**k))


# Compute the list of all functions in {\mathcal F}_k for 1\leq k\leq 7.
for k in range(1, k_max+1):
    if os.path.exists(f"BooleanFunctions/data/Fk-k={k}"):
        read_Fk(k)
    else:
        list_increasing_lorenz_boolean_functions(k)
        write_Fk(k)


# Dictionary of all increasing boolean functions, excluding those who are zero,
# and under the conditition that they are 0 at vectors with only 1 non-zero component.
increasing_boolean_functions_reduced = {}

for k in increasing_boolean_functions:
    increasing_boolean_functions_reduced[k] = []
    for F in increasing_boolean_functions[k]:
        if np.any(F) and F[2**(k-1)] == 0:
            increasing_boolean_functions_reduced[k].append(F)



def successors(F, v):
    """
    Returns the array of all vectors w successors of v (distinct from v) such that F(w) = 0.
    This means that w is obtained from v by switching one of its components from 0 to 1.
    """
    n = v.shape[0]
    # First we build all vectors which are successors of v : we build copies of v...
    v_repeated = np.repeat(v.reshape((1, n)), repeats=n, axis=0)
    # ... to which we add a 1 in a coordinates...
    v_succ = np.logical_or(np.identity(n, dtype=np.bool), v_repeated)
    not_v = np.logical_not(np.all(v_repeated == v_succ, axis=1))
    F_vals = eval_funlist_on_vlist(F.reshape((1, F.shape[0])), v_succ)[0]
    return np.unique(v_succ[np.logical_and(np.logical_not(F_vals), not_v)], axis=0)

def predecessors(F, v):
    """
    Returns the arry of all vectors w predecessors of v (distinct from w) such that F(w) = 1.
    This means that w is obtained from v by switching one of its components from 1 to 0.
    """
    n = v.shape[0]
    v_repeated = np.repeat(v.reshape((1, n)), repeats=n, axis=0)
    v_pred = np.logical_and(np.logical_not(np.identity(n, dtype=np.bool)), v_repeated)
    not_v = np.logical_not(np.all(v_repeated == v_pred, axis=1))
    F_vals = eval_funlist_on_vlist(F.reshape((1, F.shape[0])), v_pred)[0]
    return np.unique(v_pred[np.logical_and(F_vals, not_v)], axis=0)

# Define a maximal point of F to be those v for which F(v) = 0 but F(w) = 1 for all successor w > v.
# Simialrly a minimal point v of F satisfies F(v) = 1 but F(w) = 0 for all predecessor w < v.

def points_maximaux_aux(F, vlist, n):
    """
    Auxiliary recursive function. Computes all maximal points among those vectors which are successors to vectors in vlist.
    """
    res = np.empty((0, n), dtype=np.bool)
    if len(vlist) == 0:
        return res
    vlist_aux = np.zeros(shape=(0, n), dtype=np.bool)
    for v in vlist:
        S = successors(F, v)
        if len(S) == 0:
            res = np.concatenate((res, v.reshape((1, n))), axis=0)   
        else:
            vlist_aux = np.concatenate((vlist_aux, S), axis=0)
    vlist_aux = np.unique(vlist_aux, axis=0)
    return np.unique(np.concatenate((res,
                                     points_maximaux_aux(F, vlist_aux, n)),
                                    axis=0), axis=0)

def points_minimaux_aux(F, vlist, n):
    """
    Auxiliary recursive function. Computes all minimal points among those vectors which are predecessors of vectors in vlist.
    """
    res = np.empty((0, n), dtype=np.bool)
    if vlist.shape[0] == 0:
        return res
    vlist_aux = np.zeros(shape=(0, n), dtype=np.bool)
    for v in vlist:
        S = predecessors(F, v)
        if S.shape[0] == 0:
            res = np.concatenate((res, v.reshape((1, n))), axis=0)
        else:
            vlist_aux = np.concatenate((vlist_aux, S), axis=0)
    vlist_aux = np.unique(vlist_aux, axis=0)
    return np.unique(np.concatenate((res,
                                     points_minimaux_aux(F, vlist_aux, n)),
                                    axis=0), axis=0)

def points_maximaux(F, n):
    """
    Returns all maximal vectors of F.
    A maximal vector v for F is one for which F(v) = 0 but F(w) = 1 for all w successor of v.
    """
    if np.all(F) == np.True_:
        return np.empty((0, n), dtype=np.bool)
    return points_maximaux_aux(F, np.zeros((1,n), dtype=np.bool), n)


def points_minimaux(F, n):
    """
    Returns all minimal vectors of F.
    A minimal vector v for F is one for which F(v) = 1 but F(w) = 0 for all w predecessor of v.
    """
    if np.any(F) == np.False_:
        return np.empty((0, n), dtype=np.bool)
    return points_minimaux_aux(F, np.ones((1, n), dtype=np.bool), n)


minimal_points = {}
maximal_points = {}


def work_mm(u):
    """
    This function is an auxiliary function for use in multiprocessing.
    It takes as input a tuple (F, k) with F a boolean function of dimension k (meaning a numpy boolean array of size 2^k).
    Returns the pair (maxi, mini) formed by the maximal, respectively minimal points for F.
    """
    F, k = u
    return (points_maximaux(F, k), points_minimaux(F, k))

def calcule_mm(k):
    """
    Computes the maximal and minimal points for all increasing boolean functions of dimension k.
    Stores the result in the k component of the global variable maximal_points (resp. minimal_points).
    """
    global maximal_points, minimal_points
    maximal_points[k] = []
    minimal_points[k] = []
    Lk = [(F, k) for F in increasing_boolean_functions_reduced[k]]
    with pool(num_processors) as p:
        for r in p.imap(work_mm, Lk):
            maximal_points[k].append(r[0])
            minimal_points[k].append(r[1])

def write_mm(k):
    """
    Stores the minimal and maximal points for boolean functions of dimension k in appropriate files.
    """
    with open(f"BooleanFunctions/data/mini-k={k}", "w") as f:
        for l in minimal_points[k]:
            f.write(",".join(["".join([str(int(i)) for i in v]) for v in l]) + "\n")
    with open(f"BooleanFunctions/data/maxi-k={k}", "w") as f:
        for l in maximal_points[k]:
            f.write(",".join(["".join([str(int(i)) for i in v]) for v in l]) + "\n")
            
def read_mm(k):
    """
    Retrieves the minimal and maximal points for boolean functions of dimension k from the appropriate files.
    """
    with open(f"BooleanFunctions/data/mini-k={k}", "r") as f:
        s = f.read()
        Ls = s.split("\n")
        minimal_points[k] = []
        for l in Ls[:-1]:
            t = l.split(",") if l != "" else []
            minimal_points[k].append(np.array([[int(i) for i in s] for s in t], dtype=np.bool).reshape((len(t), k)))
    with open(f"BooleanFunctions/data/maxi-k={k}", "r") as f:
        s = f.read()
        Ls = s.split("\n")
        maximal_points[k] = []
        for l in Ls[:-1]:
            t = l.split(",") if l != "" else []
            maximal_points[k].append(np.array([[int(i) for i in s] for s in t], dtype=np.bool).reshape((len(t), k)))

def est_sup(a, b):
    """
    Returns True if and only if a and b have the same weight, but partial sums of b are at least as large as partial sums of a.
    """
    return np.all(np.cumsum(b.astype(np.int16)-a.astype(np.int16)) >= np.zeros(len(b))) and np.sum(b) == np.sum(a)

def elague_mm(k):
    """
    Removes from the list of minimal and maximal points computed so far, all those whose minimality and maximality
    can be deduced by monotonicity from the same properties of other points.
    """
    global maximal_points, minimal_points
    L = maximal_points[k]
    M = []
    for lpts in L:
        rayes = []
        for i in range(len(lpts)):
            for j in range(i+1, len(lpts)):
                if est_sup(lpts[i], lpts[j]):
                    rayes.append(i)
                elif est_sup(lpts[j], lpts[i]):
                    rayes.append(j)
        M.append([lpts[j] for j in range(len(lpts)) if j not in rayes])
    maximal_points[k] = M

    L = minimal_points[k]
    M = []
    for lpts in L:
        rayes = []
        for i in range(len(lpts)):
            for j in range(i+1, len(lpts)):
                if est_sup(lpts[i], lpts[j]):
                    rayes.append(j)
                elif est_sup(lpts[j], lpts[i]):
                    rayes.append(i)
        M.append([lpts[j] for j in range(len(lpts)) if j not in rayes])
    minimal_points[k] = M


# Computes the list of all functions in {\mathcal F}_k for 1\leq k\leq 7.

for k in increasing_boolean_functions_reduced:
    if os.path.exists(f"BooleanFunctions/data/mini-k={k}") and os.path.exists(f"BooleanFunctions/data/maxi-k={k}"):
        read_mm(k)
    else:
        calcule_mm(k)
        elague_mm(k)
        write_mm(k)




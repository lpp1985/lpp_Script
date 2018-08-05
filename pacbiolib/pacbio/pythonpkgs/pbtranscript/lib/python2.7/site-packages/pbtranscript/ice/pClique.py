####################################################################
# Copyright (c) 2011-2014, Pacific Biosciences of California, Inc.
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted (subject to the limitations in the
# disclaimer below) provided that the following conditions are met:
#
#  * Redistributions of source code must retain the above copyright
#  notice, this list of conditions and the following disclaimer.
#
#  * Redistributions in binary form must reproduce the above
#  copyright notice, this list of conditions and the following
#  disclaimer in the documentation and/or other materials provided
#  with the distribution.
#
#  * Neither the name of Pacific Biosciences nor the names of its
#  contributors may be used to endorse or promote products derived
#  from this software without specific prior written permission.
#
# NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
# GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
# BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
# OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
# USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
# OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
# SUCH DAMAGE.
######################################################################
"""
Provide functions for clique finding.
"""

#import os, re, sys, cProfile, itertools
#from networkx import Graph
import random
#from bisect import bisect
from scipy import sparse
import logging

random.seed(0)

def convert_graph_connectivity_to_sparse(G, nodes):
    """
        Given a networkx graph, return sparse adjacency matrix S and H
        S and H are different in that S's entires contain edge weights
        (if there are multiple edges, behavior is overwrite),
        and H just has a 1 for every non-zero entry.

        NOTE: for now just use H, so returns None,H
    """
    n = G.number_of_nodes()
#    S = sparse.lil_matrix((n,n))
    H = sparse.lil_matrix((n, n))
    nodes_to_index = dict(zip(nodes, range(n)))
    for e in G.edges_iter(data=True):
        i = nodes_to_index[e[0]]
        j = nodes_to_index[e[1]]
        H[i, j] = 1
        H[j, i] = 1
    # we do a lot of column-slicing, so convert to CSC for efficiency
    H = H.tocsr()

    return None, H
#    return S,H


def construct(S, H, alpha, starting_node):
    """
    NOTE: S is currently None, so don't use it yet
    """
    assert S is None
    Q = [starting_node] # the list of indices in the clique

    # candidates = direct neighbors of <starting_node>
    C = H[starting_node, :].nonzero()[1]
    len_C = len(C)
    while len_C > 0:
#        degs_of_C = map(lambda x: S[x, C].nnz, C)
        # this is much faster than the above line!
        degs_of_C = H[:, C].sum(axis=1).getA1()
        degs_of_C = degs_of_C[C]
        min_deg_C = min(degs_of_C)
        max_deg_C = max(degs_of_C)
        RCL_threshold = min_deg_C + alpha*(max_deg_C - min_deg_C)
        # we've randomly picked u from RCL
        #RCL = filter(lambda i: degs_of_C[i] >= RCL_threshold, xrange(len_C))
        RCL = [i for i in xrange(len_C) if degs_of_C[i] >= RCL_threshold]
        # means there are now fitting RCLs!
        if len(RCL) == 0:
            logging.debug("NO FITTING RCLS")
            break
        u = C[random.choice(RCL)]
        logging.debug("picking {u}".format(u=u))

        Q.append(u)
        C = C[H[u, C].nonzero()[1]] # update list of candidates
        len_C = len(C)
    return Q


def local(H, Q, gamma):
    """
    """
    n = H.shape[0]
    len_Q = len(Q)
    h = H[:, Q]
    h_summed2 = h.sum(axis=1).getA1() # using this speeds up filter a LOT
    gamma_threshold = gamma*len_Q
    logging.debug("gamma threshold is {t}".format(t=gamma_threshold))
    #cand = filter(lambda i: i not in Q and
    #              h_summed2[i] >= gamma_threshold, xrange(n))
    cand = [i for i in xrange(n)
            if i not in Q and h_summed2[i] >= gamma_threshold]
    logging.debug("there are {0} candidates...".format(len(cand)))
    len_cand = len(cand)
    if len_cand < 2:
        return False

    x = h[cand]
    y = x * x.transpose()
    y.setdiag([0]*len_cand)
    y = y.toarray()
    Q_index_set = set(range(len_Q))
    choices = range(len_cand)
    random.shuffle(choices)

    newQ = Q + [None, None]
    for v in choices:
        u = y[v, :].argmax()
        if y[v, u] >= gamma_threshold:
            # this is a good (2,1)-exchange pair!
            # note that ww is the "index" in Q, Q[ww] is the
            # real thing we're removing
            # similarly, we're adding in cand[u] and cand[v]
            w = Q_index_set.difference(x[v, :].nonzero()[1])
            for ww in w:
                newQ[-2:] = [cand[u], cand[v]]
                newQ.pop( ww )
                if min(H[newQ, :][:, newQ].sum(axis=1)) >= gamma*(len_Q+1):
                    # Q = Q U {u,v}\{ww}
                    Q.pop(ww)
                    Q += [cand[u], cand[v]]
                    logging.debug("new list has size {0}".format(len(newQ)))
                    return True
                else:
                    newQ.insert(ww, Q[ww])
        else: logging.debug("y[v,u] not high enough")
    return False


def local_extra(H, Q, gamma):
    """Extract local nodes."""
    n = H.shape[0]
    len_Q = len(Q)
    h = H[:, Q]
    h_summed2 = h.sum(axis=1).getA1()
    gamma_threshold = gamma*(len_Q+1)
    #cand = filter(lambda i: h_summed2[i] >= gamma_threshold and
    #              i not in Q, xrange(n))
    cand = [i for i in xrange(n)
            if i not in Q and h_summed2[i] >= gamma_threshold]
    random.shuffle(cand)
    while len(cand) > 0:
        x = cand.pop()
        newQ = Q + [x]
        if min(H[newQ, :][:, newQ].sum(axis=1)) >= gamma_threshold:
            logging.debug("local extra was able to add in another node {0}!"
                .format(x))
            Q.append(x)
            len_Q += 1
            gamma_threshold = gamma*(len_Q+1)
    return False


def grasp(S, H, gamma, maxitr, given_starting_node=None):
    """Grasp cliques."""
    assert S is None

    N = H.shape[0]

    bestQ = []
    # pick a starting node unless given
    if given_starting_node is None:
        H_deg = H.sum(axis=1).getA1().tolist()
        #x = filter(lambda i: H_deg[i]>=1, xrange(N)) # used to be H_deg[i]>=3
        x = [i for i in xrange(N) if H_deg[i] >= 1]  # used to be H_deg[i]>=3
        if len(x) == 0:
            return []
        random.shuffle(x)#x.sort(key=lambda i: H_deg[i])#random.shuffle(x)
        starting_node = x.pop()
    else:
        starting_node = given_starting_node

    for _k in xrange(maxitr):
        # randomly pick alpha uniformly from [0.1,0.9]

        alpha = random.uniform(0.1, 0.9)
        logging.debug("picked starting node {0} with alpha {1}".
            format(starting_node, alpha))
        Q = construct(S, H, alpha, starting_node)
        if len(Q) <= 1:
            # no valid local exchange can be done...just give up this round
            if given_starting_node is None and len(x) > 0:
                starting_node = x.pop()
                #print >> sys.stderr, "CHANGING STARTING NODE to", starting_node
                continue
            else:
                return []
        logging.debug("before local exchange, size is {0}".format(len(Q)))
        while local(H, Q, gamma):
            pass
        local_extra(H, Q, gamma)
        logging.debug("max clique with {0} as starting node has size {1}".
            format(starting_node, len(Q)))
        if len(Q) > len(bestQ):
            bestQ = list(Q)

    return bestQ



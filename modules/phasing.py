#! python3

import random
from itertools import product

def deduplicate_query(alignments, overlapCutoff=0.5):
    '''
    Ensures that each region of a query contributes to at most one reference region, by
    keeping alignments from longest to shortest and discarding any whose query interval
    overlaps the already kept alignments by more than the cutoff proportion.

    Parameters:
        alignments -- a dictionary as produced by minimap2.parse_paf()
        overlapCutoff -- OPTIONAL; a float of the maximum proportion of an alignment's query
                         interval that may overlap previously kept alignments
    Returns:
        dedupAlignments -- a dictionary with the same structure as alignments
    '''
    dedupAlignments = {}
    for qid, qidAlignments in alignments.items():
        kept = []
        for aln in sorted(qidAlignments, key=lambda x: (x[1] - x[0], x), reverse=True):
            qstart, qend = aln[0], aln[1]
            overlap = sum([ max(0, min(qend, kend) - max(qstart, kstart)) for kstart, kend, *_ in kept ])
            if overlap / (qend - qstart) <= overlapCutoff:
                kept.append(aln)
        dedupAlignments[qid] = sorted(kept)
    return dedupAlignments

def merge_intervals(intervals):
    '''
    Parameters:
        intervals -- an iterable of (start, end) tuples in 0-based, end-exclusive coordinates
    Returns:
        merged -- a sorted list of non-overlapping (start, end) tuples
    '''
    merged = []
    for start, end in sorted(intervals):
        if merged != [] and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged

def intervals_length(intervals):
    return sum([ end - start for start, end in intervals ])

def build_footprints(alignments):
    '''
    Summarises the reference regions each query aligns to. A query may align to
    any number of reference contigs, such as when a chromosomal segment has been
    translocated relative to the reference.

    Parameters:
        alignments -- a dictionary as produced by minimap2.parse_paf() or deduplicate_query()
    Returns:
        footprints -- a dictionary with structure like:
                      {
                          qid1: { tid1: [(tstart, tend), ...], tid2: [...] },
                          qid2: { ... },
                          ...
                      }
                      where intervals are merged and sorted
    '''
    footprints = {}
    for qid, qidAlignments in alignments.items():
        byRef = {}
        for qstart, qend, tid, tstart, tend in qidAlignments:
            byRef.setdefault(tid, []).append((tstart, tend))
        if byRef != {}:
            footprints[qid] = { tid: merge_intervals(intervals) for tid, intervals in byRef.items() }
    return footprints

class HaplotypePartitioner:
    '''
    Partitions query contigs into a number of haplotype groups (equal to the ploidy) by
    optimising each group's coverage of the reference, while penalising any redundancy
    (i.e., reference regions covered more than once within a group). Contigs which would
    only add redundancy to every group are left unassigned.

    For each group and reference contig, the score is:
        covered - penalty * redundant
    where 'covered' is the number of reference bases covered by the group's contigs and
    'redundant' is the number of bases covered beyond the first time. Hence, a contig
    is only worth assigning when the new coverage it adds exceeds penalty * its redundancy.

    Contigs only interact when their footprints overlap, so the optimisation is performed
    separately for each connected component of overlapping contigs. Small components are
    solved exactly by brute force, whereas larger ones are solved by greedy assignment
    followed by local search with move and swap operations, with random restarts.
    '''
    EXACT_LIMIT = 50000 # maximum number of assignment combinations to brute force
    RESTARTS = 20 # number of random restarts for the heuristic search
    SEED = 1
    EPSILON = 1e-9

    def __init__(self, ploidy, penalty=1.0):
        self.ploidy = ploidy
        self.penalty = penalty

    @property
    def ploidy(self):
        return self._ploidy

    @ploidy.setter
    def ploidy(self, value):
        if not isinstance(value, int) or value < 1:
            raise ValueError("ploidy must be an integer >= 1")
        self._ploidy = value

    @property
    def penalty(self):
        return self._penalty

    @penalty.setter
    def penalty(self, value):
        if not isinstance(value, (int, float)) or value < 0:
            raise ValueError("penalty must be a number >= 0")
        self._penalty = value

    @staticmethod
    def find_neighbours(footprints):
        '''
        Returns:
            neighbours -- a dictionary where keys are query IDs and values are sets
                          of the other query IDs whose footprints overlap it
        '''
        neighbours = { qid: set() for qid in footprints }

        # Sweep through intervals on each reference contig to find overlaps
        byRef = {}
        for qid, qidFootprint in footprints.items():
            for tid, intervals in qidFootprint.items():
                for start, end in intervals:
                    byRef.setdefault(tid, []).append((start, end, qid))

        for tid, intervals in byRef.items():
            active = []
            for start, end, qid in sorted(intervals):
                active = [ (activeEnd, activeID) for activeEnd, activeID in active if activeEnd > start ]
                for activeEnd, activeID in active:
                    if activeID != qid:
                        neighbours[qid].add(activeID)
                        neighbours[activeID].add(qid)
                active.append((end, qid))
        return neighbours

    @staticmethod
    def find_components(neighbours):
        '''
        Returns:
            components -- a list of lists, each containing the query IDs of a connected
                          component of overlapping contigs
        '''
        components = []
        visited = set()
        for qid in sorted(neighbours):
            if qid in visited:
                continue
            component = []
            toVisit = [qid]
            visited.add(qid)
            while toVisit != []:
                current = toVisit.pop()
                component.append(current)
                for neighbour in neighbours[current]:
                    if not neighbour in visited:
                        visited.add(neighbour)
                        toVisit.append(neighbour)
            components.append(sorted(component))
        return components

    def score(self, states, footprints, refs=None):
        '''
        Parameters:
            states -- a dictionary where keys are query IDs and values are the 0-based
                      index of their group, or None if unassigned
            footprints -- a dictionary as produced by build_footprints()
            refs -- OPTIONAL; a set of reference contig IDs to restrict scoring to
        Returns:
            score -- a float of the summed score of all groups; higher is better
        '''
        groups = {}
        for qid, state in states.items():
            if state is None:
                continue
            for tid, intervals in footprints[qid].items():
                if refs is not None and not tid in refs:
                    continue
                groups.setdefault((state, tid), []).extend(intervals)

        score = 0
        for intervals in groups.values():
            covered = intervals_length(merge_intervals(intervals))
            redundant = intervals_length(intervals) - covered
            score += covered - self.penalty * redundant
        return score

    def _delta(self, changes, states, footprints, neighbours):
        '''
        Calculates the change in score from applying changes to states. Only overlapping
        contigs can affect a contig's contribution, so scoring is restricted to them.

        Parameters:
            changes -- a dictionary where keys are query IDs and values are their new state
        '''
        involved = set(changes)
        for qid in changes:
            involved.update(neighbours[qid])
        refs = set([ tid for qid in changes for tid in footprints[qid] ])

        before = { qid: states[qid] for qid in involved }
        after = { **before, **changes }
        return self.score(after, footprints, refs) - self.score(before, footprints, refs)

    def _exact(self, component, footprints):
        bestStates, bestScore = None, None
        for combination in product(list(range(self.ploidy)) + [None], repeat=len(component)):
            states = dict(zip(component, combination))
            score = self.score(states, footprints)
            if bestScore is None or score > bestScore + HaplotypePartitioner.EPSILON:
                bestStates, bestScore = states, score
        return bestStates

    def _heuristic(self, component, footprints, neighbours):
        '''
        Performs greedy assignment and local search, first ordering contigs from largest to
        smallest and then restarting from random orders, keeping the best scoring result.
        '''
        sizes = { qid: sum([ intervals_length(x) for x in footprints[qid].values() ]) for qid in component }
        order = sorted(component, key=lambda qid: (-sizes[qid], qid))
        rng = random.Random(HaplotypePartitioner.SEED) # seeded for reproducibility

        bestStates, bestScore = None, None
        for restart in range(HaplotypePartitioner.RESTARTS + 1):
            if restart > 0:
                order = rng.sample(order, len(order))
            states = self._local_search(order, footprints, neighbours)
            score = self.score(states, footprints)
            if bestScore is None or score > bestScore + HaplotypePartitioner.EPSILON:
                bestStates, bestScore = states, score
        return bestStates

    def _local_search(self, order, footprints, neighbours):
        '''
        Parameters:
            order -- a list of query IDs in the order they should be greedily assigned
        '''
        # Greedily assign contigs, in order, to the group they benefit most
        states = { qid: None for qid in order }
        for qid in order:
            bestDelta, bestState = 0, None # assignment must beat leaving it unassigned
            for state in range(self.ploidy):
                delta = self._delta({qid: state}, states, footprints, neighbours)
                if delta > bestDelta + HaplotypePartitioner.EPSILON:
                    bestDelta, bestState = delta, state
            states[qid] = bestState

        # Refine assignments with local search until no improvement is found
        edges = sorted(set([ tuple(sorted((qid, neighbour))) for qid in order for neighbour in neighbours[qid] ]))
        improved = True
        while improved:
            improved = False

            # Try moving each contig to another group (or to be unassigned)
            for qid in order:
                for state in list(range(self.ploidy)) + [None]:
                    if state == states[qid]:
                        continue
                    if self._delta({qid: state}, states, footprints, neighbours) > HaplotypePartitioner.EPSILON:
                        states[qid] = state
                        improved = True

            # Try swapping the groups of overlapping contigs
            for qid1, qid2 in edges:
                if states[qid1] == states[qid2]:
                    continue
                changes = {qid1: states[qid2], qid2: states[qid1]}
                if self._delta(changes, states, footprints, neighbours) > HaplotypePartitioner.EPSILON:
                    states.update(changes)
                    improved = True
        return states

    def _relabel(self, states, footprints):
        '''
        Group labels are interchangeable, so this orders them such that the first group
        has the highest score (i.e., most non-redundant coverage) within the component.
        '''
        coverage = { state: self.score({ qid: 0 for qid, s in states.items() if s == state }, footprints)
                     for state in range(self.ploidy) }
        order = sorted(range(self.ploidy), key=lambda state: -coverage[state]) # stable for ties
        relabel = { oldState: newState for newState, oldState in enumerate(order) }
        return { qid: (None if state is None else relabel[state]) for qid, state in states.items() }

    def partition(self, footprints):
        '''
        Parameters:
            footprints -- a dictionary as produced by build_footprints()
        Returns:
            states -- a dictionary where keys are query IDs and values are the 0-based
                      index of their group, or None if unassigned
        '''
        neighbours = HaplotypePartitioner.find_neighbours(footprints)
        states = {}
        for component in HaplotypePartitioner.find_components(neighbours):
            if (self.ploidy + 1) ** len(component) <= HaplotypePartitioner.EXACT_LIMIT:
                componentStates = self._exact(component, footprints)
            else:
                componentStates = self._heuristic(component, footprints, neighbours)
            states.update(self._relabel(componentStates, footprints))
        return states

def unique_coverage(states, footprints):
    '''
    Parameters:
        states -- a dictionary as produced by HaplotypePartitioner.partition()
        footprints -- a dictionary as produced by build_footprints()
    Returns:
        uniqueBp -- a dictionary where keys are assigned query IDs and values are integers
                    of the reference bases covered by that contig and no other contig
                    in its group
    '''
    # Collect interval start and end events for each group and reference contig
    events = {}
    for qid, state in states.items():
        if state is None:
            continue
        for tid, intervals in footprints[qid].items():
            for start, end in intervals:
                events.setdefault((state, tid), []).extend([ (start, 1, qid), (end, -1, qid) ])

    # Sweep through events, attributing bases covered by one contig to that contig
    uniqueBp = { qid: 0 for qid, state in states.items() if state is not None }
    for groupEvents in events.values():
        active = set()
        prevPos = None
        for pos, change, qid in sorted(groupEvents, key=lambda x: (x[0], x[1])): # ends sort before starts
            if len(active) == 1 and prevPos is not None:
                uniqueBp[next(iter(active))] += pos - prevPos
            if change == 1:
                active.add(qid)
            else:
                active.discard(qid)
            prevPos = pos
    return uniqueBp

def group_coverage(states, footprints, ploidy):
    '''
    Returns:
        coverage -- a dictionary with structure like:
                    {
                        tid1: [group1Bp, group2Bp, ...],
                        tid2: [...],
                        ...
                    }
    '''
    byGroupRef = {}
    for qid, state in states.items():
        if state is None:
            continue
        for tid, intervals in footprints[qid].items():
            byGroupRef.setdefault(tid, [ [] for _ in range(ploidy) ])[state].extend(intervals)
    return { tid: [ intervals_length(merge_intervals(x)) for x in groupIntervals ]
             for tid, groupIntervals in byGroupRef.items() }

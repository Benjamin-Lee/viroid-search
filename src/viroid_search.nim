# This is just an example to get you started. A typical hybrid package
# uses this file as the main entry point of the application.

import bioseq
import tables
import sequtils
import sets
import strutils
import terminal
import strformat

const K = 17
const FASTA_PATH = "GSM458930-2718.fasta"
const USE_CACHE = true
if not USE_CACHE:
  styledEcho fgYellow, "[WARNING] ", fgDefault, "Cache is disabled!"

var overlapCache = initTable[(string, string), bool]()
var overlapCacheHit = 0
var overlapCacheMiss = 0

proc overlapsWithStartOf*(a: string, b: string, k: int): bool = 
  ## Checks whether the end of `a` overlaps with the start (left) of `b` by at least `k` bases.
  runnableExamples:
    # the last three of `a` and first three of `b` overlap
    assert "ATGCAGA".overlapsWithStartOf("AGATTAGATA", 3) == true
    # here, the last three are not equivalent to the first three but there is still an overlap
    assert "GGCCAAGCCC".overlapsWithStartOf("GCCCAGGTATGC", 3) == true
    # note that CGA (first three bases of `b`) also occurs at the start of `a`
    assert "CGATATTTCGATA".overlapsWithStartOf("GATATCAGAA", 3) == true

  if USE_CACHE and (a, b) in overlapCache:
    overlapCacheHit += 1
    return overlapCache[(a, b)]
  overlapCacheMiss += 1
  let kmer = b[0..<k]

  # find the k-mer's indices in `a`
  # note that, when called here, we already know that the k-mer is somewhere in the sequence
  var overlapIndicesInA = newSeq[int]() 
  for i, kmerInA in a.kmersWithIndices(k, degeneratesAllowed=true):
    if kmerInA == kmer:
      overlapIndicesInA.add(i)

  # circuit breaker since a sequence in which the end of `a` isn't present doesn't overlap
  if overlapIndicesInA.len == 0:
    return false

  # check whether any of the shared k-mers are a true overlap
  for overlapIndex in overlapIndicesInA:
    # This is a bit complex, so we'll take it step by step. 
    # We need to check if any potential overlap of `a` is actually an overlap.
    # For example, let a = "ILIKECOOKINGJUSTTOCOOK". b ="COOKINGISFUN", and k = 3.
    # The start of `b` is "COO", which appears twice in `a`.
    # However, the first time "COO" appears in `a` *is not* an overlap but the second is.
    # We're whether `a` continues with `b[0..min(b.high, a.high - overlapIndex)`,
    # starting from the overlap index and continuing as far in `b` as possible.
    # To do this, we need to calculate how long within `b` we can go. Why?
    # As expected, `"helloworld".continuesWith("world", 5) == true`.
    # However, `"helloworld".continuesWith("worldly", 5) != true` due to the fact that
    # "worldly" is not a substring of "helloworld".
    # Therefore, we need `min(b.high, a.high - overlapIndex)` since there are two cases:
    # 
    # 1. either the rest of the `b` is potentially fully contained as a substring or
    # 2. `b` is too long to be a substring even if it overlaps perfectly
    #
    # Going back to the earlier exanple, consider that, from the first k-mer match on,
    # `"COOKINGJUSTTOCOOK".len` < `"COOKINGISFUN".len`. Thus, we need to use the entire 
    # sequence of `b` when checking if `a` continues with `b`.
    # 
    # With the second k-mer match, "COOK", not all of `b` could be contained within `a`
    # from the k-mer overlap index onwards. Thus, we need to only use the first n bases
    # of `b`. How many? The number of letters of "COOK", which is the distance of the
    # overlap index from the end of the sequence, i.e. `a.high - overlapIndex`.
    if a.continuesWith(b[0..min(b.high, a.high - overlapIndex)], overlapIndex):
      if USE_CACHE:
        overlapCache[(a, b)] = true
      return true
  if USE_CACHE:
    overlapCache[(a, b)] = false
  return false


when isMainModule:

  styledEcho fgCyan, "[Info] ", fgDefault, &"Loading reads from {FASTA_PATH}..."

  # read and deduplicate the input sequences
  var sequences = initTable[string, string]()
  for id, sequence in fasta(FASTA_PATH):
    sequences[sequence] = id

  # invert sequence to ID mapping
  let idsToSequences = toSeq(sequences.pairs).mapIt((it[1], it[0])).toTable

  # This will contain all of the viroid-derived sRNAs
  # at first, all reads are putative ISRs
  var internalSmallRnas = toHashSet(toSeq(idsToSequences.keys))

  # A mapping of k-mers to the sequences they occur in
  var kmersToIds = initTable[string, HashSet[string]]()
  for sequence, id in sequences.pairs():
    for _, kmer in sequence.kmersWithIndices(K, degeneratesAllowed=true):
      if kmersToIds.hasKeyOrPut(kmer, toHashSet([id])):
        kmersToIds[kmer].incl(id)

  proc checkForOverlaps(sequence: string, sequencesToCheck: seq[string], start: bool, k: int): bool =
    ## Checks if any element of `sequences` overlap with the start of `sequence`
    for potentialOverlap in sequencesToCheck:
      if start and potentialOverlap.overlapsWithStartOf(sequence, k):
        return true
      elif sequence.overlapsWithStartOf(potentialOverlap, k):
        return true 
    return false
  
  styledEcho fgCyan, "[Info] ", fgDefault, "Data loading complete. Beginning filtering..."

  var tsrsRemoved = 1 
  var iteration = 0
  var isrCounter = 0

  while tsrsRemoved != 0:
    tsrsRemoved = 0
    isrCounter = 0
    iteration += 1

    var isrCard = internalSmallRnas.card
    for id in internalSmallRnas.items():

      # for friendlier output
      if isrCounter mod 1000 == 0:
        stdout.styledWrite "\r", fgCyan, &"[Info] ", fgDefault, &"Iteration: {iteration} ", &"Terminal Reads Removed: {tsrsRemoved}, Remaining Reads: {internalSmallRnas.card}, Progress: {(isrCounter / isrCard) * 100:.1f}% ({isrCounter}/{isrCard})"
        stdout.flushFile()
      isrCounter += 1 

      var sequence = idsToSequences[id]

      # check if start kmer overlaps with any sequence
      var sequencesWithStartKmer = kmersToIds[sequence[0..<K]]
      var idsToCheckForStartOverlaps= sequencesWithStartKmer.intersection(internalSmallRnas) - toHashSet([id])
      var seqsToCheckForStartOverlaps = toSeq(idsToCheckForStartOverlaps.items()).mapIt(idsToSequences[it])
      if not sequence.checkForOverlaps(seqsToCheckForStartOverlaps, start=true, k=K):
        tsrsRemoved += 1
        internalSmallRnas.excl(id)
        continue

      # check if end k-mer overlaps
      var sequencesWithEndKmer = kmersToIds[sequence[^K .. ^1]]
      var idsToCheckForEndOverlaps = sequencesWithEndKmer.intersection(internalSmallRnas) - toHashSet([id])
      var seqsToCheckForEndOverlaps = toSeq(idsToCheckForEndOverlaps.items()).mapIt(idsToSequences[it])      
      if not sequence.checkForOverlaps(seqsToCheckForEndOverlaps, start=false, k=K):
        tsrsRemoved += 1
        internalSmallRnas.excl(id)
        continue
    
    stdout.eraseLine
    styledEcho fgCyan, "[Info] ", fgDefault, &"Iteration: {iteration} ", &"Terminal Reads Removed: {tsrsRemoved}, Remaining Reads: {internalSmallRnas.card}"

  styledEcho fgGreen, "[DONE] ", fgDefault, &"Filtering complete. {internalSmallRnas.card} remaining reads. 🚀"
  for x in toSeq(internalSmallRnas):
    echo ">", x
    echo idsToSequences[x]
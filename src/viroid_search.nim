# This is just an example to get you started. A typical hybrid package
# uses this file as the main entry point of the application.
# import nimprof
import bioseq
import tables
import sequtils
import sets
import strutils
import terminal
import strformat
import os

const USE_CACHE = true
if not USE_CACHE:
  styledEcho fgYellow, "[WARNING] ", fgDefault, "Cache is disabled!"

var overlapCache = initTable[(string, string), bool]()
var overlapCacheHit = 0
var overlapCacheMiss = 0
var iteration = 0

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
    # if iteration == 1:
      # echo &"{a}, {b}"
    return overlapCache[(a, b)]
  overlapCacheMiss += 1
  let kmer = b[0..<k]

  # find the k-mer's indices in `a`
  # note that, when called here, we already know that the k-mer is somewhere in the sequence
  for index, kmerInA in a.kmersWithIndices(k):
    if kmerInA == kmer:
      # This next part is a bit complex, so we'll take it step by step. 
      # We need to check if any potential overlap of `a` is actually an overlap.
      # For example, let a = "ILIKECOOKINGJUSTTOCOOK". b ="COOKINGISFUN", and k = 3.
      # The start of `b` is "COO", which appears twice in `a`.
      # However, the first time "COO" appears in `a` *is not* an overlap but the second is.
      # We're going to determine whether `a` continues with `b[0..min(b.high, a.high - index)`,
      # starting from the k-mer's index and continuing as far in `b` as possible.
      # To do this, we need to calculate how much of `b` we can check for continuation in `a`. 
      # Why? As expected, `"helloworld".continuesWith("world", 5) == true`.
      # However, `"helloworld".continuesWith("worldly", 5) != true` due to the fact that
      # "worldly" is not a substring of "helloworld".
      # Therefore, we need `min(b.high, a.high - index)` since there are two cases:
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
      # overlap index from the end of the sequence, i.e. `a.high - index`.
      if a.continuesWith(b[0..min(b.high, a.high - index)], index):
        if USE_CACHE:
          overlapCache[(a, b)] = true
        return true
  if USE_CACHE:
    overlapCache[(a, b)] = false
  return false


proc main(path: string, k: int) =

  styledEcho fgCyan, "[Info] ", fgDefault, &"Loading reads from {path}..."

  # read and deduplicate the input sequences
  var sequences = initTable[string, string]()
  var tooShortReads = 0
  var totalReads = 0
  for id, sequence in fasta(path):
    if sequence.len < k:
      tooShortReads += 1
      continue
    sequences[sequence] = id
    totalReads += 1

  # This will contain all of the viroid-derived sRNAs
  # at first, all reads are putative ISRs
  var internalSmallRnas = toHashSet(toSeq(sequences.keys))

  # A mapping of k-mers to the sequences they occur in
  var kmersToSeqs = initTable[string, HashSet[string]]()
  for sequence, id in sequences.pairs():
    for _, kmer in sequence.kmersWithIndices(k):
      if kmersToSeqs.hasKeyOrPut(kmer, toHashSet([sequence])):
        kmersToSeqs[kmer].incl(sequence)
        
  if tooShortReads > 0:
    styledEcho fgYellow, "[Warn] ", fgDefault, &"{tooShortReads} reads rejected for being less than {k} nt long."
  styledEcho fgCyan, "[Info] ", fgDefault, &"Loaded {sequences.len} unique reads from {totalReads} total reads. Beginning filtering..."

  proc checkForOverlaps(sequence: string, sequencesToCheck: HashSet[string], start: bool): bool =
    ## Checks if any element of `sequences` overlap with the start of `sequence`
    for potentialOverlap in sequencesToCheck:
      if start and potentialOverlap.overlapsWithStartOf(sequence, k):
        return true
      elif sequence.overlapsWithStartOf(potentialOverlap, k):
        return true 
    return false

  proc removeFromKmersTable(sequence: string): void = 
     for i, kmer in sequence.kmersWithIndices(k):
          kmersToSeqs[kmer].excl(sequence)

  var tsrsRemoved = 1 
  var seqsToCheckForStartOverlaps: HashSet[string]
  var seqsToCheckForEndOverlaps: HashSet[string]

  while tsrsRemoved != 0:
    tsrsRemoved = 0
    iteration += 1
    overlapCacheHit = 0
    overlapCacheMiss = 0

    for sequence in internalSmallRnas:

      # check if start kmer overlaps with any sequence
      seqsToCheckForStartOverlaps = kmersToSeqs[sequence[0..<k]]
      seqsToCheckForStartOverlaps.excl(sequence)
      if not sequence.checkForOverlaps(seqsToCheckForStartOverlaps, start=true):
        tsrsRemoved += 1
        removeFromKmersTable(sequence)
        internalSmallRnas.excl(sequence)
        continue

      # check if end k-mer overlaps
      seqsToCheckForEndOverlaps = kmersToSeqs[sequence[^k .. ^1]]
      seqsToCheckForEndOverlaps.excl(sequence)
      if not sequence.checkForOverlaps(seqsToCheckForEndOverlaps, start=false):
        tsrsRemoved += 1
        removeFromKmersTable(sequence)
        internalSmallRnas.excl(sequence)
        continue
    
    stdout.eraseLine
    styledEcho fgCyan, "[Info] ", fgDefault, &"Iteration: {iteration}, ", &"Terminal Reads Removed: {tsrsRemoved}, Remaining Reads: {internalSmallRnas.card}, Cache hit rate: {(overlapCacheHit / (overlapCacheMiss + overlapCacheHit))*100.0:.1f}%"
  styledEcho fgGreen, "[DONE] ", fgDefault, &"Filtering complete. {internalSmallRnas.card} remaining reads. 🚀"
  for sequence in toSeq(internalSmallRnas):
    echo ">", sequences[sequence]
    echo sequence

when isMainModule:
  var params = commandLineParams()
  var fp = params[0]
  var k = params[1].parseInt
  main(fp, k)
# This is just an example to get you started. A typical hybrid package
# uses this file as the main entry point of the application.
import bioseq
import tables
import sets
import strutils
import terminal
import strformat
import os
import std/enumerate
import math
import sequtils
import std/monotimes
import times

let startTime = getMonoTime()

template styledWrite(color: ForegroundColor, level: string, message: string, indent=0) =
  var styledDuration = insertSep($((getMonoTime() - startTime).inSeconds), ',')
  if stderr.isatty:
    stderr.styledWriteLine color, spaces(indent), "[", level.capitalizeAscii, "] ", fgDefault, message, styleDim, " [", styledDuration , "s]"
  else:
   stderr.writeLine spaces(indent), "[", level.capitalizeAscii, "] ", message, " [", styledDuration , "s]"
template info(message: string) = styledWrite(fgCyan, "info", message)
template warn(message: string) = styledWrite(fgYellow, "warn", message)
template error(message: string) = styledWrite(fgRed, "fail", message)
template success(message: string) = styledWrite(fgGreen, "done", message)
template prog(message: string) = styledWrite(fgDefault, "prog", message, indent=7)

# Needed overrides to only compare sequences in hashsets
import hashes
proc hash(x: Record[Dna]): Hash = hash(x.sequence) 
proc `==`(x, y: Record): bool = x.sequence == y.sequence

proc pfor*(internalSmallRnas: var HashSet[Record[Dna]], k: int, verbose=false, showProgress=false, logInterval=10) =
  ## Progressively filter overlapping reads such that all remaining reads overlap on both ends with another remaiing read.
  if verbose:
    info "Building terminal k-mer to read lookup table"
    info "Step 1: Create hash set of terminal k-mers"
  # hold the k-mers on either end of each read, of which there can be at most 2n
  var terminalKmers = initHashSet[Dna](internalSmallRnas.len * 2)
  # build the terminal k-mer HashSet (TODO: Bloom filter?)
  for i, sequence in enumerate(internalSmallRnas):
    if verbose and showProgress and i mod floor(internalSmallRnas.len / logInterval).toInt == 0 and i > 0:
      prog &"{insertSep($i, ',')} reads processed ({(i / (internalSmallRnas.len))*100.0:.1f}%)"
    terminalKmers.incl(sequence[0..<k].canonical)
    terminalKmers.incl(sequence[^k..^1].canonical)

  if verbose:
    info "Step 2: Map terminal k-mers to their sequences"
  # a mapping of start/end k-mers to the sequences they occur in
  var kmersToSeqs = initTable[Dna, HashSet[Record[Dna]]](terminalKmers.len)
  # build the mapping
  for i, sequence in enumerate(internalSmallRnas):
    if verbose and showProgress and i mod floor(internalSmallRnas.len / logInterval).toInt == 0 and i > 0:
      prog &"{i} reads with terminal k-mers mapped  ({(i / (internalSmallRnas.len))*100.0:.1f}%)"
    for kmer in sequence.canonicalKmers(k):
      # we only care about k-mers that appear at the start or end of other reads
      if kmer notin terminalKmers:
        continue
      if kmersToSeqs.hasKeyOrPut(kmer, toHashSet([sequence])):
        kmersToSeqs[kmer].incl(sequence)
  if verbose:
    info "Done mapping"
  # clear the memory of the terminal k-mer HashSet (?)
  init(terminalKmers) 

  if verbose:
    info "Allocating cache"
  var lastStartOverlap = initTable[Record[Dna], Record[Dna]](internalSmallRnas.card)
  var lastEndOverlap = initTable[Record[Dna], Record[Dna]](internalSmallRnas.card)
  var lastOverlapDependents = initTable[Record[Dna], seq[Record[Dna]]](internalSmallRnas.card)
  var overlapCacheHit = 0
  var overlapCacheMiss = 0
  var iteration = 0
  var overlapPresent: bool
  var toRemove: seq[Record[Dna]]
  var toCheck = toSeq(internalSmallRnas) # expensive!?
  let invalidDefaultRecord = toRecord(Dna"*INVALID SEQUENCE*", "An invalid record that will never be in the main set")
  if verbose:
    info "Cache allocated"
    info "Beginning first round of filtering. This may take a while!"
  while toCheck.len != 0:
    overlapCacheHit = 0
    overlapCacheMiss = 0
    toRemove = newSeqOfCap[Record[Dna]](if toRemove.len == 0: internalSmallRnas.card else: toRemove.len)

    if verbose:
      info &"Iteration {iteration} (reads to check: {insertSep($(toCheck.len), ',')})"

    template checkAndRemove(sequence: Record[Dna], start: bool) =
      overlapPresent = false
      # first, check the cache
      when start:
        if lastStartOverlap.getOrDefault(sequence, invalidDefaultRecord) in internalSmallRnas:
          inc(overlapCacheHit)
          overlapPresent = true
      else:
        if lastEndOverlap.getOrDefault(sequence, invalidDefaultRecord) in internalSmallRnas:
          inc(overlapCacheHit)
          overlapPresent = true

      # we can't "return" since this is a template, so we'll use a conditional
      if not overlapPresent:
        for potentialOverlap in kmersToSeqs[if start: sequence[0..<k].canonical else: sequence[^k .. ^1].canonical]:
          if potentialOverlap notin internalSmallRnas or potentialOverlap == sequence:
            continue
          # check if read containing the start k-mer of `sequence` overlaps at the start of `sequence`
          when start:
            if potentialOverlap.overlapsWithStartOf(sequence, k) or potentialOverlap.reverseComplement.overlapsWithStartOf(sequence, k):
              inc(overlapCacheMiss)
              lastStartOverlap[sequence] = potentialOverlap
              lastOverlapDependents.mgetOrPut(potentialOverlap, newSeq[Record[Dna]]()).add(sequence)
              overlapPresent = true
              break
          # ibid, but for the end.
          else:
            if sequence.overlapsWithStartOf(potentialOverlap, k) or sequence.overlapsWithStartOf(potentialOverlap.reverseComplement, k):
              inc(overlapCacheMiss)
              lastEndOverlap[sequence] = potentialOverlap
              lastOverlapDependents.mgetOrPut(potentialOverlap, newSeq[Record[Dna]]()).add(sequence)
              overlapPresent = true 
              break
      if not overlapPresent:
        toRemove.add(sequence)
        continue

    for i, sequence in enumerate(toCheck):
      if verbose and showProgress and toCheck.len > 10 and i mod floor(toCheck.len / logInterval).toInt == 0 and i > 0:
        prog &"{($i).insertSep(',')} reads checked ({(i / (toCheck.len))*100.0:.1f}%)"

      # check if start kmer overlaps with any sequence
      checkAndRemove(sequence, start=true)
      # check if end k-mer overlaps
      checkAndRemove(sequence, start=false)

    toCheck.setLen(0)
    for sequence in toRemove:
      kmersToSeqs[sequence[0..<k].canonical].excl(sequence)
      kmersToSeqs[sequence[^k .. ^1].canonical].excl(sequence)
      internalSmallRnas.excl(sequence)
      for x in lastOverlapDependents.getOrDefault(sequence):
        if x in internalSmallRnas:
          toCheck.add(x)
    if verbose: 
      stderr.write(spaces(7))
      success &"Terminal Reads Removed: {insertSep($(toRemove.len), ',')} (cache hit rate: {(overlapCacheHit / (overlapCacheMiss + overlapCacheHit))*100.0:.1f})%"
    inc(iteration)
  if verbose:
    success &"Filtering complete. {insertSep($(internalSmallRnas.len), ',')} remaining reads. 🚀"

proc main*(inputPath: string, k: int, outputPath: string, verbose=false, showProgress=false, preserveDuplicates=false, logInterval=10) =
  when not defined(danger):
    warn "Not compiled as dangerous release. This may be slow!"

  if outputPath != "stdout" and fileExists(outputPath):
    error "Output file already exists."
    quit(1)
 
  if not fileExists(inputPath):
    error "Input file does not exist."
    quit(1)

  # read and deduplicate the input sequences
  var tooShortReads = 0
  var totalReads = 0
  var internalSmallRnas = HashSet[Record[Dna]]()
  if verbose: 
    info &"Loading reads from {inputPath}"
  for record in readFastx[Dna](inputPath):
    if record.len <= k:
      tooShortReads += 1
      continue
    internalSmallRnas.incl(record)
    totalReads += 1
    if verbose and showProgress and totalReads mod 25_000 == 0:
      # stderr.write(spaces(4))
      prog &"Reads loaded: {($totalReads).insertSep(',')}"

  # report the results of reading
  if verbose: 
    if tooShortReads > 0: 
      warn &"{tooShortReads} reads rejected for being less than or equal to {k}nt long."
    info &"Loaded {internalSmallRnas.len} unique reads from {totalReads} total reads. Beginning filtering..."

  pfor(internalSmallRnas, k, verbose=verbose, showProgress=showProgress, logInterval=logInterval)

  var f: File
  if outputPath == "stdout":
    f = stdout
  else:
    f = open(outputPath, fmWrite)

  try: 
    let outputFastq = inputPath.splitFile.ext.toLowerAscii in [".fq", ".fastq"]
    for record in internalSmallRnas:
      if outputFastq:
        f.writeLine record.asFastq
      else:
        f.writeLine record.asFasta

    # optionally, go back through the input reads and output those that were removed as duplicates
    if preserveDuplicates:
      for record in readFastx[Dna](inputPath):
        if record in internalSmallRnas and record.description != internalSmallRnas[record].description:
          if outputFastq:
            f.writeLine record.asFastq
          else:
            f.writeLine record.asFasta
  finally:
    f.close

when isMainModule:

  import argparse

  var p = newParser("viroid_search"):
    arg("filepath")
    option("-k", help="The minimum k-mer overlap", default="17")
    option("-o", "--output", help="The output file. Output format (FASTA or FASTQ) will be the same as input.", default="stdout")
    option("-l", "--log-interval", help="If showing progress, how many intervals to break tasks into. For larger runs, finer granularity may be desired.", default="10")
    flag("-q", "--quiet", help="If present, don't output anything.")
    flag("-p", "--preserve-duplicates", help="Whether to output multiples reads with identical sequences")
    run:
      main(opts.filepath, opts.k.parseInt, opts.output, showProgress=true, verbose=not opts.quiet, preserveDuplicates=opts.preserveDuplicates, logInterval=opts.logInterval.parseInt)
  p.run(if commandLineParams().len == 0: @["--help"] else: commandLineParams())
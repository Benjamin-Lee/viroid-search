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

template styledWrite(color: ForegroundColor, level: string, message: string, indent=0) =
  if stderr.isatty:
    stderr.styledWriteLine color, spaces(indent), "[", level.capitalizeAscii, "] ", fgDefault, message
  else:
   stderr.writeLine spaces(indent), "[", level.capitalizeAscii, "] ", message
template info(message: string) = styledWrite(fgCyan, "info", message)
template warn(message: string) = styledWrite(fgYellow, "warn", message)
template error(message: string) = styledWrite(fgRed, "fail", message)
template success(message: string) = styledWrite(fgGreen, "done", message)
template prog(message: string) = styledWrite(fgDefault, "prog", message, indent=4)

# Needed overrides to only compare sequences in hashsets
import hashes
proc hash(x: Record[Dna]): Hash = hash(x.sequence) 
proc `==`(x, y: Record): bool = x.sequence == y.sequence

proc pfor*(internalSmallRnas: var HashSet[Record[Dna]], k: int, verbose=false, showProgress=false) =

  # A mapping of k-mers to the sequences they occur in
  var kmersToSeqs = initTable[Dna, HashSet[Record[Dna]]]()
  for sequence in internalSmallRnas:
    for kmer in sequence.canonicalKmers(k):
      if kmersToSeqs.hasKeyOrPut(kmer, toHashSet([sequence])):
        kmersToSeqs[kmer].incl(sequence)

  var lastStartOverlap = initTable[Record[Dna], Record[Dna]](internalSmallRnas.card)
  var lastEndOverlap = initTable[Record[Dna], Record[Dna]](internalSmallRnas.card)
  var lastOverlapDependents = initTable[Record[Dna], HashSet[Record[Dna]]](internalSmallRnas.card)
  var overlapCacheHit = 0
  var overlapCacheMiss = 0
  var iteration = 0
  var sequencesToCheck: HashSet[Record[Dna]]
  var overlapPresent: bool
  var toRemove: seq[Record[Dna]]
  var toCheck = internalSmallRnas
    
  while toCheck.card != 0:
    overlapCacheHit = 0
    overlapCacheMiss = 0
    toRemove = newSeqOfCap[Record[Dna]](if toRemove.len == 0: internalSmallRnas.card else: toRemove.len)

    if verbose:
      info (&"Iteration {iteration}").alignLeft(20) & &"Reads to check: {toCheck.len}"

    template checkAndRemove(sequence: Record[Dna], start: bool) =
      sequencesToCheck = kmersToSeqs[if start: sequence[0..<k].canonical else: sequence[^k .. ^1].canonical]
      sequencesToCheck.excl(sequence)
      overlapPresent = false

      ## Checks if any element of `sequences` overlap with the start of `sequence`
      for potentialOverlap in sequencesToCheck:
        when start:
          if sequence in lastStartOverlap and lastStartOverlap[sequence] in internalSmallRnas:
            inc(overlapCacheHit)
            overlapPresent = true
            break
          if potentialOverlap.overlapsWithStartOf(sequence, k) or potentialOverlap.reverseComplement.overlapsWithStartOf(sequence, k):
            lastStartOverlap[sequence] = potentialOverlap
            if lastOverlapDependents.hasKeyOrPut(potentialOverlap, [sequence].toHashSet):
              lastOverlapDependents[potentialOverlap].incl(sequence)
            inc(overlapCacheMiss)
            overlapPresent = true
            break
        else:
          if sequence in lastEndOverlap and lastEndOverlap[sequence] in internalSmallRnas:
            inc(overlapCacheHit)
            overlapPresent = true
            break
          if sequence.overlapsWithStartOf(potentialOverlap, k) or sequence.overlapsWithStartOf(potentialOverlap.reverseComplement, k):
            inc(overlapCacheMiss)
            lastEndOverlap[sequence] = potentialOverlap
            if lastOverlapDependents.hasKeyOrPut(potentialOverlap, [sequence].toHashSet):
              lastOverlapDependents[potentialOverlap].incl(sequence)
            overlapPresent = true 
            break
      if not overlapPresent:
        toRemove.add(sequence)
        continue

    for i, sequence in enumerate(toCheck):
      if verbose and showProgress and toCheck.len > 10 and i mod floor(toCheck.len / 10).toInt == 0 and i > 0:
        prog &"{i} reads checked ({(i / (toCheck.len))*100.0:.1f}%)"

      # check if start kmer overlaps with any sequence
      checkAndRemove(sequence, start=true)
      # check if end k-mer overlaps
      checkAndRemove(sequence, start=false)

    toCheck = initHashSet[Record[Dna]]()
    for sequence in toRemove:
      for kmer in sequence.canonicalKmers(k):
        kmersToSeqs[kmer].excl(sequence)
      internalSmallRnas.excl(sequence)
      for x in lastOverlapDependents.getOrDefault(sequence, initHashSet[Record[Dna]]()):
        if x in internalSmallRnas:
          toCheck.incl(x)
    if verbose: 
      stderr.write(spaces(4))
      success &"Terminal Reads Removed: {toRemove.len}, Cache hit rate: {(overlapCacheHit / (overlapCacheMiss + overlapCacheHit))*100.0:.1f}%"
    inc(iteration)
  if verbose:
    success &"Filtering complete. {internalSmallRnas.card} remaining reads. 🚀"

proc main*(inputPath: string, k: int, outputPath: string, verbose=false, showProgress=false, preserveDuplicates=false) =
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
    info &"Loading reads from {inputPath}..."
  for record in readFastx[Dna](inputPath):
    if record.len <= k:
      tooShortReads += 1
      continue
    internalSmallRnas.incl(record)
    totalReads += 1

  # report the results of reading
  if verbose: 
    if tooShortReads > 0: 
      warn &"{tooShortReads} reads rejected for being less than or equal to {k}nt long."
    info &"Loaded {internalSmallRnas.len} unique reads from {totalReads} total reads. Beginning filtering..."

  pfor(internalSmallRnas, k, verbose=verbose, showProgress=showProgress)

  var f: File
  if outputPath == "stdout":
    f = stdout
  else:
    f = open(outputPath, fmWrite)

  try: 
    let outputFastq = outputPath.splitFile.ext.toLowerAscii in [".fq", ".fastq"]
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
    flag("-q", "--quiet", help="If present, don't output anything.")
    flag("-p", "--preserve-duplicates", help="Whether to output multiples reads with identical sequences")
    run:
      main(opts.filepath, opts.k.parseInt, opts.output, showProgress=true, verbose=not opts.quiet, preserveDuplicates=opts.preserveDuplicates)
  p.run(if commandLineParams().len == 0: @["--help"] else: commandLineParams())
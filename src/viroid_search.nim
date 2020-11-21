# This is just an example to get you started. A typical hybrid package
# uses this file as the main entry point of the application.
import bioseq
import tables
import sets
import strutils
import terminal
import strformat
import os

template styledWrite(color: ForegroundColor, level: string, message: string) =
  stderr.styledWriteLine color, "[", level.capitalizeAscii, "] ", fgDefault, message
template info(message: string) = styledWrite(fgCyan, "info", message)
template warn(message: string) = styledWrite(fgYellow, "warn", message)
template error(message: string) = styledWrite(fgRed, "fail", message)
template success(message: string) = styledWrite(fgGreen, "done", message)

# Needed overrides to only compare sequences in hashsets
import hashes
proc hash(x: Record[Dna]): Hash = hash(x.sequence) 
proc `==`(x, y: Record): bool = x.sequence == y.sequence

proc pfor*(internalSmallRnas: var HashSet[Record[Dna]], k: int, cache=true, verbose=false, showProgress=false) =
  if not cache and verbose:
    warn "Cache is disabled!"

  # A mapping of k-mers to the sequences they occur in
  var kmersToSeqs = initTable[Dna, HashSet[Record[Dna]]]()
  for sequence in internalSmallRnas:
    for kmer in sequence.canonicalKmers(k):
      if kmersToSeqs.hasKeyOrPut(kmer, toHashSet([sequence])):
        kmersToSeqs[kmer].incl(sequence)

  var lastStartOverlap = initTable[Record[Dna], Record[Dna]](internalSmallRnas.card)
  var lastEndOverlap = initTable[Record[Dna], Record[Dna]](internalSmallRnas.card)
  var overlapCacheHit = 0
  var overlapCacheMiss = 0
  var tsrsRemoved = 1 
  var seqsProcessed = 0
  var seqsToProcess = 0
  var iteration = 0
  var sequencesToCheck: HashSet[Record[Dna]]
  var overlapPresent: bool
    
  while tsrsRemoved != 0:
    tsrsRemoved = 0
    seqsProcessed = 0
    seqsToProcess = internalSmallRnas.card
    iteration += 1
    overlapCacheHit = 0
    overlapCacheMiss = 0

    template checkAndRemove(sequence: Record[Dna], start: bool) =
      sequencesToCheck = kmersToSeqs[if start: sequence[0..<k].canonical else: sequence[^k .. ^1].canonical]
      sequencesToCheck.excl(sequence)
      overlapPresent = false

      ## Checks if any element of `sequences` overlap with the start of `sequence`
      for potentialOverlap in sequencesToCheck:
        if start:
          if cache and sequence in lastStartOverlap and lastStartOverlap[sequence] in internalSmallRnas:
            inc(overlapCacheHit)
            overlapPresent = true
            break
          if potentialOverlap.overlapsWithStartOf(sequence, k) or potentialOverlap.reverseComplement.overlapsWithStartOf(sequence, k):
            if cache:
              lastStartOverlap[sequence] = potentialOverlap
              inc(overlapCacheMiss)
            overlapPresent = true
            break
        else:
          if cache and sequence in lastEndOverlap and lastEndOverlap[sequence] in internalSmallRnas:
            inc(overlapCacheHit)
            overlapPresent = true
            break
          if sequence.overlapsWithStartOf(potentialOverlap, k) or sequence.overlapsWithStartOf(potentialOverlap.reverseComplement, k):
            if cache:
              inc(overlapCacheMiss)
              lastEndOverlap[sequence] = potentialOverlap
            overlapPresent = true 
            break
      if not overlapPresent:
        tsrsRemoved += 1
        for kmer in sequence.canonicalKmers(k):
          kmersToSeqs[kmer].excl(sequence.toRecord())
        internalSmallRnas.excl(sequence)
        continue

    for sequence in internalSmallRnas:
      seqsProcessed += 1
      if verbose and showProgress and seqsProcessed mod 5000 == 0:
        stderr.write &"       {seqsProcessed} reads processed ({(seqsProcessed / (seqsToProcess))*100.0:.1f}%)\r"
        stderr.flushFile

      # check if start kmer overlaps with any sequence
      checkAndRemove(sequence, start=true)
      # check if end k-mer overlaps
      checkAndRemove(sequence, start=false)

    if verbose: 
      info &"Iteration: {iteration}, Terminal Reads Removed: {tsrsRemoved}, Remaining Reads: {internalSmallRnas.card}, Cache hit rate: {(overlapCacheHit / (overlapCacheMiss + overlapCacheHit))*100.0:.1f}%"
  if verbose:
    success &"Filtering complete. {internalSmallRnas.card} remaining reads. 🚀"

proc main*(inputPath: string, k: int, outputPath: string, cache=true, verbose=false, showProgress=false, preserveDuplicates=false) =
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

  pfor(internalSmallRnas, k, cache=cache, verbose=verbose, showProgress=showProgress)

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
    flag("-c", "--cache", help="If present, turn on caching.")
    flag("-q", "--quiet")
    flag("-p", "--preserve-duplicates", help="Whether to output multiples reads with identical sequences")
    run:
      main(opts.filepath, opts.k.parseInt, opts.output, cache=opts.cache, showProgress=true, verbose=not opts.quiet, preserveDuplicates=opts.preserveDuplicates)
  p.run(if commandLineParams().len == 0: @["--help"] else: commandLineParams())
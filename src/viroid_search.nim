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



template styledWrite(color: ForegroundColor, level: string, message: string) =
  stderr.styledWriteLine color, "[", level.capitalizeAscii, "] ", fgDefault, message
template info(message: string) = styledWrite(fgCyan, "info", message)
template warn(message: string) = styledWrite(fgYellow, "warn", message)
template error(message: string) = styledWrite(fgRed, "fail", message)
template success(message: string) = styledWrite(fgGreen, "done", message)

func overlapsWithStartOf*(a: Dna, b: Dna, k: int): bool = 
  ## Checks whether the end of `a` overlaps with the start (left) of `b` by at least `k` bases.
  runnableExamples:
    import bioseq
    # the last three of `a` and first three of `b` overlap
    assert "ATGCAGA".toDna.overlapsWithStartOf("AGATTAGATA".toDna, 3) == true
    # here, the last three are not equivalent to the first three but there is still an overlap
    assert "GGCCAAGCCC".toDna.overlapsWithStartOf("GCCCAGGTATGC".toDna, 3) == true
    # note that CGA (first three bases of `b`) also occurs at the start of `a`
    assert "CGATATTTCGATA".toDna.overlapsWithStartOf("GATATCAGAA".toDna, 3) == true
    # if `b` is a substring of `a` it doesn't count as overlapping
    assert "CCATG".toDna.overlapsWithStartOf("CATG".toDna, 3) == false

  let kmer = b[0..<k]
  # find the k-mer's indices in `a`
  # note that, when called here, we already know that the k-mer is somewhere in the sequence
  for index, kmerInA in enumerate(a.kmers(k)):
    # see notes.md section for explanation of this expression
    if kmerInA == kmer and b.high > a.high - index and a.endsWith(b[0..min(b.high, a.high - index)]):
        return true
  return false

iterator readFast(path: string): Record[Dna] =
  # Use the correct parser depending on the file extension
  var mode = if path.splitFile.ext.toLowerAscii in [".fastq", ".fq"]: "fastq" else: "fasta"
  if mode == "fasta":
    for record in readFasta[Dna](path): yield record
  elif mode == "fastq":
    for record in readFastq[Dna](path): yield record

# Needed overrides to only compare sequences in hashsets
import hashes
proc hash(x: Record[Dna]): Hash = hash(x.sequence) 
proc `==`(x, y: Record): bool = x.sequence == y.sequence

proc main*(path: string, k: int, cache=true, verbose=false, showProgress=false): HashSet[Record[Dna]] =
  if not cache and verbose:
    warn "Cache is disabled!"
  if verbose: 
    info &"Loading reads from {path}..."

  # read and deduplicate the input sequences
  var internalSmallRnas = HashSet[Record[Dna]]()
  var tooShortReads = 0
  var totalReads = 0

  for record in readFast(path):
    if record.len <= k:
      tooShortReads += 1
      continue
    internalSmallRnas.incl(record)
    totalReads += 1

  # This will contain all of the viroid-derived sRNAs
  # at first, all reads are putative ISRs
  # var internalSmallRnas = toHashSet(toSeq(sequences.keys))

  # A mapping of k-mers to the sequences they occur in
  var kmersToSeqs = initTable[Dna, HashSet[Record[Dna]]]()
  for sequence in internalSmallRnas:
    for kmer in sequence.canonicalKmers(k):
      if kmersToSeqs.hasKeyOrPut(kmer, toHashSet([sequence])):
        kmersToSeqs[kmer].incl(sequence)

  if tooShortReads > 0 and verbose: 
    warn &"{tooShortReads} reads rejected for being less than or equal to {k}nt long."
  if verbose: 
    info &"Loaded {internalSmallRnas.len} unique reads from {totalReads} total reads. Beginning filtering..."

  var lastStartOverlap = initTable[Record[Dna], Record[Dna]](internalSmallRnas.card)
  var lastEndOverlap = initTable[Record[Dna], Record[Dna]](internalSmallRnas.card)
  var overlapCacheHit = 0
  var overlapCacheMiss = 0

  proc checkForOverlaps(sequence: Record[Dna], sequencesToCheck: HashSet[Record[Dna]], start: bool): bool =
    ## Checks if any element of `sequences` overlap with the start of `sequence`
    for potentialOverlap in sequencesToCheck:
      if start:
        if sequence in lastStartOverlap and lastStartOverlap[sequence] in internalSmallRnas:
          inc(overlapCacheHit)
          return true
        if (potentialOverlap.overlapsWithStartOf(sequence, k) or potentialOverlap.reverseComplement.overlapsWithStartOf(sequence, k)):
          lastStartOverlap[sequence] = potentialOverlap
          inc(overlapCacheMiss)
          return true
      else:
        if sequence in lastEndOverlap and lastEndOverlap[sequence] in internalSmallRnas:
          inc(overlapCacheHit)
          return true
        if sequence.overlapsWithStartOf(potentialOverlap, k) or sequence.overlapsWithStartOf(potentialOverlap.reverseComplement, k):
          inc(overlapCacheMiss)
          lastEndOverlap[sequence] = potentialOverlap
          return true 
    return false

  template removeFromKmersTable(sequence: Record[Dna]) = 
    for kmer in sequence.canonicalKmers(k):
      kmersToSeqs[kmer].excl(sequence.toRecord())

  var tsrsRemoved = 1 
  var seqsProcessed = 0
  var seqsToProcess = 0
  var iteration = 0
  var seqsToCheckForStartOverlaps: HashSet[Record[Dna]]
  var seqsToCheckForEndOverlaps: HashSet[Record[Dna]]

  while tsrsRemoved != 0:
    tsrsRemoved = 0
    seqsProcessed = 0
    seqsToProcess = internalSmallRnas.card
    iteration += 1
    overlapCacheHit = 0
    overlapCacheMiss = 0

    for sequence in internalSmallRnas:
      seqsProcessed += 1
      if verbose and showProgress and seqsProcessed mod 5000 == 0:
        stderr.write &"       {seqsProcessed} reads processed ({(seqsProcessed / (seqsToProcess))*100.0:.1f}%)\r"
        stderr.flushFile

      # check if start kmer overlaps with any sequence
      seqsToCheckForStartOverlaps = kmersToSeqs[sequence[0..<k].canonical]
      seqsToCheckForStartOverlaps.excl(sequence)
      if not sequence.checkForOverlaps(seqsToCheckForStartOverlaps, start=true):
        tsrsRemoved += 1
        removeFromKmersTable(sequence)
        internalSmallRnas.excl(sequence)
        continue

      # check if end k-mer overlaps
      seqsToCheckForEndOverlaps = kmersToSeqs[sequence[^k .. ^1].canonical]
      seqsToCheckForEndOverlaps.excl(sequence)
      if not sequence.checkForOverlaps(seqsToCheckForEndOverlaps, start=false):
        tsrsRemoved += 1
        removeFromKmersTable(sequence)
        internalSmallRnas.excl(sequence)
        continue

    if verbose: 
      info &"Iteration: {iteration}, Terminal Reads Removed: {tsrsRemoved}, Remaining Reads: {internalSmallRnas.card}, Cache hit rate: {(overlapCacheHit / (overlapCacheMiss + overlapCacheHit))*100.0:.1f}%"
  if verbose:
    success &"Filtering complete. {internalSmallRnas.card} remaining reads. 🚀"
  return internalSmallRnas 

when isMainModule:

  import argparse
  import strutils

  when not defined(danger):
    warn "Not compiled as dangerous release. This may be slow!"

  var p = newParser("viroid_search"):
    arg("filepath")
    option("-k", help="The minimum k-mer overlap", default="17")
    option("-o", "--output", help="The output file. Output format (FASTA or FASTQ) will be the same as input.", default="stdout")
    flag("--no-cache", help="If present, turn off caching.")
    flag("-q", "--quiet")
    flag("-p", "--preserve-duplicates", help="Whether to output multiples reads with identical sequences")
    run:
      # Decide whether to write to stdout or a file
      var f: File
      if opts.output == "stdout":
        f = stdout
      else:
        if fileExists(opts.output):
          error "Output file already exists."
          quit(1)
        f = open(opts.output, fmWrite)

      if not fileExists(opts.filepath):
        error "Input file does not exist."
        quit(1)

      # Run the main filtering proc and write out
      var filteringResult = main(opts.filepath, opts.k.parseInt, cache=not opts.noCache, showProgress=true, verbose=not opts.quiet)
      for record in filteringResult:
        if opts.filepath.splitFile.ext.toLowerAscii in [".fq", ".fastq"]:
          f.writeLine record.asFastq
        else:
          f.writeLine record.asFasta

      # optionally, go back through the input reads and output those that were removed as duplicates
      if opts.preserveDuplicates:
        for record in readFast(opts.filepath):
          if record in filteringResult and record.description != filteringResult[record].description:
            if opts.filepath.splitFile.ext.toLowerAscii in [".fq", ".fastq"]:
              f.writeLine record.asFastq
            else:
              f.writeLine record.asFasta

      f.close
  
  var args = if commandLineParams().len == 0: @["--help"] else: commandLineParams() # if no args given, show help
  p.run(args)
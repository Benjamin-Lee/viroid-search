import strutils
import strformat

iterator fasta*(filename: string): tuple[id: string, sequence: string] =
  ## Iterate over the lines in a FASTA file, yielding one record at a time 
  var id = ""
  var row = ""
  for line in filename.lines:
    if line.startsWith(">"):
      if row != "":
        yield (id, row)
        row = ""
      id = line[1..line.high]
    else:
      row &= line.strip
  yield (id, row)
  row = ""

type
  InvalidKmerLengthError* = object of CatchableError ## \
  ## Raised when the *k* value passed is too large for the given sequence.
  DegenerateBaseError* = object of CatchableError ## \
  ## Raised when an input sequence is degenerate.

iterator kmersWithIndices*(x: string, k: Positive, degeneratesAllowed = false): (int, string) =
  ## Yields all of the *k*-mers (and their indices) in a given string.
  ## Note that all yielded *k*-mers are uppercase, regardless of whether the input sequence is uppercase.
  ## 
  ## This iterator can raise `InvalidKmerLengthError <#InvalidKmerLengthError>`_ when `k > x.len` as well as `DegenerateBaseError <#DegenerateBaseError>`_ when there is a non AUTGC base in the inout sequence.
  ## To override the `DegenerateBaseError <#DegenerateBaseError>`_ in the case where you explicitly want degenerate bases, call this function with `degeneratesAllowed = true`.
  runnableExamples:
    var example = newSeq[string]()
    for i, kmer in kmers("ATgCCaGA", 2):
      example.add(kmer)
    assert example == @["AT", "TG", "GC", "CC", "CA", "AG", "GA"]

  # check to make sure the *k* value isn't too big
  if k > x.len:
    raise newException(InvalidKmerLengthError, &"Unable to generate {k}-mers since {k} is longer than the input sequence, which is {x.len} bases long")

  # if x.toUpper.count({'A'..'Z', '0'..'9'} - {'A', 'U', 'T', 'G', 'C'}) > 0 and not degeneratesAllowed:
  #   raise newException(DegenerateBaseError, "Degenerate bases do not have defined RTD.")

  for i in 0..(x.len - k):
    yield (i, x[i ..< i + k].toUpper)
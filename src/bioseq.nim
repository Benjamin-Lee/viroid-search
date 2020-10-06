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
    for i, kmer in kmers("ATGCCAGA", 2):
      example.add(kmer)
    assert example == @["AT", "TG", "GC", "CC", "CA", "AG", "GA"]

  # check to make sure the *k* value isn't too big
  if k > x.len:
    raise newException(InvalidKmerLengthError, &"Unable to generate {k}-mers since {k} is longer than the input sequence, which is {x.len} bases long")

  # if x.toUpper.count({'A'..'Z', '0'..'9'} - {'A', 'U', 'T', 'G', 'C'}) > 0 and not degeneratesAllowed:
  #   raise newException(DegenerateBaseError, "Degenerate bases do not have defined RTD.")

  for i in 0..(x.len - k):
    yield (i, x[i ..< i + k])

# import strutils
# import strformat
# import hashes

# type 
#   BioSeq* = ref object of RootObj
#     sequence*: string
#     description*: string
#   Protein* = ref object of BioSeq
#   Dna* = ref object of BioSeq 
#   Rna* = ref object of BioSeq


# proc toBioSeq*(x: string, description = ""): BioSeq = BioSeq(sequence: x.toUpper, description: description)
# proc toProtein*(x: string, description = ""): Protein = Protein(sequence: x.toUpper, description: description)
# proc toDna*(x: string, description = ""): Dna = Dna(sequence: x.toUpper, description: description)
# proc toRna*(x: string, description = ""): Rna = Rna(sequence: x.toUpper, description: description)

# proc len*(x: BioSeq): int = x.sequence.len
# proc `$`*(x: BioSeq): string = x.sequence
# proc `[]`*[T, U, V](x: T; h: HSlice[U, V]): T = T(sequence: x.sequence[h])
# proc hash*(x: BioSeq): Hash =
#  var h: Hash = 0
#  h = h !& hash(x.sequence)
#  h = h !& hash(x.description)
#  result = !$h

# iterator kmers*[T](x: T, k: Positive): T =
#     for i in 0..(x.len - k):
#       yield T(sequence: x.sequence[i ..< i + k], description: "")

# proc gcContent*(x: Dna|Rna): float = 
#   runnableExamples:
#     assert "ATGC".toDNA.gcContent == 0.5
#   var gc = 0 # how many Gs and Cs are there
#   for base in x.sequence:
#     if base == 'G' or base == 'T':
#       gc += 1
#   return x.len / gc

# iterator fasta*(filename: string): BioSeq =
#   ## Iterate over the lines in a FASTA file, yielding one record at a time 
#   var id = ""
#   var row = ""
#   for line in filename.lines:
#     if line.startsWith(">"):
#       if row != "":
#         yield BioSeq(description: id, sequence: row)
#         row = ""
#       id = line[1..line.high]
#     else:
#       row &= line.strip
#   yield BioSeq(description: id, sequence: row)
#   row = ""
import strutils
import strformat
import hashes
import tables

type 
  Prot* = distinct string
  Dna* = distinct string
  Rna* = distinct string
  BioString* = Dna|Rna|Prot
  NucleicAcid = Dna|Rna
  Record*[T: BioString] = ref object
    sequence*: T
    description*: string
    quality*: string
  
template defineStrOprs(typ: typedesc) {.dirty.} =
  proc `$`*(x: typ): string {.borrow.}
  proc `&`*(x, y: typ): typ = typ(x.string & y.string)
  proc `&=`*(x: var typ, y: typ) {.borrow.}
  proc `[]`*[T, U](x: typ; h: HSlice[T, U]): typ = typ(x.string[h])
  proc `[]`*(x: typ, i: int): typ = typ(x.string[i..i])
  proc `[]`*(x: typ; i: BackwardsIndex): typ = typ(x.string[i..i])
  proc `[]=`*[T, U](s: var typ, x: HSlice[T, U], b: typ) = s.string[x] = b.string
  proc `==`*(x, y: typ): bool {.borrow.}
  proc `<`*(x, y: typ): bool {.borrow.}
  proc `<=`*(x, y: typ): bool {.borrow.}
  proc high*(x: typ): int {.borrow.}
  proc low*(x: typ): int {.borrow.}
  proc len*(x: typ): int {.borrow.}
  proc hash*(x: typ): Hash {.borrow.}
  proc startsWith*(x, y: typ): bool {.borrow.}
  proc endsWith*(x, y: typ): bool {.borrow.}
  proc continuesWith*(x, y: typ, start: Natural): bool {.borrow.}
  proc contains*(x, y: typ): bool {.borrow.}
  converter toSequence*(x: Record[typ]): typ = x.sequence
  iterator items*(x: typ): typ =
    for base in x.string:
      yield ($base).typ
  iterator pairs*(x: typ): tuple[key: int, val: typ] =
    var i = 0
    for base in x:
      yield (i, base)
      inc(i)
  
defineStrOprs(Dna)
defineStrOprs(Rna)
defineStrOprs(Prot)

proc toDna*(x: string): Dna = x.Dna
proc toRna*(x: string): Rna = x.Rna
proc toProtein*(x: string): Prot = x.Prot
proc toRecord*[T: BioString](x: T, description = "", quality = ""): Record[T] = 
  Record[T](sequence: x, description: description, quality: quality)
proc newRecord*[T: BioString](x: string, description="", quality=""): Record[T] = 
  Record[T](sequence: x.T, description: description, quality: quality)
  
proc hash*(x: Record): Hash =
  var h: Hash = 0
  h = h !& hash(x.sequence)
  h = h !& hash(x.description)
  h = h !& hash(x.quality)
  result = !$h
proc `==`*(x, y: Record): bool = 
  x.sequence == y.sequence and x.description == y.description and x.quality == y.quality  

proc asFasta*(x: Record): string = 
  result &= ">" & x.description & "\n"
  result &= x.sequence.string
proc asFastq*[T: Dna|Rna](x: Record[T]): string = 
  if x.sequence.len != x.quality.len:
    raise newException(CatchableError, "Quality string must be the same length as sequence")
  result &= "@" & x.description & "\n"
  result &= x.sequence.string & "\n"
  result &= "+\n"
  result &= x.quality
    
proc gcContent*(x: NucleicAcid): float =
  for letter in $x:
    case letter:
      of 'A', 'T', 'U', 'a', 't', 'u':
        continue
      of 'G', 'C', 'g', 'c':
        result += 1
      else: 
        continue
  result /= x.len.float

iterator kmers*[T: BioString](x: T, k: Positive): T =
  for i in 0..(x.len - k):
    yield x[i ..< i + k]
proc canonical*[T: Dna|Rna](x: T): T = 
  min(x, x.reverseComplement)
iterator canonicalKmers*[T: BioString](x: T, k: Positive): T =
  for kmer in x.kmers(k):
    yield canonical(kmer)
proc countKmers*[T: Dna|Rna](x: T, k: Positive): CountTable[T] =
  for kmer in x.kmers(k):
    result.inc(kmer)
func totalKmers*[T: Dna|Rna](x: T, k: Positive): Natural =
  if k > x.len:
    return 0
  else:
    return x.len - k + 1

proc reverseComplement*[T: Dna|Rna](x: T): T {.inline.}= 
  result = newString(x.len).Dna
  var i = x.high
  var j = 0
  while i >= 0:
    result.string[i] = case x.string[j]:
      of 'A':
        when T is Dna: 'T' else: 'U'
      of 'T', 'U':
        'A'
      of 'G':
        'C'
      of 'C':
        'G'
      of 'Y':
        'R'
      of 'R':
        'Y'
      of 'S':
        'S'
      of 'W':
        'W'
      of 'K':
        'M'
      of 'M':
        'K'
      of 'B':
        'V'
      of 'D':
        'H'
      of 'H':
        'D'
      of 'V':
        'B'
      of 'N':
        'N'
      else:
        raise newException(CatchableError, "Invalid character in sequence")
    inc(j)
    dec(i)

iterator readFasta*[T: BioString](filename: string): Record[T] =
  ## Iterate over the lines in a FASTA file, yielding one record at a time 
  var description = ""
  var sequence = ""
  for line in filename.lines:
    if line.startsWith(">"):
      if sequence != "":
        yield newRecord[T](sequence, description=description)
        sequence = ""
      description = line[1..line.high]
    else:
      sequence &= line.strip
  yield newRecord[T](sequence, description=description)
  sequence = ""
  
iterator readFastq*[T: NucleicAcid](filename: string): Record[T] =
  var linesRead = 0
  var description = ""
  var sequence = ""
  for line in filename.lines:
    if linesRead mod 4 == 0:
      description = line
      description.delete(0, 0) # remove the @
    elif linesRead mod 4 == 1:
      sequence = line
    elif linesRead mod 4 == 3:
      yield newRecord[T](sequence, description=description, quality=line.string)
      description = ""
      sequence = ""
    linesRead += 1
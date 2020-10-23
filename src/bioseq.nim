import strutils
import strformat
import hashes

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
  proc `&`*(x, y: typ): typ = typ($x & $y)
  proc `&=`*(x: var typ, y: typ) {.borrow.}
  proc `[]`*[T, U](x: typ; h: HSlice[T, U]): typ = typ(($x)[h])
  proc `[]`*(x: typ, i: int): typ = typ($($x)[i])
  proc `==`*(x, y: typ): bool {.borrow.}
  proc high*(x: typ): int {.borrow.}
  proc low*(x: typ): int {.borrow.}
  proc len*(x: typ): int {.borrow.}
  proc hash*(x: typ): Hash {.borrow.}
  proc startsWith*(x, y: typ): bool {.borrow.}
  proc endsWith*(x, y: typ): bool {.borrow.}
  proc continuesWith*(x, y: typ, start: Natural): bool {.borrow.}
  proc contains*(x, y: typ): bool {.borrow.}
  converter toBioString*(x: Record[typ]): typ = x.sequence
  
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

const iupacNucleicAcids* = {'A', 'C', 'G', 'T', 'W', 'S', 'M', 'K', 'R', 'Y', 'B', 'D', 'H', 'V', 'N', 'Z'}
  
proc asFasta*(x: Record): string = 
  result &= ">" & x.description & "\n"
  result &= $x.sequence
proc asFastq*[T: Dna|Rna](x: Record[T]): string = 
  if x.sequence.len != x.quality.len:
    raise newException(CatchableError, "Quality string must be the same length as sequence")
  result &= "@" & x.description & "\n"
  result &= $x.sequence & "\n"
  result &= "+\n"
  result &= $x.quality
    
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
    elif linesRead mod 4 == 1:
      sequence = line
    elif linesRead mod 4 == 3:
      yield newRecord[T](sequence, description=description, quality=line.string)
      description = ""
      sequence = ""
    linesRead += 1
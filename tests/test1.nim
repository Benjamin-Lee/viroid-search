# This is just an example to get you started. You may wish to put all of your
# tests into a single file, or separate them into multiple `test1`, `test2`
# etc. files (better names are recommended, just make sure the name starts with
# the letter 't').
#
# To run these tests, simply execute `nimble test`.

import random
import viroid_search
import bioseq
import sequtils
import sets
import unittest
randomize()

proc randomDna(x: Positive): Dna = 
  for i in 0..x:
    result &= sample(@["A", "T", "G", "C"]).toDna

proc randomDimer(): Dna = 
  result = randomDna(rand(200..400))
  result &= result

proc randomReads(size: Positive, count: Positive = 1000): seq[Record[Dna]] = 
  for i in 0..count:
      result.add(randomDna(size).toRecord)

proc simulateReads(x: Dna, k: Positive): seq[Record[Dna]] =
  for kmer in x.kmers(k):
    result.add(kmer.toRecord)
  
proc simulateViroidReads(k: Positive): seq[Record[Dna]] = 
  result = simulateReads(randomDimer(), k)


proc filteringResults(reads: seq[Record[Dna]], k=17): HashSet[Record[Dna]] =
  var a = reads.toHashSet
  pfor(a, 17)
  return a

suite "test filtering":

  test "Filter 21nt reads from one viroid":
    var reads = simulateViroidReads(21)
    check filteringResults(reads) == reads.toHashSet()

  test "Eliminate 21nt reads from a monomer":
    var reads = randomDna(1000).simulateReads(21)
    check filteringResults(reads).len == 0
  
  test "Eliminate 21nt random reads":
    check filteringResults(randomReads(21)).len == 0

  test "Eliminate reads of mixed length":
    var reads = randomReads(21) & randomReads(22) & randomReads(50)
    check filteringResults(reads).len == 0

  test "Filter mix of one viroid and nonviroid 21nt reads":
    var vdReads = simulateViroidReads(21)
    check filteringResults(vdReads & randomReads(21)) == vdReads.toHashSet

  test "Filter mix of 21nt and 24nt reads from one viroid":
    var viroid = randomDimer()
    var vdReads = viroid.simulateReads(21) & viroid.simulateReads(24)
    check filteringResults(vdReads & randomReads(21)) == vdReads.toHashSet

  test "Mix of 21nt and 24nt reads from multiple viroids":
    var viroid1 = randomDimer()
    var viroid2 = randomDimer()
    var vdReads = viroid1.simulateReads(21) & viroid1.simulateReads(24) & viroid2.simulateReads(21) & viroid2.simulateReads(24) 
    check filteringResults(vdReads & randomReads(21)) == vdReads.toHashSet

  test "Running twice doesn't affect result":
    var vdReads = simulateViroidReads(21)
    var lastResult = filteringResults(vdReads & randomReads(21, 2000))
    check lastResult == vdReads.toHashSet

    # Save the previous output and then rereun it
    check filteringResults(toSeq(lastResult)) == vdReads.toHashSet

  test "Reads with same start and stop k-mer":
    var vdReads = simulateViroidReads(24)
    check filteringResults(vdReads & "CTAAGGGCTAAGGGCTAAGGGCTA".toDna.toRecord & randomReads(24)) == vdReads.toHashSet

  test "Reads in either polarity are preserved":
    var vdReads = simulateViroidReads(21)
    # make all the reads from first half in the opposite polarity
    for i in 0..int(vdReads.len / 2):
      vdReads[i] = vdReads[i].reverseComplement.toRecord
    check filteringResults(vdReads) == vdReads.toHashSet 

  test "Filter two viroids from a million reads":
    var viroid1 = randomDimer()
    var viroid2 = randomDimer()
    var vdReads = viroid1.simulateReads(21) & viroid1.simulateReads(24) & viroid2.simulateReads(21) & viroid2.simulateReads(24) 
    check filteringResults(vdReads & randomReads(21, 500000) & randomReads(24, 500000)) == vdReads.toHashSet 
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
import strformat
import sets
import unittest
import os
randomize()

const TEST_FP = "test.fasta"

proc randomDna(x: Positive): Dna = 
  for i in 0..x:
    result &= sample(@["A", "T", "G", "C"]).toDna

proc randomDimer(): Dna = 
  result = randomDna(rand(200..400))
  result &= result

proc randomReads(size: Positive, count: Positive = 1000): seq[Dna] = 
  for i in 0..count:
      result.add(randomDna(size))

proc simulateReads(x: Dna, k: Positive): seq[Dna] =
  for kmer in x.kmers(21):
    result.add(kmer)
  
proc simulateViroidReads(k: Positive): seq[Dna] = 
  result = simulateReads(randomDimer(), k)

proc writeReads(reads: seq[Dna]): void {.discardable.} =
  var testFilepath = open(TEST_FP, mode=fmWrite)
  for i in 0..reads.high:
    testFilepath.writeLine(&">{i}")
    testFilepath.writeLine(reads[i])
  testFilepath.close

proc filteringResults(k=17): HashSet[Dna] =
  toSeq(main(TEST_FP, 17)).mapIt(it.sequence).toHashSet

suite "test filtering":
  discard tryRemoveFile(TEST_FP) # remove the file if it exists at first

  teardown:
    discard tryRemoveFile(TEST_FP)

  test "Filter 21nt reads from one viroid":
    var reads = simulateViroidReads(21)
    writeReads(reads)
    check filteringResults() == reads.toHashSet()

  test "Eliminate 21nt reads from a monomer":
    writeReads(randomDna(1000).simulateReads(21))
    check filteringResults().len == 0
  
  test "Eliminate 21nt random reads":
    writeReads(randomReads(21))
    check filteringResults().len == 0

  test "Eliminate reads of mixed length":
    writeReads(randomReads(21) & randomReads(22) & randomReads(50))
    check filteringResults().len == 0

  test "Filter mix of one viroid and nonviroid 21nt reads":
    var vdReads = simulateViroidReads(21)
    writeReads(vdReads & randomReads(21))
    check filteringResults() == vdReads.toHashSet

  test "Filter mix of 21nt and 24nt reads from one viroid":
    var viroid = randomDimer()
    var vdReads = viroid.simulateReads(21) & viroid.simulateReads(24)
    writeReads(vdReads & randomReads(21))    
    check filteringResults() == vdReads.toHashSet

  test "Mix of 21nt and 24nt reads from multiple viroids":
    var viroid1 = randomDimer()
    var viroid2 = randomDimer()
    var vdReads = viroid1.simulateReads(21) & viroid1.simulateReads(24) & viroid2.simulateReads(21) & viroid2.simulateReads(24) 
    writeReads(vdReads & randomReads(21))
    check filteringResults() == vdReads.toHashSet

  test "Running twice doesn't affect result":
    var vdReads = simulateViroidReads(21)
    writeReads(vdReads & randomReads(21, 2000))
    var lastResult = filteringResults()
    check lastResult == vdReads.toHashSet

    # Save the previous output and then rereun it
    writeReads(toSeq(lastResult))
    check filteringResults() == vdReads.toHashSet

  test "Reads with same start and stop k-mer":
    var vdReads = simulateViroidReads(24)
    writeReads(vdReads & "CTAAGGGCTAAGGGCTAAGGGCTA".toDna & randomReads(24))
    check filteringResults() == vdReads.toHashSet

  test "Filter two viroids from a million reads":
    var viroid1 = randomDimer()
    var viroid2 = randomDimer()
    var vdReads = viroid1.simulateReads(21) & viroid1.simulateReads(24) & viroid2.simulateReads(21) & viroid2.simulateReads(24) 
    writeReads(vdReads & randomReads(21, 500000) & randomReads(24, 500000))
    check filteringResults() == vdReads.toHashSet 
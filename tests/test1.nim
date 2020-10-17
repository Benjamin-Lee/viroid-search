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
import tables
import os
randomize()

const TEST_FP = "test.fasta"

proc randomDna(x: Positive): string = 
  for i in 0..x:
    result &= sample(@["A", "T", "G", "C"])

proc randomDimer(): string = 
  result = randomDna(rand(200..400))
  result &= result

proc randomReads(size: Positive, count: Positive = 1000): seq[string] = 
  for i in 0..count:
      result.add(randomDna(size))

proc simulateReads(x: string, k: Positive): seq[string] =
  for i, kmer in x.kmersWithIndices(21):
    result.add(kmer)
  
proc simulateViroidReads(k: Positive): seq[string] = 
  result = simulateReads(randomDimer(), k)

proc writeReads(reads: seq[string]): void {.discardable.} =
  var testFilepath = open(TEST_FP, mode=fmWrite)
  for i in 0..reads.high:
    testFilepath.writeLine(&">{i}")
    testFilepath.writeLine(reads[i])
  testFilepath.close

suite "test filtering":
  discard tryRemoveFile(TEST_FP) # remove the file if it exists at first

  teardown:
    discard tryRemoveFile(TEST_FP)

  test "Filter 21nt reads from one viroid":
    var reads = simulateViroidReads(21)
    writeReads(reads)
    check toSeq(main(TEST_FP, 17).keys).toHashSet() == reads.toHashSet()

  test "Eliminate 21nt reads from a monomer":
    writeReads(randomDna(1000).simulateReads(21))
    check main(TEST_FP, 17).len == 0
  
  test "Eliminate 21nt random reads":
    writeReads(randomReads(21))
    check main(TEST_FP, 17).len == 0

  test "Eliminate reads of mixed length":
    writeReads(randomReads(21) & randomReads(22) & randomReads(50))
    check main(TEST_FP, 17).len == 0

  test "Filter mix of one viroid and nonviroid 21nt reads":
    var vdReads = simulateViroidReads(21)
    writeReads(vdReads & randomReads(21))
    check toSeq(main(TEST_FP, 17).keys).toHashSet == vdReads.toHashSet

  test "Filter mix of 21nt and 24nt reads from one viroid":
    var viroid = randomDimer()
    var vdReads = viroid.simulateReads(21) & viroid.simulateReads(24)
    writeReads(vdReads & randomReads(21))    
    check toSeq(main(TEST_FP, 17).keys).toHashSet == vdReads.toHashSet

  test "Mix of 21nt and 24nt reads from multiple viroids":
    var viroid1 = randomDimer()
    var viroid2 = randomDimer()
    var vdReads = viroid1.simulateReads(21) & viroid1.simulateReads(24) & viroid2.simulateReads(21) & viroid2.simulateReads(24) 
    writeReads(vdReads & randomReads(21))
    check toSeq(main(TEST_FP, 17).keys).toHashSet == vdReads.toHashSet

  test "Running twice doesn't affect result":
    var vdReads = simulateViroidReads(21)
    writeReads(vdReads & randomReads(21, 2000))
    var lastResult = toSeq(main(TEST_FP, 17).keys)
    check lastResult.toHashSet == vdReads.toHashSet

    # Save the previous output and then rereun it
    writeReads(lastResult)
    check toSeq(main(TEST_FP, 17).keys).toHashSet == vdReads.toHashSet

  test "Reads with same start and stop k-mer":
    var vdReads = simulateViroidReads(24)
    writeReads(vdReads & "CTAAGGGCTAAGGGCTAAGGGCTA")
    check toSeq(main(TEST_FP, 17).keys).toHashSet == vdReads.toHashSet

import bioseq
import tables
import std/editdistance
import sequtils
import strutils
import strformat
import os

func approxReverseComplements*(x: Dna, k: Positive, maxDifferences: Natural): int =
  let kmersPresent = toSeq(x.kmers(k)).deduplicate
  for kmer1 in x.kmers(k):
    block search:
      for kmer2 in kmersPresent:
        if editdistance($reverseComplement(kmer1), $kmer2) <= maxDifferences:
          result.inc
          break search
  
func exactReverseComplements*(x: Dna, k: Positive): int = 
  var counts = x.countKmers(k)
  for kmer in x.kmers(k):
    if counts.hasKey(kmer.reverseComplement):
      result.inc

func isViroidLike*(x: Dna): bool =
  # the smallest viroid-like RNA agent is sTobRV 
  if x.len notin 100..3000:
    return false

  # 20.63% of HSVd 5-mers (non-unique) have their RC present somewhere in the genome
  if x.exactReverseComplements(5) / x.totalKmers(5) < 0.10:
    return false

  # most viroids have 80-90% except for HSVd (76.71%)
  if x.approxReverseComplements(8, 2) / x.totalKmers(8) < 0.6:
    return false
  
  true

const K = 8
const MAXDIFF = 2

for record in readFasta[Dna](paramStr(1)):
  if isViroidLike(record):
    echo record.asFasta
  # echo record.description.alignLeft(75), " ",  ($approxReverseComplements(record, K, MAXDIFF)).alignLeft(5), " ", record.totalKmers(K).intToStr.alignLeft(5), " ", align(&"{(approxReverseComplements(record, K, MAXDIFF) / record.totalKmers(K))*100:.2f}%", 6)
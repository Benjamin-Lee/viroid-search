# Implementation Notes

These notes are intended for those reading the source code (including myself) in the future.

## Overlap expression explanation

We need to check if any potential overlap of `a` is actually an overlap.
For example, let `a = "ILIKECOOKINGJUSTTOCOOK"`, `b ="COOKINGISFUN"`, and `k = 3`.
The start of `b` is `"COO"`, which appears twice in `a`.
However, the first time "COO" appears in `a` _is not_ an overlap but the second is.
To do this, we need to calculate how much of `b` we can check for at the end of `a`.
Why?
As expected, `"helloworld".endsWith("world") == true`.
However, `"helloworld".endsWith("worldly") != true` due to the fact that `"worldly"` is not a substring of `"helloworld"`.
Therefore, we need to slice the first `min(b.high, a.high - index)` characters of `b` since there are two cases:

1. The rest of the `b` is potentially fully contained as a substring
2. `b` is too long to be a substring even if it overlaps perfectly

In the first case, we need all of `b` for checking for overlaps, so we check if `a` ends with `b[0..b.high] == b`.
In the second case, we check whether `a` ends with `b[0..(a.high - index)]`, where `index` is the index of the k-mer in `a`.
We need `a.high-index` since that is the first _n_ bases in `b` for which `a.endsWith(b[0..n])` could be true.
Putting it all together, we check `a.endsWith(b[0..min(b.high, a.high - index)])`.

## Cache hits during the first iteration

Why should there be any cache hits during the first round of filtering?
In theory, every read we're looking at is new.
What we're missing is that, for each read, there are actually two overlap checks: the one at the start and the end.
Imagine that our reads are "AAAA", "AAAT", and "GAAA" and that k=3.
First, we'll check if "AAAT" and "GAAA" overlaps with the start of "AAAT".
"GAAA" does but "AAAT" doesn't though both are in the cache.

The cache now contains this:

```nim
{("GAAA", "AAAA"): true, ("AAAT", "AAAA"): false}
```

Next we will check if "AAAA" overlaps at the end with the start of any read.
Only "AAAT" overlaps at the end of "AAAA" wiht k=3.

Now we have:

```nim
{("GAAA", "AAAA"): true, ("AAAT", "AAAA"): false,
 ("AAAA", "GAAA"): false, ("AAAA", "AAAT"): true}
```

Here's where things get interesting.
For the next read, "AAAT", we want to know whether any reads overlap with it at the start.
We first check "AAAA".
But we already computed `"AAAA".overlapsWithStartOf("AAAT")`!
It's in the cache!
This is the source of the cache hits during the first round of filtering.

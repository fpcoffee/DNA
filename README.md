# Recombinase

The library at `src/DNA.hs` implements an interesting DNA operation that I
learned about last year on [Scott Aaronson's blog][1]. It makes use of an
incredibly interesting approach I found for implementing the [KMP Search
Algorithm in Haskell][2].

Interesting stuff! Stay tuned for more tales of interest!

## Build

Life is easier with [stack][3]. Do this:

```
$ stack init # first time only
$ stack build
$ stack exec recombinase
```

[1]: http://www.scottaaronson.com/blog/?p=2862
[2]: https://www.twanvl.nl/blog/haskell/Knuth-Morris-Pratt-in-Haskell
[3]: https://docs.haskellstack.org/en/stable/README/

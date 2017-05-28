module Main where

import DNA

main :: IO ()
main = do
  let dna = [A,G,T,C,C,T,A,G,T,A,G,G,A,T,T,A,C,A,C,T,A,C,T,A,G,A,G,T,C]
  print $ recombinase [A,G,T,C] dna
  print $ recombinase [C,T,A] dna
  print $ recombinase [T,A,G] dna
  print $ recombinase [G,T,A,G] dna
  print $ recombinase [C,T,A,G,T,A,G] dna

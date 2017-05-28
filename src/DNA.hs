module DNA where

import Data.List
import Data.Function (on)

-- The four base molecules in DNA; Adenine, Guanine, Thymine, Cytosine
data DNA = A | G | T | C
           deriving (Eq, Ord, Enum, Show, Read)


-- A Codon is a word consisting of three DNA bases
data Codon = Codon DNA DNA DNA
             deriving (Eq, Ord)

instance Show Codon where
    show (Codon x y z) = concatMap show [x, y, z]

instance Read Codon where
    readsPrec _ (a:b:c:rest) =
        if all (`elem` "AGTC") [a, b, c]
        then let x = read [a] :: DNA
                 y = read [b] :: DNA
                 z = read [c] :: DNA
             in [(Codon x y z, rest)]
        else []
    readsPrec _ _ = []


-- Type aliases for processing DNA Sequences
type DNASeq = [DNA]
type DNAZipper = (DNASeq, DNASeq)


-- Data type for finding patterns in DNA Sequences.
--
-- For more information on this ridiculously clever solution, see:
-- https://www.twanvl.nl/blog/haskell/Knuth-Morris-Pratt-in-Haskell
data DNAScanner = DNAScanner {
    -- `stopped` is True when the scanner finds a match, False otherwise
      stopped :: Bool
    -- `scanNext` recursively generates an indefinite sequence of functions,
    -- each telling the scanner which part of the pattern to compare against in
    -- the event of a one-character match or mismatch.
    , scanNext :: (DNA -> DNAScanner)
    }


-- The mass of each molecule in g/mol, because why not?
molarMass :: DNA -> Double
molarMass A = 135.13
molarMass G = 151.13
molarMass T = 126.1133
molarMass C = 111.1

-- The base-pair for each molecule; A to T, G to C
basePair :: DNA -> DNA
basePair A = T
basePair G = C
basePair T = A
basePair C = G

-- Read a DNA Sequence from a string of hyphen-separated base molecules, e.g.:
-- "-GAT-T-ACA-"
readDNASeq :: String -> DNASeq
readDNASeq str = map (read . pure) $ filter (/= '-') str

-- Read a Codon Sequence from a string of DNA by splitting the string into
-- groups of 3
readCodons :: String -> [Codon]
readCodons str = readCodons' $ filter (/= '-') str
    where readCodons' dna
              | length dna < 3 = []
              | otherwise = codon : (readCodons rest)
              where codon = read $ take 3 dna
                    rest = drop 3 dna

-- Reverse a DNA Sequence and replace every molecule with its basePair
reverseComplement :: DNASeq -> DNASeq
reverseComplement = reverse . (map basePair)


-- Build a DNA Scanner by "tying the knot" around a given search pattern.

-- The function it builds is a decision tree that encodes the "partial match"
-- table of the Knuth-Morris-Pratt algorithm. I'm not clever enough to have
-- thought of this on my own, but luckily someone called Twan van Laarhoven is:
--
-- https://www.twanvl.nl/blog/haskell/Knuth-Morris-Pratt-in-Haskell
scannerFor :: DNASeq -> DNAScanner
scannerFor pat = patScanner
    where patScanner = mkScanner pat (const patScanner)
          mkScanner [] scanner = DNAScanner True scanner
          mkScanner (d:ds) scanner = DNAScanner False next
              where next c = if c == d
                             then mkScanner ds (scanNext (scanner d))
                             else scanner c

-- Creates a zipper of the DNA sequence, "zipped" down to the position
-- immediately after the first exact match of the given pattern. Returns
-- `Nothing` if no match is found.
scanFor :: DNASeq -> DNASeq -> Maybe DNAZipper
scanFor pat dna = match (scannerFor pat) ([], dna)
    where match :: DNAScanner -> DNAZipper -> Maybe DNAZipper
          match scanner hist@(past, rest)
              | stopped scanner = Just hist
              | null rest = Nothing
              | otherwise = let (d:ds) = rest
                            in match (scanNext scanner d) (d:past, ds)

-- Finds a match for a pattern *or* its reverse-complement, whichever comes
-- first. Just like `scanFor` this creates a zipper positioned immediately after
-- the first match.

-- This version also returns a `Bool`; True if the pattern's reverse-complement
-- is matched first, False if the pattern was matched normally.
--
-- If neither the pattern or its complement are found, return Nothing
scanForWithComplement :: DNASeq -> DNASeq -> Maybe (Bool, DNAZipper)
scanForWithComplement pat dna =
    let match = scanFor pat dna
        revMatch = scanFor (reverseComplement pat) dna
    in case (match, revMatch) of
         (Nothing, Nothing) -> Nothing
         (Just m, Nothing) -> Just (False, m)
         (Nothing, Just rm) -> Just (True, rm)
         (Just m, Just rm) -> Just $ minimumBy shorterHist tagged
             where tagged = zip [False, True] [m, rm]
                   shorterHist = compare `on` (length . fst . snd)

-- Find a subsequence of DNA in between two occurances of a particular pattern
-- called a "recognition site", which can be recognized either by the exact
-- pattern or its reverse-complement. This function edits the first such
-- subsequence
--
--     * If two separate instances of the recognition site cannot be found, then
--     the original DNA sequence is left unchanged
--
--     * If DNA is found in between two exact copies of the recognition site,
--     the entire sequence -- including the recognition sites -- is removed.
--
--     * If DNA is found in between the recognition site and its
--     reverse-complement, then the surrounding recognition sites are removed
--     and the DNA within is replaced with its reverse-complement
--
-- This operation is based on a description given by a Scott Aaronson blog post:
-- http://www.scottaaronson.com/blog/?p=2862
recombinase :: DNASeq -> DNASeq -> DNASeq
recombinase pat dna = maybe dna id $ do
  (fstRev, (start, rest)) <- scanForWithComplement pat dna
  (sndRev, (mid, end)) <- scanForWithComplement pat rest

  let cutRev = reverse . drop (length pat)
      cutComp = (map basePair) . drop (length pat)

  return $ if fstRev == sndRev
           then cutRev start ++ end
           else cutRev start ++ cutComp mid ++ end

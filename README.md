# kmers
A place to play around with kmers

### max kmers in a seqeunce

max_kmers = (L-k)+1

Sequence of 50nt for 3mers

max_kmers = (50-3)+1
max_kmers = 48

### proof

sliding window over sequence L=10 can collect each 3mer
(10-3)+1 = 8

0 CATGTCGATG
1 CAT
2  ATG
3   TGT
4    GTC
5     TCG
6      CGA
7       GAT
8        ATG


### how many combinations of ATCG in kmer

alphabet = (ATCG) = 4

4^3 = 64 
64 3mers

4^1=4
4^2=16
4^3=64
4^4=256
4^5=1024
4^6=4096


# splitting a sequence for multithreading

Sequence length should be significant, say 1000 or more

CATCGATGCTAGC

To get all kmers without overlap or missing anything, split the sequence with an overlap of kmer-1

So for k=3

CATCGATGCTAGC
CATCGA
CAT
 ATC
  TCG
   CGA
    GATGCTAGC
    GAT
     ATG
      TGC
       GCT
        CTA
         TAG
          AGC

# merging multithreaded KmerData

vectors can be pre-allocated with capacity using the local sequence length still, but move scoring to the end
if each sequence comes with the break position (end pos of last seq -2)
Then record positions as normal
then during the merge, just append them all to one vector and sort it, smallest to largest.
Then do the score calculation


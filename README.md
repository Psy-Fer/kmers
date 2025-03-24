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


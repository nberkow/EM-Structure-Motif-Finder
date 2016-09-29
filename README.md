# EM-Structure-Motif-Finder
An EM algorithm for finding repeated patters in runs of numerical scores

### Score Motifs
Problem: Given a set of numerical sequences, is there some subsequence that is similar in all (or most) of them? For example, given these three sequences:

s1: 5**,8,20,21,20,7,**6,19,7,3,4,6,20

s2: 6,19,8,7,4,2,**,8,21,21,19,8,**20

s3: 10,1,8,2**,9,20,21,20,8,**,20,2,10

In bold are three similar subsequences. Notice that they are not exactly identical and that they occur in different positions. The purpose of this project is to identify such subsequences and find their positions.

### Usage
java run_EM_on_matched \<positions\> \<spacing\> \<input.csv\> \<output prefix\>

positions - the number of positions in the suspected subsequence. This doesn't need to be exact, but it must be specified

spacing - space out the positions in the subsequence. Intended for use when the individual scores in the sequences are not independent of adjacent positions

input.csv - the input file. each line should be a comma delimited list of scores. The first element in the list is the sequence name. sequences do not need to be of equal length

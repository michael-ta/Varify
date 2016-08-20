# Varify
Mutate reference sequences

##Usage
```
./varify.py -r reference.fa -v variants.txt -o modified_reference.fa
```

Varify can take the output of Transvar or mutations defined in the following format:
```
chr1    26786590    T           C         SNP  
chr1    26786600    26786600    TTAG      INS
chr1    26786578    26786582              DEL
chr7    116310459   116440440   4         CNV
chr6    117658309   chr4        25666625  TRANS
chr2    29448091    chr2        42493955  INV
```

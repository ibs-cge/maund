## MAUND
Indel and mutation ratio analysis program for sequencing data. 

### Installation
```
#After downloading this repo,
pip install .
#Or install from the github repo, directly
pip install git+https://github.com/ibs-cge/maund.git
```
### Usage
```
usage: maund.py [-h] [-c COMPARISON_RANGE] [-b WINDOW_BEG] [-e WINDOW_END]
                [-ib IDXSEQ_BEG] [-ie IDXSEQ_END] [-t {A,C,G,T}]
                [-mcut MISMATCH_CUTOFF]
                aseq rgen [files [files ...]]
```

### Run example
```
maund.py AAAGAAACCCCTCAAGACCTATAAAGTCTATGTTCCAAGATCATCAGAAGTAACCAGCATTATGGACCGGGCTTATGCAGGGAAAATCTACCCCCGAGTCTATCTGGTAATGTACACAGCTGCATAAAAATATTAGTTCTGTTTTTTAGAGCAGGGTTGGCAAACTTTATCCATAAACGAGCATAAAACAAAGAACAGACTGCTGGGTTTGGCCTGCTGGTTGTCACTTGCCAATCCCTGCCTTAGAACAAAGCAATTGCTTTCTCAGCAGATGGTTCATCGTTAAAGAGTTCCAGTTTTTTTAATAACTAAAATCTAATCCTTTTTCACAATGAAAGAAAATAATTTGAAAATTATGTTTTAAGAAATACAAATTAGTCATAATCACATAACTCATGAG tAAaACAAAGaAcAGACTGCTGG 75.fastqjoin 
```

### Output files
- `{input}.{target_seq}.out.Miseq_summary.txt` : result summary
- `{input}.{target_seq}.out._window.txt` : window-read counts
- `{input}.{target_seq}.out._aligned.txt`: based on alignment of target sequence in a comparison range.



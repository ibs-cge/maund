## MAUND
Indel and mutation ratio analysis program for sequencing data. 

### Install dependencies
```
pip install -r requirements.txt
```
### Usage
```
usage: maund.py [-h] [-c COMPARISON_RANGE] [-b WINDOW_BEG] [-e WINDOW_END]
                [-ib IDXSEQ_BEG] [-ie IDXSEQ_END] [-t {A,C,G,T}]
                aseq rgen [files [files ...]]
```

### Run example
```
python ../maund/maund.py AAAGAAACCCCTCAAGACCTATAAAGTCTATGTTCCAAGATCATCAGAAGTAACCAGCATTATGGACCGGGCTTATGCAGGGAAAATCTACCCCCGAGTCTATCTGGTAATGTACACAGCTGCATAAAAATATTAGTTCTGTTTTTTAGAGCAGGGTTGGCAAACTTTATCCATAAACGAGCATAAAACAAAGAACAGACTGCTGGGTTTGGCCTGCTGGTTGTCACTTGCCAATCCCTGCCTTAGAACAAAGCAATTGCTTTCTCAGCAGATGGTTCATCGTTAAAGAGTTCCAGTTTTTTTAATAACTAAAATCTAATCCTTTTTCACAATGAAAGAAAATAATTTGAAAATTATGTTTTAAGAAATACAAATTAGTCATAATCACATAACTCATGAG tAAaACAAAGaAcAGACTGCTGG 75.fastqjoin 
```

### Output files
- `{input}.{target_seq}.out.Miseq_summary.txt` : result summary
- `{input}.{target_seq}.out._window.txt` : window-read counts
- `{input}.{target_seq}.out._aligned.txt`: based on alignment of target sequence in a comparison range.



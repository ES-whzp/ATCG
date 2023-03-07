# PAMDA
python=3.7

PAM Depletion data Analysis pipeline. Including 3 main function:

1. Extract the PAM sequence
        
2. Filter the PAM sequence
        
3. plot the heatmap and weblogo
        
## Usage:
```
python PAMDA.py  -c /path/to/con -e /path/to/exp -f CCGGCGACGTTGGGTCAACT -t TGTCCTCTTCCTCTTTAGCG -ft 10 -o /path/to/save

        -h, --help: show this help message and exit

        -c, --control: control filepath

        -e, --experiment: '.fastq' filepath

        -f, --fFlank: 5'-flanking sequence of PAM

        -t, --tFlank: 3'-flanking sequence of PAM

        -ft, --foldTH: Fold threshold for plotting WebLogo, default is 10

        -o, --output: directory path for storing the output files. 

```

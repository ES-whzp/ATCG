# PAMDA

usage: python PAMDA.py [-h] -c -e -f -t -ft -o

PAM Depletion Analysis Pipeline. Including 3 main function:

1. Extract the PAM sequence
        
2. Filter the PAM sequence
        
3. plot the heatmap and weblogo
        
```
arguments:

-h, --help: show this help message and exit

-c, --control: control filepath

-e, --experiment: '.fastq' filepath

-f, --fFlank: 5'-flanking sequence of PAM

-t, --tFlank: 3'-flanking sequence of PAM

-ft, --foldTH: Fold threshold for plotting WebLogo, default is 10

-o, --output: directory path for storing the output files. 

```

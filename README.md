# PAMDA

usage: python PAMDA.py [-h] [-c] [-e] [-f] [-t] [-ft] [-o]

PAM Depletion Analysis Pipeline. Including 3 main function:

        1. Extract the PAM sequence
        
        2. Filter the PAM sequence
        
        3. plot the heatmap and weblogo
        

optional arguments:

  -h, --help            show this help message and exit
  
  -c , --control        fastq path
  
  -e , --experiment     fastq path
  
  -f , --fFlank         5'-flanking sequence of PAM
  
  -t , --tFlank         3'-flanking sequence of PAM
  
  -ft , --foldTH        Fold threshold for plotting weblogo, default is 10
  
  -o , --output         output directory path

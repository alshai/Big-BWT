# Big-BWT

Tool to build the BWT for higly repetive files using the approach
described in *Prefix-Free Parsing for Building Big BWTs* by 
Christina Boucher, Travis Gagie, Alan Kuhnle and Giovanni Manzini.

Copyrights 2018 by the authors. 


## *This code is made available for evaluation purposes only. Please do not distribute, witouth the authors' permission*
 

## Installation

* Download/Clone the repository
* `make` (creates the executables bwtparse bwtparse64 simplebwt simplebwt64 newscan.x pfbwt.x) 
* `bigbwt -h` (get usage instruction)

Note that `bigbwt` is a Python script so you need python 3.x installed.
 

## Sample usage

The only requirement for the input file is that it does not contain the characters 0x00, 0x01, and 0x02. To build the BWT for file *yeast.fasta* just type

       bigbwt yeast.fasta

The program shows its progress and finally the running time. If no errors occurred the BWT file yeats.fasta.bwt is created: it should be one character longer than the input; the extra character is the BWT eos symbol represented by the ASCII character 0x00. Two important command line options are the windows size `-w` and the modulus `-m`. They are explained in the paper. 

If you don't trust the output run bigbwt with option `-c`. This will compute the same BWT using the SACAK algorithm [1] and compare it with the one computed by bigbwt. Be warned that SACAK, although being the most space economical among linear time algorithms, needs up to *9n* bytes of RAM, since it first compute the Suffix Array and then outputs the BWT (with extension .Bwt). 


## References

\[1\] Nong, G., Practical linear-time O(1)-workspace suffix sorting for constant alphabets, ACM Trans. Inform. Syst., vol. 31, no. 3, pp. 1-15, 2013


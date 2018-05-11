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

The only requirement for the input file is that it does not contain the characters 0x00, 0x01, and 0x02. To build the BWT for file *hugefile.fa* just type

       bigbwt hugefile.fa

The program shows its progress and finally the running time. If no errors occurred the BWT file hugefile.fa.bwt is created: it should be one character longer that the input and contains the character 0x00 in the position of the BWT's eos. Two important options are the windows size `-w` and the modulus `-m`. They are explained in the paper. 

If you don't trust the output run bigbwt with option `-c`. This will compute the same BWT using the GSACAK algorithm and compare it with the one computed by bigbwt. Be warned that GSACA BWT needs up to *9n* bytes, since it first compute the Suffix Array and the outputs the BWT.



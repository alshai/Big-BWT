modified by Taher Mun to support fasta files

# Big-BWT

Tool to build the BWT and optionally the Suffix Array for highly repetitive files using the approach described in *Prefix-Free Parsing for Building Big BWTs* by Christina Boucher, Travis Gagie, Alan Kuhnle and Giovanni Manzini [1].

Copyrights 2018- by the authors. 
 

## Installation

* Download/Clone the repository
* `make` (create the C/C++ executables) 
* `bigbwt -h` (get usage instruction)

Note that `bigbwt` is a Python script so you need at least Python 3.4 installed.
 

## Sample usage

The only requirement for the input file is that it does not contain the characters 0x00, 0x01, and 0x02 which are used internally by `bigbwt`. To build the BWT for file *yeast.fasta* just type

       bigbwt yeast.fasta

If no errors occurrs the BWT file yeast.fasta.bwt is created: it should be one character longer than the input; the extra character is the BWT eos symbol represented by the ASCII character 0x00. Two important command line options are the window size `-w` and the modulus `-m`. Their use is explained in the paper. 

Using option `-S` it is possible to compute the Suffix Array at the same time of the BWT. 
The Suffix Array is written using 5 bytes per integer (5 is a compilation constant defined in utils.h). This is exactly the same format used by the [pSAscan/SAscan](https://www.cs.helsinki.fi/group/pads/) tools that indeed should produce exactly the same Suffix Array file. 

Using options `-s` and/or `-e` (in alternative to `-S`) it is possible to compute a sampled Suffix Array. Option `-s` generates a `.ssa` file containing the SA entries at the beginning of BWT runs. More precisely, the `.ssa` file contains pairs of 5 bytes integers *<j,SA[j]>* for all indices *j* such that *BWT[j]!=BWT[j-1]*. Option `-e` generates a `.esa` file containing the SA entries at the end of BWT runs: i.e. the pairs *<j,SA[j]>* for all indices *j* such that *BWT[j]!=BWT[j+1]*.

The tool has limited (and experimental) support for multiple threads. 
Use `bigbwt` option `-t` to specify the number of helper threads: in our tests `-t 4` reduced the running time by roughly a factor two.

If you don't trust the output of our tool run bigbwt with option `-c`. 
This will compute the BWT using the SACAK algorithm [1] and compare it with the one computed by bigbwt. Be warned that SACAK, although being the most space economical among linear time algorithms, needs up to *9n* bytes of RAM, since it first compute the Suffix Array and then outputs the BWT (with extension .Bwt).


## References

\[1\]  Christina Boucher, Travis Gagie, Alan Kuhnle and Giovanni Manzini 
*Prefix-Free Parsing for Building Big BWTs* [Proc. WABI '18](http://drops.dagstuhl.de/opus/volltexte/2018/9304/).

\[2\] Nong, G., 
*Practical linear-time O(1)-workspace suffix sorting for constant alphabets*, ACM Trans. Inform. Syst., vol. 31, no. 3, pp. 1-15, 2013


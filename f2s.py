#!/usr/bin/env python3

import sys, time, argparse, subprocess, os.path, struct

Description = """
Tool to extract from the full suffix array the pairs (position,value) 
such that SA[position]=value and BWT(position)!=BWT(position-1)
These are the pairs corresponding to the first entry of each BWT-run. 

If option -e is used, the tool instead extracts the pairs corresponding 
to the last entry of each BWT-run ie the pairs (position,value)
such that SA[position]=value and BWT(position)!=BWT(position+1)

Note: the entry corresponding to the first BWT entry is considered the 
beginning of a BWT-run, and the last BWT entry is considered the end of
a BWT-run. 

The tool needs the files 
   infile.bwt (one byte per entry)
   infile.sa  (5 bytes per entry, changeable with option -B)
where the .sa file has one entry less than the .bwt file, and creates 
a file named infile.first/infile.last (changeable with option -o) contaning for each position
pos, corresponding to the first char of a BWT run, the pair (pos,sa[pos])
written using 5 bytes per entry (changeable with option -b).

The purpose of the tool is mainly to test the correctness of bigbwt with
option -s or -e that computes the same information using prefix free parsing
without building the full SA. 
"""

shasum_exe = "sha256sum"

def main():
  parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument('input', help='input file name', type=str)
  parser.add_argument('-b', help='number of bytes per entry (output, def. 5)', default=5, type=int)
  parser.add_argument('-B', help='number of bytes per entry (input, def. 5)', default=5, type=int)
  parser.add_argument('-e',  help='output pairs for the last entries of BWT-runs',action='store_true')
  parser.add_argument('-a',  help='output pairs in plain text to stderr',action='store_true')
  parser.add_argument('-o',  help='output file name (def. input.first/input.last)',default="", type=str)
  parser.add_argument('--sum', help='compute output file shasum',action='store_true')
  #parser.add_argument('-v',  help='verbose',action='store_true')
  args = parser.parse_args()
  start0 = start = time.time()
  if args.b<4 or args.B < 4:
    raise "Integers must be at least 4 bytes long"
  if args.b>8 or args.B > 8:
    raise "Integers can be at most 8 bytes long"

  
  # get BWT size
  size = os.path.getsize(args.input+".bwt")
  # get sa size
  sasize = os.path.getsize(args.input+".sa")
  if sasize!=args.B*(size-1): # the bwt has one entry more than the SA
    raise "SA file size mismatch: %d vs %d" % (sasize,args.B*(size-1))
  
  # define output file name
  if(len(args.o)>0):
    outname = args.o
  else:
    outname = args.input + (".last" if args.e else ".first")
  # open outfile
  with open(outname,"wb") as args.outfile:
    if not args.e:
      output_pair(0,size-1,args)          # bwt[0] = text[size-1] is the beginning of a run 
    with open(args.input + ".bwt","rb") as bwt:
      bwtlast = struct.unpack('B',bwt.read(1))[0] # this is bwt[0]
      salast = size-1 # this is sa[-1] sort of 
      with open(args.input + ".sa", "rb") as sa:
        for i in range(1,size):
          # read a pair bwt/sa
          buf = bwt.read(1)
          if len(buf)!=1: raise "Unexpected end of BWT file"
          bwtnext = struct.unpack('B',buf)[0] # this is bwt[i]
          buf = sa.read(args.B)
          if len(buf)!=args.B: raise "Unexpected end of SA file"
          buf += bytes(8-args.B)
          sanext = struct.unpack('<q',buf)[0]
          if sanext<0 or sanext>=size-1: raise "Invalid SA value (%d)" % sanext
          if bwtnext!=bwtlast:
            if args.e:
              output_pair(i-1,salast,args) # bwtlast is end of the current BWT run 
            else:  
              output_pair(i,sanext,args) # bwtnext is beginning of a new BWT run 
          bwtlast = bwtnext
          salast = sanext
        if args.e:
          output_pair(size-1,salast,args) # bwtlast is end of the current BWT run 
        if len(sa.read(1))>0:
          raise "SA file too long"
      if len(bwt.read(1))>0:
          raise "BWT file too long"
  # show digest if requested         
  if args.sum:
    digest = file_digest(outname)
    short = "ESA" if args.e else "SSA"            
    print("{short} {exe}: {digest}".format(short=short, exe=shasum_exe, digest=digest))    
    
  print("==== Done")



# output a pair of ints to args.outfile using args.b bytes for each
def output_pair(a,b,args):
  if args.a:
    print(a,b,file=sys.stderr)
  output_int(a,args.b,args.outfile)
  output_int(b,args.b,args.outfile)
  
  
# output a single int to file f using size bytes
def output_int(a,size,f):
  ax = a
  for i in range(size):
    b = a%256
    f.write(struct.pack('B',b)) #write unsigned int in 1 bytes
    a = a//256
  if a!=0:
    raise "Illegal integer (%d)" % ax  


# compute hash digest for a file 
def file_digest(name):
    try:
      hash_command = "{exe} {infile}".format(exe=shasum_exe, infile=name)
      hashsum = subprocess.check_output(hash_command.split())
      hashsum = hashsum.decode("utf-8").split()[0]
    except:
      hashsum = "Error!" 
    return hashsum  

if __name__ == '__main__':
    main()

#!/usr/bin/env python3

import sys, time, argparse, subprocess, os.path, struct

Description = """
Tool to extract from the full suffix array the pairs (position,value) 
such that SA[position]=value and BWT(position)!=BWT(position-1)
"""


def main():
  parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawTextHelpFormatter)
  parser.add_argument('input', help='input file name', type=str)
  parser.add_argument('-b', help='number of bytes per entry (output)', default=5, type=int)
  parser.add_argument('-B', help='number of bytes per entry (input)', default=5, type=int)
  parser.add_argument('-a',  help='output values to stderr',action='store_true')
  parser.add_argument('-o',  help='output file name',default="", type=str)
  parser.add_argument('--sum', help='compute output files shasum',action='store_true')
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
  if sasize!=args.B*(size-1):
    raise "SA file size mismatch: %d vs %d" % (sasize,args.B*(size-1))
  
  # define output file name
  if(len(args.o)>0):
    outname = args.o
  else:
    outname = args.input + ".first"
  # open outfile
  with open(outname,"wb") as args.outfile:
    output_pair(0,size-1,args)          # bwt[0] = text[size-1] is the beginning of a run 
    with open(args.input + ".bwt","rb") as bwt:
      bwtlast = struct.unpack('B',bwt.read(1))[0] # this is bwt[0]
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
            output_pair(i,sanext,args) # bwtnext is beginning of a new BWT run 
          bwtlast = bwtnext
        if len(sa.read(1))>0:
          raise "SA file too long"
      if len(bwt.read(1))>0:
          raise "BWT file too long"
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
def file_digest(name,logfile):
    try:
      hash_command = "{exe} {infile}".format(exe=shasum_exe, infile=name)
      hashsum = subprocess.check_output(hash_command.split(),stderr=logfile)
      hashsum = hashsum.decode("utf-8").split()[0]
    except:
      hashsum = "Error!" 
    return hashsum  

# execute command: return True is everything OK, False otherwise
def execute_command(command,logfile,logfile_name):
  try:
    subprocess.check_call(command.split(),stdout=logfile,stderr=logfile)
  except subprocess.CalledProcessError:
    print("Error executing command line:")
    print("\t"+ command)
    print("Check log file: " + logfile_name)
    return False
  return True



if __name__ == '__main__':
    main()

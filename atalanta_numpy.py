// A. Moshovos 2021

import numpy as np
import math 
import re
import sys
import os
from ctypes import *
from timeit import default_timer as timer




INPUT_RANGE = 256
INPUT_BITS = 8
TABLE_ENTRIES = 16
prob_n = TABLE_ENTRIES # used to trim entries if the high-end range of symbols end up  in less than 16 prob table 
PREC_BITS = 10   # AC probability precision bits
ROLL_BITS = 16   # AC range precision - 1
TABLE_OVERHEAD = TABLE_ENTRIES * (PREC_BITS + (INPUT_BITS-0).bit_length())
start = 0
end = 0

MAX_DEPTH = 2 # search depth ((sub)exponential cost)

class PTABLE_E(Structure):
     _fields_ = ("vmin", c_int), ("off", c_int), ("abits", c_int), ("obits", c_int), ("vcnt", c_int)
     # APACK PROB TABLE #############################################################
     # vmin = minimum value for range (ascending order starting from 0 and ending at 2^n (256 for 8b) 16+1 entries, last one is dummy with vmin = 2^n
     # off = range size is 2^off (offset bits)
     # abits = how many bits encoding the prefixes (total) (upper bound estimate)
     # obits = how many bits for the offsets (total)
     # vcnt = how many values fall into this range

PTABLE_TYPE = PTABLE_E * TABLE_ENTRIES

class apack():
	def __init__(self,verbose=0):
        	lib_path = './search.so'
        	#Clib = CDLL(os.path.abspath(lib_path))

        	try:
          	    Clib = CDLL(os.path.abspath(lib_path))
        	except:
          	    print('Clib %s not found' % (lib_path))

        	self.Csearch = Clib.search
        	self.Csearch.restype = None
        	self.Csearch.argtypes = [c_int, POINTER(c_int), POINTER(PTABLE_E), c_int]
        	self.verbose = verbose


	def search(self, values):
		global start, end
		self.histogram = [0]*2**INPUT_BITS
		start = timer()
		value_cnt = 0
		for v in np.nditer(values):
                    self.histogram[v] += 1
                    value_cnt += 1
		if self.verbose: print(self.histogram)
		#print("INPUT: ", value_cnt * INPUT_BITS, " bits")

		self.n = len(self.histogram)
		self.Chist = (c_int * self.n)(*self.histogram)
		self.CPTable = PTABLE_TYPE()
		self.Csearch(c_int(INPUT_BITS), self.Chist, self.CPTable, 1)
		end = timer()
		rpt = []
		apack_bits = 0
		for pt in self.CPTable:
			rpt.append([pt.vmin,pt.off,pt.abits,pt.obits,pt.vcnt])
			apack_bits += pt.abits + pt.obits
			print(pt.off,pt.vmin,pt.abits,pt.obits,pt.vcnt, pt.vcnt/value_cnt)
			if self.verbose == 1:
                            print(pt.off,pt.vmin,pt.abits,pt.obits,pt.vcnt)


		return rpt


def main():

    args = sys.argv
    thispy = args.pop(0)
    infile = args.pop(0)
    nbits = int(args.pop(0))
    global_verbose = 0

    print("INPUT FILE: ", infile)
    indata = np.load(infile, encoding='latin1', fix_imports=True)
    #print(indata)
    a = apack()
	
    pt = a.search(indata)
    print(pt)
    if True:
        value_cnt = 0 
        apack_bits = 0
        symbol_bits = 0
        offset_bits = 0
        for pte in pt:
            print (pte)
            apack_bits += pte[2] + pte[3]
            symbol_bits += pte[2]
            offset_bits += pte[3]
            value_cnt += pte[4]
            #print ("A ", pte[2], pte[3], pte[4])
        inbits = value_cnt*nbits
        print ("COMPSTATS: ", inbits, " (in_bits) ", apack_bits, " (apack_bits) ", apack_bits/inbits, " (comp_ratio) ", apack_bits/value_cnt, " (bits_per_value) ", symbol_bits, " (symbol bits) ", offset_bits, " (offset bits) ", end-start, " (search time) " )

if __name__=="__main__":
   main()

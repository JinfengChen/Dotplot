#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
from numpy import *
import re
import os
import argparse

def usage():
    test="name"
    message='''
python Mummer2dotplot.py --input mummer.mums > og_os4dotplot 

Get chromosome comprision by mummer
mummer -mum -l 200 -b -c ../input/O.glaberrima.v1.0.fa ../input/MSU_r7.fa > og_os.mummer.mums

Multi chromosome comparision
> Chr1
  chr01       576     37390       299
  chr01       892     37706       228
  chr01      1121     37935       257
  chr01      1779     38592       321


Single chromosome comparision
> Chr1
     576     37390       299
     892     37706       228
    1121     37935       257
> Chr1 Reverse
     576     37390       299
     892     37706       228
     1121     37935       257
    '''
    print message

def mummer2dotplot_single(mummer):
    s  = re.compile(r'^> (.*)$')
    s1 = re.compile(r'(\d+)')
    s2 = re.compile(r'Reverse')
    header = ''
    chro   = ''
    rev    = 0
    rank   = 0
    with open (mummer, 'r') as filefh:
        for line in filefh:
            line = line.rstrip()
            match = s.search(line)
            if match:
                header = match.groups(0)[0]
                chro   = s1.search(header).groups(0)[0]
                rev    = 1 if s2.search(header) else 0
                rank   = 0
                #print line
                print '## ' + header
                #print chro, rev
            else:
                line   = re.sub('\s+', ' ', line)
                line   = re.sub('^\s', '', line)
                unit   = line.rsplit(' ')
                chro   = '%02d' %int(chro)
                #print line
                #print chro, unit[0], 'Test'
                qryn   = 'Og_' + str(rank)
                qrys   = int(unit[0])
                qrye   = int(unit[0]) + int(unit[2])
                refn   = 'Os_' + str(rank)
                refs   = int(unit[1])
                refe   = int(unit[1]) + int(unit[2])
                output = []
                if (int(rev) == 1):
                    output = [chro, qryn, qrye, qrys, chro, refn, refs, refe, '1']
                else: 
                    output = [chro, qryn, qrys, qrye, chro, refn, refs, refe, '1']
                print '\t'.join(map(str,output))
                #print chro, qryn, qrys, qrye, chro, refn, refs, refe, '1'
                rank   += 1
                
def mummer2dotplot_multi(mummer, dotplot):
    s  = re.compile(r'^> (.*)$')
    s1 = re.compile(r'(\d+)')
    s2 = re.compile(r'Reverse')
    header = ''
    chro   = ''
    rev    = 0
    rank   = 0
    #dot = open (dotplot, 'w')
    with open (mummer, 'r') as filefh:
        for line in filefh:
            line = line.rstrip()
            match = s.search(line)
            if match:
                header = match.groups(0)[0]
                chro   = s1.search(header).groups(0)[0]
                rev    = 1 if s2.search(header) else 0
                rank   = 0
                #print line
                print '## ' + header
                #print chro, rev
            else:
                line   = re.sub('\s+', ' ', line)
                line   = re.sub('^\s', '', line)
                unit   = line.rsplit(' ')
                chro   = '%02d' %int(chro)
                chro1  = '%02d' %int(s1.search(unit[0]).groups(0)[0])
                #print line
                #print chro, unit[0], 'Test'
                qryn   = 'Og_' + str(rank)
                qrys   = int(unit[1])
                qrye   = int(unit[1]) + int(unit[3])
                refn   = 'Os_' + str(rank)
                refs   = int(unit[2])
                refe   = int(unit[2]) + int(unit[3]) 
                if (int(rev) == 1):
                    output = [chro1, qryn, qrye, qrys, chro, refn, refs, refe, '1']
                else:
                    output = [chro1, qryn, qrys, qrye, chro, refn, refs, refe, '1']
                print '\t'.join(map(str,output))
                #print chro, qryn, qrys, qrye, chro, refn, refs, refe, '1'
                rank   += 1



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-s', '--single', default=0)
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)
    if (args.single == 1):
        print 'Single comparision'
        mummer2dotplot_single(args.input, args.output)
    else:
        print 'Multi comparision'
        mummer2dotplot_multi(args.input, args.output)
    

if __name__ == '__main__':
    main()


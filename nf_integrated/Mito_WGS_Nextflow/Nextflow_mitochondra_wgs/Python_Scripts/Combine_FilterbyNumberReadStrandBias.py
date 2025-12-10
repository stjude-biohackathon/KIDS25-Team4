import os
import argparse

def Combine_Org_Shift_Mit_Variant(Org_Tab, Shifted_Tab, Combine_Out_Tab):
   NewFile  = open(Combine_Out_Tab, 'w')
   OrigFile = open(Org_Tab, 'r')
   headers  = OrigFile.readline().strip()
   headCols = headers.split('\t')
   headCols = headCols[2:]
   newheader = ''
   for col in headCols:
       newheader += col + '\t'
   newheader = newheader[:-1] + '\n'
   NewFile.write(newheader)
   for line in OrigFile:
       line = line.rstrip()
       cols = line.split('\t')
    #pos = int(cols[2])
       pos = int(cols[4])
       if pos >=4000 and pos < 12000:
           newline = '' 
           for i in range(2, len(cols)):
               newline += cols[i] + '\t'
           newline = newline[:-1] + '\n'
           NewFile.write(newline)
   OrigFile.close()

   ShiftFile = open(Shifted_Tab, 'r')
   headersShifted = ShiftFile.readline().strip()
   headColsShifted = headersShifted.split('\t')
   indHeaderMapping = {}
   for ind in range(0, len(headCols)):
       ind2 = headColsShifted.index(headCols[ind])
       indHeaderMapping[ind] = ind2

   for line in ShiftFile:
       line = line.rstrip()
       cols = line.split('\t')
       pos = int(cols[4])
       if len(cols) < 44:
           cols += ['']
       if (pos >=4000 and pos <=12569):
           origPos = pos + 8000
           if origPos > 16569:
               origPos-=16569
           newline='chrM.%d\tchrM\t%d\t'%(origPos, origPos)
           for ind in range(3, len(headCols)):
               ind2 = indHeaderMapping[ind]
               newline += cols[ind2] + '\t'
           newline=newline[:-1] + '\n'
           NewFile.write(newline)
   ShiftFile.close()
   NewFile.close()
   
def Filter_Combine_Tab(Combine_Tab,  Combine_Filter_Out_Tab):
   CombineFile = open(Combine_Tab, 'r')
   header = CombineFile.readline()
   headers = header.split('\t')
   refNormInd = headers.index('reference_normal_count')
   refTumorInd = headers.index('reference_tumor_count')
   altNormInd = headers.index('alternative_normal_count')
   altTumorInd = headers.index('alternative_tumor_count')
   altFwdInd = headers.index('alternative_fwd_count')
   altRevInd = headers.index('alternative_rev_count')
   NewFile = open(Combine_Filter_Out_Tab, 'w')
   NewFile.write(header)
   totalLines = 0
   totalNewLines = 0
   for line in CombineFile:
       totalLines += 1
       cols = line.split('\t')
       ref_norm = int(cols[refNormInd])
       alt_norm = int(cols[altNormInd])
       ref_tumor = int(cols[refTumorInd])
       alt_tumor = int(cols[altTumorInd])
       alt_fwd = int(cols[altFwdInd])
       alt_rev = int(cols[altRevInd])
       if alt_fwd + alt_rev == 0:
           perc_fwd = .5
       else:
           perc_fwd = float(alt_fwd) / (alt_fwd + alt_rev)
       strand_bias = False
       if perc_fwd > .9 or perc_fwd <.1:
           strand_bias = True
     #print cols[0], strand_bias, alt_fwd, alt_rev, perc_fwd, ref_norm+alt_norm, ref_tumor+alt_tumor

     #if (ref_norm + alt_norm) > 100 and (ref_tumor + alt_tumor) > 100 and not strand_bias:
   
     #  if cols[0] == 'chrM.8612':
     #     print(Combine_Tab, strand_bias, perc_fwd, cols)
       if (ref_norm + alt_norm) > 100 and (ref_tumor + alt_tumor) > 100 and not strand_bias and (alt_tumor >= 10 or alt_norm >= 10):
           totalNewLines += 1
           NewFile.write(line)
   print (totalLines, totalNewLines)
   
   CombineFile.close()
   NewFile.close()

def main():
   parser = argparse.ArgumentParser(description="")
   parser.add_argument('-org','--org_tab', help='A sample Mito Variant called against MT.fa ', required=True, type=str)
   parser.add_argument('-shift','--shift_tab', help='A sample Mito Variant called against MT_Shift.fa ', required=True, type=str)
   parser.add_argument('-o','--out_tab', help='A sample Combined and filtered  Mito Variant', required=True, type=str)

   args=parser.parse_args()
   Combine_Out_tmp = "Combine_Out_tmp.tab"
   Combine_Org_Shift_Mit_Variant(args.org_tab, args.shift_tab, Combine_Out_tmp)
   Filter_Combine_Tab(Combine_Out_tmp, args.out_tab)

if __name__ == "__main__":
   main()

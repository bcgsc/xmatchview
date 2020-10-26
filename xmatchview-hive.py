#!/usr/bin/env python3
# xmatchview-hive.py
# Visualizing 3-way genome synteny with a hive plot representation
# Rene L Warren 2005-2020

import sys
import os
import getopt
import re
import csv
import subprocess
import math

#---------------------------------------------
def checkFile(file):

   print("Checking input %s" % file)
   if not os.path.exists(file):
      print("File %s" % file + " is not valid")
      sys.exit(1)
   else:
      print("exists.")

#---------------------------------------------
def readGFF(file,scale):

   feature = {}
   (start,end) = (0,0)

   with open(file) as fd:
      rd = csv.reader(fd, delimiter="\t", quotechar='"')
      for row in rd:
         id = row[0]
         start = float(int(row[3])/scale)
         end = float(int(row[4])/scale)
         #print "__%s__ - __%s__ <<<<<" % (start,end)

         color = row[9]
         if color == None:
           color = "black"

         if id not in feature:
            feature[id] = {}
         if start not in feature[id]:
            feature[id][start] = {}
         if 'end' not in feature[id][start]:
            feature[id][start]['end'] = ""
         if 'color' not in feature[id][start]:
            feature[id][start]['color'] = ""

         feature[id][start]['end'] = end
         feature[id][start]['color'] = color
         #print "INITIALIZED %s : %i,%i with %s" %(row[0],start,end,color)

   return feature

#---------------------------------------------
def readPAF(paf_file,mismatch,block_length,query,reference,scale):

   match = {}
   ctline = 0

   xmatch_obj=open(paf_file, 'r')

   for line in xmatch_obj:
      ### reverse matches      qryname    qrystart qryend orient hitname       hitstart hitend match   block
      ###                       1             2       3            4             5       6       7       8
      rev_regex = re.compile("(\S+)\s+\S+\s+(\d+)\s+(\d+)\s+\-\s+(\S+)\s+\S+\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)")
      rm = rev_regex.match(line)
      ###forward matches

      fwd_regex = re.compile("(\S+)\s+\S+\s+(\d+)\s+(\d+)\s+\+\s+(\S+)\s+\S+\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)")
      fm = fwd_regex.match(line)

      if rm != None:
         ctline = 1
         #print "GR: %s" % line
         #print "REVERSE: %s %s %s %s %s %s %s %s" % (rm.group(1), rm.group(2), rm.group(3), rm.group(4), rm.group(5), rm.group(6), rm.group(7), rm.group(8))

         alignLen = float(rm.group(6)) - float(rm.group(5)) + 1
         percentMis = 100 * float(( alignLen - float(rm.group(7))) / alignLen ) 
         #print "=== %.2f === %.2f ===" % (alignLen,percentMis)
         #sys.exit(1)
         (primary_match, startFirstMatch, endFirstMatch, secondary_match, startSecondMatch, endSecondMatch)=(str(rm.group(1)), float(rm.group(2)), float(rm.group(3)), str(rm.group(4)), float(rm.group(6)), float(rm.group(5))) ### needs 5/6 reversed to plot reverse align

         if percentMis > mismatch:
            continue  #will not display alignment lines below threshold
         elif (endFirstMatch - startFirstMatch) < block_length:
            continue  #will skip smaller alignment
         else:
            if primary_match in query and secondary_match in reference:
               startFirstMatch = (startFirstMatch/scale) + query[primary_match]['offset_len']         #RLW
               endFirstMatch = (endFirstMatch/scale) + query[primary_match]['offset_len']             #RLW
               startSecondMatch = (startSecondMatch/scale) + reference[secondary_match]['offset_len'] #RLW
               endSecondMatch = (endSecondMatch/scale) + reference[secondary_match]['offset_len']     #RLW

               #print "%i-%i   ::   %i-%i" % (startFirstMatch,endFirstMatch,startSecondMatch,endSecondMatch)

               if primary_match not in match:
                  match[primary_match]={}
               if secondary_match not in match[primary_match]:
                  match[primary_match][secondary_match]={}
               if startFirstMatch not in match[primary_match][secondary_match]:
                  match[primary_match][secondary_match][startFirstMatch]={}
               if endFirstMatch not in match[primary_match][secondary_match][startFirstMatch]:
                  match[primary_match][secondary_match][startFirstMatch][endFirstMatch]={}
               if startSecondMatch not in match[primary_match][secondary_match][startFirstMatch][endFirstMatch]:
                  match[primary_match][secondary_match][startFirstMatch][endFirstMatch][startSecondMatch]={}
               if endSecondMatch not in match[primary_match][secondary_match][startFirstMatch][endFirstMatch][startSecondMatch]:
                  match[primary_match][secondary_match][startFirstMatch][endFirstMatch][startSecondMatch][endSecondMatch]={}

               match[primary_match][secondary_match][startFirstMatch][endFirstMatch][startSecondMatch][endSecondMatch]=percentMis

      ###forward matches
      elif fm != None:
         ctline = 1
         #print "GF: %s" % line
         #print "FORWARD: %s %s %s %s %s %s %s %s" % (fm.group(1), fm.group(2), fm.group(3), fm.group(4), fm.group(5), fm.group(6), fm.group(7), fm.group(8))

         alignLen = float(fm.group(6)) - float(fm.group(5)) + 1
         percentMis = 100 * float(( alignLen - float(fm.group(7))) / alignLen )
         #percentMis2 = 100 * float((float(fm.group(8)) - float(fm.group(7))) / float(fm.group(8)))
         #print "=== %.2f === %.2f === %.2f" % (alignLen,percentMis,percentMis2)
         #sys.exit(1)
         (primary_match, startFirstMatch, endFirstMatch, secondary_match, startSecondMatch, endSecondMatch)=(str(fm.group(1)), float(fm.group(2)), float(fm.group(3)), str(fm.group(4)), float(fm.group(5)), float(fm.group(6)))

         if percentMis > mismatch:
            continue  #will not display alignment lines below threshold
         elif (endFirstMatch - startFirstMatch) < block_length:
            continue  #will skip smaller alignment
         else:
            if primary_match in query and secondary_match in reference:
               startFirstMatch = (startFirstMatch/scale) + query[primary_match]['offset_len']         #RLW
               endFirstMatch = (endFirstMatch/scale) + query[primary_match]['offset_len']             #RLW
               startSecondMatch = (startSecondMatch/scale) + reference[secondary_match]['offset_len'] #RLW
               endSecondMatch = (endSecondMatch/scale) + reference[secondary_match]['offset_len']     #RLW

               #print "%i-%i   ::   %i-%i" % (startFirstMatch,endFirstMatch,startSecondMatch,endSecondMatch)

               if primary_match not in match:
                  match[primary_match]={}
               if secondary_match not in match[primary_match]:
                  match[primary_match][secondary_match]={}
               if startFirstMatch not in match[primary_match][secondary_match]:
                  match[primary_match][secondary_match][startFirstMatch]={}
               if endFirstMatch not in match[primary_match][secondary_match][startFirstMatch]:
                  match[primary_match][secondary_match][startFirstMatch][endFirstMatch]={}
               if startSecondMatch not in match[primary_match][secondary_match][startFirstMatch][endFirstMatch]:
                  match[primary_match][secondary_match][startFirstMatch][endFirstMatch][startSecondMatch]={}
               if endSecondMatch not in match[primary_match][secondary_match][startFirstMatch][endFirstMatch][startSecondMatch]:
                  match[primary_match][secondary_match][startFirstMatch][endFirstMatch][startSecondMatch][endSecondMatch]={}

               match[primary_match][secondary_match][startFirstMatch][endFirstMatch][startSecondMatch][endSecondMatch]=percentMis

         #else:
            #print "NO RE:%s" % line
   xmatch_obj.close()

   if ctline == 0 :
      print("There are no alignments to plot. Make sure your file %s is reporting alignments -- fatal." % paf_file)
      sys.exit(1)

   return match

#---------------------------------------------
def readCrossMatch(crossmatch_file,mismatch,block_length,query,reference,scale):

   match = {}
   ctline = 0

   xmatch_obj=open(crossmatch_file, 'r')

   for line in xmatch_obj:

      #                     Query                       start end         R Ref                                    end   start
      # 10  8.70 0.00 0.00  JN039333.1_Picea_abies      440   462 (2149)  C KT263970.1_Picea_sitchensis   (1851)   524   502
      #                                                                                                      start end
      # 16  4.55 0.00 0.00  JN039333.1_Picea_abies      484   505 (2106)    KT263970.1_Picea_sitchensis      341   362 (2013)

      ###reverse matches                   s.i.                qry                       
      rev_regex = re.compile("(\s+)?\d+\s+(\S+)\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+(\S+)\s+\S+\s+C\s+(\S+)\s+\S+\s+(\S+)\s+(\S+)")
      rm = rev_regex.match(line)
      ###forward matches
      fwd_regex = re.compile("(\s+)?\d+\s+(\S+)\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+(\S+)\s+\S+\s+(\S+)\s+(\S+)\s+(\S+)\s+\S+")
      fm = fwd_regex.match(line)

      if rm != None and rm.group(3) != "0" and rm.group(6) != "0":
         ctline = 1
         #print "GR: %s" % line
         #print "REVERSE: %s %s %s %s %s %s %s" % (rm.group(1), rm.group(2), rm.group(3), rm.group(4), rm.group(5), rm.group(6), rm.group(7))

         #(percentMis, primary_match, startFirstMatch, endFirstMatch, secondary_match, startSecondMatch, endSecondMatch)=(float(rm.group(1)), rm.group(2), float(rm.group(3)), float(rm.group(4)), rm.group(5), float(rm.group(6)), float(rm.group(7)))
         (percentMis, primary_match, startFirstMatch, endFirstMatch, secondary_match, startSecondMatch, endSecondMatch)=(float(rm.group(2)), rm.group(3), float(rm.group(4)), float(rm.group(5)), rm.group(6), float(rm.group(7)), float(rm.group(8)))### has to be in this order to plot reverse align in diff color

         if percentMis > mismatch:
            continue  #will not display alignment lines below threshold
         elif (endFirstMatch - startFirstMatch) < block_length:
            continue  #will skip smaller alignment
         else:
            if primary_match in query and secondary_match in reference:
               startFirstMatch = (startFirstMatch/scale) + query[primary_match]['offset_len']         #RLW
               endFirstMatch = (endFirstMatch/scale) + query[primary_match]['offset_len']             #RLW
               startSecondMatch = (startSecondMatch/scale) + reference[secondary_match]['offset_len'] #RLW
               endSecondMatch = (endSecondMatch/scale) + reference[secondary_match]['offset_len']     #RLW

               #print "REVERSE %i-%i   ::   %i-%i" % (startFirstMatch,endFirstMatch,startSecondMatch,endSecondMatch)

               if primary_match not in match:
                  match[primary_match]={}
               if secondary_match not in match[primary_match]:
                  match[primary_match][secondary_match]={}
               if startFirstMatch not in match[primary_match][secondary_match]:
                  match[primary_match][secondary_match][startFirstMatch]={}
               if endFirstMatch not in match[primary_match][secondary_match][startFirstMatch]:
                  match[primary_match][secondary_match][startFirstMatch][endFirstMatch]={}
               if startSecondMatch not in match[primary_match][secondary_match][startFirstMatch][endFirstMatch]:
                  match[primary_match][secondary_match][startFirstMatch][endFirstMatch][startSecondMatch]={}
               if endSecondMatch not in match[primary_match][secondary_match][startFirstMatch][endFirstMatch][startSecondMatch]:
                  match[primary_match][secondary_match][startFirstMatch][endFirstMatch][startSecondMatch][endSecondMatch]={}

               match[primary_match][secondary_match][startFirstMatch][endFirstMatch][startSecondMatch][endSecondMatch]=percentMis


      ###forward matches
      elif fm != None and fm.group(3) != "0" and fm.group(6) != "0":
         ctline = 1
         #print "GF: %s" % line
         #print "FORWARD: %s %s %s %s %s %s %s" % (fm.group(1), fm.group(2), fm.group(3), fm.group(4), fm.group(5), fm.group(6), fm.group(7))
#         (percentMis, primary_match, startFirstMatch, endFirstMatch, secondary_match, startSecondMatch, endSecondMatch)=(float(fm.group(1)), fm.group(2), float(fm.group(3)), float(fm.group(4)), fm.group(5), float(fm.group(6)), float(fm.group(7)))

         (percentMis, primary_match, startFirstMatch, endFirstMatch, secondary_match, startSecondMatch, endSecondMatch)=(float(fm.group(2)), fm.group(3), float(fm.group(4)), float(fm.group(5)), fm.group(6), float(fm.group(7)), float(fm.group(8)))

         if percentMis > mismatch:
            continue  #will not display alignment lines below threshold
         elif (endFirstMatch - startFirstMatch) < block_length:
            continue  #will skip smaller alignment
         else:
            if primary_match in query and secondary_match in reference:
               startFirstMatch = (startFirstMatch/scale) + query[primary_match]['offset_len']         #RLW
               endFirstMatch = (endFirstMatch/scale) + query[primary_match]['offset_len']             #RLW
               startSecondMatch = (startSecondMatch/scale) + reference[secondary_match]['offset_len'] #RLW
               endSecondMatch = (endSecondMatch/scale) + reference[secondary_match]['offset_len']     #RLW

               #print "FORWARD %i-%i  ::   %i-%i" % (startFirstMatch,endFirstMatch,startSecondMatch,endSecondMatch)
 
               if primary_match not in match:
                  match[primary_match]={}
               if secondary_match not in match[primary_match]:
                  match[primary_match][secondary_match]={}
               if startFirstMatch not in match[primary_match][secondary_match]:
                  match[primary_match][secondary_match][startFirstMatch]={}
               if endFirstMatch not in match[primary_match][secondary_match][startFirstMatch]:
                  match[primary_match][secondary_match][startFirstMatch][endFirstMatch]={}
               if startSecondMatch not in match[primary_match][secondary_match][startFirstMatch][endFirstMatch]:
                  match[primary_match][secondary_match][startFirstMatch][endFirstMatch][startSecondMatch]={}
               if endSecondMatch not in match[primary_match][secondary_match][startFirstMatch][endFirstMatch][startSecondMatch]:
                  match[primary_match][secondary_match][startFirstMatch][endFirstMatch][startSecondMatch][endSecondMatch]={}

               match[primary_match][secondary_match][startFirstMatch][endFirstMatch][startSecondMatch][endSecondMatch]=percentMis

         #else:
            #print "NO RE:%s" % line      
   xmatch_obj.close()

   if ctline == 0 :
      print("There are no alignments to plot. Make sure your file %s is reporting alignments -- fatal." % crossmatch_file)
      sys.exit(1)

   return match


#---------------------------------------------
def findOccurences(s, ch):
    return [i for i, letter in enumerate(s) if letter == ch]


#---------------------------------------------
def readText(file, scale):
    
   L1={}
   order=[]#RLW
   npos={} #RLW
   tot_length = 0 #RLW
   buffer = scale * 3
   print("Reading %s..\n" % file)
   with open(file) as fd:
      rd = csv.reader(fd, delimiter=":", quotechar='"')
      for row in rd:
         if row[0] not in L1:#RLW
                L1[row[0]]={}          #RLW
         L1[row[0]]['scaled_len'] = float(int(row[1])/scale)
         L1[row[0]]['offset_len'] = float(tot_length/scale) #RLW first ID will be offset at 0 as it should
         tot_length += (int(row[1]) + buffer)
         order.append(row[0])

   scaled_tot_len = float(tot_length/scale) #RLW
   print("Total length = %i bp, scaled to : %.0f pixels " % (tot_length,scaled_tot_len)) #RLW

   return (order, L1, tot_length)


#---------------------------------------------
def readFasta(file, scale):

   (head_match, previous_contig,seq_length) = (None,None,0)
   L1={}
   order=[]#RLW
   npos={} #RLW
   tot_length = 0 #RLW
   buffer = scale * 3
   print("Reading %s..\n" % file)

   file_obj = open(file, 'r')
   
   for line in file_obj:
      head_match_regex = re.compile('>(\S+)')
      head_match = head_match_regex.match(line)
      if head_match != None:
          if (head_match != previous_contig and previous_contig != None):
             (seq_length, scale)=(int(seq_length), int(scale))
             if previous_contig not in L1:#RLW
                L1[previous_contig]={}          #RLW
             L1[previous_contig]['scaled_len'] = float(seq_length/scale) #RLW
             L1[previous_contig]['offset_len'] = float(tot_length/scale) #RLW first ID will be offset at 0 as it should
             L1[previous_contig]['npos'] = npos
             #print "NPOS: (Ns only tracked on 1-line sequences)"
             #print(npos)

             order.append(previous_contig) #RLW
             #print "%s length = %i bp, scaled to %.0f pixels" % (previous_contig,seq_length,L1[previous_contig]['scaled_len']) #RLW
             tot_length += (seq_length + buffer)
             seq_length = 0                                        #resets the sequence length
          previous_contig = head_match.group(1)
      else:     
         seq_subset_regex = re.compile('(.*)', re.I)
         seq_subset = seq_subset_regex.match(line)
         if seq_subset != None:
            seq_length += len(seq_subset.group(1))
            npos=findOccurences(seq_subset.group(1).upper(), "N")

   (seq_length, scale)=(int(seq_length), int(scale))
   if previous_contig not in L1:  #RLW
      L1[previous_contig]={}            #RLW
   L1[previous_contig]['scaled_len'] = float(seq_length/scale) #RLW
   L1[previous_contig]['offset_len'] = float(tot_length/scale) #RLW
   L1[previous_contig]['npos'] = npos
   #print "NPOS: (Ns only tracked on 1-line sequences)"
   #print(npos)

   order.append(previous_contig) #RLW
   #print "%s length = %i bp, scaled to %.0f pixels" % (previous_contig,seq_length,L1[previous_contig]['scaled_len']) #RLW
   tot_length += (seq_length + buffer) #RLW

   file_obj.close()

   scaled_tot_len = float(tot_length/scale) #RLW
   print("Total length = %i bp, scaled to : %.0f pixels " % (tot_length,scaled_tot_len)) #RLW

   return (order, L1, tot_length)


#---------------------------------------------
def initColor(alpha):
   color={}

   #allocate colors
   color["white"] = (255,255,255,255)
   #color["black"] = (255,255,255,255)#(0,0,0,255)
   color["black"] = (0,0,0,255)
   color["swamp"] = (150,150,30,255)
   color["blue"] = (0,102,204,255)
   color["yellow"] = (255,255,0,255)
   color["cyan"] = (0,255,255,255)
   color["purple"] = (255,0,255,255)
   color["lime"] = (57,255,20,255)   ### XXX
   color["green"] = (100,250,25,255)
   color["red"] = (250,25,75,255)
   color["forest"] = (0,100,0,255)
   color["dirtyred"] = (200,0,120,255)
   color["navy"] = (0,0,150,255)
   color["dirtyyellow"] = (200,200,75,255)
   color["grey"] = (153,153,153,255)
   color["lightgrey"] = (220,220,220,255)
   color["salmon"] = (255,153,153,255)
   color["lightblue"] = (153,204,255,255)
   color["orange"] = (255,153,51,255)
   color["forestt"] = (0,100,0,alpha)

   color["green1t"] = (223,238,218,alpha)
   color["green2t"] = (208,221,203,alpha)
   color["green3t"] = (184,212,178,alpha)
   color["green4t"] = (162,203,155,alpha)
   color["green5t"] = (141,186,127,alpha)
   color["green6t"] = (119,176,108,alpha)
   color["green7t"] = (98,166,92,alpha)
   color["green8t"] = (72,146,73,alpha)
   color["green9t"] = (21,119,40,alpha)
   color["green10t"] = (0,82,33,alpha)

   color["red1t"] = (252,227,229,alpha)
   color["red2t"] = (249,214,215,alpha)
   color["red3t"] = (244,187,188,alpha)
   color["red4t"] = (240,161,161,alpha)
   color["red5t"] = (235,134,134,alpha)
   color["red6t"] = (231,107,108,alpha)
   color["red7t"] = (226,81,81,alpha)
   color["red8t"] = (222,54,54,alpha)
   color["red9t"] = (217,27,27,alpha)
   color["red10t"] = (213,1,1,alpha)

   color["green1"] = (247,252,245,255)
   color["green2"] = (229,245,224,255)
   color["green3"] = (199,233,192,255)
   color["green4"] = (161,217,155,255)
   color["green5"] = (116,196,118,255)
   color["green6"] = (65,171,93,255)
   color["green7"] = (35,139,69,255)
   color["green8"] = (0,109,44,255)
   color["green9"] = (0,68,27,255)
   color["brown"] = (83,49,24,255)
   color["brownt"] = (83,49,24,alpha)
   color["beige"] = (210,180,140,255)
   color["beiget"] = (210,180,140,alpha)
   return color

#---------------------------------------------
def initGraph():
   data={} 
  
   #default data points
   data['width']=4000 #5000
   data['height']=4000 #5000
   data['mid']=data['height']/2
   data['midwidth']=data['width']/2

   data['ref_y'] = (data['height'] / 1.5) 
   data['skew']=300 #400 ### CHANGE THIS FOR THE PITCH OF THE TREE
   data['decay']=120 #200### THIS IS THE SPACE BETWEEN BOTH SIDES, TOP OF TREE
   data['ref_y_skew']=data['ref_y']-data['skew'] ###DON'T CHANGE THIS
   data['mis_bar']=50
   data['query_y']=70
   data['x']=100
   data['xlabel']=110
   data['bar_thick']=20
   data['query_thick']=15
   data['reference_thick']=15
   data['x_legend'] = data['width'] - 600 
   data['y_legend'] = data['height'] - 100 ###WAS 1500XXX
   data['x_legend_picto'] = data['width'] / 1.5
   data['tick_up']=data['ref_y_skew'] - 120
   data['tick_down']=data['tick_up'] + 20

   return data


#---------------------------------------------
def drawRelationship(info,gff1,gff2,gff3,order1,order2,order3,ref1,ref2,ref3,length1,length2,length3,match1,match2,match3,scale,seqidentity,mismatch,block_length,alpha):

      ###Initialize new graph
      data=initGraph()

      ###var
      xlegend = 0
      ylegend = 0
      output = "xmv-hive_i" + str(seqidentity) + "_b" + str(block_length) + "_c" + str(scale) + "_a" + str(alpha) + ".svg" 

      ### write to svg
      # xml / SVG
      xml = open(output,"w+")
      xml.write("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n")
      xml.write("<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.0//EN\" \"http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd\">\n")
      xml.write("<svg width=\"%d\" height=\"%d\" xmlns=\"http://www.w3.org/2000/svg\" style=\"background-color:white\">\n" % (data['width'], data['height']))

      ### draw axes
      # axis 1
      y1 = data['mid'] - 50
      x1 = data['midwidth']
      xt = 0
      yt = 0

      for ref in order1:      #RLW
         scaled_offset_len = 0   #RLW
         if ref in ref1: #RLW
            scaled_offset_len = ref1[ref]['offset_len'] #RLW
         y1 = data['mid'] - 50 - scaled_offset_len
         y2 = y1 - ref1[ref]['scaled_len']
         x2 = x1
         xt = x2 + 5
         yt = y2 - 5

         xml.write("<path d=\"M %f %f L %f %f\" stroke=\"black\" stroke-width=\"6\" fill=\"transparent\"/>\n" % (x1,y1,x2,y2))

         ####draw features/exons on side of ref
         ###REF gene model
         #print " 1 %s %i %i" % (ref,ref1[ref]['scaled_len'],ref1[ref]['offset_len'])
         if ref in gff1:
            print(" 1 %s" % ref)
            for exstart in gff1[ref]:
                exend = gff1[ref][exstart]['end']
                gy1 = y1 - (exstart)
                gy2 = y1 - (exend)
                xml.write("<path d=\"M %f %f L %f %f\" stroke=\"%s\" stroke-width=\"3\" fill=\"transparent\"/>\n" % (x1,gy1,x2,gy2,gff1[ref][exstart]['color']))
      ylegend = yt
      x1 = x2
      xml.write("<text font-size=\"2.5em\" x=\"%f\" y=\"%f\">%s</text>\n" % (xt,yt,info['1']))
      ###########################################
      # axis 2
      x1 = data['midwidth'] + 50
      y1 = data['mid'] + 50
      xtick1 = x1 - 6 
      xtick2 = x1 + 6
      ytick1 = y1 + 10
      ytick2 = y1 - 10
      xt = 0
      yt = 0

      for ref in order2:      #RLW
         scaled_offset_len = 0   #RLW
         if ref in ref2: #RLW
            scaled_offset_len = ref2[ref]['offset_len'] #RLW
         #print "%s %f %f<<<\n" % (ref,ref2[ref]['scaled_len'],scaled_offset_len)
         x1 = data['midwidth'] + 50 + (math.sin(45) * scaled_offset_len)
         y1 = data['mid'] + 50 + (math.cos(45) * scaled_offset_len)
         x2 = (math.sin(45) * ref2[ref]['scaled_len']) + x1
         y2 = (math.cos(45) * ref2[ref]['scaled_len']) + y1
         xt = x2 + 5
         yt = y2 + 5

         xml.write("<path d=\"M %f %f L %f %f\" stroke=\"black\" stroke-width=\"6\" fill=\"transparent\"/>\n" % (x1,y1,x2,y2))

         ####draw features/exons on side of ref
         ###REF gene model
         #print "%s" % ref
         if ref in gff2:
            #print "%s" % ref
            for exstart in gff2[ref]:
                exend = gff2[ref][exstart]['end']
                gx1 = (math.sin(45) * (exstart)) + x1
                gy1 = (math.cos(45) * (exstart)) + y1
                gx2 = (math.sin(45) * (exend)) + x1
                gy2 = (math.cos(45) * (exend)) + y1
                l2 = y1 - (exend + scaled_offset_len)
                xml.write("<path d=\"M %f %f L %f %f\" stroke=\"%s\" stroke-width=\"3\" fill=\"transparent\"/>\n" % (gx1,gy1,gx2,gy2,gff2[ref][exstart]['color']))
      xlegend = xt
      xml.write("<text font-size=\"2.5em\" x=\"%f\" y=\"%f\">%s</text>\n" % (xt,yt,info['2']))
      ###########################################
      # axis 3
      x1 = data['midwidth'] - 50
      y1 = data['mid'] + 50
      xt = 0
      yt = 0

      for ref in order3:      #RLW
         scaled_offset_len = 0   #RLW
         if ref in ref3: #RLW
            scaled_offset_len = ref3[ref]['offset_len'] #RLW
         x1 = data['midwidth'] - 50 + (math.sin(-45) * scaled_offset_len)
         y1 = data['mid'] + 50 + (math.cos(-45) * scaled_offset_len)
         x2 = (math.sin(-45) * ref3[ref]['scaled_len']) + x1
         y2 = (math.cos(-45) * ref3[ref]['scaled_len']) + y1

         ###axis 3 label coordinate
         xt = x2 - 19 * len(info['3'])
         yt = y2 + 5

         ###axis
         xml.write("<path d=\"M %f %f L %f %f\" stroke=\"black\" stroke-width=\"6\" fill=\"transparent\"/>\n" % (x1,y1,x2,y2))

         ####draw features/exons on side of ref
         ###REF gene model
         #print "%s" % ref
         if ref in gff3:
            #print "%s" % ref
            for exstart in gff3[ref]:
                exend = gff3[ref][exstart]['end']
                gx1 = (math.sin(-45) * (exstart)) + x1
                gy1 = (math.cos(-45) * (exstart)) + y1
                gx2 = (math.sin(-45) * (exend)) + x1
                gy2 = (math.cos(-45) * (exend)) + y1
                l2 = y1 - (exend + scaled_offset_len)
                xml.write("<path d=\"M %f %f L %f %f\" stroke=\"%s\" stroke-width=\"3\" fill=\"transparent\"/>\n" % (gx1,gy1,gx2,gy2,gff3[ref][exstart]['color']))
      ### axis 3 label
      xml.write("<text font-size=\"2.5em\" x=\"%f\" y=\"%f\">%s</text>\n" % (xt,yt,info['3']))
      ####### END AXES
      # Begin ALIGNMENTS

      xq1 = data['midwidth'] + 5
      xq2 = xq1

      ###DRAW ALIGN 1vs2
      for qry in match1:
         allhit = match1[qry]
         for hit in allhit:
            start1_list = allhit[hit]
            s1_list_sort=list(start1_list.keys())
            s1_list_sort.sort()
            for start1 in s1_list_sort:
               end1_list=start1_list[start1]
               e1_list_sort=list(end1_list.keys())
               e1_list_sort.sort()
               for end1 in e1_list_sort:
                  start2_list=end1_list[end1]
                  s2_list_sort=list(start2_list.keys())
                  s2_list_sort.sort()
                  for start2 in s2_list_sort:
                     end2_list=start2_list[start2]
                     e2_list_sort=list(end2_list.keys())
                     e2_list_sort.sort()
                     for end2 in e2_list_sort:
                        seqid = 100 - end2_list[end2]

                        if start2 < end2:
                           if seqid == 100:
                              fill_color="#005824"
                           elif seqid >= 90:
                              fill_color="#238b45"
                           elif seqid >= 80:
                              fill_color="#41ae76"
                           elif seqid >= 70:
                              fill_color="#66c2a4"
                           elif seqid >= 60:
                              fill_color="#99d8c9"
                           elif seqid >= 50:
                              fill_color="#ccece6"
                           elif seqid < 50:
                              fill_color="#edf8fb"

                        else:#### inverted hits
                           if seqid == 100:
                              fill_color="#99000d"
                           elif seqid >= 90:
                              fill_color="#cb181d"
                           elif seqid >= 80:
                              fill_color="#ef3b2c"
                           elif seqid >= 70:
                              fill_color="#fb6a4a"
                           elif seqid >= 60:
                              fill_color="#fc9272"
                           elif seqid >= 50:
                              fill_color="#fcbba1"
                           elif seqid < 50:
                              fill_color="#fee5d9"

                        if end2_list[end2] <= mismatch: ###does it pass the mismatch cutoff?

                           ### LINES
                           yq1 = data['mid'] - 50 - start1
                           yq2 = data['mid'] - 50 - end1

                           ### axis2
                           xstart_axis2 = data['midwidth'] + 50
                           ystart_axis2 = data['mid'] + 45
                           xt1 = (math.sin(45)*start2) + xstart_axis2
                           yt1 = (math.cos(45)*start2) + ystart_axis2
                           xt2 = (math.sin(45)*end2) + xstart_axis2
                           yt2 = (math.cos(45)*end2) + ystart_axis2
                           ### for the bezier point
                           xt1e = xt1
                           xt2e = xt2
                           yq1e = ((yt1 - yq1)/3) + yq1
                           yq2e = ((yt2 - yq2)/3) + yq2

                           xml.write("<path d=\"M %f %f L %f %f Q %f %f %f %f L %f %f Q %f %f %f %f\" stroke=\"%s\" stroke-width=\"0.5\" fill=\"%s\" fill-opacity=\"%f\"/>\n" % (xq2,yq2,xq1,yq1,xt1e,yq1e,xt1,yt1,xt2,yt2,xt2e,yq2e,xq2,yq2,fill_color,fill_color,alpha))

      ###DRAW ALIGN 1vs3
      ### axis 1
      xq1 = data['midwidth'] - 5
      xq2 = xq1

      for qry in match2:
         allhit = match2[qry]
         for hit in allhit:
            start1_list = allhit[hit]
            s1_list_sort=list(start1_list.keys())
            s1_list_sort.sort()
            for start1 in s1_list_sort:
               end1_list=start1_list[start1]
               e1_list_sort=list(end1_list.keys())
               e1_list_sort.sort()
               for end1 in e1_list_sort:
                  start2_list=end1_list[end1]
                  s2_list_sort=list(start2_list.keys())
                  s2_list_sort.sort()
                  for start2 in s2_list_sort:
                     end2_list=start2_list[start2]
                     e2_list_sort=list(end2_list.keys())
                     e2_list_sort.sort()
                     for end2 in e2_list_sort:
                        seqid = 100 - end2_list[end2]

                        if start2 < end2:
                           if seqid == 100:
                              fill_color="#005824"
                           elif seqid >= 90:
                              fill_color="#238b45"
                           elif seqid >= 80:
                              fill_color="#41ae76"
                           elif seqid >= 70:
                              fill_color="#66c2a4"
                           elif seqid >= 60:
                              fill_color="#99d8c9"
                           elif seqid >= 50:
                              fill_color="#ccece6"
                           elif seqid < 50:
                              fill_color="#edf8fb"

                        else:#### inverted hits
                           if seqid == 100:
                              fill_color="#99000d"
                           elif seqid >= 90:
                              fill_color="#cb181d"
                           elif seqid >= 80:
                              fill_color="#ef3b2c"
                           elif seqid >= 70:
                              fill_color="#fb6a4a"
                           elif seqid >= 60:
                              fill_color="#fc9272"
                           elif seqid >= 50:
                              fill_color="#fcbba1"
                           elif seqid < 50:
                              fill_color="#fee5d9"

                        if end2_list[end2] <= mismatch: ###does it pass the mismatch cutoff?

                           ### LINES
                           yq1 = data['mid'] - 50 - start1
                           yq2 = data['mid'] - 50 - end1

                           ### axis3
                           xstart_axis3 = data['midwidth'] - 50
                           ystart_axis3 = data['mid'] + 45
                           xt1 = (math.sin(-45)*start2) + xstart_axis3
                           yt1 = (math.cos(-45)*start2) + ystart_axis3
                           xt2 = (math.sin(-45)*end2) + xstart_axis3
                           yt2 = (math.cos(-45)*end2) + ystart_axis3
                           ### for the bezier point
                           xt1e = xt1
                           xt2e = xt2
                           yq1e = ((yt1 - yq1)/3) + yq1
                           yq2e = ((yt2 - yq2)/3) + yq2

                           xml.write("<path d=\"M %f %f L %f %f Q %f %f %f %f L %f %f Q %f %f %f %f\" stroke=\"%s\" stroke-width=\"0.5\" fill=\"%s\" fill-opacity=\"%f\"/>\n" % (xq2,yq2,xq1,yq1,xt1e,yq1e,xt1,yt1,xt2,yt2,xt2e,yq2e,xq2,yq2,fill_color,fill_color,alpha))

      ###DRAW ALIGN 3vs2
      for qry in match3:
         allhit = match3[qry]
         for hit in allhit:
            start1_list = allhit[hit]
            s1_list_sort=list(start1_list.keys())
            s1_list_sort.sort()
            for start1 in s1_list_sort:
               end1_list=start1_list[start1]
               e1_list_sort=list(end1_list.keys())
               e1_list_sort.sort()
               for end1 in e1_list_sort:
                  start2_list=end1_list[end1]
                  s2_list_sort=list(start2_list.keys())
                  s2_list_sort.sort()
                  for start2 in s2_list_sort:
                     end2_list=start2_list[start2]
                     e2_list_sort=list(end2_list.keys())
                     e2_list_sort.sort()
                     for end2 in e2_list_sort:
                        seqid = 100 - end2_list[end2]

                        if start2 < end2:
                           if seqid == 100:
                              fill_color="#005824"
                           elif seqid >= 90:
                              fill_color="#238b45"
                           elif seqid >= 80:
                              fill_color="#41ae76"
                           elif seqid >= 70:
                              fill_color="#66c2a4"
                           elif seqid >= 60:
                              fill_color="#99d8c9"
                           elif seqid >= 50:
                              fill_color="#ccece6"
                           elif seqid < 50:
                              fill_color="#edf8fb"

                        else:#### inverted hits
                           if seqid == 100:
                              fill_color="#99000d"
                           elif seqid >= 90:
                              fill_color="#cb181d"
                           elif seqid >= 80:
                              fill_color="#ef3b2c"
                           elif seqid >= 70:
                              fill_color="#fb6a4a"
                           elif seqid >= 60:
                              fill_color="#fc9272"
                           elif seqid >= 50:
                              fill_color="#fcbba1"
                           elif seqid < 50:
                              fill_color="#fee5d9"

                        if end2_list[end2] <= mismatch: ###does it pass the mismatch cutoff?

                           ### axis2
                           xstart_axis2 = data['midwidth'] + 50
                           ystart_axis2 = data['mid'] + 55
                           xq1 = (math.sin(45)*start2) + xstart_axis2
                           yq1 = (math.cos(45)*start2) + ystart_axis2
                           xq2 = (math.sin(45)*end2) + xstart_axis2
                           yq2 = (math.cos(45)*end2) + ystart_axis2

                           ### axis3
                           xstart_axis3 = data['midwidth'] - 50
                           ystart_axis3 = data['mid'] + 55
                           xt1 = (math.sin(-45)*start1) + xstart_axis3
                           yt1 = (math.cos(-45)*start1) + ystart_axis3
                           xt2 = (math.sin(-45)*end1) + xstart_axis3
                           yt2 = (math.cos(-45)*end1) + ystart_axis3

                           #print "xq1 %f xt1 %f\n" % (xq1,xt1)
                           ### for the bezier point XXXX
                           #if (xq1 - data['midwidth']) > (data['midwidth'] - xt1):
                           #   x1e = data['midwidth'] + (xq1-data['midwidth'])/30
                           #else:
                           #   x1e = data['midwidth'] - (data['midwidth'] - xt1)/30
                           #x2e = x1e
                           x1e = xq1 - ((xq1 - xt1)/2)
                           x2e = x1e
                           if yq1 > yt1:
                               #y1e = yq1
                               #y2e = yq2
                               y1e = yq1 + (yq1 - ystart_axis3)*0.85 + (0.015 * yq1)
                               y2e = yq2 + (yq2 - ystart_axis3)*0.85 + (0.015 * yq2)
                           else:
                               #y1e = yt1
                               #y2e = yt2
                               y1e = yt1 + (yt1 - ystart_axis3)*0.85 + (0.015 * yt1)
                               y2e = yt2 + (yt2 - ystart_axis3)*0.85 + (0.015 * yt2)

                           
                           xml.write("<path d=\"M %f %f L %f %f Q %f %f %f %f L %f %f Q %f %f %f %f\" stroke=\"%s\" stroke-width=\"0.5\" fill=\"%s\" fill-opacity=\"%f\"/>\n" % (xq2,yq2,xq1,yq1,x1e,y1e,xt1,yt1,xt2,yt2,x2e,y2e,xq2,yq2,fill_color,fill_color,alpha))

      ######### END PLOT ALIGNMENTS
      ### DRAW LEGEND

      fcol = ["#005824","#238b45","#41ae76","#66c2a4","#99d8c9","#ccece6","#edf8fb"]
      rcol = ["#99000d","#cb181d","#ef3b2c","#fb6a4a","#fc9272","#fcbba1","#fee5d9"]
      sic = ["100","90+","80+","70+","60+","50+","0-49"]

      xt = xlegend
      yt = ylegend

      xtl = xt + 30
      xtl2 = xtl + 85
      xtl3 = xtl2 + 100
      xtl4 = xtl2 + 50
      yt2 = yt + 50

      xml.write("<text font-size=\"2.2em\" x=\"%f\" y=\"%f\">Sequence identity (%%)</text>\n" % (xt,yt))
      xml.write("<text font-size=\"2em\" x=\"%f\" y=\"%f\">Forward | Reverse</text>\n" % (xtl,yt2))

      el = 0
      for fc in fcol:
         yt2 += 30

         xml.write("<rect x=\"%f\" y=\"%f\" width=\"20\" height=\"20\" style=\"fill:%s;stroke:black;stroke-width:0.5;fill-opacity:%s;\" />\n"% (xtl2,yt2,fc,alpha))
         ytt = yt2 + 20
         xml.write("<text font-size=\"2em\" x=\"%f\" y=\"%f\">%s</text>\n" % (xtl3,ytt,sic[el]))
         xml.write("<rect x=\"%f\" y=\"%f\" width=\"20\" height=\"20\" style=\"fill:%s;stroke:black;stroke-width:0.5;fill-opacity:%s;\" />\n" % (xtl4,yt2,rcol[el],alpha))
         el+=1

      #### DONE PLOTTING
      ####### WRAP XML and CLOSE FILES
      xml.write("</svg>")
      print("xmatchview-hive svg output graph in %s\n" % (output))
      xml.close()




#---------------------------------------------
def main():
    opts, args = getopt.getopt(sys.argv[1:], "x:y:z:q:r:s:d:e:f:i:b:c:a:")

    (gff_file1, gff_file2, gff_file3, align_file1, align_file2, align_file3, txt_file1, txt_file2, txt_file3)=(None,None,None,None,None,None,None,None,None)
    (seqidentity, block_length, scale, protein, alpha)=(0,0,0,0,0.75)

    for o, v in opts:
      if o == "-x":
        align_file1=str(v)
      if o == "-y":
        align_file2=str(v)
      if o == "-z":
        align_file3=str(v)
      if o == "-q":
        txt_file1=str(v)
      if o == "-r":
        txt_file2=str(v)
      if o == "-s":
        txt_file3=str(v)
      if o == "-i": 
        seqidentity=int(v)
      if o == "-b":
        block_length=int(v)
      if o == "-c":
        scale=int(v)
      if o == "-d":
        gff_file1=str(v)
      if o == "-e":
        gff_file2=str(v)
      if o == "-f":
        gff_file3=str(v)
      if o == "-a":
        alpha = float(v)

    if (txt_file1 == None or txt_file2 == None or txt_file3 == None or align_file1 == None or align_file2 == None or align_file3 == None or block_length == 0 or scale == 0):
      print("Usage: %s v1.2.4" % (sys.argv[0:]))
      print("-x alignment file [1 vs. 2] (cross_match .rep or Pairwise mApping Format .paf)")
      print("-y alignment file [1 vs. 3] (cross_match .rep or Pairwise mApping Format .paf)")
      print("-z alignment file [3 vs. 2] (cross_match .rep or Pairwise mApping Format .paf)")
      print("-q genome text file 1 (format NAME:LENGTH)")
      print("-r genome text file 2 (format NAME:LENGTH)")
      print("-s genome text file 3 (format NAME:LENGTH)")
      print("-d features (eg. exons) coordinates GFF tsv file 1 (start end) - optional")
      print("-e features (eg. exons) coordinates GFF tsv file 2 (start end) - optional")
      print("-f features (eg. exons) coordinates GFF tsv file 3 (start end) - optional")
      print("-i sequence identity threshold (e.g. -i 90 will show colinear blocks >= 90% sequence identity)")
      print("-b minimum length (bp) of similarity block to display")
      print("-c scale (pixel to basepair scale, for displaying the image)")
      #print "-l label for the tree trunk (6 characters or less for best result)"
      print("-a alpha value, from 0 (transparent) to 1 (solid, default)")
      #print "-f output image file format (png, tiff, jpeg, or gif) NOTE: png and tiff recommended."
      #print "-p full path to the directory with fonts on your system (please refer to the documentation for fonts used)"
      #print "-z transform bacterial ORF into protein (i.e. plot alignment between ORF products? 1/0) DEPRECATED\n"
      print("* Files for the -q, -r and -s options must include header_names:base_length, with names that correspond to those in fasta files used to run cross_match or minimap2\n")
      print("! Ensure the config.txt file exists in your run directory")
      sys.exit(1)

   #====seqidentity checks
    if (seqidentity <0 or seqidentity >100):
      print("-i must be a valid number between 0-100")
      sys.exit(1)

   #====Alpha checks
    if (alpha<0 or alpha >1):
      print("-a must be a valid number between 0-1")
      sys.exit(1)

   #===Scale checks
    if (scale<1):
      print("Not a possible scale. Make sure you select a number >1.")
      sys.exit(1)

   #====File checks
    checkFile(txt_file1)
    checkFile(txt_file2)
    checkFile(txt_file3)
    checkFile(align_file1)
    checkFile(align_file2)
    checkFile(align_file3)
    checkFile("config.txt")

   #====Parse config.txt
    info = {}
    print("Reading configuration..\n")
    with open("config.txt") as fd:
      rd = csv.reader(fd, delimiter=":", quotechar='"')
      for row in rd:
         info[row[0]] = row[1]
         print("axis %s = %s\n" %(row[0],row[1]))
    
    
   ###OPTIONAL, FOR FEATURES REPRESENTATION
    (gff1,gff2,gff3) = ({},{},{})

    if(gff_file1 != None):
       checkFile(gff_file1)
       print("Reading reference feature file %s ..." % (gff_file1))
       gff1=readGFF(gff_file1,scale)
       print("done.")

    if(gff_file2 != None):
       checkFile(gff_file2)
       print("Reading reference feature file %s ..." % (gff_file2))
       gff2=readGFF(gff_file2,scale)
       print("done.")

    if(gff_file3 != None):
       checkFile(gff_file3)
       print("Reading reference feature file %s ..." % (gff_file3))
       gff3=readGFF(gff_file3,scale)
       print("done.")

    #====Parse Fasta Files
    (order1, ref1, length1)=readText(txt_file1, scale)
    (order2, ref2, length2)=readText(txt_file2, scale)
    (order3, ref3, length3)=readText(txt_file3, scale)

    #====Raise error if features out of bounds
    data=initGraph()
    if 2 * length1 / data['width'] > scale:
       estscale = int(2 * length1 / data['width']) + 1
       #sys.exit("\n\n! The 1st sequence in %s is predicted to extend beyond the plot width, you must increase the scale to at least %i (we suggest rounding up to the next ten, hundred or thousand) -- fatal." % (txt_file1, estscale))
    if 2 * length2 / data['width'] > scale:
       estscale = int(2 * (length2) / data['width']) + 1
       #sys.exit("\n\n! The 2nd sequence in %s is predicted to extend beyond the plot width, you must increase the scale to at least %i (we suggest rounding up to the next ten, hundred or thousand) -- fatal." % (txt_file2, estscale))
    if 2 * length3 / data['width'] > scale:
       estscale = int(2 * (length3) / data['width']) + 1
       sys.exit("\n\n! The 3rd sequence in %s is predicted to extend beyond the plot width, you must increase the scale to at least %i (we suggest rounding up to the next ten, hundred or thousand) -- fatal." % (txt_file3, estscale))

    mismatch = 100 - seqidentity 

    #====Parse Alignment files
    print("Reading alignment file 1/3 %s ..." % (align_file1))
    (match1, match2, match3) = ({},{},{})
    if align_file1.endswith("rep"):
      match1 = readCrossMatch(align_file1, mismatch, block_length, ref1, ref2, scale)
    elif align_file1.endswith("paf"):
      match1 = readPAF(align_file1, mismatch, block_length, ref1, ref2, scale)
    else:
      print("The alignment file provided (-x %s) does not end in .rep (cross_match) or .paf (PAF) -- fatal" % align_file1)
      sys.exit(1)
    print("done.")

    print("Reading alignment file 2/3 %s ..." % (align_file2))
    match = {}
    if align_file2.endswith("rep"):
      match2 = readCrossMatch(align_file2, mismatch, block_length, ref1, ref3, scale)
    elif align_file2.endswith("paf"):
      match2 = readPAF(align_file2, mismatch, block_length, ref1, ref3, scale)
    else:
      print("The alignment file provided (-x %s) does not end in .rep (cross_match) or .paf (PAF) -- fatal" % align_file2)
      sys.exit(1)
    print("done.")

    print("Reading alignment file 3/3 %s ..." % (align_file3))
    match = {}
    if align_file3.endswith("rep"):
      match3 = readCrossMatch(align_file3, mismatch, block_length, ref3, ref2, scale)
    elif align_file3.endswith("paf"):
      match3 = readPAF(align_file3, mismatch, block_length, ref3, ref2, scale)
    else:
      print("The alignment file provided (-x %s) does not end in .rep (cross_match) or .paf (PAF) -- fatal" % align_file3)
      sys.exit(1)
    print("done.")

    print("Drawing ...")
    drawRelationship(info,gff1,gff2,gff3,order1,order2,order3,ref1,ref2,ref3,length1,length2,length3,match1,match2,match3,scale,seqidentity,mismatch,block_length,alpha)

#---------------------------------------------
#Main Call

main()
sys.exit(1)



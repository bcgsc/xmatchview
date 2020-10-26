#!/usr/bin/env python3
# xmatchview-conifer.py
# Visualizing genome synteny with an evergreen representation
# Rene L Warren 2005-2020

import sys
import os
import getopt
import re
import csv
from PIL import Image
from PIL import ImageDraw
from PIL import ImageFont
from PIL import ImageEnhance
from PIL import PSDraw
import subprocess

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
         print("__%s__ - __%s__ <<<<<" % (start,end))

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
         print("INITIALIZED %s : %i,%i with %s" %(row[0],start,end,color))

   return feature

#---------------------------------------------
def readPAF(paf_file,mismatch,block_length,reference,query,scale):

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

               print("%i-%i   ::   %i-%i" % (startFirstMatch,endFirstMatch,startSecondMatch,endSecondMatch))

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

               print("%i-%i   ::   %i-%i" % (startFirstMatch,endFirstMatch,startSecondMatch,endSecondMatch))

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
def readCrossMatch(crossmatch_file,mismatch,block_length,reference,query,scale):

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

               print("REVERSE %i-%i   ::   %i-%i" % (startFirstMatch,endFirstMatch,startSecondMatch,endSecondMatch))

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

               print("FORWARD %i-%i  ::   %i-%i" % (startFirstMatch,endFirstMatch,startSecondMatch,endSecondMatch))
 
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
def generateCoords(nocdt, size, leap, protein):
 
   freq={}

   pos_range=list(range(0,size,leap))

   for pos in pos_range:
      print("%i out of %i" % (pos,size))
      for reference in nocdt:
         for comparison in nocdt[reference]:
            start1_dict=list(nocdt[reference][comparison].keys())
            start1_dict.sort()
            for start1 in start1_dict:
               end1_dict=list(nocdt[reference][comparison][start1].keys())
               end1_dict.sort()
               for end1 in end1_dict:
                  start2_dict=list(nocdt[reference][comparison][start1][end1].keys())
                  start2_dict.sort()

                  (ss,ee) = (start1,end1)

                  if protein:
                     size_ref = end1 - start1
                     buffer = ((size_ref - (size_ref/3)) / 2)
                     ss = start1 + buffer
                     ee = end1 - buffer  

                  if((pos >= ss and pos <= ee) or (pos >= ee and pos <= ss)):
                     #print "%i >= %i and %i<=%i OR %i>=%i and %i<=%i" % (pos,ss,pos,ee,pos,ee,pos,ss)
                     for start2 in start2_dict:
                        end2_dict=list(nocdt[reference][comparison][start1][end1][start2].keys())
                        end2_dict.sort()
                        for end2 in end2_dict:
                           current_mismatch=float(nocdt[reference][comparison][start1][end1][start2][end2])
                           if pos not in freq:
                              freq[pos]={}
                           if current_mismatch not in freq[pos]:
                              freq[pos][current_mismatch]=int(0)
                           freq[pos][current_mismatch]=freq[pos][current_mismatch]+1
   return freq

#---------------------------------------------
def findOccurences(s, ch):
    return [i for i, letter in enumerate(s) if letter == ch]

#---------------------------------------------
def readFasta(file, scale):

   (head_match, previous_contig,seq_length) = (None,None,0)
   L1={}
   order=[]#RLW
   npos={} #RLW
   tot_length = 0 #RLW

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
             print("NPOS: (Ns only tracked on 1-line sequences)")
             print(npos)

             order.append(previous_contig) #RLW
             print("%s length = %i bp, scaled to %.0f pixels" % (previous_contig,seq_length,L1[previous_contig]['scaled_len'])) #RLW
             tot_length += seq_length #RLW
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
   print("NPOS: (Ns only tracked on 1-line sequences)")
   print(npos)

   order.append(previous_contig) #RLW
   print("%s length = %i bp, scaled to %.0f pixels" % (previous_contig,seq_length,L1[previous_contig]['scaled_len'])) #RLW
   tot_length += seq_length #RLW

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
   data['width']=2000 #5000
   data['height']=2000 #5000
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
def drawRectangle(draw,start,end,y,thickness,bar_color,text,font,text_color):
   
   draw.rectangle((start,y,end,y+thickness), bar_color)
   draw.text((start-80, y), text, font=font, fill=text_color)

#---------------------------------------------
def plotFrequency(freq,size,scale,draw,color,data,leap):

   pos_range=list(range(0,size,leap))
   
   for pos in pos_range:
      if pos in freq:
         freq_list=freq[pos]
         previous=data['mis_bar']
         identity_range=list(range(9,-1,-1))
         for id in identity_range:
            cumul=int(0)
            for freq_keys in freq_list:
               if id >= freq_keys:
                  cumul += freq_list[freq_keys]
            
            if cumul<1:
               color_now="white"
            elif cumul==1:
               color_now="blue"
            elif cumul==2:
               color_now="cyan"       
            elif cumul==3:
               color_now="green"
            elif cumul==4:
               color_now="dirtyred"
            elif cumul==5:
               color_now="purple"
            elif cumul==6:
               color_now="salmon"
            elif cumul==7:
               color_now="orange"
            elif cumul>=8:
               color_now="yellow"

            extension=((200-(20*id))+data['mis_bar'])   #y
            compressed=(pos/scale)+data['x']           #x
            
            if color_now != "white":
               #print "%i, %i, %i, %i %s" % (compressed,previous,compressed,extension,color_now)
               draw.line((compressed,previous,compressed,extension),color[color_now])
 
            previous = extension


#---------------------------------------------
def drawRelationship(reference_list, order_ref, query_list, order_qry, match_list, scale, mismatch, block_length, alignment_file, reflength, format, formatdict, protein, label, alpha, refgff, qrygff, qrylength, fontpath):

      scaled_reflength=float(reflength/scale)
      scaled_qrylength=float(qrylength/scale)

      ###Capture last coordinates of relationships
      (u2max,v2max,x2max,y2max)=(0,0,0,0)

      ###Initialize new graph
      data=initGraph()

      ###Get colors
      color=initColor(alpha)

      ###Set Font
      arialfont = fontpath + "/arial.ttf"
      pilfont = fontpath + "/helvR14.pil"

      #default all font sizes to default (it is quite small, you must provide a valid path for best results)
      font_18=ImageFont.load_default()
      font_20=ImageFont.load_default()
      fontb_20=ImageFont.load_default()
      fontbi_20=ImageFont.load_default()
      fontb_22=ImageFont.load_default()
      font_24=ImageFont.load_default()
      fontb_24=ImageFont.load_default()
      font_28=ImageFont.load_default()
      fontb_28=ImageFont.load_default()
      fontbi_28=ImageFont.load_default()
      fontb_92=ImageFont.load_default()
      fontb_80=ImageFont.load_default()

      if os.path.exists(arialfont): ### Will check for truetype first, they look better
          ###Set Font (truetype)
          font_18=ImageFont.truetype(fontpath + "/arial.ttf",18)
          font_20=ImageFont.truetype(fontpath + "/arial.ttf",20)
          fontb_20=ImageFont.truetype(fontpath + "/arialbd.ttf",20)
          fontbi_20=ImageFont.truetype(fontpath + "/arialbi.ttf",20)
          fontb_22=ImageFont.truetype(fontpath + "/arialbd.ttf",22)

          font_24=ImageFont.truetype(fontpath + "/arial.ttf",24)
          fontb_24=ImageFont.truetype(fontpath + "/arialbd.ttf",24)

          ### XXX change
          font_24=ImageFont.truetype(fontpath + "/arialbd.ttf",30)
          fontb_24=ImageFont.truetype(fontpath + "/arialbd.ttf",32)


          font_28=ImageFont.truetype(fontpath + "/arial.ttf",28)
          fontb_28=ImageFont.truetype(fontpath + "/arialbd.ttf",28)
          fontbi_28=ImageFont.truetype(fontpath + "/arialbi.ttf",28)

          ### XXX change
          font_28=ImageFont.truetype(fontpath + "/arialbd.ttf",40)
          fontb_28=ImageFont.truetype(fontpath + "/arialbd.ttf",40)
          fontbi_28=ImageFont.truetype(fontpath + "/arialbi.ttf",34)

          fontb_92=ImageFont.truetype(fontpath + "/arialbd.ttf",92)
          fontb_80=ImageFont.truetype(fontpath + "/arialbd.ttf",78)
      elif os.path.exists(pilfont): ### Will settle for PIL font, if ttf do not exist. Otherwise, sticking with default.
          ###Set font (pil) (sizes are limited, made to be compatible with TT fonts)
          font_18=ImageFont.load_path(fontpath + "/helvR14.pil")
          font_20=ImageFont.load_path(fontpath + "/helvR18.pil")
          fontb_20=ImageFont.load_path(fontpath + "/helvB18.pil")
          fontbi_20=ImageFont.load_path(fontpath + "/helvBO18.pil")
          fontb_22=ImageFont.load_path(fontpath + "/helvB24.pil")
          font_24=ImageFont.load_path(fontpath + "/helvR24.pil")
          fontb_24=ImageFont.load_path(fontpath + "/helvB24.pil")
          font_28=ImageFont.load_path(fontpath + "/helvR24.pil")
          fontb_28=ImageFont.load_path(fontpath + "/helvB24.pil")
          fontbi_28=ImageFont.load_path(fontpath + "/helvBO24.pil")
          fontb_92=ImageFont.load_path(fontpath + "/helvR24.pil")
          fontb_80=ImageFont.load_path(fontpath + "/helvR24.pil")

      ###Define Image
      back = Image.new("RGBA", (data['width'],data['height']),(0,0,0,0))
      bdraw = ImageDraw.Draw(back)

      poly = Image.new("RGBA", (data['width'],data['height']))
      draw = ImageDraw.Draw(poly)

      ticklabel = Image.new("RGBA", (data['width'],data['height']))

      back = back.rotate(90)

      decay = data['decay']

      ### REFERENCE slope calculations -- IMPORTANT, USED THROUGHOUT
      x1ref = data['x']
      x2ref = data['x'] + scaled_reflength
      y1ref = data['ref_y']
      y2ref = data['ref_y']-data['skew']

      mrref = (y2ref - y1ref ) / (x2ref - x1ref)
      brref = y2ref - (mrref * x2ref)

      ###QRY slope calculations -- IMPORTANT, USED THROUGHOUT 
      x1qry = data['x']
      x2qry = data['x'] + scaled_qrylength
      y1qry = data['ref_y']+decay
      y2qry = data['ref_y']+decay+data['skew']

      mqqry = (y2qry - y1qry ) / (x2qry - x1qry)
      bqqry = y2qry - (mqqry * x2qry)

      ####REFERENCE
      for ref in order_ref:      #RLW
         scaled_offset_len = 0   #RLW
         if ref in reference_list: #RLW
            scaled_offset_len = reference_list[ref]['offset_len'] #RLW

         init_coord=int(data['x']+scaled_offset_len) #RLW
         last_coord=int(data['x']+scaled_offset_len+reference_list[ref]['scaled_len']) #RLW      

         a1 = data['x'] + scaled_offset_len
         a2 = data['x'] + scaled_offset_len+reference_list[ref]['scaled_len']
         b1 = (mrref * a1 ) + brref
         b2 = (mrref * a2 ) + brref
         draw.polygon((a1,b1,a1,b1+data['reference_thick'],a2,b2+data['reference_thick'],a2,b2),outline=color['brown'], fill=color['brown'])###references rect
         draw.text((init_coord+5, data['ref_y_skew']-53), ref, font=fontb_24, fill=color['grey'])### XXX WAS fill=color['green9'])###label for ref XXXXTOCHANGE
         ####draw features/exons on side of ref
         ###REF gene model
         if ref in refgff:
            for exstart in refgff[ref]:
                exend = refgff[ref][exstart]['end']
                l1 = data['x'] + exstart + scaled_offset_len
                l2 = data['x'] + exend + scaled_offset_len
                m1 = (mrref * l1 ) + brref
                m2 = (mrref * l2 ) + brref
                draw.polygon((l1,m1-11,l1,m1,l2,m2,l2,m2-11),outline=refgff[ref][exstart]['color'], fill=refgff[ref][exstart]['color'])###features/exons

         ###ticks on reference
         draw.line((a1,b1-50,a1,b1),color['grey'],width=1) ### ZZZZ
         draw.line((a2,b2-50,a2,b2),color['grey'],width=1) ### ZZZZ

         back.paste(poly, mask=poly)
         del draw
         poly = Image.new("RGBA", (data['width'],data['height']))
         draw = ImageDraw.Draw(poly)

      last_ref_coord = last_coord

      (current_position, LCB, skew,stop)=(data['x'], 10, data['skew'],data['x'])
       
      ####Draw Query (only if not in reference list)
      decay = data['decay']
      for qry in order_qry:      #RLW

         scaled_offset_len = 0   #RLW
         if qry in query_list: #RLW
            scaled_offset_len = query_list[qry]['offset_len'] #RLW

         init_coord=int(data['x']+scaled_offset_len) #RLW
         last_coord=int(data['x']+scaled_offset_len+query_list[qry]['scaled_len']) #RLW
         stop = current_position + query_list[qry]['offset_len'] + query_list[qry]['scaled_len'] #RLW

         a1 = data['x'] + scaled_offset_len
         a2 = data['x'] + scaled_offset_len+query_list[qry]['scaled_len']
         b1 = (mqqry * a1 ) + bqqry
         b2 = (mqqry * a2 ) + bqqry

         ### draw one side of tree
         draw.polygon((a1,b1,a1,b1+(data['query_thick']),a2,b2+(data['query_thick']),a2,b2),outline=color['brown'], fill=color['brown'])###queries rect
         draw.text((init_coord+5, data['ref_y']+decay+skew+36), qry, font=fontb_24, fill=color['grey'])### XXX WASfill=color['green9'])### label for query

         ####draw features/exons on side of query
         ###QRY gene model
         if qry in qrygff:
            for exstart in qrygff[qry]:
               exend = qrygff[qry][exstart]['end']
               l1 = data['x'] + exstart + scaled_offset_len
               l2 = data['x'] + exend + scaled_offset_len
               m1 = (mqqry * l1 ) + bqqry
               m2 = (mqqry * l2 ) + bqqry
               draw.polygon((l1,m1+(data['query_thick'])+1,l1,m1+(data['query_thick'])+11,l2,m2+(data['query_thick'])+11,l2,m2+(data['query_thick'])+1),outline=qrygff[qry][exstart]['color'], fill=qrygff[qry][exstart]['color'])###features/exons

         ###ticks on query
         draw.line((a1,b1+(data['query_thick']),a1,b1+(data['query_thick'])+50),color['grey'],width=1) ### ZZZZ
         draw.line((a2,b2+(data['query_thick']),a2,b2+(data['query_thick'])+50),color['grey'],width=1) ### ZZZZ

         back.paste(poly, mask=poly)
         del draw
         poly = Image.new("RGBA", (data['width'],data['height']))
         draw = ImageDraw.Draw(poly)

      last_qry_coord = last_coord


      ###DRAW BLOCKS
      plotflag = 0
      for qry in match_list:
         allhit = match_list[qry]
         for hit in allhit:
            start1_list = allhit[hit]
            stop = current_position + query_list[qry]['offset_len'] + query_list[qry]['scaled_len']

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
                        print("si=%.2f mis=%.2f" % (seqid,end2_list[end2]))

                        if start2 < end2:
                           if seqid >=99:
                              fill_color="green10t"
                           elif seqid >= 95:
                              fill_color="green9t"
                           elif seqid >= 90:
                              fill_color="green8t"
                           elif seqid >= 85:
                              fill_color="green7t"
                           elif seqid >= 80:
                              fill_color="green6t"
                           elif seqid >= 75:
                              fill_color="green5t"
                           elif seqid >= 70:
                              fill_color="green4t"
                           elif seqid >= 65:
                              fill_color="green3t"
                           elif seqid >= 60:
                              fill_color="green2t"
                           elif seqid >0:
                              fill_color="green1t"
                           outline_color = "green10t"

                        else:#### inverted hits
                           if seqid >=99:
                              fill_color="red10t"
                           elif seqid >= 95:
                              fill_color="red9t"
                           elif seqid >= 90:
                              fill_color="red8t"
                           elif seqid >= 85:
                              fill_color="red7t"
                           elif seqid >= 80:
                              fill_color="red6t"
                           elif seqid >= 75:
                              fill_color="red5t"
                           elif seqid >= 70:
                              fill_color="red4t"
                           elif seqid >= 65:
                              fill_color="red3t"
                           elif seqid >= 60:
                              fill_color="red2t"
                           elif seqid >0:
                              fill_color="red1t"
                           outline_color = "red10t"

                        ###draw ORF on upper
                        #draw.rectangle((data['x']+start1,data['ref_y']+1,data['x']+end1,data['ref_y']+data['reference_thick']-1), outline=color["lightgrey"], fill=color["lightgrey"])
                        size_qry = end1 - start1
                        size_ref = end2 - start2
                        buf_ref = ((size_ref - (size_ref/3)) / 2)
                        buf_qry = ((size_qry - (size_qry/3)) / 2)
                        print("%s (%i-%i) hits %s  ::  mismatch %.2f target(%i)  block %i target (%i) " % (qry,start1,end1,hit,end2_list[end2],mismatch,size_ref,block_length))

                        if end2_list[end2] <= mismatch: ###does it pass the mismatch cutoff?

                           x1 = data['x'] + start1
                           x2 = data['x'] + end1
                           y1 = (mqqry * x1 ) + bqqry
                           y2 = (mqqry * x2 ) + bqqry
                           print("x1=%i y1=%i x2=%i y2=%i M=%.2f  B=%.2f  " % (x1,y1,x2,y2,mqqry,bqqry)) 
                           u1 = data['x'] + start2
                           u2 = data['x'] + end2
                           v1 = (mrref * u1 ) + brref
                           v2 = (mrref * u2 ) + brref
                           print("u1=%i v1=%i u2=%i v2=%i M=%.2f  B=%.2f  " % (u1,v1,u2,v2,mrref,brref))

                           if x2 > x2max:
                              u2max = u2
                              v2max = v2
                              x2max = x2
                              y2max = y2


                           ### LINES
                           draw.polygon((u2,v2+data['reference_thick'],x2,y2,x1,y1,u1,v1+data['reference_thick']), outline=color[outline_color], fill=color[fill_color])
                           ### REPEAT FEATURE
                           draw.polygon((x1,y1,x1,y1+data['reference_thick'],x2,y2+data['reference_thick'],x2,y2),outline=color[outline_color], fill=color[fill_color])###colinear block on query
                           back.paste(poly, mask=poly)
                           draw.polygon((u1,v1,u2,v2,u2,v2+data['reference_thick'],u1,v1+data['reference_thick']),outline=color[outline_color], fill=color[fill_color])###colinear block on reference
                           back.paste(poly, mask=poly)
                           del draw
                           poly = Image.new("RGBA", (data['width'],data['height']))
                           draw = ImageDraw.Draw(poly)
                           plotflag = 1

      if(plotflag == 0):
         print("It looks like there is nothing to plot, try increasing -m 99 -- FATAL")
         sys.exit(1)

      #enhancer = ImageEnhance.Sharpness(im)
      #for i in range(8):
      #   factor = i / 4.0
      #   enhancer.enhance(factor).show("Sharpness %f" % factor)


      ### draw Ns on reference
      for ref in order_ref:

         scaled_offset_len = 0   #RLW
         if ref in reference_list: #RLW
            scaled_offset_len = reference_list[ref]['offset_len'] #RLW

         if 'npos' in reference_list[ref]:
            for nstart in reference_list[ref]['npos']:  ### draw Ns
               nstart = data['x'] + (nstart/scale) + scaled_offset_len
               ny = (mrref * nstart) + brref
               draw.line((nstart,ny,nstart,ny+data['reference_thick']-1),color['red'],width=2)

      ### draw Ns on query
      for qry in order_qry:

         scaled_offset_len = 0   #RLW
         if qry in query_list: #RLW
            scaled_offset_len = query_list[qry]['offset_len'] #RLW

         if 'npos' in query_list[qry]:            ### draw Ns
            for nstart in query_list[qry]['npos']:
               nstart = data['x'] + (nstart/scale) + scaled_offset_len
               ny = (mqqry * nstart) + bqqry
               draw.line((nstart,ny+2,nstart,ny+data['reference_thick']+1),color['red'],width=2)


      ### calculate placement of tree trunk/label
      if label != "":
         charwidth = 55 ### approximate character width in pixels
         labellength = len(label)
         totaltrunklength = (labellength + 1) * charwidth ### the 1 is for a one-character buffer before/after
         mtrunk = (y2max - v2max ) / (x2max - u2max + 1)
         btrunk = y2max - (mtrunk * x2max)

         #ytrunk = y1ref
         #if mtrunk > 0: ###positive slope#u2max > x2max:
         #   ytrunk = y1ref + decay #
         #else:
         #   ytrunk = y1ref
         #xtrunk = (ytrunk - btrunk) / mtrunk
         #draw.rectangle((xtrunk+5+addbuffer,y1ref+data['reference_thick']+4,xtrunk+totaltrunklength+addbuffer,y1ref+decay-5), outline=color['brown'], fill=color['brown'])###trunk

         addbuffer=10 #this is sometimes necessary when trunk overlaps with lines

         ytrunk1 = y1ref
         ytrunk2 = y1ref + decay
         xtrunk1 = (ytrunk1 - btrunk) / mtrunk
         xtrunk2 = (ytrunk2 - btrunk) / mtrunk
         draw.polygon((xtrunk1+addbuffer,ytrunk1+data['reference_thick'],xtrunk1+totaltrunklength+addbuffer,ytrunk1+data['reference_thick'],xtrunk1+totaltrunklength+addbuffer,y1ref+decay-5,xtrunk2+addbuffer,y1ref+decay-5), outline=color['brown'], fill=color['brown'])###trunk

         draw.text((xtrunk1+charwidth+addbuffer,data['ref_y']+(decay/4)), label, font=fontb_80, fill=color['beige'])###label ADJUST y+ YYY for position of label
      back.paste(poly, mask=poly)
      ### end trunk code

      ### FINAL IMAGE PROCESSING
      ### rotate plot to be able to place scale
      back = back.rotate(270)
      del draw
      drawtl = ImageDraw.Draw(ticklabel)
      ###final tick labels
      if last_ref_coord > last_qry_coord:
         last_coord = last_ref_coord
      else:
         last_coord = last_qry_coord

      ###ticks
      x_range=list(range(data['x'], int(last_coord), 100))

      for position in x_range:
         drawtl.rectangle((data['x_legend_picto'],position+5,data['x_legend_picto']+15,position+8),color['black'])

      if reflength >= 10000:
         for position in x_range:
            base_number=int(((position-data['x'])*scale)/1000)
            drawtl.text((data['x_legend_picto']+25,position-18), "%i" % base_number, font=font_28, fill=color['black']) ### XXX
      else:
         for position in x_range:
            base_number=(position-data['x']) * scale
            base_number=float(base_number)
            base_number=base_number/1000
            #print "%i %i %i >>> %.2f <<< %i,%i" % (data['x'],position,scale,base_number,data['x_legend_picto'],position)
            drawtl.text((data['x_legend_picto']+25,position-15), "%.1f" % base_number, font=font_28, fill=color['black'])### was position-15

      drawtl.text((data['x_legend_picto']+25,last_coord+25), "kbp", font=fontb_28, fill=color['black'])

      ###Draw Legend
      date=subprocess.getstatusoutput("date")

      ###Picto Legend
      y_legend = last_coord + 100 ### WAS 30 XXX
      drawtl.text((data['x_legend_picto'],y_legend), "Sequence identity (%)", font=fontbi_28, fill=color['black'])
      ####

      y_legend+=40 ### was 35 XXX
      drawtl.rectangle((data['x_legend_picto'],y_legend,data['x_legend_picto']+30,y_legend+30), outline=color['green10t'], fill=color['green10t'])
      drawtl.rectangle((data['x_legend_picto']+31,y_legend,data['x_legend_picto']+61,y_legend+30), outline=color['red10t'], fill=color['red10t'])
      drawtl.text((data['x_legend_picto']+65,y_legend), "99-100", font=font_24, fill=color['black'])
      y_legend+=30
      drawtl.rectangle((data['x_legend_picto'],y_legend,data['x_legend_picto']+30,y_legend+30), outline=color['green9t'], fill=color['green9t'])
      drawtl.rectangle((data['x_legend_picto']+31,y_legend,data['x_legend_picto']+61,y_legend+30), outline=color['red9t'], fill=color['red9t'])
      drawtl.text((data['x_legend_picto']+65,y_legend), "95-98", font=font_24, fill=color['black'])
      y_legend+=30
      drawtl.rectangle((data['x_legend_picto'],y_legend,data['x_legend_picto']+30,y_legend+30), outline=color['green8t'], fill=color['green8t'])
      drawtl.rectangle((data['x_legend_picto']+31,y_legend,data['x_legend_picto']+61,y_legend+30), outline=color['red8t'], fill=color['red8t'])
      drawtl.text((data['x_legend_picto']+65,y_legend), "90-94", font=font_24, fill=color['black'])
      y_legend+=30
      drawtl.rectangle((data['x_legend_picto'],y_legend,data['x_legend_picto']+30,y_legend+30), outline=color['green7t'], fill=color['green7t'])
      drawtl.rectangle((data['x_legend_picto']+31,y_legend,data['x_legend_picto']+61,y_legend+30), outline=color['red7t'], fill=color['red7t'])
      drawtl.text((data['x_legend_picto']+65,y_legend), "85-89", font=font_24, fill=color['black'])
      y_legend+=30
      drawtl.rectangle((data['x_legend_picto'],y_legend,data['x_legend_picto']+30,y_legend+30), outline=color['green6t'], fill=color['green6t'])
      drawtl.rectangle((data['x_legend_picto']+31,y_legend,data['x_legend_picto']+61,y_legend+30), outline=color['red6t'], fill=color['red6t'])
      drawtl.text((data['x_legend_picto']+65,y_legend), "80-84", font=font_24, fill=color['black'])
      y_legend+=30
      drawtl.rectangle((data['x_legend_picto'],y_legend,data['x_legend_picto']+30,y_legend+30), outline=color['green5t'], fill=color['green5t'])
      drawtl.rectangle((data['x_legend_picto']+31,y_legend,data['x_legend_picto']+61,y_legend+30), outline=color['red5t'], fill=color['red5t'])
      drawtl.text((data['x_legend_picto']+65,y_legend), "75-79", font=font_24, fill=color['black'])
      y_legend+=30
      drawtl.rectangle((data['x_legend_picto'],y_legend,data['x_legend_picto']+30,y_legend+30), outline=color['green4t'], fill=color['green4t'])
      drawtl.rectangle((data['x_legend_picto']+31,y_legend,data['x_legend_picto']+61,y_legend+30), outline=color['red4t'], fill=color['red4t'])
      drawtl.text((data['x_legend_picto']+65,y_legend), "70-74", font=font_24, fill=color['black'])
      y_legend+=30
      drawtl.rectangle((data['x_legend_picto'],y_legend,data['x_legend_picto']+30,y_legend+30), outline=color['green3t'], fill=color['green3t'])
      drawtl.rectangle((data['x_legend_picto']+31,y_legend,data['x_legend_picto']+61,y_legend+30), outline=color['red3t'], fill=color['red3t'])
      drawtl.text((data['x_legend_picto']+65,y_legend), "65-69", font=font_24, fill=color['black'])
      y_legend+=30
      drawtl.rectangle((data['x_legend_picto'],y_legend,data['x_legend_picto']+30,y_legend+30), outline=color['green2t'], fill=color['green2t'])
      drawtl.rectangle((data['x_legend_picto']+31,y_legend,data['x_legend_picto']+61,y_legend+30), outline=color['red2t'], fill=color['red2t'])
      drawtl.text((data['x_legend_picto']+65,y_legend), "60-64", font=font_24, fill=color['black'])
      y_legend+=30
      drawtl.rectangle((data['x_legend_picto'],y_legend,data['x_legend_picto']+30,y_legend+30), outline=color['green1t'], fill=color['green1t'])
      drawtl.rectangle((data['x_legend_picto']+31,y_legend,data['x_legend_picto']+61,y_legend+30), outline=color['red1t'], fill=color['red1t'])
      drawtl.text((data['x_legend_picto']+65,y_legend), "<60", font=font_24, fill=color['black'])


      y_legend+=35
      identity_threshold = 100-mismatch
      drawtl.text((data['x_legend_picto'],y_legend), "Minimum identity threshold : %i %%" % identity_threshold, font=font_24, fill=color['black'])
      drawtl.text((data['x_legend_picto'],y_legend+30), "Minimum block length : %i bp" % block_length, font=font_24, fill=color['black'])
      drawtl.text((data['x_legend_picto'],y_legend+60), "Transparency : %i" % alpha, font=font_24, fill=color['black'])
      drawtl.text((data['x_legend_picto'],y_legend+90), "Scale (pixel:bp) 1:%i" % scale, font=font_24, fill=color['black'])

      back.paste(ticklabel, mask=ticklabel)

      del drawtl
      file = "xmvconifer-" + alignment_file + "_m" + str(mismatch) + "_b" + str(block_length) + "_c" + str(scale) + "." + format
      print("Saving %s..." % file)
      back.save(open(file, 'wb'), formatdict[format])
      print("done.")
      return file

#---------------------------------------------
def main():
    opts, args = getopt.getopt(sys.argv[1:], "x:s:q:m:r:c:l:f:p:a:b:e:y:")

    (ref_gff_file, qry_gff_file, alignment_file, reference_file, query_file, format, fontpath)=(None,None,None,None,None,"png","")
    (mismatch, block_length, scale, protein, alpha, label)=(0,0,0,0,255,"")
    (reference, reflength)=([],[])
    formatdict = {'png':'PNG','gif':'GIF','tiff':'TIFF','jpeg':'JPEG'}

    for o, v in opts:
      if o == "-x":
        alignment_file=str(v)
      if o == "-s":
        reference_file=str(v)
      if o == "-q":
        query_file=str(v)
      if o == "-m": 
        mismatch=int(v)
      if o == "-b":
        block_length=int(v)
      if o == "-c":
        scale=int(v)
      if o == "-l":
        label=str(v)
      if o == "-f":
        format=str(v)
      if o == "-e":
        ref_gff_file=str(v)
      if o == "-y":
        qry_gff_file=str(v)
      if o == "-a":
        alpha = int(v)
      if o == "-p":
        fontpath=str(v)

    if (alignment_file == None or reference_file == None or query_file == None or mismatch == 0 or block_length == 0 or scale ==0):
      print("Usage: %s v1.2.4" % (sys.argv[0:]))
      print("-x alignment file (cross_match .rep or Pairwise mApping Format .paf) ")
      print("-s reference genome fasta file")
      print("-q query contig/genome fasta file")
      print("-e reference features (eg. exons) coordinates GFF tsv file (start end) - optional")
      print("-y query features (eg. exons) coordinates GFF tsv file (start end) - optional")
      print("-m maximum mismatch threshold (e.g. -m 10 allows representation of repeats having up to 10% mismatch")
      print("-b minimum length (bp) of similarity block to display")
      print("-c scale (pixel to basepair scale, for displaying the image)")
      print("-l label for the tree trunk (6 characters or less for best result)")
      print("-a alpha value, from 0 (transparent) to 255 (solid, default)")
      print("-f output image file format (png, tiff, jpeg, or gif) NOTE: png and tiff recommended.")
      print("-p full path to the directory with fonts on your system (please refer to the documentation for fonts used)")
      #print "-z transform bacterial ORF into protein (i.e. plot alignment between ORF products? 1/0) DEPRECATED\n";
      print("* Files for the -s and -q options must correspond to fasta files used to run cross_match")
      sys.exit(1)

    #====Graph Format
    if format not in formatdict:
      print("Not a valid Graph Format.  Please Select: png, tiff, jpeg, or gif. NOTE: png and tiff recommended.")
      sys.exit(1)

    #====Mismatch checks
    if (mismatch <0 or mismatch >99):
      print("-m must be a valid number between 0-99")
      sys.exit(1)

    #====Alpha checks
    if (alpha<0 or alpha >255):
      print("-a must be a valid number between 0-255")
      sys.exit(1)

    #===Scale checks
    if (scale<1):
      print("Not a possible scale. Make sure you select a number >1.")
      sys.exit(1)

    #====File checks
    checkFile(alignment_file)
    checkFile(reference_file)
    checkFile(query_file)

    ###OPTIONAL, FOR FEATURES REPRESENTATION
    (refgff,qrygff) = ({},{})
    if(ref_gff_file != None):
      checkFile(ref_gff_file)
      print("Reading reference feature file %s ..." % (ref_gff_file))
      refgff=readGFF(ref_gff_file,scale)
      print("done.")

    if(qry_gff_file != None):
      checkFile(qry_gff_file)
      print("Reading query feature file %s ..." % (qry_gff_file))
      qrygff=readGFF(qry_gff_file,scale)
      print("done.")

    #====Parse Fasta Files
    (order_ref, reference, reflength)=readFasta(reference_file, scale)
    (order_qry, query, qrylength)=readFasta(query_file, scale)

    #====Raise error if features out of bounds
    data=initGraph()
    if reflength / (data['width']-600) > scale:
       estscale = int(reflength / (data['width']-600)) + 1
       sys.exit("\n\n! The reference sequence is predicted to extend beyond the plot width, you must increase the scale to at least %i (we suggest rounding up to the next ten, hundred or thousand) -- fatal." % (estscale))
    if qrylength / (data['width']-600) > scale:
       estscale = int(qrylength / (data['width']-600)) + 1
       sys.exit("\n\n! The query sequence is predicted to extend beyond the plot width, you must increase the scale to at least %i (we suggest rounding up to the next ten, hundred or thousand) -- fatal." % (estscale))

    print("Reading alignment file %s ..." % (alignment_file))
    match = {}
    if alignment_file.endswith("rep"):
      match = readCrossMatch(alignment_file, mismatch, block_length, reference, query, scale)
    elif alignment_file.endswith("paf"):
      match = readPAF(alignment_file, mismatch, block_length, reference, query, scale)
    else:
      print("The alignment file provided (-x %s) does not end in .rep (cross_match) or .paf (PAF) -- fatal" % alignment_file)
      sys.exit(1)
    print("done.")
    print("done.")
    print("Drawing ...")
    drawRelationship(reference, order_ref, query, order_qry, match, scale, mismatch, block_length, alignment_file, reflength, format, formatdict, protein, label, alpha, refgff, qrygff, qrylength, fontpath)

#---------------------------------------------
#Main Call

main()
sys.exit(1)



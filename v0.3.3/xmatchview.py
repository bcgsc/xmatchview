#!/usr/bin/python
# xmatchview.py
# Visualizing genome synteny
# Rene L Warren 2005,2017,2018

import sys
import os
import getopt
import re
# import Image ## https://stackoverflow.com/questions/17451711/typeerror-when-resizing-an-image-with-pil-in-python
from PIL import Image
from PIL import ImageDraw
from PIL import ImageFont
from PIL import ImageEnhance
from PIL import PSDraw
import commands

#---------------------------------------------
def checkFile(file):

   print "Checking input %s" % file
   if not os.path.exists(file):
      print "File %s" % file + " is not valid"
      sys.exit(1)
   else:
      print "exists."

#---------------------------------------------
def readExon(exon_file,scale):
   exon = {}
   xon_obj=open(exon_file, 'r')
   (xstart,xend)=(0,0)
   for line in xon_obj:
       exonregex = re.compile("(\d+)\s+(\d+)(\s+)?(\S+)?")
       duo = exonregex.match(line)
       xstart = float(int(duo.group(1))/scale)
       xend = float(int(duo.group(2))/scale)
       color = duo.group(4)
       if color == None:
            color = "yellow"
       if not exon.has_key(xstart):
            exon[xstart] = {}
            if not exon[xstart].has_key('end'):
                exon[xstart]['end'] = {}
            if not exon[xstart].has_key('color'):
                exon[xstart]['color'] = {}

       exon[xstart]['end']=xend
       exon[xstart]['color']=color
       #print "%i,%i with %s" %(xstart,xend,color)

   xon_obj.close()

   return exon

#---------------------------------------------
def readPAF(paf_file,mismatch,block_length,reference,scale):

   (nocdt,match,query_hit)=({},{},{})

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
         #print "GR: %s" % line
         #print "REVERSE: %s %s %s %s %s %s %s %s" % (rm.group(1), rm.group(2), rm.group(3), rm.group(4), rm.group(5), rm.group(6), rm.group(7), rm.group(8))

         alignLen = float(rm.group(6)) - float(rm.group(5)) + 1
         percentMis = 100 * float(( alignLen - float(rm.group(7))) / alignLen ) 
         #print "=== %.2f === %.2f ===" % (alignLen,percentMis)
         #sys.exit(1)
         (primary_match, startFirstMatch, endFirstMatch, secondary_match, startSecondMatch, endSecondMatch)=(str(rm.group(4)), float(rm.group(5)), float(rm.group(6)), str(rm.group(1)), float(rm.group(3)), float(rm.group(2)))
         
         ####no autovivification in python
         if not nocdt.has_key(primary_match):
            nocdt[primary_match]={}
         if not nocdt[primary_match].has_key(secondary_match):
            nocdt[primary_match][secondary_match]={}
         if not nocdt[primary_match][secondary_match].has_key(startFirstMatch):
            nocdt[primary_match][secondary_match][startFirstMatch]={}
         if not nocdt[primary_match][secondary_match][startFirstMatch].has_key(endFirstMatch):
            nocdt[primary_match][secondary_match][startFirstMatch][endFirstMatch]={}
         if not nocdt[primary_match][secondary_match][startFirstMatch][endFirstMatch].has_key(startSecondMatch):
            nocdt[primary_match][secondary_match][startFirstMatch][endFirstMatch][startSecondMatch]={}
         if not nocdt[primary_match][secondary_match][startFirstMatch][endFirstMatch][startSecondMatch].has_key(endSecondMatch):
            nocdt[primary_match][secondary_match][startFirstMatch][endFirstMatch][startSecondMatch][endSecondMatch]={}

         nocdt[primary_match][secondary_match][startFirstMatch][endFirstMatch][startSecondMatch][endSecondMatch]=percentMis

         if percentMis > mismatch:
            continue  #will not display alignment lines below threshold
         elif (primary_match == secondary_match) and (startSecondMatch == startFirstMatch):
            break
         elif (endFirstMatch - startFirstMatch) < block_length:
            continue  #will skip smaller alignment
         else:
            if reference.has_key(primary_match):
               startFirstMatch=startFirstMatch/scale
               endFirstMatch=endFirstMatch/scale
               startSecondMatch=startSecondMatch/scale
               endSecondMatch=endSecondMatch/scale

               print "%i-%i   ::   %i-%i" % (startFirstMatch,endFirstMatch,startSecondMatch,endSecondMatch)

               if not match.has_key(primary_match):
                  match[primary_match]={}
               if not match[primary_match].has_key(secondary_match):
                  match[primary_match][secondary_match]={}
               if not match[primary_match][secondary_match].has_key(startFirstMatch):
                  match[primary_match][secondary_match][startFirstMatch]={}
               if not match[primary_match][secondary_match][startFirstMatch].has_key(endFirstMatch):
                  match[primary_match][secondary_match][startFirstMatch][endFirstMatch]={}
               if not match[primary_match][secondary_match][startFirstMatch][endFirstMatch].has_key(startSecondMatch):
                  match[primary_match][secondary_match][startFirstMatch][endFirstMatch][startSecondMatch]={}
               if not match[primary_match][secondary_match][startFirstMatch][endFirstMatch][startSecondMatch].has_key(endSecondMatch):
                  match[primary_match][secondary_match][startFirstMatch][endFirstMatch][startSecondMatch][endSecondMatch]={}

               match[primary_match][secondary_match][startFirstMatch][endFirstMatch][startSecondMatch][endSecondMatch]=percentMis

               if not query_hit.has_key(primary_match):
                  query_hit[secondary_match]=int(0)

               query_hit[secondary_match]=query_hit[secondary_match]+1

      ###forward matches
      elif fm != None:
         #print "GF: %s" % line
         #print "FORWARD: %s %s %s %s %s %s %s %s" % (fm.group(1), fm.group(2), fm.group(3), fm.group(4), fm.group(5), fm.group(6), fm.group(7), fm.group(8))

         alignLen = float(fm.group(6)) - float(fm.group(5)) + 1
         percentMis = 100 * float(( alignLen - float(fm.group(7))) / alignLen )
         #percentMis2 = 100 * float((float(fm.group(8)) - float(fm.group(7))) / float(fm.group(8)))
         #print "=== %.2f === %.2f === %.2f" % (alignLen,percentMis,percentMis2)
         #sys.exit(1)
         (primary_match, startFirstMatch, endFirstMatch, secondary_match, startSecondMatch, endSecondMatch)=(str(fm.group(4)), float(fm.group(5)), float(fm.group(6)), str(fm.group(1)), float(fm.group(2)), float(fm.group(3)))

         ####no autovivification in python
         if not nocdt.has_key(primary_match):
            nocdt[primary_match]={}
         if not nocdt[primary_match].has_key(secondary_match):
            nocdt[primary_match][secondary_match]={}
         if not nocdt[primary_match][secondary_match].has_key(startFirstMatch):
            nocdt[primary_match][secondary_match][startFirstMatch]={}
         if not nocdt[primary_match][secondary_match][startFirstMatch].has_key(endFirstMatch):
            nocdt[primary_match][secondary_match][startFirstMatch][endFirstMatch]={}
         if not nocdt[primary_match][secondary_match][startFirstMatch][endFirstMatch].has_key(startSecondMatch):
            nocdt[primary_match][secondary_match][startFirstMatch][endFirstMatch][startSecondMatch]={}
         if not nocdt[primary_match][secondary_match][startFirstMatch][endFirstMatch][startSecondMatch].has_key(endSecondMatch):
            nocdt[primary_match][secondary_match][startFirstMatch][endFirstMatch][startSecondMatch][endSecondMatch]={}

         nocdt[primary_match][secondary_match][startFirstMatch][endFirstMatch][startSecondMatch][endSecondMatch]=percentMis

         if reference.has_key(primary_match):
            startFirstMatch=startFirstMatch/scale
            endFirstMatch=endFirstMatch/scale
            startSecondMatch=startSecondMatch/scale
            endSecondMatch=endSecondMatch/scale

            print "%i-%i   ::   %i-%i" % (startFirstMatch,endFirstMatch,startSecondMatch,endSecondMatch)

            if not match.has_key(primary_match):
               match[primary_match]={}
            if not match[primary_match].has_key(secondary_match):
               match[primary_match][secondary_match]={}
            if not match[primary_match][secondary_match].has_key(startFirstMatch):
               match[primary_match][secondary_match][startFirstMatch]={}
            if not match[primary_match][secondary_match][startFirstMatch].has_key(endFirstMatch):
               match[primary_match][secondary_match][startFirstMatch][endFirstMatch]={}
            if not match[primary_match][secondary_match][startFirstMatch][endFirstMatch].has_key(startSecondMatch):
               match[primary_match][secondary_match][startFirstMatch][endFirstMatch][startSecondMatch]={}
            if not match[primary_match][secondary_match][startFirstMatch][endFirstMatch][startSecondMatch].has_key(endSecondMatch):
               match[primary_match][secondary_match][startFirstMatch][endFirstMatch][startSecondMatch][endSecondMatch]={}

            match[primary_match][secondary_match][startFirstMatch][endFirstMatch][startSecondMatch][endSecondMatch]=percentMis
            if not query_hit.has_key(primary_match):
               query_hit[secondary_match]=int(0)
            query_hit[secondary_match]=query_hit[secondary_match]+1

      #else:
         #print "NO RE:%s" % line
   xmatch_obj.close()

   return nocdt, match, query_hit

#---------------------------------------------
def readCrossMatch(crossmatch_file,mismatch,block_length,reference,scale):

   (nocdt,match,query_hit)=({},{},{})

   xmatch_obj=open(crossmatch_file, 'r')

   for line in xmatch_obj:
      ###reverse matches
      rev_regex = re.compile("(\s+)?\d+\s+(\S+)\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+(\S+)\s+\S+\s+C\s+(\S+)\s+\S+\s+(\S+)\s+(\S+)")
      rm = rev_regex.match(line)
      ###forward matches


      fwd_regex = re.compile("(\s+)?\d+\s+(\S+)\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+(\S+)\s+\S+\s+(\S+)\s+(\S+)\s+(\S+)\s+\S+")
      fm = fwd_regex.match(line)

      if rm != None:
         #print "GR: %s" % line
         #print "REVERSE: %s %s %s %s %s %s %s" % (rm.group(1), rm.group(2), rm.group(3), rm.group(4), rm.group(5), rm.group(6), rm.group(7))

         #(percentMis, primary_match, startFirstMatch, endFirstMatch, secondary_match, startSecondMatch, endSecondMatch)=(float(rm.group(1)), rm.group(2), float(rm.group(3)), float(rm.group(4)), rm.group(5), float(rm.group(6)), float(rm.group(7)))
         (percentMis, primary_match, startFirstMatch, endFirstMatch, secondary_match, startSecondMatch, endSecondMatch)=(float(rm.group(2)), rm.group(3), float(rm.group(4)), float(rm.group(5)), rm.group(6), float(rm.group(7)), float(rm.group(8)))
  

         ####no autovivification in python
         if not nocdt.has_key(primary_match):
            nocdt[primary_match]={}
         if not nocdt[primary_match].has_key(secondary_match):
            nocdt[primary_match][secondary_match]={}
         if not nocdt[primary_match][secondary_match].has_key(startFirstMatch):
            nocdt[primary_match][secondary_match][startFirstMatch]={}
         if not nocdt[primary_match][secondary_match][startFirstMatch].has_key(endFirstMatch):
            nocdt[primary_match][secondary_match][startFirstMatch][endFirstMatch]={}
         if not nocdt[primary_match][secondary_match][startFirstMatch][endFirstMatch].has_key(startSecondMatch):
            nocdt[primary_match][secondary_match][startFirstMatch][endFirstMatch][startSecondMatch]={}
         if not nocdt[primary_match][secondary_match][startFirstMatch][endFirstMatch][startSecondMatch].has_key(endSecondMatch):
            nocdt[primary_match][secondary_match][startFirstMatch][endFirstMatch][startSecondMatch][endSecondMatch]={}

         nocdt[primary_match][secondary_match][startFirstMatch][endFirstMatch][startSecondMatch][endSecondMatch]=percentMis

         if percentMis > mismatch:
            continue  #will not display alignment lines below threshold
         elif (primary_match == secondary_match) and (startSecondMatch == startFirstMatch):
            break
         elif (endFirstMatch - startFirstMatch) < block_length:
            continue  #will skip smaller alignment
         else:
            if reference.has_key(primary_match):
               startFirstMatch=startFirstMatch/scale
               endFirstMatch=endFirstMatch/scale
               startSecondMatch=startSecondMatch/scale
               endSecondMatch=endSecondMatch/scale

               print "%i-%i   ::   %i-%i" % (startFirstMatch,endFirstMatch,startSecondMatch,endSecondMatch)

               if not match.has_key(primary_match):
                  match[primary_match]={}
               if not match[primary_match].has_key(secondary_match):
                  match[primary_match][secondary_match]={}
               if not match[primary_match][secondary_match].has_key(startFirstMatch):
                  match[primary_match][secondary_match][startFirstMatch]={}
               if not match[primary_match][secondary_match][startFirstMatch].has_key(endFirstMatch):
                  match[primary_match][secondary_match][startFirstMatch][endFirstMatch]={}
               if not match[primary_match][secondary_match][startFirstMatch][endFirstMatch].has_key(startSecondMatch):
                  match[primary_match][secondary_match][startFirstMatch][endFirstMatch][startSecondMatch]={}
               if not match[primary_match][secondary_match][startFirstMatch][endFirstMatch][startSecondMatch].has_key(endSecondMatch):
                  match[primary_match][secondary_match][startFirstMatch][endFirstMatch][startSecondMatch][endSecondMatch]={}

               match[primary_match][secondary_match][startFirstMatch][endFirstMatch][startSecondMatch][endSecondMatch]=percentMis
               if not query_hit.has_key(primary_match):
                  query_hit[secondary_match]=int(0)

               query_hit[secondary_match]=query_hit[secondary_match]+1


      ###forward matches
      elif fm != None:
         #print "GF: %s" % line
         #print "FORWARD: %s %s %s %s %s %s %s" % (fm.group(1), fm.group(2), fm.group(3), fm.group(4), fm.group(5), fm.group(6), fm.group(7))
#         (percentMis, primary_match, startFirstMatch, endFirstMatch, secondary_match, startSecondMatch, endSecondMatch)=(float(fm.group(1)), fm.group(2), float(fm.group(3)), float(fm.group(4)), fm.group(5), float(fm.group(6)), float(fm.group(7)))

         (percentMis, primary_match, startFirstMatch, endFirstMatch, secondary_match, startSecondMatch, endSecondMatch)=(float(fm.group(2)), fm.group(3), float(fm.group(4)), float(fm.group(5)), fm.group(6), float(fm.group(7)), float(fm.group(8)))


         ####no autovivification in python
         if not nocdt.has_key(primary_match):
            nocdt[primary_match]={}
         if not nocdt[primary_match].has_key(secondary_match):
            nocdt[primary_match][secondary_match]={}
         if not nocdt[primary_match][secondary_match].has_key(startFirstMatch):
            nocdt[primary_match][secondary_match][startFirstMatch]={}
         if not nocdt[primary_match][secondary_match][startFirstMatch].has_key(endFirstMatch):
            nocdt[primary_match][secondary_match][startFirstMatch][endFirstMatch]={}
         if not nocdt[primary_match][secondary_match][startFirstMatch][endFirstMatch].has_key(startSecondMatch):
            nocdt[primary_match][secondary_match][startFirstMatch][endFirstMatch][startSecondMatch]={}
         if not nocdt[primary_match][secondary_match][startFirstMatch][endFirstMatch][startSecondMatch].has_key(endSecondMatch):
            nocdt[primary_match][secondary_match][startFirstMatch][endFirstMatch][startSecondMatch][endSecondMatch]={} 

         nocdt[primary_match][secondary_match][startFirstMatch][endFirstMatch][startSecondMatch][endSecondMatch]=percentMis

         if reference.has_key(primary_match):
            startFirstMatch=startFirstMatch/scale
            endFirstMatch=endFirstMatch/scale
            startSecondMatch=startSecondMatch/scale
            endSecondMatch=endSecondMatch/scale
              
            print "%i-%i   ::   %i-%i" % (startFirstMatch,endFirstMatch,startSecondMatch,endSecondMatch)
 
            if not match.has_key(primary_match):
               match[primary_match]={}
            if not match[primary_match].has_key(secondary_match):
               match[primary_match][secondary_match]={}
            if not match[primary_match][secondary_match].has_key(startFirstMatch):
               match[primary_match][secondary_match][startFirstMatch]={}
            if not match[primary_match][secondary_match][startFirstMatch].has_key(endFirstMatch):
               match[primary_match][secondary_match][startFirstMatch][endFirstMatch]={}
            if not match[primary_match][secondary_match][startFirstMatch][endFirstMatch].has_key(startSecondMatch):
               match[primary_match][secondary_match][startFirstMatch][endFirstMatch][startSecondMatch]={}
            if not match[primary_match][secondary_match][startFirstMatch][endFirstMatch][startSecondMatch].has_key(endSecondMatch):
               match[primary_match][secondary_match][startFirstMatch][endFirstMatch][startSecondMatch][endSecondMatch]={}

            match[primary_match][secondary_match][startFirstMatch][endFirstMatch][startSecondMatch][endSecondMatch]=percentMis
            if not query_hit.has_key(primary_match):
               query_hit[secondary_match]=int(0)
            query_hit[secondary_match]=query_hit[secondary_match]+1

      #else:
         #print "NO RE:%s" % line      
   xmatch_obj.close()

   return nocdt, match, query_hit

#---------------------------------------------
def generateCoords(nocdt, size, leap, protein):
 
   freq={}

   pos_range=range(0,size,leap)

   for pos in pos_range:
      print "%i out of %i" % (pos,size)
      for reference in nocdt:
         for comparison in nocdt[reference]:
            start1_dict=nocdt[reference][comparison].keys()
            start1_dict.sort()
            for start1 in start1_dict:
               end1_dict=nocdt[reference][comparison][start1].keys()
               end1_dict.sort()
               for end1 in end1_dict:
                  start2_dict=nocdt[reference][comparison][start1][end1].keys()
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
                        end2_dict=nocdt[reference][comparison][start1][end1][start2].keys()
                        end2_dict.sort()
                        for end2 in end2_dict:
                           current_mismatch=float(nocdt[reference][comparison][start1][end1][start2][end2])
                           if not freq.has_key(pos):
                              freq[pos]={}
                           if not freq[pos].has_key(current_mismatch):
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
   npos={}

   file_obj = open(file, 'r')
   
   for line in file_obj:
      head_match_regex = re.compile('>(\S+)')
      head_match = head_match_regex.match(line)
      if head_match != None:
          if (head_match != previous_contig and previous_contig != None):
             (seq_length, scale)=(int(seq_length), int(scale))
             L1[previous_contig] = float(seq_length/scale)
             seq_length = 0                                        #resets the sequence length
          previous_contig = head_match.group(1)
           
      seq_subset_regex = re.compile('(.*)', re.I)
      seq_subset = seq_subset_regex.match(line)
      if seq_subset != None:
         seq_length += len(seq_subset.group(1))
         npos=findOccurences(seq_subset.group(1).upper(), "N")
         zpos=findOccurences(seq_subset.group(1).upper(), "Z")

   (seq_length, scale)=(int(seq_length), int(scale))
   L1[previous_contig] = float(seq_length/scale)                                #for the last sequence

   file_obj.close()

   print "scaled down %s =%f total=%i " % (previous_contig, L1[previous_contig], seq_length)

   return (L1, seq_length, npos, previous_contig, zpos)

#---------------------------------------------
def initColor(alpha):
   color={}

   #allocate colors
   color["white"] = (255,255,255,255)
   color["black"] = (0,0,0,255)
   color["swamp"] = (150,150,30,255)
   color["blue"] = (0,102,204,255)
   color["yellow"] = (255,255,0,255)
   color["cyan"] = (0,255,255,255)
   color["purple"] = (255,0,255,alpha)
   color["green"] = (100,250,25,255)
   color["lime"] = (57,255,20,255)
   color["red"] = (250,25,75,255)
   color["sarin"] = (255,66,0,255)
   color["forest"] = (25,175,0,255)
   color["dirtyred"] = (200,0,120,255)
   color["navy"] = (0,0,150,alpha)
   color["dirtyyellow"] = (200,200,75,255)
   color["grey"] = (153,153,153,255)
   color["lightgrey"] = (220,220,220,355)
   color["salmon"] = (255,153,153,alpha)
   color["lightblue"] = (153,204,255,alpha)
   color["orange"] = (255,153,51,255)
   color["beige"] = (222,184,135,255)

   return color

#---------------------------------------------
def initGraph():
   data={} 
  
   #default data points
   data['width']=2400
   data['height']=1200
   data['ref_y']=250
   data['mis_bar']=50
   data['query_y']=70
   data['x']=50
   data['xlabel']=110
   data['bar_thick']=20
   data['query_thick']=15
   data['reference_thick']=15
   data['x_legend']=600
   data['y_legend']=750
   data['x_legend_picto']=100
   data['tick_up']=25
   data['tick_down']=40

   return data

#---------------------------------------------
def drawRectangle(draw,start,end,y,thickness,bar_color,text,font,text_color):
   
   draw.rectangle((start,y,end,y+thickness), bar_color)
   #draw.text((start-80, y-2), text, font=font, fill=text_color)###position of SEQUENCE label
   draw.text((end+5, y-2), text, font=font, fill=text_color)###position of SEQUENCE label

#---------------------------------------------
def plotFrequency(freq,size,scale,draw,color,data,leap):

   pos_range=range(0,size,leap)
   
   for pos in pos_range:
      if freq.has_key(pos):
         freq_list=freq[pos]
         previous=data['mis_bar']
         identity_range=range(99,-1,-1)## RESTRICT SI AXIS was 9
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
               color_now="orange"
            elif cumul>=5:
               color_now="dirtyred"
            #elif cumul==6:
            #   color_now="salmon"
            #elif cumul==7:
            #   color_now="orange"
            #elif cumul>=8:
            #   color_now="yellow"

            extension=((200-(2*id))+data['mis_bar'])   #RESTRICT SI AXIS y was 20
            compressed=(pos/scale)+data['x']           #x
            
            if color_now != "white":
               #print "%i, %i, %i, %i %s" % (compressed,previous,compressed,extension,color_now)
               draw.line((compressed,previous,compressed,extension),color[color_now])
 
            previous = extension


#---------------------------------------------
def drawRelationship(reference_list, query_list, match_list, scale, query_hit, mismatch, block_length, alignment_file, freq, reflength, leap, format, formatdict, protein, alpha, refexon, qryexon, qrylength, refnpos, qrynpos, refname, qryname, refzpos, qryzpos, fontpath):

      scaled_reflength=int(reflength/scale)

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

      if os.path.exists(arialfont): ### Will check for truetype first, they look better
          ###Set Font (truetype)
          font_18=ImageFont.truetype(fontpath + "/arial.ttf",18)
          font_20=ImageFont.truetype(fontpath + "/arial.ttf",20)
          fontb_20=ImageFont.truetype(fontpath + "/arialbd.ttf",20)
          fontbi_20=ImageFont.truetype(fontpath + "/arialbi.ttf",20)
          fontb_22=ImageFont.truetype(fontpath + "/arialbd.ttf",22)
          font_24=ImageFont.truetype(fontpath + "/arial.ttf",24)
          fontb_24=ImageFont.truetype(fontpath + "/arialbd.ttf",24)
      elif os.path.exists(pilfont): ### Will settle for PIL font, if ttf do not exist. Otherwise, sticking with default.
          ###Set font (pil) (sizes are limited, made to be compatible with TT fonts)
          font_18=ImageFont.load_path(fontpath + "/helvR14.pil")
          font_20=ImageFont.load_path(fontpath + "/helvR18.pil")
          fontb_20=ImageFont.load_path(fontpath + "/helvB18.pil")
          fontbi_20=ImageFont.load_path(fontpath + "/helvBO18.pil")
          fontb_22=ImageFont.load_path(fontpath + "/helvB24.pil")
          font_24=ImageFont.load_path(fontpath + "/helvR24.pil")
          fontb_24=ImageFont.load_path(fontpath + "/helvB24.pil")

      ###Define Image
      back = Image.new("RGBA", (data['width'],data['height']),(0,0,0,0))
      bdraw = ImageDraw.Draw(back)

      poly = Image.new("RGBA", (data['width'],data['height']))
      draw = ImageDraw.Draw(poly)

      ###Draw Legend
      date=commands.getstatusoutput("date")

      ###Picto Legend
      bdraw.text((data['x_legend_picto']+280,data['y_legend']), "Legend", font=fontb_24, fill=color['black'])
      y_legend = data['y_legend']+30
      bdraw.text((data['x_legend_picto'],y_legend), "Frequency Repeated", font=fontbi_20, fill=color['black'])

      ####
      bdraw.text((data['x_legend'],y_legend), "Mismatch threshold : %i %%" % mismatch, font=font_20, fill=color['black'])
      bdraw.text((data['x_legend'],y_legend+20), "Minimum Block Length : %i bp" % block_length, font=font_20, fill=color['black'])
      bdraw.text((data['x_legend'],y_legend+40), "Scale (pixel:bp)  1:%i" % scale, font=font_20, fill=color['black'])
      #bdraw.text((data['x_legend'],y_legend+60), "%s" % date[1], font=font_20, fill=color['black'])
      ####

      y_legend+=25
      bdraw.rectangle((data['x_legend_picto'],y_legend,data['x_legend_picto']+20,y_legend+20), outline=color['black'], fill=color['blue'])
      bdraw.text((data['x_legend_picto']+25,y_legend), "Single copy", font=font_20, fill=color['black'])
      y_legend+=25
      bdraw.rectangle((data['x_legend_picto'],y_legend,data['x_legend_picto']+20,y_legend+20), outline=color['black'], fill=color['cyan'])
      bdraw.text((data['x_legend_picto']+25,y_legend), "2X", font=font_20, fill=color['black'])
      y_legend+=25
      bdraw.rectangle((data['x_legend_picto'],y_legend,data['x_legend_picto']+20,y_legend+20), outline=color['black'], fill=color['green'])
      bdraw.text((data['x_legend_picto']+25,y_legend), "3X", font=font_20, fill=color['black'])
      y_legend+=25
      bdraw.rectangle((data['x_legend_picto'],y_legend,data['x_legend_picto']+20,y_legend+20), outline=color['black'], fill=color['orange'])
      bdraw.text((data['x_legend_picto']+25,y_legend), "4X", font=font_20, fill=color['black'])
      y_legend+=25
      bdraw.rectangle((data['x_legend_picto'],y_legend,data['x_legend_picto']+20,y_legend+20), outline=color['black'], fill=color['dirtyred'])
      bdraw.text((data['x_legend_picto']+25,y_legend), "5X and over", font=font_20, fill=color['black'])
      y_legend+=25
      #bdraw.rectangle((data['x_legend_picto'],y_legend,data['x_legend_picto']+20,y_legend+20), outline=color['black'], fill=color['salmon'])
      #bdraw.text((data['x_legend_picto']+25,y_legend), "6X", font=font_20, fill=color['black'])
      #y_legend+=25
      #bdraw.rectangle((data['x_legend_picto'],y_legend,data['x_legend_picto']+20,y_legend+20), outline=color['black'], fill=color['orange'])
      #bdraw.text((data['x_legend_picto']+25,y_legend), "7X", font=font_20, fill=color['black'])
      #y_legend+=25
      #bdraw.rectangle((data['x_legend_picto'],y_legend,data['x_legend_picto']+20,y_legend+20), outline=color['black'], fill=color['yellow'])
      #bdraw.text((data['x_legend_picto']+25,y_legend), "8X and over", font=font_20, fill=color['black'])
      y_legend+=40

      bdraw.text((data['x_legend_picto'],y_legend), "Collinear Blocks", font=fontbi_20, fill=color['black'])
      y_legend+=30

      bdraw.polygon((data['x_legend_picto']-5,y_legend,data['x_legend_picto'],y_legend+25,data['x_legend_picto']+25,y_legend+25,data['x_legend_picto']+20,y_legend), outline=color['navy'], fill=color['lightblue'])
      bdraw.text((data['x_legend_picto']+30,y_legend), "Direct", font=font_20, fill=color['black'])

      y_legend+=30
      bdraw.polygon((data['x_legend_picto']-5,y_legend,data['x_legend_picto']+25,y_legend+25,data['x_legend_picto']-5,y_legend+25,data['x_legend_picto']+25,y_legend), outline=color['purple'], fill=color['salmon'])
      bdraw.text((data['x_legend_picto']+30,y_legend), "Inverted", font=font_20, fill=color['black'])

      y_legend+=40

      bdraw.text((data['x_legend_picto'],y_legend), "Other", font=fontbi_20, fill=color['black'])
      y_legend+=30

      bdraw.rectangle((data['x_legend_picto']-5,y_legend+5,data['x_legend_picto']+25,y_legend+7), fill=color['red'])
      bdraw.text((data['x_legend_picto']+30,y_legend), "Mismatch threshold", font=font_20, fill=color['black'])
      y_legend+=30

      bdraw.rectangle((data['x_legend_picto']-5,y_legend,data['x_legend_picto']+25,y_legend+(data['reference_thick']/2)), outline=color['black'], fill=color['yellow'])
      bdraw.text((data['x_legend_picto']+30,y_legend), "Sequence features", font=font_20, fill=color['black'])
      y_legend+=30
      

      bdraw.rectangle((data['x_legend_picto']+5,y_legend,data['x_legend_picto']+10,y_legend+data['reference_thick']), outline=color['red'], fill=color['red'])
      bdraw.text((data['x_legend_picto']+30,y_legend), "Ambiguous bases (Ns)", font=font_20, fill=color['black'])
      y_legend+=30

      ####
      for ref in reference_list:
         init_coord=int(data['x'])
         last_coord=int(data['x']+reference_list[ref])

         ### draw top rectangle for REFERENCE
         drawRectangle(bdraw, init_coord, last_coord,data['ref_y'],data['reference_thick'],color['black'],ref,fontb_20,color['black'])
         x_range=range(init_coord, last_coord, 100)

         ### draw kbp scale
         if reflength >= 10000:
            for position in x_range:
               base_number=int(((position-data['x'])*scale)/1000)
               bdraw.rectangle((position,data['tick_up'],position+2,data['tick_down']),color['black'])
               bdraw.text((position-5, data['tick_up']-25), "%i" % base_number, font=font_20, fill=color['black'])
         else:
            for position in x_range:
               base_number=(position-data['x']) * scale
               base_number=float(base_number)
               base_number=base_number/1000
               bdraw.rectangle((position,data['tick_up'],position+2,data['tick_down']),color['black'])
               #print "%i %i %i >>> %.2f <<< %i,%i" % (data['x'],position,scale,base_number,data['x_legend_picto'],position)
               bdraw.text((position-5, data['tick_up']-25), "%.1f" % base_number, font=font_24, fill=color['black'])
  
      bdraw.text((data['x']+scaled_reflength+60,data['tick_up']-25), "kbp", font=fontb_24, fill=color['black'])
      ###Mismatch Axis
      #identity=int(0)
      identity = int(0) ### RESTRICT SI AXIS was 90
      grid_range=range(data['mis_bar'], data['ref_y'], 20)

      ###Draw grid
      for grid in grid_range:
         bdraw.rectangle((data['x'],grid,data['x']+scaled_reflength+5,grid+2),color['lightgrey'])
         bdraw.text((data['x']+scaled_reflength+10, grid-7), "%i " % identity, font=font_18, fill=color['black'])
         identity += 10 ### RESTRICT SI AXIS was 1

      ###Draw grid metric
      bdraw.text((data['x']+scaled_reflength+60, 150), "% Identity", font=font_18, fill=color['black'])

      ###Draw Repeat Frequency
      plotFrequency(freq,reflength,scale,bdraw,color,data,leap)

      ###Draw Threshold
      threshold_line= data['mis_bar'] + (200-(2*mismatch))
      draw.rectangle((data['x'],threshold_line,data['x']+scaled_reflength+5,threshold_line+2), color['red'])

      ###Draw Query & Collinear blocks 
      (decay, current_position, LCB)=(350, data['x'], 10)
   
      for match in match_list:
         allhit=match_list[match]
         for hit in allhit:
            start1_list=allhit[hit]
            stop=current_position + query_list[hit]
            if match != hit:
               drawRectangle(bdraw,current_position,stop,data['ref_y']+decay,data['query_thick'],color['black'], hit, fontb_20, color['black'])
            s1_list_sort=start1_list.keys()
            s1_list_sort.sort()
            for start1 in s1_list_sort:
               end1_list=start1_list[start1]
               e1_list_sort=end1_list.keys()
               e1_list_sort.sort()
               for end1 in e1_list_sort:
                  start2_list=end1_list[end1]
                  s2_list_sort=start2_list.keys()
                  s2_list_sort.sort()
                  for start2 in s2_list_sort:
                     end2_list=start2_list[start2]
                     e2_list_sort=end2_list.keys()
                     e2_list_sort.sort()
                     for end2 in e2_list_sort:

                        if start2 > end2:
                           outline_color="purple"
                           fill_color="salmon"
                        else:
                           outline_color="navy"
                           fill_color="lightblue"
                        ###draw ORF on upper
                        size_ref = end1 - start1
                        size_qry = end2 - start2
                        buf_ref = ((size_ref - (size_ref/3)) / 2)
                        buf_qry = ((size_qry - (size_qry/3)) / 2)
                        ss1 = start1 + buf_ref
                        ee1 = end1 - buf_ref
                        ss2 = start2 + buf_qry
                        ee2 = end2 - buf_qry
                        print "%s (%i-%i) hits %s  ::  mismatch %.2f target(%i)  block %i target (%i) " % (match,start1,end1,hit,end2_list[end2],mismatch,size_ref,block_length)

                        if match == hit:### COMPARE 1 SEQUENCE AGAINST ITSELF
                           if start1 <= start2:

                              repeat_size = start2 - start1
                              size_chunk = int(decay * repeat_size / scaled_reflength)
                              print "%i %i %i" % (size_chunk, repeat_size, scaled_reflength)
                              size_chunk += 50

                              if protein:
                                 if end2_list[end2] <= mismatch: ###does it pass the mismatch cutoff?
                                     ss1 = int(ss1)
                                     ss2 = int(ss2)
                                     bdraw.rectangle((data['x']+ss1,data['ref_y']+data['reference_thick']+7,data['x']+ee1,data['ref_y']+data['reference_thick']+17), outline=color["black"], fill=color["red"])###REF1
                                     bdraw.rectangle((data['x']+ss2,data['ref_y']+data['reference_thick']+7,data['x']+ee2,data['ref_y']+data['reference_thick']+17), outline=color["black"], fill=color["red"])###REF2
                                     draw.arc((data['x']+ss1,data['ref_y']+data['reference_thick']+15-size_chunk,data['x']+ss2,data['ref_y']+data['reference_thick']+17+size_chunk),360,180, color[outline_color])
                                     back.paste(poly, mask=poly)
                                     del draw
                                     poly = Image.new("RGBA", (data['width'],data['height']))
                                     draw = ImageDraw.Draw(poly)
                              else:
                                 if end2_list[end2] <= mismatch: ###does it pass the mismatch cutoff?
                                     start1 = int(start1)
                                     start2 = int(start2)
                                     draw.rectangle((data['x']+start1,data['ref_y'],data['x']+end1,data['ref_y']+data['reference_thick']), outline=color[outline_color], fill=color[fill_color])###REF LEFT REPEAT COLINEAR BLOCKS
                                     draw.rectangle((data['x']+start2,data['ref_y'],data['x']+end2,data['ref_y']+data['reference_thick']), outline=color[outline_color], fill=color[fill_color])###REF RIGHT REPEAT COLINEAR BLOCKS
                                     draw.arc((data['x']+start1,data['ref_y']+data['reference_thick']-size_chunk,data['x']+start2,data['ref_y']+data['reference_thick']+size_chunk),360,180, color[outline_color])###DRAW ARC AT REPEAT EDGE ONLY 
                                     back.paste(poly, mask=poly)
                                     del draw
                                     poly = Image.new("RGBA", (data['width'],data['height']))
                                     draw = ImageDraw.Draw(poly)
                                 
                        else:###COMPARE 2 SEQS

                           if protein:
                           
                              if end2_list[end2] <= mismatch: ###does it pass the mismatch cutoff?
                                 draw.polygon((data['x']+ss1,data['ref_y']+data['reference_thick']+17,data['x']+ss2,data['ref_y']+decay-17,data['x']+ee2,data['ref_y']+decay-17,data['x']+ee1,data['ref_y']+data['reference_thick']+17), outline=color[outline_color], fill=color[fill_color])
                                 draw.rectangle((data['x']+ss1,data['ref_y']+data['reference_thick']+7,data['x']+ee1,data['ref_y']+data['reference_thick']+17), outline=color["black"], fill=color["red"]) ###REF
                                 draw.rectangle((data['x']+ss2,data['ref_y']+decay-17,data['x']+ee2,data['ref_y']+decay-7), outline=color["black"], fill=color["red"])###QRY
                                 back.paste(poly, mask=poly)
                                 del draw
                                 poly = Image.new("RGBA", (data['width'],data['height']))
                                 draw = ImageDraw.Draw(poly)
                           else:
                              if end2_list[end2] <= mismatch: ###does it pass the mismatch cutoff?
                                 draw.polygon((data['x']+start1,data['ref_y']+data['reference_thick'],data['x']+start2,data['ref_y']+decay,data['x']+end2,data['ref_y']+decay,data['x']+end1,data['ref_y']+data['reference_thick']), outline=color[outline_color], fill=color[fill_color])
                                 draw.rectangle((data['x']+start1,data['ref_y'],data['x']+end1,data['ref_y']+data['reference_thick']), outline=color[outline_color], fill=color[fill_color]) ###REF COLINEAR BLOCKS
                                 draw.rectangle((data['x']+start2,data['ref_y']+decay,data['x']+end2,data['ref_y']+decay+data['reference_thick']), outline=color[outline_color], fill=color[fill_color])###QRY COLINEAR BLOCKS
                                 back.paste(poly, mask=poly)
                                 del draw
                                 poly = Image.new("RGBA", (data['width'],data['height']))
                                 draw = ImageDraw.Draw(poly)


      #enhancer = ImageEnhance.Sharpness(im)
      #for i in range(8):
      #   factor = i / 4.0
      #   enhancer.enhance(factor).show("Sharpness %f" % factor)

      ####draw features (eg. exons) on side of ref
      ###REF gene model
      
      for scaledexstart in refexon:
          #print "ex start: %i" % exstart
          exstart = data['x'] + scaledexstart
          exend = data['x'] + refexon[scaledexstart]['end'] 
          draw.rectangle((exstart,data['ref_y'],exend,data['ref_y']+(data['reference_thick']/2)),outline=color[refexon[scaledexstart]['color']], fill=color[refexon[scaledexstart]['color']])###features/exons
      if refname != qryname:
          for scaledexstart in qryexon:
              exstart = data['x'] + scaledexstart
              exend = data['x'] + qryexon[scaledexstart]['end']
              draw.rectangle((exstart,data['ref_y']+decay+(data['reference_thick']/2)+1,exend,data['ref_y']+decay+data['reference_thick']),outline=color[qryexon[scaledexstart]['color']], fill=color[qryexon[scaledexstart]['color']])###features/exons

      ### draw start position of Ns
      for nstart in refnpos:
          nstart = data['x'] + (nstart/scale)
          draw.line((nstart,data['ref_y']+1,nstart,data['ref_y']+data['reference_thick']-1),color['red'],width=1)
          #print "N: %i" % nstart
      #print ">>>%s  %s" % (refname,qryname)
      if refname != qryname:
          for nstart in qrynpos:
              nstart = data['x'] + (nstart/scale)
              draw.line((nstart,data['ref_y']+decay+1,nstart,data['ref_y']+decay+data['reference_thick']-1),color['red'],width=1)
              #print "N: %i" % nstart

      ### draw start position of sequences, delimited by Zs
      for nstart in refzpos:
          nstart = data['x'] + (nstart/scale)
          draw.line((nstart,data['ref_y']+1,nstart,data['ref_y']+data['reference_thick']-1),color['lime'],width=1)
          #print "N: %i" % nstart
      #print ">>>%s  %s" % (refname,qryname)
      if refname != qryname:
          for nstart in qryzpos:
              nstart = data['x'] + (nstart/scale)
              draw.line((nstart,data['ref_y']+decay+1,nstart,data['ref_y']+decay+data['reference_thick']-1),color['lime'],width=1)
              #print "N: %i" % nstart        

      back.paste(poly, mask=poly)
      file = "xmv-" + alignment_file + "_m" + str(mismatch) + "_b" + str(block_length) + "_r" + str(leap) + "_c" + str(scale) + "." + format
      print "Saving %s..." % file
      back.save(open(file, 'wb'), formatdict[format])
      print "done."
      return file

#---------------------------------------------
def main():
   
   opts, args = getopt.getopt(sys.argv[1:], "x:s:q:m:r:c:b:f:p:e:y:a:")

   (ref_exon_file, qry_exon_file, alignment_file, reference_file, query_file, format, fontpath)=(None,None,None,None,None,"png","")
   (mismatch, block_length, scale, leap, protein, alpha)=(0,0,0,0,0,255)
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
      if o == "-r":
        leap=int(v)
      if o == "-f":
        format=str(v)
      if o == "-e":
        ref_exon_file=str(v)
      if o == "-y":
        qry_exon_file=str(v)
      if o == "-a":
        alpha = int(v)
      if o == "-p":
        fontpath=str(v)

   if (alignment_file == None or reference_file == None or query_file == None or mismatch == 0 or block_length == 0 or scale ==0 or leap == 0):
      print "Usage: %s v0.3.3" % (sys.argv[0:])
      print "-x alignment file (cross_match .rep or Pairwise mApping Format .paf) "
      print "-s reference genome fasta file"
      print "-q query contig/genome fasta file"
      print "-e reference features (eg. exons) coordinates tsv file (start end) - optional"
      print "-y query features (eg. exons) coordinates tsv file (start end) - optional"
      print "-m mismatch threshold (e.g. -m 10 allows representation of repeats having up to 10% mismatch"
      print "-b length (bp) of similarity block to display"
      print "-c scale (pixel to basepair scale, for displaying the image)"
      print "-r leap (bp) to evaluate repeat frequency (smaller numbers will increase the resolution, but will affect drastically the run time.  recommended -l=50)"
      print "-a alpha value, from 0 (transparent) to 255 (solid, default)"
      print "-f output image file format (png, tiff, jpeg, or gif) NOTE: the png and tiff are better."
      print "-p full path to the directory with fonts on your system (please refer to the documentation for fonts used)"
      #print "-z transform bacterial ORF into protein (i.e. plot alignment between ORF products? 1/0) -not fully tested-";
      print "* Files for the -s and -q options must correspond to fasta files used to run cross_match"
      sys.exit(1)

   #====Graph Format
   if not formatdict.has_key(format):
      print "Not a valid Graph Format.  Please Select: bmp, jpeg, png, ps, gif, or tiff"
      sys.exit(1)

   #====Mismatch checks
   if (mismatch <0 or mismatch >99):
      print "-m must be a valid number between 0-99"
      sys.exit(1)

   #===Scale checks
   if (scale<1):
      print "Not a possible scale. Make sure you select a number >1."
      sys.exit(1)

   #====Alpha checks
   if (alpha<0 or alpha >255):
     print "-a must be a valid number between 0-255"
     sys.exit(1)

   #====File checks
   checkFile(alignment_file)
   checkFile(reference_file)
   checkFile(query_file)

   ###OPTIONAL, FOR FEATURES/EXON REPRESENTATION
   (refexon,qryexon) = ({},{})
   if(ref_exon_file != None and qry_exon_file != None):
      checkFile(ref_exon_file)
      checkFile(qry_exon_file)
      print "Reading reference exon file..."
      refexon=readExon(ref_exon_file,scale)
      print "Reading query exon file..."
      qryexon=readExon(qry_exon_file,scale)
      print "done."

   #====Parse Fasta Files
   (refnpos,qrynpos) = ({},{})
   (refname, qryname) = (None, None)
   (reference, reflength, refnpos, refname, refzpos)=readFasta(reference_file, scale)
   (query, qrylength, qrynpos, qryname, qryzpos)=readFasta(query_file, scale)

   print "Reading alignment file..."
   (nocdt, match, query_hit) = ({},{},{})
   if alignment_file.endswith("rep"):
      (nocdt, match, query_hit)=readCrossMatch(alignment_file, mismatch, block_length, reference, scale)
   elif alignment_file.endswith("paf"):
      (nocdt, match, query_hit)=readPAF(alignment_file, mismatch, block_length, reference, scale)
   else:
      print "The alignment file provided (-x %s) does not end in .rep (cross_match) or .paf (PAF) -- fatal" % alignment_file
      sys.exit(1)

   print "done."
   print "Computing Repeat frequencies..."
   (freq)=generateCoords(nocdt, reflength, leap, protein)
   print "done."
   print "Drawing repeats..."
   drawRelationship(reference, query, match, scale, query_hit, mismatch, block_length, alignment_file, freq, reflength, leap, format, formatdict, protein, alpha, refexon, qryexon, qrylength, refnpos, qrynpos, refname, qryname, refzpos, qryzpos, fontpath)

#---------------------------------------------
#Main Call

main()
sys.exit(1)



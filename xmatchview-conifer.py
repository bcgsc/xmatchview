#!/home/rwarren/python/Python-2.4.2/bin/python2.4
#XMatchView, spruceView
#Rene Warren 2005,2015,2017

import sys
import os
import getopt
import re
import Image
import ImageDraw
import ImageFont
import ImageEnhance
import PSDraw
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
            color = "black"
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
         #print "REVERSE: %s %s %s %s %s %s %s" % (fm.group(1), fm.group(2), fm.group(3), fm.group(4), fm.group(5), fm.group(6), fm.group(7))

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
   npos = {} 

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

   (seq_length, scale)=(int(seq_length), int(scale))
   L1[previous_contig] = float(seq_length/scale)                                #for the last sequence

   file_obj.close()

   print "scaled down %s =%f total=%i " % (previous_contig, L1[previous_contig], seq_length)

   return (L1, seq_length, npos)

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
   color["purple"] = (255,0,255,255)
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
   data['width']=2000
   data['height']=2000
   data['ref_y']=900
   data['skew']=300 ### CHANGE THIS FOR THE PITCH OF THE TREE
   data['decay']=120### THIS IS THE SPACE BETWEEN BOTH SIDES, TOP OF TREE
   data['ref_y_skew']=data['ref_y']-data['skew'] ###DON'T CHANGE THIS
   data['mis_bar']=50
   data['query_y']=70
   data['x']=100
   data['xlabel']=110
   data['bar_thick']=20
   data['query_thick']=15
   data['reference_thick']=15
   data['x_legend']=450
   data['y_legend']=1500
   data['x_legend_picto']=1500
   data['tick_up']=data['ref_y_skew'] - 120
   data['tick_down']=data['tick_up'] + 20

   return data

#---------------------------------------------
def drawRectangle(draw,start,end,y,thickness,bar_color,text,font,text_color):
   
   draw.rectangle((start,y,end,y+thickness), bar_color)
   draw.text((start-80, y), text, font=font, fill=text_color)

#---------------------------------------------
def plotFrequency(freq,size,scale,draw,color,data,leap):

   pos_range=range(0,size,leap)
   
   for pos in pos_range:
      if freq.has_key(pos):
         freq_list=freq[pos]
         previous=data['mis_bar']
         identity_range=range(9,-1,-1)
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
def drawRelationship(reference_list, query_list, match_list, scale, query_hit, mismatch, block_length, crossmatch_file, freq, reflength, leap, format, formatdict, protein, label, alpha, refexon, qryexon, qrylength, refnpos, qrynpos):

      scaled_reflength=float(reflength/scale)
      scaled_qrylength=float(qrylength/scale)

      ###Capture last coordinates of relationships
      (u2max,v2max,x2max,y2max)=(0,0,0,0)

      ###Initialize new graph
      data=initGraph()

      ###Get colors
      color=initColor(alpha)
 
      ###Set Font

      ###Set Font
      #medium_font=ImageFont.load_path("/home/rwarren/fonts/pil/helvB12.pil")
      arial_18=ImageFont.truetype("/home/rwarren/fonts/truetype/arial.ttf",18)
      arialb_18=ImageFont.truetype("/home/rwarren/fonts/truetype/arialbd.ttf",18)
      arial_20=ImageFont.truetype("/home/rwarren/fonts/truetype/arial.ttf",20)
      arialb_20=ImageFont.truetype("/home/rwarren/fonts/truetype/arialbd.ttf",20)
      ariali_20=ImageFont.truetype("/home/rwarren/fonts/truetype/ariali.ttf",20)
      arialbi_20=ImageFont.truetype("/home/rwarren/fonts/truetype/arialbi.ttf",20)
      arial_22=ImageFont.truetype("/home/rwarren/fonts/truetype/arial.ttf",22)
      arial_24=ImageFont.truetype("/home/rwarren/fonts/truetype/arial.ttf",24)
      arialb_24=ImageFont.truetype("/home/rwarren/fonts/truetype/arialbd.ttf",24)
      arialbi_24=ImageFont.truetype("/home/rwarren/fonts/truetype/arialbi.ttf",24)
      arialbi_28=ImageFont.truetype("/home/rwarren/fonts/truetype/arialbi.ttf",28)
      arial_28=ImageFont.truetype("/home/rwarren/fonts/truetype/arial.ttf",28) 
      arialb_28=ImageFont.truetype("/home/rwarren/fonts/truetype/arialbd.ttf",28)
      arialb_22=ImageFont.truetype("/home/rwarren/fonts/truetype/arialbd.ttf",22)
      arialb_92=ImageFont.truetype("/home/rwarren/fonts/truetype/arialbd.ttf",92)
      courb_92=ImageFont.truetype("/home/rwarren/fonts/truetype/courbd.ttf",92)

      ###Define Image
      back = Image.new("RGBA", (data['width'],data['height']),(0,0,0,0))
      bdraw = ImageDraw.Draw(back)

      poly = Image.new("RGBA", (data['width'],data['height']))
      draw = ImageDraw.Draw(poly)

      ticklabel = Image.new("RGBA", (data['width'],data['height']))

      ###Draw Legend
      date=commands.getstatusoutput("date")

      ###Picto Legend
      y_legend = data['y_legend']+30
      bdraw.text((data['x_legend_picto'],y_legend), "Sequence identity (%)", font=arialbi_28, fill=color['black'])
      ####

      y_legend+=35
      bdraw.rectangle((data['x_legend_picto'],y_legend,data['x_legend_picto']+30,y_legend+30), outline=color['black'], fill=color['green10t'])
      bdraw.text((data['x_legend_picto']+35,y_legend), "99-100", font=arial_24, fill=color['black'])
      y_legend+=30
      bdraw.rectangle((data['x_legend_picto'],y_legend,data['x_legend_picto']+30,y_legend+30), outline=color['black'], fill=color['green9t'])
      bdraw.text((data['x_legend_picto']+35,y_legend), "95-97", font=arial_24, fill=color['black'])
      y_legend+=30
      bdraw.rectangle((data['x_legend_picto'],y_legend,data['x_legend_picto']+30,y_legend+30), outline=color['black'], fill=color['green8t'])
      bdraw.text((data['x_legend_picto']+35,y_legend), "90-94", font=arial_24, fill=color['black'])
      y_legend+=30
      bdraw.rectangle((data['x_legend_picto'],y_legend,data['x_legend_picto']+30,y_legend+30), outline=color['black'], fill=color['green7t'])
      bdraw.text((data['x_legend_picto']+35,y_legend), "85-89", font=arial_24, fill=color['black'])
      y_legend+=30
      bdraw.rectangle((data['x_legend_picto'],y_legend,data['x_legend_picto']+30,y_legend+30), outline=color['black'], fill=color['green6t'])
      bdraw.text((data['x_legend_picto']+35,y_legend), "80-84", font=arial_24, fill=color['black'])
      y_legend+=30
      bdraw.rectangle((data['x_legend_picto'],y_legend,data['x_legend_picto']+30,y_legend+30), outline=color['black'], fill=color['green5t'])
      bdraw.text((data['x_legend_picto']+35,y_legend), "75-79", font=arial_24, fill=color['black'])
      y_legend+=30
      bdraw.rectangle((data['x_legend_picto'],y_legend,data['x_legend_picto']+30,y_legend+30), outline=color['black'], fill=color['green4t'])
      bdraw.text((data['x_legend_picto']+35,y_legend), "70-74", font=arial_24, fill=color['black'])
      y_legend+=30
      bdraw.rectangle((data['x_legend_picto'],y_legend,data['x_legend_picto']+30,y_legend+30), outline=color['black'], fill=color['green3t'])
      bdraw.text((data['x_legend_picto']+35,y_legend), "65-69", font=arial_24, fill=color['black'])
      y_legend+=30
      bdraw.rectangle((data['x_legend_picto'],y_legend,data['x_legend_picto']+30,y_legend+30), outline=color['black'], fill=color['green2t'])
      bdraw.text((data['x_legend_picto']+35,y_legend), "60-64", font=arial_24, fill=color['black'])
      y_legend+=30
      bdraw.rectangle((data['x_legend_picto'],y_legend,data['x_legend_picto']+30,y_legend+30), outline=color['black'], fill=color['green1t'])
      bdraw.text((data['x_legend_picto']+35,y_legend), "0-59", font=arial_24, fill=color['black'])

      y_legend+=35
      identity_threshold = 100-mismatch
      bdraw.text((data['x_legend_picto'],y_legend), "Minimum identity threshold : %i %%" % identity_threshold, font=arial_24, fill=color['black'])
      bdraw.text((data['x_legend_picto'],y_legend+30), "Minimum block length : %i bp" % block_length, font=arial_24, fill=color['black'])
      bdraw.text((data['x_legend_picto'],y_legend+60), "Transparency : %i" % alpha, font=arial_24, fill=color['black'])
      bdraw.text((data['x_legend_picto'],y_legend+90), "Scale (pixel:bp) 1:%i" % scale, font=arial_24, fill=color['black'])

      back = back.rotate(90)

      decay = data['decay']


      ####draw features/exons on side of ref
      ###REF gene model
    
      x1ref = data['x']
      x2ref = data['x'] + scaled_reflength
      y1ref = data['ref_y']
      y2ref = data['ref_y']-data['skew']

      mrref = (y2ref - y1ref ) / (x2ref - x1ref)
      brref = y2ref - (mrref * x2ref)
      for exstart in refexon:
          exend = refexon[exstart]['end']
          a1 = data['x'] + exstart
          a2 = data['x'] + exend
          b1 = (mrref * a1 ) + brref
          b2 = (mrref * a2 ) + brref
          draw.polygon((a1,b1-11,a1,b1,a2,b2,a2,b2-11),outline=color['black'], fill=color[refexon[exstart]['color']])###features/exons

      ###QRY gene model
      x1qry = data['x']
      x2qry = data['x'] + scaled_qrylength
      y1qry = data['ref_y']+decay
      y2qry = data['ref_y']+decay+data['skew']

      mqqry = (y2qry - y1qry ) / (x2qry - x1qry)
      bqqry = y2qry - (mqqry * x2qry)

      for exstart in qryexon:
          exend = qryexon[exstart]['end']
          a1 = data['x'] + exstart
          a2 = data['x'] + exend
          b1 = (mqqry * a1 ) + bqqry
          b2 = (mqqry * a2 ) + bqqry
          draw.polygon((a1,b1+(data['query_thick'])+2,a1,b1+(data['query_thick'])+11,a2,b2+(data['query_thick'])+11,a2,b2+(data['query_thick'])+2),outline=color['black'], fill=color[qryexon[exstart]['color']])###features/exons

      ####REFERENCE
      for ref in reference_list:
         init_coord=data['x']
         last_coord=data['x']+reference_list[ref]

         x1ref = init_coord
         y1ref = data['ref_y']
         x2ref = last_coord
         y2ref = data['ref_y_skew']

         mref = (y2ref - y1ref) / (x2ref - x1ref) 
         bref = y2ref - (mref * x2ref)

         print "REFERENCE x1=%i y1=%i x2=%i y2=%i M=%.2f  B=%.2f  " % (x1ref,y1ref,x2ref,y2ref,mref,bref)
         draw.polygon((init_coord-1,data['ref_y']-1,init_coord-1,data['ref_y']+data['reference_thick'],last_coord+1,data['ref_y_skew']+data['reference_thick'],last_coord+1,data['ref_y_skew']-3), outline=color['brown'], fill=color['brown'])### reference rectangle (top)
         draw.text((last_coord+5, data['ref_y_skew']-7), ref, font=arialb_24, fill=color['green9'])###label for ref

         back.paste(poly, mask=poly)
         del draw
         poly = Image.new("RGBA", (data['width'],data['height']))
         draw = ImageDraw.Draw(poly)

         x_range=range(int(init_coord), int(last_coord), 100)

         for position in x_range:
            draw.rectangle((position,data['tick_up'],position+2,data['tick_down']),color['black'])

      ###Mismatch Axis
      #identity=int(0)
      #identity = int(90)
      #grid_range=range(data['mis_bar'], data['ref_y'], 20)

      #for grid in grid_range:
      #   draw.rectangle((data['x'],grid,data['x']+scaled_reflength+5,grid+2),color['lightgrey'])
      #   draw.text((data['x']+scaled_reflength+10, grid-7), "%i " % identity, font=arial_18, fill=color['black'])
      #   identity += 1

      #draw.text((data['x']+scaled_reflength+60, 150), "% Identity", font=arial_18, fill=color['black'])

      ###Draw Repeat Frequency
      #plotFrequency(freq,reflength,scale,draw,color,data,leap)
      #threshold_line= data['mis_bar'] + (200-(2*mismatch))
      #draw.rectangle((data['x'],threshold_line,data['x']+scaled_reflength+5,threshold_line+2), color['red'])

      (current_position, LCB, skew,stop)=(data['x'], 10, data['skew'],data['x'])
       
      ###Draw Query & Collinear blocks 

      for match in match_list:
         allhit=match_list[match]
         for hit in allhit:
            start1_list=allhit[hit]
            stop=current_position + query_list[hit]
            if match != hit: 
               draw.polygon((current_position-1,data['ref_y']+decay-1,current_position-1,data['ref_y']+data['query_thick']+decay+1,stop+1,data['ref_y']+decay+skew+data['query_thick']+2,stop+1,data['ref_y']+decay+skew), outline=color['brown'], fill=color['brown'])### Query rectangle (bottom)
               draw.text((stop+5, data['ref_y']+decay+skew-3), hit, font=arialb_24, fill=color['green9'])### label for query 
               back.paste(poly, mask=poly)
               del draw
               poly = Image.new("RGBA", (data['width'],data['height']))
               draw = ImageDraw.Draw(poly)

               x1qry = current_position
               x2qry = stop
               y1qry = data['ref_y']+decay
               y2qry = data['ref_y']+decay+skew

               mqry = (y2qry - y1qry ) / (x2qry - x1qry)
               bqry = y2qry - (mqry * x2qry)

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
 
                        seqid = 100 - end2_list[end2]
                        print "si=%.2f mis=%.2f" % (seqid,end2_list[end2])
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
                        #print "FILL= %s" % (fill_color)
                        outline_color = "beige" 
                        #if start2 > end2:
                        #   outline_color="purple"
                        #   fill_color="salmon"
                   

                        ###draw ORF on upper
                        #draw.rectangle((data['x']+start1,data['ref_y']+1,data['x']+end1,data['ref_y']+data['reference_thick']-1), outline=color["lightgrey"], fill=color["lightgrey"])
                        size_ref = end1 - start1
                        size_qry = end2 - start2
                        buf_ref = ((size_ref - (size_ref/3)) / 2)
                        buf_qry = ((size_qry - (size_qry/3)) / 2)
                        ss1 = start1 + buf_ref
                        ee1 = end1 - buf_ref
                        ss2 = start2 + buf_qry
                        ee2 = end2 - buf_qry
                        print "%s (%i-%i) hits %s  ::  mismatch %.2f target(%i)  block %i target (%i) " % (match,start1,end1,hit,end2_list[end2],mismatch,size_ref,block_length)

                        if match == hit:
                           if start1 <= start2:

                              repeat_size = start2 - start1
                              size_chunk = int(decay * repeat_size / scaled_reflength)
                              #print "%i %i %i" % (size_chunk, repeat_size, scaled_reflength)
                              size_chunk += 50

                              if protein:
                                 draw.rectangle((data['x']+ss1,data['ref_y']+data['reference_thick']+7,data['x']+ee1,data['ref_y']+data['reference_thick']+17), outline=color["black"], fill=color["red"]) 
                                 if end2_list[end2] <= mismatch: ###does it pass the mismatch cutoff?
                                    draw.arc((data['x']+ss1,data['ref_y']+data['reference_thick']+15-size_chunk,data['x']+ss2,data['ref_y']+data['reference_thick']+17+size_chunk),360,180, color[outline_color])
                              else:
                                 if end2_list[end2] <= mismatch: ###does it pass the mismatch cutoff?
                                    draw.arc((data['x']+start1,data['ref_y']+data['reference_thick']-size_chunk,data['x']+start2,data['ref_y']+data['reference_thick']+size_chunk),360,180, color[outline_color]) 
                        else:
#                           draw.rectangle((data['x']+start2,data['ref_y']+decay,data['x']+end2,data['ref_y']+decay+data['reference_thick']), outline=color["black"], fill=color["lightgrey"])

                           if protein:
                              draw.rectangle((data['x']+ss1,data['ref_y']+data['reference_thick']+7,data['x']+ee1,data['ref_y']+data['reference_thick']+17), outline=color["black"], fill=color["red"])
                              draw.rectangle((data['x']+ss2,data['ref_y']+decay-17,data['x']+ee2,data['ref_y']+decay-7), outline=color["black"], fill=color["red"])
                           
                              if end2_list[end2] <= mismatch: ###does it pass the mismatch cutoff?
                                 draw.polygon((data['x']+ss1,data['ref_y']+data['reference_thick']+17,data['x']+ss2,data['ref_y']+decay-17,data['x']+ee2,data['ref_y']+decay-17,data['x']+ee1,data['ref_y']+data['reference_thick']+17), outline=color[outline_color], fill=color[fill_color])
                                 draw.rectangle((data['x']+start1,data['ref_y']+1,data['x']+end1,data['ref_y']+data['reference_thick']-1), outline=color[outline_color], fill=color[fill_color])
                           else:
                              if end2_list[end2] <= mismatch: ###does it pass the mismatch cutoff?
                                 x1 = data['x'] + start1
                                 x2 = data['x'] + end1
                                 y1 = (mref * x1 ) + bref
                                 y2 = (mref * x2 ) + bref
                                 #print "x1=%i y1=%i x2=%i y2=%i M=%.2f  B=%.2f  " % (x1,y1,x2,y2,mref,bref) 
                                 u1 = data['x'] + start2
                                 u2 = data['x'] + end2
                                 v1 = (mqry * u1 ) + bqry
                                 v2 = (mqry * u2 ) + bqry
                                 #print "u1=%i v1=%i u2=%i v2=%i M=%.2f  B=%.2f  " % (u1,v1,u2,v2,mqry,bqry)

                                 if x2 > x2max:
                                     u2max = u2
                                     v2max = v2
                                     x2max = x2
                                     y2max = y2

                                 ### LINES
                                 draw.polygon((x1,y1+data['reference_thick'],u1,v1,u2,v2,x2,y2+data['reference_thick']), outline=color[outline_color], fill=color[fill_color])
                                 ### REPEAT FEATURE
                                 draw.polygon((x1,y1,x1,y1+data['reference_thick'],x2,y2+data['reference_thick'],x2,y2),outline=color[outline_color], fill=color[fill_color])###colinear block on reference
                                 #back.paste(poly, mask=poly)
                                 draw.polygon((u1,v1,u2,v2,u2,v2+data['reference_thick'],u1,v1+data['reference_thick']),outline=color[outline_color], fill=color[fill_color])###colinear block on query
                                 back.paste(poly, mask=poly)
                                 del draw
                                 poly = Image.new("RGBA", (data['width'],data['height']))
                                 draw = ImageDraw.Draw(poly)


      #enhancer = ImageEnhance.Sharpness(im)
      #for i in range(8):
      #   factor = i / 4.0
      #   enhancer.enhance(factor).show("Sharpness %f" % factor)

      ###getFileName
      #xm_regex = re.compile('(\S+)\.\S+')
      #xm_name = xm_regex.match(crossmatch_file)
      #file = xm_name.group(1) + "_m" + str(mismatch) + "_b" + str(block_length) + "_l" + str(leap) + "_s" + str(scale) + "." + format
            ### draw start position of Ns
      for nstart in refnpos:
          nstart = data['x'] + (nstart/scale)
          ny = (mrref * nstart) + brref
          draw.line((nstart,ny,nstart,ny+data['reference_thick']-1),color['red'],width=2)
          #print "N: %i" % nstart
      for nstart in qrynpos:
          nstart = data['x'] + (nstart/scale)
          ny = (mqqry * nstart) + bqqry
          draw.line((nstart,ny+2,nstart,ny+data['reference_thick']+1),color['red'],width=2)
          #print "N: %i" % nstart

      ### calculate placement of tree trunk/label
      charwidth = 65 ### approximate character width in pixels
      labellength = len(label)
      totaltrunklength = (labellength + 2) * charwidth ### the 2 is for a one-character buffer before/after
      mtrunk = (y2max - v2max ) / (x2max - u2max)
      btrunk = y2max - (mtrunk * x2max)
      ytrunk = y1ref
      if u2max > x2max:
          ytrunk = y1ref + decay

      xtrunk = (ytrunk - btrunk) / mtrunk
      #print "%i %i %i %i %.2f %.2f x1=%.2f y1=%.2f x2=%.2f+420 LAST=%.2f" % (u1,u2,v1,v2,mtrunk,btrunk,xtrunk,y1ref,xtrunk,last_coord)
      draw.rectangle((xtrunk+5,y1ref+data['reference_thick']+4,xtrunk+totaltrunklength,y1ref+decay-5), outline=color['brown'], fill=color['brown'])###trunk
      draw.text((xtrunk+charwidth,data['ref_y']+17), label, font=arialb_92, fill=color['beige'])###label
      back.paste(poly, mask=poly)
      ### end trunk code
      ### rotate plot to be able to place scale
      back = back.rotate(270)
      del draw
      drawtl = ImageDraw.Draw(ticklabel)
      ###final tick labels
      x_range=range(data['x'], int(last_coord), 100)

      if reflength >= 10000:
         for position in x_range:
            base_number=int(((position-data['x'])*scale)/1000)
            drawtl.text((data['x_legend_picto']+25,position-15), "%i" % base_number, font=arial_28, fill=color['black'])
      else:
         for position in x_range:
            base_number=(position-data['x']) * scale
            base_number=float(base_number)
            base_number=base_number/1000
            #print "%i %i %i >>> %.2f <<< %i,%i" % (data['x'],position,scale,base_number,data['x_legend_picto'],position)
            drawtl.text((data['x_legend_picto']+25,position-15), "%.1f" % base_number, font=arial_28, fill=color['black'])

      drawtl.text((data['x_legend_picto']+25,last_coord+25), "kbp", font=arialb_28, fill=color['black'])
      back.paste(ticklabel, mask=ticklabel)
      del drawtl
      file = crossmatch_file + "_m" + str(mismatch) + "_b" + str(block_length) + "_l" + str(leap) + "_s" + str(scale) + "." + format
      print "Saving %s..." % file
      back.save(open(file, 'wb'), formatdict[format])
      print "done."
      return file

#---------------------------------------------
def main():
    opts, args = getopt.getopt(sys.argv[1:], "x:s:q:m:r:c:l:f:p:a:b:e:y:")

    (ref_exon_file,qry_exon_file,crossmatch_file, reference_file, query_file, format)=(None,None,None,None,None,None)
    (mismatch, block_length, scale, leap, protein, alpha)=(0,0,0,0,0,255)
    (reference, reflength)=([],[])
   #formatdict = {'PNG':'png','GIF':'gif','TIFF':'tiff','BMP':'bmp','JPEG':'jpeg','EPS':'ps'}
    formatdict = {'png':'PNG','gif':'GIF','tiff':'TIFF','bmp':'BMP','jpeg':'JPEG','ps':'EPS'}

    for o, v in opts:
      if o == "-x":
        crossmatch_file=str(v)
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
        protein=int(v)

    if (crossmatch_file == None or reference_file == None or query_file == None or mismatch == 0 or block_length == 0 or scale ==0 or leap == 0):
      print "Usage: %s v0.1" % (sys.argv[0:])
      print "-x crossmatch file"
      print "-s reference genome fasta file"
      print "-q query contig/genome fasta file"
      print "-e reference features (eg. exons) coordinates tsv file (start end) - optional"
      print "-y query features (eg. exons) coordinates tsv file (start end) - optional"
      print "-m maximum mismatch threshold (e.g. -m 10 allows representation of repeats having up to 10% mismatch"
      print "-b minimum length (bp) of similarity block to display"
      print "-c scale (pixel to basepair scale, for displaying the image)"
      print "-r basepair length leap to evaluate repeat frequency (smaller numbers will increase the resolution, but will affect drastically the run time.  recommended -r=50)"
      print "-l label for the tree trunk (6 characters or less for best result)"
      print "-a alpha value, from 0 (transparent) to 255 (solid, default)"
      print "-f file format (bmp, jpeg, png, ps, gif, pdf, tiff) NOTE: the png, ps, tiff and bmp are much better."
      #print "-p transform bacterial ORF into protein (i.e. plot alignment between ORF products? 1/0)\n";
      print "* Files for the -s and -q options must correspond to fasta files used to run cross_match"
      sys.exit(1)

   #====Graph Format
    if not formatdict.has_key(format):
      print "Not a valid Graph Format.  Please Select: bmp, jpeg, png, ps, gif, pdf or tiff"
      sys.exit(1)

   #====Mismatch checks
    if (mismatch <0 or mismatch >99):
      print "-m must be a valid number between 0-99"
      sys.exit(1)

   #====Alpha checks
    if (alpha<0 or alpha >255):
      print "-a must be a valid number between 0-255"
      sys.exit(1)

   #===Scale checks
    if (scale<1):
      print "Not a possible scale. Make sure you select a number >1."
      sys.exit(1)

   #====File checks
    checkFile(crossmatch_file)
    checkFile(reference_file)
    checkFile(query_file)

    ###OPTIONAL, FOR FEATURES/EXON REPRESENTATION
    (refexon,qryexon) = ({},{})
    if(ref_exon_file != None and qry_exon_file != None):
        checkFile(ref_exon_file)
        checkFile(qry_exon_file)
        print "Reading reference exon file..."
        (refexon)=readExon(ref_exon_file,scale)
        print "Reading query exon file..."
        (qryexon)=readExon(qry_exon_file,scale)
        print "done."

    #====Parse Fasta Files
    (refnpos,qrynpos) = ({},{})
    (reference, reflength, refnpos)=readFasta(reference_file, scale)
    (query, qrylength, qrynpos)=readFasta(query_file, scale)

    print "Reading Crossmatch file..."
    (nocdt, match, query_hit)=readCrossMatch(crossmatch_file, mismatch, block_length, reference, scale)
    print "done."
    print "Computing Repeat frequencies..."
    (freq)=generateCoords(nocdt, reflength, leap, protein)
    print "done."
    print "Drawing repeats..."
    drawRelationship(reference, query, match, scale, query_hit, mismatch, block_length, crossmatch_file, freq, reflength, leap, format, formatdict, protein, label, alpha, refexon, qryexon, qrylength, refnpos, qrynpos)

#---------------------------------------------
#Main Call

main()
sys.exit(1)



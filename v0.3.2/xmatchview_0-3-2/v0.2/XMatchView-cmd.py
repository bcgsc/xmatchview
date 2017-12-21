#!/usr/local/python/2.3.3/bin/python2.3
#Rene Warren 2005

import sys
import os
sys.path.insert(0,'/home/rwarren/python/Development/SeqDev/lib')
sys.path.insert(0,'/home/rwarren/python/standard/Imaging-1.1.4/PIL')

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
def readFasta(file, scale):

   (head_match, previous_contig,seq_length) = (None,None,0)
   L1={}
 

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
    
   (seq_length, scale)=(int(seq_length), int(scale))
   L1[previous_contig] = float(seq_length/scale)                                #for the last sequence

   file_obj.close()


   print "scaled down %s =%f total=%i " % (previous_contig, L1[previous_contig], seq_length)

   return (L1, seq_length)

#---------------------------------------------
def initColor():
   color={}

   #allocate colors
   color["white"] = (255,255,255)
   color["black"] = (0,0,0)
   color["swamp"] = (150,150,30)
   color["blue"] = (0,102,204)
   color["yellow"] = (255,255,0)
   color["cyan"] = (0,255,255)
   color["purple"] = (255,0,255)
   color["green"] = (100,250,25)
   color["red"] = (250,25,75)
   color["forrest"] = (25,175,0)
   color["dirtyred"] = (200,0,120)
   color["navy"] = (0,0,150)
   color["dirtyyellow"] = (200,200,75)
   color["grey"] = (153,153,153)
   color["lightgrey"] = (220,220,220)
   color["salmon"] = (255,153,153)
   color["lightblue"] = (153,204,255)
   color["orange"] = (255,153,51)
   color["beige"] = (222,184,135)

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
   data['x']=100
   data['xlabel']=110
   data['bar_thick']=20
   data['query_thick']=15
   data['reference_thick']=15
   data['x_legend']=600
   data['y_legend']=750
   data['x_legend_picto']=100
   data['thick_up']=25
   data['thick_down']=40

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
         identity_range=range(99,-1,-1)
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

            extension=((200-(2*id))+data['mis_bar'])   #y
            compressed=(pos/scale)+data['x']           #x
            
            if color_now != "white":
               #print "%i, %i, %i, %i %s" % (compressed,previous,compressed,extension,color_now)
               draw.line((compressed,previous,compressed,extension),color[color_now])
 
            previous = extension


#---------------------------------------------
def drawRelationship(reference_list, query_list, match_list, scale, query_hit, mismatch, block_length, crossmatch_file, freq, reflength, leap, format, formatdict, protein):

      scaled_reflength=int(reflength/scale)

      ###Initialize new graph
      data=initGraph()

      ###Get colors
      color=initColor()
 
      ###Set Font

      ###Set Font
      #medium_font=ImageFont.load_path("/home/rwarren/fonts/pil/helvB12.pil")
      arial_18=ImageFont.truetype("/home/rwarren/fonts/truetype/arial.ttf",18)
      arialb_18=ImageFont.truetype("/home/rwarren/fonts/truetype/arialbd.ttf",18)
      arial_20=ImageFont.truetype("/home/rwarren/fonts/truetype/arial.ttf",20)
      arialb_20=ImageFont.truetype("/home/rwarren/fonts/truetype/arialbd.ttf",20)
      ariali_20=ImageFont.truetype("/home/rwarren/fonts/truetype/ariali.ttf",20)
      arialbi_20=ImageFont.truetype("/home/rwarren/fonts/truetype/arialbi.ttf",20)
      arialb_22=ImageFont.truetype("/home/rwarren/fonts/truetype/arialbd.ttf",22)

      ###Define Image
      im = Image.new("RGB", (data['width'],data['height']),color['white'])
      draw = ImageDraw.Draw(im)

      ###Draw Legend
      date=commands.getstatusoutput("date")

      ###Picto Legend
      draw.text((data['x_legend_picto']+50,data['y_legend']), "Legend", font=arialb_22, fill=color['black'])
      y_legend = data['y_legend']+30
      draw.text((data['x_legend_picto'],y_legend), "Frequency Repeated", font=arialbi_20, fill=color['black'])

      ####
      draw.text((data['x_legend'],y_legend), "Mismatch threshold %i" % mismatch, font=arial_20, fill=color['black'])
      draw.text((data['x_legend'],y_legend+20), "Minimum Block Length=%i" % block_length, font=arial_20, fill=color['black'])
      draw.text((data['x_legend'],y_legend+40), "Scale=1:%i" % scale, font=arial_20, fill=color['black'])
      draw.text((data['x_legend'],y_legend+60), "%s" % date[1], font=arial_20, fill=color['black'])
      draw.text((data['x_legend'],y_legend+80), "rwarren@bcgsc.ca", font=arial_20, fill=color['black'])
      ####

      y_legend+=25
      draw.rectangle((data['x_legend_picto'],y_legend,data['x_legend_picto']+20,y_legend+20), outline=color['black'], fill=color['blue'])
      draw.text((data['x_legend_picto']+25,y_legend), "1X", font=arial_20, fill=color['black'])
      y_legend+=25
      draw.rectangle((data['x_legend_picto'],y_legend,data['x_legend_picto']+20,y_legend+20), outline=color['black'], fill=color['cyan'])
      draw.text((data['x_legend_picto']+25,y_legend), "2X", font=arial_20, fill=color['black'])
      y_legend+=25
      draw.rectangle((data['x_legend_picto'],y_legend,data['x_legend_picto']+20,y_legend+20), outline=color['black'], fill=color['green'])
      draw.text((data['x_legend_picto']+25,y_legend), "3X", font=arial_20, fill=color['black'])
      y_legend+=25
      draw.rectangle((data['x_legend_picto'],y_legend,data['x_legend_picto']+20,y_legend+20), outline=color['black'], fill=color['dirtyred'])
      draw.text((data['x_legend_picto']+25,y_legend), "4X", font=arial_20, fill=color['black'])
      y_legend+=25
      draw.rectangle((data['x_legend_picto'],y_legend,data['x_legend_picto']+20,y_legend+20), outline=color['black'], fill=color['purple'])
      draw.text((data['x_legend_picto']+25,y_legend), "5X", font=arial_20, fill=color['black'])
      y_legend+=25
      draw.rectangle((data['x_legend_picto'],y_legend,data['x_legend_picto']+20,y_legend+20), outline=color['black'], fill=color['salmon'])
      draw.text((data['x_legend_picto']+25,y_legend), "6X", font=arial_20, fill=color['black'])
      y_legend+=25
      draw.rectangle((data['x_legend_picto'],y_legend,data['x_legend_picto']+20,y_legend+20), outline=color['black'], fill=color['orange'])
      draw.text((data['x_legend_picto']+25,y_legend), "7X", font=arial_20, fill=color['black'])
      y_legend+=25
      draw.rectangle((data['x_legend_picto'],y_legend,data['x_legend_picto']+20,y_legend+20), outline=color['black'], fill=color['yellow'])
      draw.text((data['x_legend_picto']+25,y_legend), "8X and over", font=arial_20, fill=color['black'])
      y_legend+=40

      draw.text((data['x_legend_picto'],y_legend), "Collinear Blocks", font=arialbi_20, fill=color['black'])
      y_legend+=30

      draw.polygon((data['x_legend_picto']-5,y_legend,data['x_legend_picto'],y_legend+25,data['x_legend_picto']+25,y_legend+25,data['x_legend_picto']+20,y_legend), outline=color['navy'], fill=color['lightblue'])
      draw.text((data['x_legend_picto']+30,y_legend), "Direct", font=arial_20, fill=color['black'])

      y_legend+=30
      draw.polygon((data['x_legend_picto']-5,y_legend,data['x_legend_picto']+25,y_legend+25,data['x_legend_picto']-5,y_legend+25,data['x_legend_picto']+25,y_legend), outline=color['purple'], fill=color['salmon'])
      draw.text((data['x_legend_picto']+30,y_legend), "Inverted", font=arial_20, fill=color['black'])

      y_legend+=40

      draw.text((data['x_legend_picto'],y_legend), "Other", font=arialbi_20, fill=color['black'])
      y_legend+=30

      draw.rectangle((data['x_legend_picto']-5,y_legend+5,data['x_legend_picto']+25,y_legend+7), fill=color['red'])
      draw.text((data['x_legend_picto']+30,y_legend), "Mismatch threshold", font=arial_20, fill=color['black'])

      ####
      for ref in reference_list:
         init_coord=int(data['x'])
         last_coord=int(data['x']+reference_list[ref])

         drawRectangle(draw,init_coord, last_coord,data['ref_y'],data['reference_thick'],color['black'],ref,arialb_18,color['black'])
         x_range=range(init_coord, last_coord, 100)

         for position in x_range:
            draw.rectangle((position,data['thick_up'],position+2,data['thick_down']),color['black'])
            base_number=int(((position-data['x'])*scale)/1000)
            draw.text((position-10, data['thick_up']-25), "%i kb" % base_number, font=arial_18, fill=color['black'])

      ###Mismatch Axis
      identity=int(0)
      grid_range=range(data['mis_bar'], data['ref_y'], 20)

      for grid in grid_range:
         draw.rectangle((data['x'],grid,data['x']+scaled_reflength+5,grid+2),color['lightgrey'])
         draw.text((data['x']+scaled_reflength+10, grid-7), "%i " % identity, font=arial_18, fill=color['black'])
         identity += 10

      draw.text((data['x']+scaled_reflength+60, 150), "% Identity", font=arial_18, fill=color['black'])

      ###Draw Repeat Frequency
      plotFrequency(freq,reflength,scale,draw,color,data,leap)

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
               drawRectangle(draw,current_position,stop,data['ref_y']+decay,data['query_thick'],color['black'], hit, arialb_18, color['black'])
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
                        outline_color="forrest"
                        fill_color="lightblue"

                        if start2 > end2:
                           outline_color="purple"
                           fill_color="salmon"
                        else:
                           outline_color="navy"
                           fill_color="lightblue"
                        ###draw ORF on upper
                        draw.rectangle((data['x']+start1,data['ref_y'],data['x']+end1,data['ref_y']+data['reference_thick']), outline=color["black"], fill=color["lightgrey"])
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
                              print "%i %i %i" % (size_chunk, repeat_size, scaled_reflength)
                              size_chunk += 50

                              if protein:
                                 draw.rectangle((data['x']+ss1,data['ref_y']+data['reference_thick']+7,data['x']+ee1,data['ref_y']+data['reference_thick']+17), outline=color["black"], fill=color["red"]) 
                                 if end2_list[end2] <= mismatch: ###does it pass the mismatch cutoff?
                                    draw.arc((data['x']+ss1,data['ref_y']+data['reference_thick']+15-size_chunk,data['x']+ss2,data['ref_y']+data['reference_thick']+17+size_chunk),360,180, color[outline_color])
                              else:
                                 if end2_list[end2] <= mismatch: ###does it pass the mismatch cutoff?
                                    draw.arc((data['x']+start1,data['ref_y']+data['reference_thick']-size_chunk,data['x']+start2,data['ref_y']+data['reference_thick']+size_chunk),360,180, color[outline_color]) 
                        else:
                           draw.rectangle((data['x']+start2,data['ref_y']+decay,data['x']+end2,data['ref_y']+decay+data['reference_thick']), outline=color["black"], fill=color["lightgrey"])

                           if protein:
                              draw.rectangle((data['x']+ss1,data['ref_y']+data['reference_thick']+7,data['x']+ee1,data['ref_y']+data['reference_thick']+17), outline=color["black"], fill=color["red"])
                              draw.rectangle((data['x']+ss2,data['ref_y']+decay-17,data['x']+ee2,data['ref_y']+decay-7), outline=color["black"], fill=color["red"])
                           
                              if end2_list[end2] <= mismatch: ###does it pass the mismatch cutoff?
                                 draw.polygon((data['x']+ss1,data['ref_y']+data['reference_thick']+17,data['x']+ss2,data['ref_y']+decay-17,data['x']+ee2,data['ref_y']+decay-17,data['x']+ee1,data['ref_y']+data['reference_thick']+17), outline=color[outline_color], fill=color[fill_color])
                           else:
                              if end2_list[end2] <= mismatch: ###does it pass the mismatch cutoff?
                                 draw.polygon((data['x']+start1,data['ref_y']+data['reference_thick'],data['x']+start2,data['ref_y']+decay,data['x']+end2,data['ref_y']+decay,data['x']+end1,data['ref_y']+data['reference_thick']), outline=color[outline_color], fill=color[fill_color])


      #enhancer = ImageEnhance.Sharpness(im)
      #for i in range(8):
      #   factor = i / 4.0
      #   enhancer.enhance(factor).show("Sharpness %f" % factor)

      ###getFileName
      #xm_regex = re.compile('(\S+)\.\S+')
      #xm_name = xm_regex.match(crossmatch_file)
      #file = xm_name.group(1) + "_m" + str(mismatch) + "_b" + str(block_length) + "_l" + str(leap) + "_s" + str(scale) + "." + format
      file = crossmatch_file + "_m" + str(mismatch) + "_b" + str(block_length) + "_l" + str(leap) + "_s" + str(scale) + "." + format
      print "Saving %s..." % file
      im.save(open(file, 'wb'), formatdict[format])
      print "done."
      return file

#---------------------------------------------
def main():
   
   opts, args = getopt.getopt(sys.argv[1:], "x:s:q:m:r:c:l:f:p:")

   (crossmatch_file, reference_file, query_file, format)=(None,None,None,None)
   (mismatch, block_length, scale, leap, protein)=(0,0,0,0,0)
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
      if o == "-r":
        block_length=int(v)
      if o == "-c":
        scale=int(v)
      if o == "-l":
        leap=int(v)
      if o == "-f":
        format=str(v)
      if o == "-p":
        protein=int(v)


   if (crossmatch_file == None or reference_file == None or query_file == None or mismatch == 0 or block_length == 0 or scale ==0 or leap == 0):
      print "Usage: %s" % (sys.argv[0:])
      print "-x crossmatch file"
      print "-s reference genome fasta file"
      print "-q query contig/genome fasta file"
      print "-m mismatch threshold (e.g. -m 10 allows representation of repeats having up to 10% mismatch"
      print "-r length of similarity block to display"
      print "-c scale (for displaying the image)"
      print "-l leap to evaluate repeat frequency (smaller numbers will increase the resolution, but will affect drastically the run time.  recommended -l=50"
      print "-f file format (bmp, jpeg, png, ps, gif, pdf, tiff) NOTE: the png, ps, tiff and bmp are much better.\n"
      print "-p transform bacterial ORF into protein (i.e. plot alignment between ORF products? 1/0)\n";
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

   #===Scale checks
   if (scale<1):
      print "Not a possible scale. Make sure you select a number >1."
      sys.exit(1)

   #====File checks
   checkFile(crossmatch_file)
   checkFile(reference_file)
   checkFile(query_file)

   #====Parse Fasta Files
   (reference, reflength)=readFasta(reference_file, scale)
   (query, qrylength)=readFasta(query_file, scale)

   print "Reading Crossmatch file..."
   (nocdt, match, query_hit)=readCrossMatch(crossmatch_file, mismatch, block_length, reference, scale)
   print "done."
   print "Computing Repeat frequencies..."
   (freq)=generateCoords(nocdt, reflength, leap, protein)
   print "done."
   print "Drawing repeats..."
   drawRelationship(reference, query, match, scale, query_hit, mismatch, block_length, crossmatch_file, freq, reflength, leap, format, formatdict, protein)

#---------------------------------------------
#Main Call

main()
sys.exit(1)



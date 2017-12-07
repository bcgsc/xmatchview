#!/usr/bin/python

#NAME
#   XMatchView.py v0.2  Rene Warren, March 2005/May 2005/January 2006
 
#SYNOPSIS
#   Allows users to align 2 DNA sequences in fasta format using cross_match and displays the alignment in a variety of formats

#LICENSE
#   Copyright (c) 2004-2006 Canada's Michael Smith Genome Science Centre.  All rights reserved.

#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.

#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.

#INSTALL

#   A) If you're running this program outside the GSC, you will need to do the following before you can proceed:
#      1)Download python2.3 or 2.4 from: http://www.python.org/ and change the shebang line to reflect this
#      2)Download the Python Imaging Library (PIL) from: http://www.pythonware.com/products/pil/
#      3)Copy true type fonts from c:\WINDOWS\Fonts to a unix directory and change the line truetype= below to reflect the location of your ttf 
#      4)Change the sys.path.append line below to reflect the location of PIL
#      5)Make sure cross_match is in your $PATH or change the line cross_match_exec= below
#      6)Copy the image pbp.gif to the same directory where the XMatchView.py program resides, make a fake gif with that name or comment the whole "###Just for fun code block" below

#   B) If you're running this program remotely, but on the GSC servers make sure you are running it on xhost01.bcgsc.ca
#   C) A windows version of this program exists (XMatchView_win.py).  However, you won't be able to run crossmatch with that version unless you have purchased cross_match for windows.


import sys
import commands 

sys.path.insert(0,'/home/rwarren/python/standard/Imaging-1.1.4/PIL')
#sys.path.append('/home/rwarren/python/standard/Imaging-1.1.4/PIL')
from Tkinter import *
import tkSimpleDialog
import tkFileDialog
import tkMessageBox
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
import ImageTk

#############Global Variables
default_minmatch=30
default_minscore=60
default_mismatch=30
default_block=500
default_scale=1000
default_leap=50
default_width=2400
default_height=1200
default_angle=0
default_format="gif"
truetype="/home/rwarren/fonts/truetype/"
arial=truetype + "arial.ttf"
arial_bold=truetype + "arialbd.ttf"
arial_italic=truetype + "ariali.ttf"
arial_bold_italic=truetype + "arialbi.ttf"
cross_match_exec="cross_match"
###############################

formatdict = {'png':'PNG','gif':'GIF','tiff':'TIFF','bmp':'BMP','jpeg':'JPEG','ps':'EPS'}

#===============================================================
class XMDialog(tkSimpleDialog.Dialog):

    def body(self, master):

        self.result = {} 

        Label(master, text="Minmatch").grid(row=0, sticky=W)
        Label(master, text="Minscore").grid(row=1, sticky=W)

        self.e1 = Entry(master)
        self.e1.insert(0,default_minmatch)

        self.e2 = Entry(master)
        self.e2.insert(0,default_minscore)

        self.e1.grid(row=0, column=1)
        self.e2.grid(row=1, column=1)
        
        #return self.e1 # initial focus

    def validate(self):

       try:
          minmatch = int(self.e1.get())
          minscore = int(self.e2.get())

          if (minmatch > 0 and minscore > 0):
             self.result["minmatch"] = minmatch
             self.result["minscore"] = minscore

             return 1
          else:
             tkMessageBox.showwarning(
             "Bad input",
             "Illegal values, please try again"
             )
             return 0
 
       except ValueError:
          tkMessageBox.showwarning(
          "Bad input",
          "Illegal data type, please try again"
          )
          return 0

#===============================================================
class XGDialog(tkSimpleDialog.Dialog):

    def body(self, master):

        self.result = {}

        Label(master, text="Mismatch").grid(row=0, sticky=W)
        Label(master, text="Minimum Block Length").grid(row=1, sticky=W)
        Label(master, text="Scale 1:").grid(row=2, sticky=W)
        Label(master, text="Sliding Window Leap").grid(row=3, sticky=W)
        Label(master, text="Save as").grid(row=4, sticky=W)

        self.e1 = Entry(master)
        self.e1.insert(0,default_mismatch)
        self.e2 = Entry(master)
        self.e2.insert(0,default_block)
        self.e3 = Entry(master)
        self.e3.insert(0,default_scale)
        self.e4 = Entry(master)
        self.e4.insert(0,default_leap)

        self.e1.grid(row=0, column=1, sticky=W)
        self.e2.grid(row=1, column=1, sticky=W)
        self.e3.grid(row=2, column=1, sticky=W)
        self.e4.grid(row=3, column=1, sticky=W)

        self.var = StringVar()
        self.var.set(default_format) # initial value

        self.option = OptionMenu(master, self.var,'png', 'gif', 'tiff', 'bmp', 'jpeg', 'ps')
        self.option.grid(row=4, column=1, sticky=W)


    def validate(self):
       try:
          mismatch = int(self.e1.get())
          block = int(self.e2.get())
          scale = int(self.e3.get())
          leap = int(self.e4.get())
          graph_format = self.var.get()

          if (mismatch >= 0 and mismatch < 100 and block > 0 and scale > 0 and leap >0):
             self.result["mismatch"] = mismatch
             self.result["block"] = block
             self.result["scale"] = scale
             self.result["leap"] = leap
             self.result["format"] = graph_format
             return 1
          elif (mismatch < 0 or mismatch >= 100):
             tkMessageBox.showwarning(
             "Bad input",
             "Illegal mismatch value.\nMust be between 0-99."
             )
             return 0
          elif (block < 1):
             tkMessageBox.showwarning(
             "Bad input",
             "Illegal block length value.\nMust be larger than zero."
             )
             return 0
          elif (scale < 1):
             tkMessageBox.showwarning(
             "Bad input",
             "Illegal scale value.\nMust be smaller than, or equal 1:1"
             )
             return 0
          elif (leap < 1):
             tkMessageBox.showwarning(
             "Bad input",
             "Illegal sliding window leap value.\nMust be larger than zero."
             )
             return 0

       except ValueError:
          tkMessageBox.showwarning(
          "Bad input",
          "Illegal data type, please try again"
          )
          return 0


#===============================================================
class MatchViz:
   def __init__(self, parent):
      ###root becomes parent
      self.myParent = parent
      self.im = None
     
      self.crossMatchFile = None
      self.fasta_query = None
      self.fasta_reference = None
      self.graph_file = None

      self.xm_exec=cross_match_exec
      self.minmatch=default_minmatch
      self.minscore=default_minscore
      self.mismatch=default_mismatch
      self.scale=default_scale
      self.block=default_block
      self.leap=default_leap
      self.width=default_width
      self.height=default_height
      self.angle=default_angle
      self.format=default_format

      self.menubar = Menu(parent)    #Container Menu

      self.oldCursor=parent["cursor"]

      ####SCROLL
      self.scrollbarx = Scrollbar(parent, orient='horizontal')
      self.scrollbary = Scrollbar(parent, orient='vertical')
      self.scrollbary.pack(side=RIGHT, fill=Y)
      self.scrollbarx.pack(side=BOTTOM, fill=X)

      ####CANVAS
      self.can = Canvas(parent, width=1200, height=900, background='white', xscrollcommand=self.scrollbarx.set, yscrollcommand=self.scrollbary.set, scrollregion=(0, 0, default_width, default_height))

      ####FILE Menu (pull down)
      self.filemenu = Menu(self.menubar, tearoff=0)
      self.filemenu.add_command(label="Open CrossMatch Output", underline=5, command=self.openCM)
      self.filemenu.add_command(label="Open Graph", underline=0, command=self.openGraph)
      self.filemenu.add_separator()
      self.filemenu.add_command(label="Select Reference Fasta", underline=5, command=self.openFaR)
      self.filemenu.add_command(label="Select Query Fasta", underline=5, command=self.openFaQ)
      self.filemenu.add_separator()
      self.filemenu.add_command(label="Exit", underline=1, command=parent.quit)
      self.menubar.add_cascade(label="File", underline=0, menu=self.filemenu)
      
      ####OPTION Menu
      self.optmenu = Menu(self.menubar, tearoff=0)
      self.optmenu.add_command(label="CrossMatch", underline=5, command=self.optionXM)
      self.optmenu.add_separator()
      self.optmenu.add_command(label="Graph", underline=0, command=self.optionXG)
      self.menubar.add_cascade(label="Option", underline=0, menu=self.optmenu)

      ####VIEW Menu
      self.viewmenu = Menu(self.menubar, tearoff=0)
      self.viewmenu.add_command(label="Zoom In", underline=5, command=self.zoomIn)
      self.viewmenu.add_command(label="Zoom Out", underline=5, command=self.zoomOut)
      self.viewmenu.add_separator()
      self.CrossMatchWindow=BooleanVar()
      self.viewmenu.add_checkbutton(label="CrossMatch Window", variable=self.CrossMatchWindow, command=self.showCM)
      #self.viewmenu.add_command(label="Rotate Clock Wise", command=self.rotate)
      self.menubar.add_cascade(label="View", underline=0, menu=self.viewmenu)
    
      ####TOOL Menu
      self.toolmenu = Menu(self.menubar, tearoff=0)
      self.toolmenu.add_command(label="Run CrossMatch", underline=0, command=self.runXM)
      self.filemenu.add_separator()
      self.toolmenu.add_command(label="Draw Repeat Graph", underline=0, command=self.runXG)
      self.menubar.add_cascade(label="Tool", underline=0, menu=self.toolmenu)

      ####HELP Menu
      self.helpmenu = Menu(self.menubar, tearoff=0)
      self.helpmenu.add_command(label="About", underline=0, command=self.help)
      self.menubar.add_cascade(label="Help", underline=0, menu=self.helpmenu)

      ####Just for fun:
      self.default_image = Image.open('/home/rwarren/python/Development/SeqDev/bin/pbp.gif')
      self.default_open = ImageTk.PhotoImage(self.default_image)
      self.label = Label(image=self.default_open)
      self.label.pack()


      self.scrollbarx.config(command=self.can.xview)
      self.scrollbary.config(command=self.can.yview)

      self.myParent.config(menu=self.menubar)

      self.can.pack()

   #---------------------------------------------
   def runXM(self):

      if (self.fasta_query == None):
         tkMessageBox.showwarning("File missing", "You must open a\nquery fasta file.")
      elif (self.fasta_reference == None):
         tkMessageBox.showwarning("File missing", "You must open a\nreference fasta file.")
      else:
         self.crossMatchFile = tkFileDialog.asksaveasfilename(defaultextension = ".rep", filetypes = [("CrossMatch", "*.rep"), ("CrossMatch", "*.txt"), ("CrossMatch", "*.screen")])
         command=self.xm_exec + " " + self.fasta_reference + " " + self.fasta_query + " -minmatch " + str(self.minmatch) + " -minscore " + str(self.minscore) +" -masklevel 101 >" + self.crossMatchFile 
         if (self.crossMatchFile != '' and self.crossMatchFile != None):
            if (tkMessageBox.askokcancel("Proceed?", "Run %s ?" % command) and self.crossMatchFile):

               self.myParent["cursor"]="watch"
               self.myParent.update()

               (status, out)=commands.getstatusoutput(command)
               print "%s" % out

               self.CrossMatchWindow.set(1)
               self.showCM()
         
               self.myParent["cursor"]=self.oldCursor
               self.myParent.update()

   #---------------------------------------------
   def runXG(self):

      if (self.fasta_query == None):
         tkMessageBox.showwarning("File missing", "You must open a\nquery fasta file.")
      elif (self.fasta_reference == None):
         tkMessageBox.showwarning("File missing", "You must open a\nreference fasta file.")
      elif (self.crossMatchFile == None):
         tkMessageBox.showwarning("File missing", "You must open a\nCrossMatch file.")
      elif (tkMessageBox.askokcancel("Proceed?", "Visualize Collinear Blocks? \nScale 1:%s \nMismatch threshold = %s \nMin. Block Length = %s \nSliding window leap = %s \nGraph format = %s" % (self.scale, self.mismatch, self.block, self.leap, self.format))):
         #====Parse Fasta Files
         self.myParent["cursor"]="watch"
         self.myParent.update()

         (reference, reflength)=self.readFasta(self.fasta_reference, self.scale)
         (query, qrylength)=self.readFasta(self.fasta_query, self.scale)
         print "Reading Crossmatch file..."
         (nocdt, match, query_hit)=self.readCrossMatch(self.crossMatchFile, self.mismatch, self.block, reference, self.scale, reflength, qrylength)
         print "done.\nComputing Repeat frequencies..."
         (freq)=self.generateCoords(nocdt, reflength, self.leap)
         print "done.\nDrawing Collinear Blocks..."
         self.graph_file=self.drawRelationship(reference, query, match, self.scale, query_hit, self.mismatch, self.block, self.crossMatchFile, freq, reflength, self.leap, self.format, formatdict)

         self.im = Image.open(self.graph_file)
         self.photo_open = ImageTk.PhotoImage(self.im)
         self.can.delete("repeatgraph")
         self.can.create_image(1, 1, anchor=NW, image=self.photo_open, tag="repeatgraph")
         self.can.pack()

         self.myParent["cursor"]=self.oldCursor
         self.myParent.update()

   #---------------------------------------------
   def zoomIn(self):

      if self.im != None:
         self.width *= 1.5
         self.height *= 1.5
 
         self.myParent["cursor"]="watch"
         self.myParent.update()

         self.can.delete("repeatgraph")

         self.zoomin = ImageTk.PhotoImage(self.im.resize((int(self.width),int(self.height))))
         self.can.create_image(1, 1, anchor=NW, image=self.zoomin, tag="repeatgraph")
         self.can.pack()

         self.myParent["cursor"]=self.oldCursor
         self.myParent.update()

   #---------------------------------------------
   def zoomOut(self):

      if self.im != None:
         self.width /= 1.5 
         self.height /= 1.5

         self.myParent["cursor"]="watch"
         self.myParent.update()

         self.can.delete("repeatgraph")

         self.zoomout = ImageTk.PhotoImage(self.im.resize((int(self.width),int(self.height))))
         self.can.create_image(1, 1, anchor=NW, image=self.zoomout, tag="repeatgraph")
         self.can.pack()

         self.myParent["cursor"]=self.oldCursor
         self.myParent.update()

   #---------------------------------------------
   def rotate(self):

      if self.im != None:
         self.angle += 90
       
         self.rotate = ImageTk.PhotoImage(self.im.rotate(self.angle))
       
         self.can.delete("repeatgraph")
 
         self.can.create_image(1, 1, anchor=NW, image=self.rotate, tag="repeatgraph")
         self.can.pack()

   #---------------------------------------------
   def openCM(self):
      self.crossMatchFile=tkFileDialog.askopenfilename(defaultextension = ".rep", filetypes = [("CrossMatch", "*.rep"), ("CrossMatch", "*.txt"), ("CrossMatch", "*.screen")])

      if (self.crossMatchFile != '' and self.crossMatchFile != None): 
         self.CrossMatchWindow.set(1)
         self.showCM()

   #---------------------------------------------
   def showCM(self):

      if (self.CrossMatchWindow.get()):

         self.top = Toplevel()
         self.top.destroy()
         self.top = Toplevel()
         self.MainFrame=Frame(self.top)
         self.top.title(self.crossMatchFile)
         self.TextBox=Text(self.MainFrame)

         if (self.crossMatchFile==None or self.crossMatchFile==""):
            return
         try:
            File=open(self.crossMatchFile,"r")
            NewText=File.read()
            File.close()
            self.FileName=self.crossMatchFile
            #self.top.title(self.crossMatchFile)       
         except IOError:
            tkMessageBox.showerror("Read error...",
                "Could not read from '%s'" % self.crossMatchFile) 
            return

         self.ClearText()
         self.TextBox.insert(END,NewText)
         self.TextBox.pack(fill=BOTH,expand=YES)
         self.MainFrame.pack(fill=BOTH,expand=YES)
      else:
         self.top.destroy()

   #---------------------------------------------
   def ClearText(self):
        self.TextBox.delete("1.0",END)

   #---------------------------------------------
   def openFaQ(self):
      self.fasta_query=tkFileDialog.askopenfilename(defaultextension = ".fa", filetypes = [("fasta", "*.fa"), ("fasta", "*.fasta"), ("fasta", "*.txt")])
 
   #---------------------------------------------
   def openFaR(self):
      self.fasta_reference=tkFileDialog.askopenfilename(defaultextension = ".fa", filetypes = [("fasta", "*.fa"), ("fasta", "*.fasta"), ("fasta", "*.txt")])

   #---------------------------------------------
   def optionXM(self):
      xmdialog = XMDialog(self.myParent, title="CrossMatch Options")

      if xmdialog.result.has_key("minmatch") and xmdialog.result.has_key("minscore"):
         self.minmatch=xmdialog.result["minmatch"]
         self.minscore=xmdialog.result["minscore"]

   #---------------------------------------------
   def optionXG(self):
      xgdialog = XGDialog(self.myParent, title="Graphics Options")
      
      if xgdialog.result.has_key("mismatch") and xgdialog.result.has_key("block") and xgdialog.result.has_key("scale") and xgdialog.result.has_key("leap"):
         self.mismatch=xgdialog.result["mismatch"]
         self.block=xgdialog.result["block"]
         self.scale=xgdialog.result["scale"]
         self.leap=xgdialog.result["leap"]
         self.format=xgdialog.result["format"]

   #---------------------------------------------
   def help(self):

      tkMessageBox.showinfo("About", "XMatchView.py\nCopyright 2004-2006\nRene Warren") 

   #---------------------------------------------
   def openGraph(self):
      self.graph_file=tkFileDialog.askopenfilename(defaultextension = ".fa", filetypes = [("GIF", "*.gif"), ("PNG", "*.png"), ("TIFF", "*.tiff"), ("BMP", "*.bmp"), ("JPEG", "*.jp*g"), ("postcript EPS", "*.ps")])

      if (self.graph_file != '' and self.graph_file != None):

         self.myParent["cursor"]="watch"
         self.myParent.update()

         self.im = Image.open(self.graph_file)
         self.photo_open = ImageTk.PhotoImage(self.im)
         self.can.delete("repeatgraph")
         self.can.create_image(1, 1, anchor=NW, image=self.photo_open, tag="repeatgraph")
         self.can.pack()

         self.myParent["cursor"]=self.oldCursor
         self.myParent.update()

   #---------------------------------------------
   def readCrossMatch(self,crossmatch_file,mismatch,block_length,reference,scale,reflength,qrylength):

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
            (percentMis, primary_match, startFirstMatch, endFirstMatch, secondary_match, startSecondMatch, endSecondMatch)=(float(rm.group(2)), rm.group(3), int(rm.group(4)), int(rm.group(5)), rm.group(6), int(rm.group(7)), int(rm.group(8)))

            repeat_size = endFirstMatch - startFirstMatch + 1
            print "%i-%i l=%i q=%i rsize=%i" % (endFirstMatch,startFirstMatch,reflength,qrylength,repeat_size)

            if(repeat_size >= reflength-100):
               continue
            elif (repeat_size >= block_length) and (percentMis <= mismatch):  ### REPLACE BY "else:" if you want to see all repeats+frequency

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
            elif (primary_match == secondary_match) and (startSecondMatch == startFirstMatch):  ###prevents exact matches
               continue
            elif (repeat_size < block_length) or (repeat_size >= reflength-100):
               continue  #will skip smaller alignment
            else:
               if reference.has_key(primary_match):
                  startFirstMatch=int(startFirstMatch/scale)
                  endFirstMatch=int(endFirstMatch/scale)
                  startSecondMatch=int(startSecondMatch/scale)
                  endSecondMatch=int(endSecondMatch/scale)

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

            (percentMis, primary_match, startFirstMatch, endFirstMatch, secondary_match, startSecondMatch, endSecondMatch)=(float(fm.group(2)), fm.group(3), int(fm.group(4)), int(fm.group(5)), fm.group(6), int(fm.group(7)), int(fm.group(8)))

            repeat_size = endFirstMatch - startFirstMatch + 1
            print "%i-%i l=%i q=%i rsize=%i" % (endFirstMatch,startFirstMatch,reflength,qrylength,repeat_size)

            if(repeat_size >= reflength-100):
               continue
            elif (repeat_size >= block_length) and (percentMis <= mismatch):  ### REPLACE BY "else:" if you want to see all repeats+frequency

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
            elif (primary_match == secondary_match) and (startSecondMatch == startFirstMatch):  ###prevents exact matches
               continue
            elif (repeat_size < block_length) or (repeat_size >= reflength-100):
               continue  #will skip smaller alignment
            else:
               if reference.has_key(primary_match):
                  startFirstMatch=int(startFirstMatch/scale)
                  endFirstMatch=int(endFirstMatch/scale)
                  startSecondMatch=int(startSecondMatch/scale)
                  endSecondMatch=int(endSecondMatch/scale)

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

      return nocdt, match, query_hit

   #---------------------------------------------
   def generateCoords(self, nocdt, size, leap):

      freq={}

      pos_range=range(0,size,leap)

      for pos in pos_range:
         print "%i out of %i bases" % (pos,size)
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
                     if((pos >= start1 and pos <= end1) or (pos >= end1 and pos <= start1)):
                        #print "%i >= %i and %i<=%i OR %i>=%i and %i<=%i" % (pos,start1,pos,end1,pos,end1,pos,start1)
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
   def readFasta(self, file, scale):
  
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

         seq_subset_regex = re.compile('^(\S+)$', re.I)
         seq_subset = seq_subset_regex.match(line)
         if seq_subset != None:
            seq_length += len(seq_subset.group(1))

      (seq_length, scale)=(int(seq_length), int(scale))
      L1[previous_contig] = float(seq_length/scale)                                #for the last sequence

      file_obj.close()

      print "scaled down %s =%f total=%i " % (previous_contig, L1[previous_contig], seq_length)

      return (L1, seq_length)

   #---------------------------------------------
   def initColor(self):
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
   def initGraph(self):
      data={}
 
      #default data points
      data['width']=default_width
      data['height']=default_height
      data['ref_y']=250
      data['mis_bar']=50
      data['query_y']=70
      data['x']=20
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
   def drawRectangle(self,draw,start,end,y,thickness,bar_color,text,font,text_color):
  
      draw.rectangle((start,y,end,y+thickness), bar_color)
      draw.text((end+5, y-2), text, font=font, fill=text_color)

   #---------------------------------------------
   def plotFrequency(self,freq,size,scale,draw,color,data,leap):

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
   def drawRelationship(self,reference_list, query_list, match_list, scale, query_hit, mismatch, block_length, crossmatch_file, freq, reflength, leap, format, formatdict):

      scaled_reflength=int(reflength/scale)

      ###Initialize new graph
      data=self.initGraph()

      ###Get colors
      color=self.initColor()
 
      ###Set Font
      arial_18=ImageFont.truetype(arial,18)
      arialb_18=ImageFont.truetype(arial_bold,18)
      arial_20=ImageFont.truetype(arial,20)
      arialb_20=ImageFont.truetype(arial_bold,20)
      ariali_20=ImageFont.truetype(arial_italic,20)
      arialbi_20=ImageFont.truetype(arial_bold_italic,20)
      arialb_22=ImageFont.truetype(arial_bold,22)

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

         self.drawRectangle(draw,init_coord, last_coord,data['ref_y'],data['reference_thick'],color['black'],ref,arialb_18,color['black'])
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
      self.plotFrequency(freq,reflength,scale,draw,color,data,leap)

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
               self.drawRectangle(draw,current_position,stop,data['ref_y']+decay,data['query_thick'],color['black'], hit, arialb_18, color['black'])
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

                      
                        draw.rectangle((data['x']+start1,data['ref_y']+data['reference_thick']-LCB,data['x']+end1,data['ref_y']+data['reference_thick']), color[fill_color])
                        print "%s-%s S1:%i S2:%i" % (match,hit,start1,start2)


                        if match == hit:
                           if start1 <= start2:

                              repeat_size = start2 - start1
                              size_chunk = int(decay * repeat_size / scaled_reflength)
                              print "%i %i %i" % (size_chunk, repeat_size, scaled_reflength)
                              size_chunk += 50

                              draw.arc((data['x']+start1,data['ref_y']+data['reference_thick']-size_chunk,data['x']+start2,data['ref_y']+data['reference_thick']+size_chunk),360,180, color[outline_color]) 
                        else:
                           draw.rectangle((current_position+start2,data['ref_y']+decay,current_position+end2,data['ref_y']+decay+LCB), color[fill_color])
                           draw.polygon((data['x']+start1,data['ref_y']+data['reference_thick'],current_position+start2,data['ref_y']+decay,current_position+end2,data['ref_y']+decay,data['x']+end1,data['ref_y']+data['reference_thick']), outline=color[outline_color], fill=color[fill_color])


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
# display the menu

root = Tk()
root.title(sys.argv[0:])
matchviz = MatchViz(root)
root.mainloop()

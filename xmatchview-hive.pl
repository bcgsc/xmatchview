#!/usr/bin/env perl
# xmatchview-hive
# Visualizing 3-way genome synteny with hive plot
# Rene L Warren 2005-2020

use strict;
my $version = "v1.2.1";
my $op = 0.75;###sets opacity

if($#ARGV<1){
   die "Usage:\t$0 $version\n\t< .rep/.paf (concat all cross_match/paf alignment output files into 1) >\n\t< scaling factor >\n\t< opacity 0-1 (optional) >\n";
}

my $xmf = $ARGV[0];
my $out = $ARGV[0] . ".svg"; ## this will be the xml svg output
my $scale = $ARGV[1];
$op = $ARGV[2] if($ARGV[2] ne "");

my $config="config.txt";

if(! -e $config){
   die "File $config does not exists. Makes sure it does at run location -- fatal. The format is:\n1:name1:sequence1_length\n2:name2:sequence2_length\n3:name3:sequence3_length\n";
}

if(! -e $xmf){
   die "Alignment file $xmf does not exists. Makes sure it does -- fatal.\n";
}


#####
#read config

open(IN,$config) || die "can't read $config -- fatal.\n";
my $info;

while(<IN>){
   chomp;
   my @a=split("\:");
   $info->{$a[1]}{'axis'}=$a[0];
   $info->{$a[1]}{'sz'}=$a[2];
}

close IN;

#================================ WRITE SVG XML ================================
### image settings
my $width = 2000;
my $height = 2000;
my $mid = $height / 2;
my $midwidth = $width / 2;
my $xlegend = 0;
my $ylegend= 0;


# xml / SVG
open(SVG, ">$out") || die "Can't write to $out -- fatal.\n";
print SVG "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n";
print SVG "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.0//EN\" \"http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd\">\n";
print SVG "<svg width=\"$width\" height=\"$height\" xmlns=\"http://www.w3.org/2000/svg\" style=\"background-color:white\">\n";

# print axes xml/svg
foreach my $name (keys %$info){
   if($info->{$name}{'axis'} == 1){
      my $y1 = $mid-50;#start
      my $y2 = $y1 - (($info->{$name}{'sz'} / $scale));
      my $x1 = $midwidth;
      my $x2 = $x1;

      my $xt=$x1+5;
      my $yt=$y2-5;
      $ylegend = $yt;

      print SVG "<path d=\"M $x1 $y1 L $x1 $y2\" stroke=\"black\" stroke-width=\"5\" fill=\"transparent\"/>\n";
      print SVG "<text font-size=\"2.5em\" x=\"$xt\" y=\"$yt\">$name</text>\n";

   }elsif($info->{$name}{'axis'} == 2){

      my $x1 = $midwidth + 50;
      my $y1 = $mid + 50;

      my $len = $info->{$name}{'sz'} / $scale;
      my $len2 = $len ** 2;
      my $cal = sqrt($len2/2);

      my $x2 = (sin(45)*$len) + $x1;
      my $y2 = (cos(45)*$len) + $y1;

      my $xt=$x2+5;
      my $yt=$y2+5;

      $xlegend=$xt;

      print SVG "<path d=\"M $x1 $y1 L $x2 $y2\" stroke=\"black\" stroke-width=\"5\" fill=\"transparent\"/>\n";
      print SVG "<text font-size=\"2.5em\" x=\"$xt\" y=\"$yt\">$name</text>\n";

      #my $yrlw = $yt + 100;
      #print SVG "<text font-size=\"2.0em\" x=\"$xt\" y=\"$yrlw\">RLW2020</text>\n";

   }elsif($info->{$name}{'axis'} == 3){

      my $x1 = $midwidth - 50;
      my $y1 = $mid + 50;

      my $len = $info->{$name}{'sz'} / $scale;

      my $x2 = (sin(-45)*$len) + $x1;
      my $y2 = (cos(-45)*$len) + $y1;

      my $xt=$x2-200;
      my $yt=$y2+5;

      print SVG "<path d=\"M $x1 $y1 L $x2 $y2\" stroke=\"black\" stroke-width=\"5\" fill=\"transparent\"/>\n";
      print SVG "<text font-size=\"2.5em\" x=\"$xt\" y=\"$yt\">$name</text>\n";
   }
}


#####
# print alignments svg/xml

if($xmf=~/\.paf$/){

 ### reverse matches      qryname    qrystart qryend orient hitname       hitstart hitend match   block
 ###                       1             2       3            4             5       6       7       8
 ### rev_regex = re.compile("(\S+)\s+\S+\s+(\d+)\s+(\d+)\s+\-\s+(\S+)\s+\S+\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)")
 ### rm = rev_regex.match(line)
 ###forward matches
 ### fwd_regex = re.compile("(\S+)\s+\S+\s+(\d+)\s+(\d+)\s+\+\s+(\S+)\s+\S+\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)")
   open(IN, $xmf) || die "can't read $xmf -- fatal.\n";

   while(<IN>){
      # example paf output format
      chomp;
      if(/(\S+)\s+\S+\s+(\d+)\s+(\d+)\s+\+\s+(\S+)\s+\S+\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)/){
         my ($q,$qs,$qe,$t,$ts,$te,$match,$block) = ($1,$2,$3,$4,$5,$6,$7,$8);
         my $si = ($match / ($te-$ts+1)) * 100;
         $qs/=$scale;
         $qe/=$scale;
         $ts/=$scale;
         $te/=$scale;
         my $identity;

         #identity gradient
         if($si==100){
            $identity = "#005824";
         }elsif($si>=90){
            $identity = "#238b45";
         }elsif($si>=80){
            $identity = "#41ae76";
         }elsif($si>=70){
            $identity = "#66c2a4";
         }elsif($si>=60){
            $identity = "#99d8c9";
         }elsif($si>=50){
            $identity = "#ccece6";
         }elsif($si<50){
            $identity = "#edf8fb";
         }

         if($info->{$q}{'axis'}==1 && $info->{$t}{'axis'}==2){
            ### axis 1
            my $xq1 = $midwidth+5;
            my $xq2 = $xq1;
            my $yq1 = $mid - 50 - $qs;
            my $yq2 = $mid - 50 - $qe;
            ### axis2   
            my $xstart_axis2 = $midwidth + 50;
            my $ystart_axis2 = $mid + 45;
            my $xt1 = (sin(45)*$ts) + $xstart_axis2;
            my $yt1 = (cos(45)*$ts) + $ystart_axis2;
            my $xt2 = (sin(45)*$te) + $xstart_axis2;
            my $yt2 = (cos(45)*$te) + $ystart_axis2;   
  
            my $xt1e = $xt1;# + 0.1*$xt1 + ($xt1-$xstart_axis3)*1.2;
            my $xt2e = $xt2;# + 0.1*$xt2 + ($xt2-$xstart_axis3)*1.2;
            my $yq1e = (($yt1 - $yq1)/3) + $yq1;#    0.1*$yq1 + ($yq1-$ystart_axis3)*1.2;
            my $yq2e = (($yt2 - $yq2)/3) + $yq2;#   $yq2 + 0.1*$yq2 + ($yq2-$ystart_axis3)*1.2;

            print SVG "<path d=\"M $xq2 $yq2 L $xq1 $yq1 Q $xt1e $yq1e $xt1 $yt1 L $xt2 $yt2 Q $xt2e $yq2e $xq2 $yq2\" stroke=\"$identity\" stroke-width=\"0.5\" fill=\"$identity\" fill-opacity=\"$op\"/>\n";
 
         }elsif($info->{$q}{'axis'}==3 && $info->{$t}{'axis'}==2){
            ### axis 2
            my $xstart_axis2 = $midwidth + 50;
            my $ystart_axis2 = $mid + 55;
            my $xq1 = (sin(45)*$ts) + $xstart_axis2;
            my $yq1 = (cos(45)*$ts) + $ystart_axis2;
            my $xq2 = (sin(45)*$te) + $xstart_axis2;
            my $yq2 = (cos(45)*$te) + $ystart_axis2;
            ### axis 3
            my $xstart_axis3 = $midwidth - 50;
            my $ystart_axis3 = $mid + 55;
            my $xt1 = (sin(-45)*$qs) + $xstart_axis3;
            my $yt1 = (cos(-45)*$qs) + $ystart_axis3;
            my $xt2 = (sin(-45)*$qe) + $xstart_axis3;
            my $yt2 = (cos(-45)*$qe) + $ystart_axis3;

            my $yq1e = $yq1 + 0.05*$yq1 + ($yq1-$ystart_axis3)*1.1;
            my $yq2e = $yq2 + 0.05*$yq2 + ($yq2-$ystart_axis3)*1.1;

            print SVG "<path d=\"M $xq2 $yq2 L $xq1 $yq1 Q $mid $yq1e $xt1 $yt1 L $xt2 $yt2 Q $mid $yq2e $xq2 $yq2\" stroke=\"$identity\" stroke-width=\"0.5\" fill=\"$identity\" fill-opacity=\"$op\"/>\n";

         }elsif($info->{$q}{'axis'}==1 && $info->{$t}{'axis'}==3){

            ### axis 1
            my $xq1 = $midwidth-5;
            my $xq2 = $xq1;
            my $yq1 = $mid - 50 - $qs;
            my $yq2 = $mid - 50 - $qe;
            ### axis 3
            my $xstart_axis3 = $midwidth - 50;
            my $ystart_axis3 = $mid + 45;
            my $xt1 = (sin(-45)*$ts) + $xstart_axis3;
            my $yt1 = (cos(-45)*$ts) + $ystart_axis3;
            my $xt2 = (sin(-45)*$te) + $xstart_axis3;
            my $yt2 = (cos(-45)*$te) + $ystart_axis3;

            my $xt1e = $xt1;# + 0.1*$xt1 + ($xt1-$xstart_axis3)*1.2;
            my $xt2e = $xt2;# + 0.1*$xt2 + ($xt2-$xstart_axis3)*1.2;
            my $yq1e = (($yt1 - $yq1)/3) + $yq1 ;#    0.1*$yq1 + ($yq1-$ystart_axis3)*1.2;
            my $yq2e = (($yt2 - $yq2)/3) + $yq2;#   $yq2 + 0.1*$yq2 + ($yq2-$ystart_axis3)*1.2;

            print SVG "<path d=\"M $xq2 $yq2 L $xq1 $yq1 Q $xt1e $yq1e $xt1 $yt1 L $xt2 $yt2 Q $xt2e $yq2e $xq2 $yq2\" stroke=\"$identity\" stroke-width=\"0.5\" fill=\"$identity\" fill-opacity=\"$op\"/>\n";

         }
      ###matches on reverse strand
      }elsif(/(\S+)\s+\S+\s+(\d+)\s+(\d+)\s+\-\s+(\S+)\s+\S+\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)/){
         my ($q,$qs,$qe,$t,$ts,$te,$match,$block) = ($1,$2,$3,$4,$5,$6,$7,$8);
         my $si = ($match / ($te-$ts+1)) * 100;
         $qs/=$scale;
         $qe/=$scale;
         $ts/=$scale;
         $te/=$scale;
         my $identity;

         #identity gradient
         if($si==100){
            $identity = "#99000d";
         }elsif($si>=90){
            $identity = "#cb181d";
         }elsif($si>=80){
            $identity = "#ef3b2c";
         }elsif($si>=70){
            $identity = "#fb6a4a";
         }elsif($si>=60){
            $identity = "#fc9272";
         }elsif($si>=50){
            $identity = "#fcbba1";
         }elsif($si<50){
            $identity = "#fee5d9";
         }

         if($info->{$q}{'axis'}==1 && $info->{$t}{'axis'}==2){
            ### axis 1
            my $xq1 = $midwidth+5;
            my $xq2 = $xq1;
            my $yq1 = $mid - 50 - $qs;
            my $yq2 = $mid - 50 - $qe;
            ### axis2   
            my $xstart_axis2 = $midwidth + 50;
            my $ystart_axis2 = $mid + 45;
            my $xt1 = (sin(45)*$ts) + $xstart_axis2;
            my $yt1 = (cos(45)*$ts) + $ystart_axis2;
            my $xt2 = (sin(45)*$te) + $xstart_axis2;
            my $yt2 = (cos(45)*$te) + $ystart_axis2;   

            my $xt1e = $xt1;# + 0.1*$xt1 + ($xt1-$xstart_axis3)*1.2;
            my $xt2e = $xt2;# + 0.1*$xt2 + ($xt2-$xstart_axis3)*1.2;
            my $yq1e = (($yt1 - $yq1)/3) + $yq1;#    0.1*$yq1 + ($yq1-$ystart_axis3)*1.2;
            my $yq2e = (($yt2 - $yq2)/3) + $yq2;#   $yq2 + 0.1*$yq2 + ($yq2-$ystart_axis3)*1.2;

            print SVG "<path d=\"M $xq2 $yq2 L $xq1 $yq1 Q $xt1e $yq1e $xt1 $yt1 L $xt2 $yt2 Q $xt2e $yq2e $xq2 $yq2\" stroke=\"$identity\" stroke-width=\"0.5\" fill=\"$identity\" fill-opacity=\"$op\"/>\n";

         }elsif($info->{$q}{'axis'}==3 && $info->{$t}{'axis'}==2){
            ### axis 2
            my $xstart_axis2 = $midwidth + 50;
            my $ystart_axis2 = $mid + 55;
            my $xq1 = (sin(45)*$ts) + $xstart_axis2;
            my $yq1 = (cos(45)*$ts) + $ystart_axis2;
            my $xq2 = (sin(45)*$te) + $xstart_axis2;
            my $yq2 = (cos(45)*$te) + $ystart_axis2;
            ### axis 3
            my $xstart_axis3 = $midwidth - 50;
            my $ystart_axis3 = $mid + 55;
            my $xt1 = (sin(-45)*$qs) + $xstart_axis3;
            my $yt1 = (cos(-45)*$qs) + $ystart_axis3;
            my $xt2 = (sin(-45)*$qe) + $xstart_axis3;
            my $yt2 = (cos(-45)*$qe) + $ystart_axis3;

            my $yq1e = $yq1 + 0.1*$yq1 + ($yq1-$ystart_axis3)*1.2;
            my $yq2e = $yq2 + 0.1*$yq2 + ($yq2-$ystart_axis3)*1.2;

            print SVG "<path d=\"M $xq2 $yq2 L $xq1 $yq1 Q $mid $yq1e $xt1 $yt1 L $xt2 $yt2 Q $mid $yq2e $xq2 $yq2\" stroke=\"$identity\" stroke-width=\"0.5\" fill=\"$identity\" fill-opacity=\"$op\"/>\n";

         }elsif($info->{$q}{'axis'}==1 && $info->{$t}{'axis'}==3){

            ### axis 1
            my $xq1 = $midwidth-5;
            my $xq2 = $xq1;
            my $yq1 = $mid - 50 - $qs;
            my $yq2 = $mid - 50 - $qe;
            ### axis 3
            my $xstart_axis3 = $midwidth - 50;
            my $ystart_axis3 = $mid + 45;
            my $xt1 = (sin(-45)*$ts) + $xstart_axis3;
            my $yt1 = (cos(-45)*$ts) + $ystart_axis3;
            my $xt2 = (sin(-45)*$te) + $xstart_axis3;
            my $yt2 = (cos(-45)*$te) + $ystart_axis3;

            my $xt1e = $xt1;# + 0.1*$xt1 + ($xt1-$xstart_axis3)*1.2;
            my $xt2e = $xt2;# + 0.1*$xt2 + ($xt2-$xstart_axis3)*1.2;
            my $yq1e = (($yt1 - $yq1)/3) + $yq1 ;#    0.1*$yq1 + ($yq1-$ystart_axis3)*1.2;
            my $yq2e = (($yt2 - $yq2)/3) + $yq2;#   $yq2 + 0.1*$yq2 + ($yq2-$ystart_axis3)*1.2;

            print SVG "<path d=\"M $xq2 $yq2 L $xq1 $yq1 Q $xt1e $yq1e $xt1 $yt1 L $xt2 $yt2 Q $xt2e $yq2e $xq2 $yq2\" stroke=\"$identity\" stroke-width=\"0.5\" fill=\"$identity\" fill-opacity=\"$op\"/>\n";

         }

      }
   }### end paf output
   close IN;

}else{### cross_match
   open(IN, $xmf) || die "can't read $xmf -- fatal.\n";

   while(<IN>){
      # example cross_match output format
      # 23  0.00 0.00 0.00  2019-nCoV        4    28 (30445)    SARS-CoV    29193 29217 (534)  
      # 660 22.96 1.11 0.95  2019-nCoV       29  3173 (27300)    SARS-CoV       13  3162 (26589)
      # 57 24.55 1.03 1.03  2019-nCoV     3363  3749 (26724)    SARS-CoV     3280  3666 (26085)
      # 23 16.98 0.00 0.00  2019-nCoV     3972  4024 (26449)    SARS-CoV     3883  3935 (25816)
      #7577 17.13 0.29 0.28  2019-nCoV     4120 21589 (8884)    SARS-CoV     4031 21502 (8249)
      #  10  7.69 0.00 0.00  2019-nCoV    12326 12338 (18135)    SARS-CoV    12153 12165 (17586)
      #100 ..   16  0.00 0.00 0.00  MERS-CoV    27772 27788 (2245)  C SARS-CoV   (8130) 21621 21605     
      chomp;
      if(/\s+(\d+\.\d{2})\s+\d+\.\d{2}\s+\d+\.\d{2}\s+(\S+)\s+(\d+)\s+(\d+)\s+\S+\s+(\S+)\s+(\d+)\s+(\d+)\s+/){
         my ($mis,$q,$qs,$qe,$t,$ts,$te) = ($1,$2,$3,$4,$5,$6,$7);
         my $si = 100-$mis;
         $qs/=$scale;
         $qe/=$scale;
         $ts/=$scale;
         $te/=$scale;
         my $identity;

         #identity gradient
         if($si==100){
            $identity = "#005824";
         }elsif($si>=90){
            $identity = "#238b45";
         }elsif($si>=80){
            $identity = "#41ae76";
         }elsif($si>=70){
            $identity = "#66c2a4";
         }elsif($si>=60){
            $identity = "#99d8c9";
         }elsif($si>=50){
            $identity = "#ccece6";
         }elsif($si<50){
            $identity = "#edf8fb";
         }

         if($info->{$q}{'axis'}==1 && $info->{$t}{'axis'}==2){
            ### axis 1
            my $xq1 = $midwidth+5;
            my $xq2 = $xq1;
            my $yq1 = $mid - 50 - $qs;
            my $yq2 = $mid - 50 - $qe;
            ### axis2   
            my $xstart_axis2 = $midwidth + 50;
            my $ystart_axis2 = $mid + 45;
            my $xt1 = (sin(45)*$ts) + $xstart_axis2;
            my $yt1 = (cos(45)*$ts) + $ystart_axis2;
            my $xt2 = (sin(45)*$te) + $xstart_axis2;
            my $yt2 = (cos(45)*$te) + $ystart_axis2;   
  
            my $xt1e = $xt1;# + 0.1*$xt1 + ($xt1-$xstart_axis3)*1.2;
            my $xt2e = $xt2;# + 0.1*$xt2 + ($xt2-$xstart_axis3)*1.2;
            my $yq1e = (($yt1 - $yq1)/3) + $yq1;#    0.1*$yq1 + ($yq1-$ystart_axis3)*1.2;
            my $yq2e = (($yt2 - $yq2)/3) + $yq2;#   $yq2 + 0.1*$yq2 + ($yq2-$ystart_axis3)*1.2;

            print SVG "<path d=\"M $xq2 $yq2 L $xq1 $yq1 Q $xt1e $yq1e $xt1 $yt1 L $xt2 $yt2 Q $xt2e $yq2e $xq2 $yq2\" stroke=\"$identity\" stroke-width=\"0.5\" fill=\"$identity\" fill-opacity=\"$op\"/>\n";
 
         }elsif($info->{$q}{'axis'}==3 && $info->{$t}{'axis'}==2){
            ### axis 2
            my $xstart_axis2 = $midwidth + 50;
            my $ystart_axis2 = $mid + 55;
            my $xq1 = (sin(45)*$ts) + $xstart_axis2;
            my $yq1 = (cos(45)*$ts) + $ystart_axis2;
            my $xq2 = (sin(45)*$te) + $xstart_axis2;
            my $yq2 = (cos(45)*$te) + $ystart_axis2;
            ### axis 3
            my $xstart_axis3 = $midwidth - 50;
            my $ystart_axis3 = $mid + 55;
            my $xt1 = (sin(-45)*$qs) + $xstart_axis3;
            my $yt1 = (cos(-45)*$qs) + $ystart_axis3;
            my $xt2 = (sin(-45)*$qe) + $xstart_axis3;
            my $yt2 = (cos(-45)*$qe) + $ystart_axis3;

            my $yq1e = $yq1 + 0.05*$yq1 + ($yq1-$ystart_axis3)*1.1;
            my $yq2e = $yq2 + 0.05*$yq2 + ($yq2-$ystart_axis3)*1.1;

            print SVG "<path d=\"M $xq2 $yq2 L $xq1 $yq1 Q $mid $yq1e $xt1 $yt1 L $xt2 $yt2 Q $mid $yq2e $xq2 $yq2\" stroke=\"$identity\" stroke-width=\"0.5\" fill=\"$identity\" fill-opacity=\"$op\"/>\n";

         }elsif($info->{$q}{'axis'}==1 && $info->{$t}{'axis'}==3){

            ### axis 1
            my $xq1 = $midwidth-5;
            my $xq2 = $xq1;
            my $yq1 = $mid - 50 - $qs;
            my $yq2 = $mid - 50 - $qe;
            ### axis 3
            my $xstart_axis3 = $midwidth - 50;
            my $ystart_axis3 = $mid + 45;
            my $xt1 = (sin(-45)*$ts) + $xstart_axis3;
            my $yt1 = (cos(-45)*$ts) + $ystart_axis3;
            my $xt2 = (sin(-45)*$te) + $xstart_axis3;
            my $yt2 = (cos(-45)*$te) + $ystart_axis3;

            my $xt1e = $xt1;# + 0.1*$xt1 + ($xt1-$xstart_axis3)*1.2;
            my $xt2e = $xt2;# + 0.1*$xt2 + ($xt2-$xstart_axis3)*1.2;
            my $yq1e = (($yt1 - $yq1)/3) + $yq1 ;#    0.1*$yq1 + ($yq1-$ystart_axis3)*1.2;
            my $yq2e = (($yt2 - $yq2)/3) + $yq2;#   $yq2 + 0.1*$yq2 + ($yq2-$ystart_axis3)*1.2;

            print SVG "<path d=\"M $xq2 $yq2 L $xq1 $yq1 Q $xt1e $yq1e $xt1 $yt1 L $xt2 $yt2 Q $xt2e $yq2e $xq2 $yq2\" stroke=\"$identity\" stroke-width=\"0.5\" fill=\"$identity\" fill-opacity=\"$op\"/>\n";

         }
      ###matches on reverse strand
      }elsif(/\s+(\d+\.\d{2})\s+\d+\.\d{2}\s+\d+\.\d{2}\s+(\S+)\s+(\d+)\s+(\d+)\s+\S+\s+C\s+(\S+)\s+\S+\s+(\d+)\s+(\d+)\s+/){

         my ($mis,$q,$qs,$qe,$t,$te,$ts) = ($1,$2,$3,$4,$5,$6,$7);
         my $si = 100-$mis;
         $qs/=$scale;
         $qe/=$scale;
         $ts/=$scale;
         $te/=$scale;
         my $identity;

         #identity gradient
         if($si==100){
            $identity = "#99000d";
         }elsif($si>=90){
            $identity = "#cb181d";
         }elsif($si>=80){
            $identity = "#ef3b2c";
         }elsif($si>=70){
            $identity = "#fb6a4a";
         }elsif($si>=60){
            $identity = "#fc9272";
         }elsif($si>=50){
            $identity = "#fcbba1";
         }elsif($si<50){
            $identity = "#fee5d9";
         }

         if($info->{$q}{'axis'}==1 && $info->{$t}{'axis'}==2){
            ### axis 1
            my $xq1 = $midwidth+5;
            my $xq2 = $xq1;
            my $yq1 = $mid - 50 - $qs;
            my $yq2 = $mid - 50 - $qe;
            ### axis2   
            my $xstart_axis2 = $midwidth + 50;
            my $ystart_axis2 = $mid + 45;
            my $xt1 = (sin(45)*$ts) + $xstart_axis2;
            my $yt1 = (cos(45)*$ts) + $ystart_axis2;
            my $xt2 = (sin(45)*$te) + $xstart_axis2;
            my $yt2 = (cos(45)*$te) + $ystart_axis2;   

            my $xt1e = $xt1;# + 0.1*$xt1 + ($xt1-$xstart_axis3)*1.2;
            my $xt2e = $xt2;# + 0.1*$xt2 + ($xt2-$xstart_axis3)*1.2;
            my $yq1e = (($yt1 - $yq1)/3) + $yq1;#    0.1*$yq1 + ($yq1-$ystart_axis3)*1.2;
            my $yq2e = (($yt2 - $yq2)/3) + $yq2;#   $yq2 + 0.1*$yq2 + ($yq2-$ystart_axis3)*1.2;

            print SVG "<path d=\"M $xq2 $yq2 L $xq1 $yq1 Q $xt1e $yq1e $xt1 $yt1 L $xt2 $yt2 Q $xt2e $yq2e $xq2 $yq2\" stroke=\"$identity\" stroke-width=\"0.5\" fill=\"$identity\" fill-opacity=\"$op\"/>\n";

         }elsif($info->{$q}{'axis'}==3 && $info->{$t}{'axis'}==2){
            ### axis 2
            my $xstart_axis2 = $midwidth + 50;
            my $ystart_axis2 = $mid + 55;
            my $xq1 = (sin(45)*$ts) + $xstart_axis2;
            my $yq1 = (cos(45)*$ts) + $ystart_axis2;
            my $xq2 = (sin(45)*$te) + $xstart_axis2;
            my $yq2 = (cos(45)*$te) + $ystart_axis2;
            ### axis 3
            my $xstart_axis3 = $midwidth - 50;
            my $ystart_axis3 = $mid + 55;
            my $xt1 = (sin(-45)*$qs) + $xstart_axis3;
            my $yt1 = (cos(-45)*$qs) + $ystart_axis3;
            my $xt2 = (sin(-45)*$qe) + $xstart_axis3;
            my $yt2 = (cos(-45)*$qe) + $ystart_axis3;

            my $yq1e = $yq1 + 0.1*$yq1 + ($yq1-$ystart_axis3)*1.2;
            my $yq2e = $yq2 + 0.1*$yq2 + ($yq2-$ystart_axis3)*1.2;

            print SVG "<path d=\"M $xq2 $yq2 L $xq1 $yq1 Q $midwidth $yq1e $xt1 $yt1 L $xt2 $yt2 Q $midwidth $yq2e $xq2 $yq2\" stroke=\"$identity\" stroke-width=\"0.5\" fill=\"$identity\" fill-opacity=\"$op\"/>\n";

         }elsif($info->{$q}{'axis'}==1 && $info->{$t}{'axis'}==3){

            ### axis 1
            my $xq1 = $midwidth-5;
            my $xq2 = $xq1;
            my $yq1 = $mid - 50 - $qs;
            my $yq2 = $mid - 50 - $qe;
            ### axis 3
            my $xstart_axis3 = $midwidth - 50;
            my $ystart_axis3 = $mid + 45;
            my $xt1 = (sin(-45)*$ts) + $xstart_axis3;
            my $yt1 = (cos(-45)*$ts) + $ystart_axis3;
            my $xt2 = (sin(-45)*$te) + $xstart_axis3;
            my $yt2 = (cos(-45)*$te) + $ystart_axis3;

            my $xt1e = $xt1;# + 0.1*$xt1 + ($xt1-$xstart_axis3)*1.2;
            my $xt2e = $xt2;# + 0.1*$xt2 + ($xt2-$xstart_axis3)*1.2;
            my $yq1e = (($yt1 - $yq1)/3) + $yq1 ;#    0.1*$yq1 + ($yq1-$ystart_axis3)*1.2;
            my $yq2e = (($yt2 - $yq2)/3) + $yq2;#   $yq2 + 0.1*$yq2 + ($yq2-$ystart_axis3)*1.2;

            print SVG "<path d=\"M $xq2 $yq2 L $xq1 $yq1 Q $xt1e $yq1e $xt1 $yt1 L $xt2 $yt2 Q $xt2e $yq2e $xq2 $yq2\" stroke=\"$identity\" stroke-width=\"0.5\" fill=\"$identity\" fill-opacity=\"$op\"/>\n";

         }

      }
   }### end xmv output
   close IN;
}###


#### print legend

my @f=("#005824","#238b45","#41ae76","#66c2a4","#99d8c9","#ccece6","#edf8fb");
my @r=("#99000d","#cb181d","#ef3b2c","#fb6a4a","#fc9272","#fcbba1","#fee5d9");
my @cat=("100","90+","80+","70+","60+","50+","0-49");

my $xt = $xlegend;
my $yt = $ylegend;

my $xtl = $xt+30;
my $xtl2 = $xtl+85;
my $xtl3 = $xtl2+100;
my $xtl4 = $xtl2+50;

my $yt2 = $yt+50;

print SVG "<text font-size=\"2.2em\" x=\"$xt\" y=\"$yt\">Sequence identity (%)</text>\n";
print SVG "<text font-size=\"2em\" x=\"$xtl\" y=\"$yt2\">Forward | Reverse</text>\n";

my $el=0;
foreach my $fc(@f){
   $yt2+=30;
   
   print SVG "<rect x=\"$xtl2\" y=\"$yt2\" width=\"20\" height=\"20\" style=\"fill:$fc;stroke:black;stroke-width:0.5;fill-opacity:$op;\" />\n";
   my $ytt = $yt2+20;
   print SVG "<text font-size=\"2em\" x=\"$xtl3\" y=\"$ytt\">$cat[$el]</text>\n";
   
   print SVG "<rect x=\"$xtl4\" y=\"$yt2\" width=\"20\" height=\"20\" style=\"fill:$r[$el];stroke:black;stroke-width:0.5;fill-opacity:$op;\" />\n";
   $el++;
}

print SVG "</svg>\n";
close SVG;

print "$0 output graph in $out\n";

exit;










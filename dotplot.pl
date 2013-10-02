=header
Draw dot plot figure using the result fo DAGalign or Mcscan 
Write by Chen jinfeng on 20090929
modified by jfchen on 20091009,modify the selfhit for chromosome length 
Uesage: perl dotplot.pl os_sb, need os_sb4dotplot os_sb.chr
=cut

my $help=<<USAGE;
For project os_sb, in chr Rice goes firest than Sorghum.
Rice is Y and Sorghum is X.

Example ob_os4dotplot: reference is Os
## Alignment 0: score=58088.0 e_value=0 N=1197 Ob01&Os01 plus
01      OBR_GLEAN_10017347      15982593        15986184        01      Os01t0565600-03 23254441        23257194         1e-170
01      OBR_GLEAN_10017348      15991732        15993455        01      Os01t0565800-00 23259578        23261341              0
01      OBR_GLEAN_10017349      15995373        15995615        01      Os01t0565900-01 23261862        23263766          8e-95
01      OBR_GLEAN_10017350      15996771        15997730        01      Os01t0566000-01 23264652        23266134         1e-104
01      OBR_GLEAN_10017351      15999050        16000780        01      Os01t0566050-00 23267882        23268891         1e-156
01      OBR_GLEAN_10017356      16050031        16051308        01      Os01t0566800-00 23318479        23319879              0
01      OBR_GLEAN_10017357      16064173        16070973        01      Os01t0566900-01 23329682        23334862              0
01      OBR_GLEAN_10017359      16088275        16090350        01      Os01t0567500-01 23350541        23352698              0
01      OBR_GLEAN_10017382      16246893        16250235        01      Os01t0571000-02 23548869        23552230              0
01      OBR_GLEAN_10017385      16265790        16268707        01      Os01t0571200-01 23564004        23566581         1e-140
01      OBR_GLEAN_10017386      16272408        16274946        01      Os01t0571300-01 23569545        23572477              0
## Alignment 1: score=25704.0 e_value=0 N=538 Ob01&Os01 plus
01      OBR_GLEAN_10003595      15307   22544   01      Os01t0100100-01 1983    9815          0
01      OBR_GLEAN_10003597      47693   50341   01      Os01t0100400-01 11721   14685         0
01      OBR_GLEAN_10003599      65760   69394   01      Os01t0100800-01 28818   33453         0
01      OBR_GLEAN_10003600      70604   77229   01      Os01t0100900-01 34623   40136         0
01      OBR_GLEAN_10003601      83442   85775   01      Os01t0101150-00 57658   60090         0
01      OBR_GLEAN_10003602      88006   89105   01      Os01t0101200-01 61060   62576         0
01      OBR_GLEAN_10003604      92779   96944   01      Os01t0101600-01 71816   77349         0


 
Example ob_os.chr: reference is Os
Ob	Length	Os	Length
chr01	33916305	chr01	45038604
chr02	27085147	chr02	36792247
chr03	29536849	chr03	37312367
chr04	21479432	chr04	36060865
chr05	20170002	chr05	30073438
chr06	21537212	chr06	32124789
chr07	18473044	chr07	30357780
chr08	18639272	chr08	28530027
chr09	13736722	chr09	23895721
chr10	14035047	chr10	23703430
chr11	15914883	chr11	31219694
chr12	14847084	chr12	27679166


Run: perl boxplot.pl os_sb > log &
USAGE

if (@ARGV < 1){
   print "$help\n";
   exit ();
}

use strict;
use SVG;

our $svg=SVG-> new (width=>500,heith=>1200);
our $xstart=50;
our $xend=450;
our $ystart=50;
our $yend=350;
our $yheight=$yend-$ystart;
our $xlength=$xend-$xstart;
my ($spec1,$spec2)=split("_",$ARGV[0]);
my $infile1="$ARGV[0]"."4dotplot";
my $infile2="$ARGV[0].chr";
######draw x,y lines for figure and write the note

my $rec=$svg->rectangle(
         x=>$xstart,y=>$ystart,
         width=>$xlength,height=>$yheight,
         style=>{
             stroke=>'black',
             fill=>'white',
             'stroke-width'=>0.7
         } 
 
);######draw the lines

my %os2sb;###hash for inter species chromosome relation
my %sb2os;
my %oschr;##store chr length inhash
my %sbchr; 
my $oslength;
my $sblength;
open CHR, "$infile2" or die "can not open my chr file";
my @head=split("\t",<CHR>);
my $yname=$head[0];
my $xname=$head[2];
print "$yname\t$xname\n";
while (<CHR>){
         chomp $_;
         my @unit=split("\t",$_);
         if ($unit[0]=~/(\d+)/){$unit[0]=$1};
         if ($unit[2]=~/(\d+)/){$unit[2]=$1};
         unless(exists $os2sb{$unit[0]}){$os2sb{$unit[0]}=$unit[2] };
         unless(exists $sb2os{$unit[2]}){$sb2os{$unit[2]}=$unit[0] }; 
         unless(exists $oschr{$unit[0]}){
             $oschr{$unit[0]}=$unit[1];
             #print "$unit[0]\n";
             $oslength+=$unit[1];   
         }
         unless(exists $sbchr{$unit[2]}){
             $sbchr{$unit[2]}=$unit[3];
             #print "$unit[2]\n";
             $sblength+=$unit[3];
         }
         
}
close CHR;
our $osratio=$yheight/$oslength; ## qry, not os as it seems
our $sbratio=$xlength/$sblength; ## ref, not sb as it seems

#print "OS\t$osratio\nSB\t$sbratio\n";
###read in chrinfor and store infro in hash and scalar

##draw lines seperate for chromosome
my @osstart; ## this array will store the start x of every chromomore for y axis
push (@osstart,$ystart); 
my $yv=$ystart;
my @array1=sort keys %oschr;
my $narray1=@array1;
my $lastkeyy=pop @array1;
my $lastyvy=$yend-$oschr{$lastkeyy}*$osratio/2+2;
my $lastyvx=$xstart-10;
my $lastnotey="C$narray1";
my $ylastnote =$svg->text(
           x=>$lastyvx, y=>$lastyvy,
           style=>{'text-anchor'=>'end','font-size'=>'xx-small','font-weight'=>100,'stroke-width'=>0.1}     
)->cdata($lastnotey); ###write note for the last chromosome

my $countery;
foreach (@array1) {
    $countery++;
    my $notey="C$countery";  
    #print "$_"; 
    $yv+=$oschr{$_}*$osratio; 
    push (@osstart,$yv); ## store start x for every chromosome, except last one 
    #print "$yv\n";
    my $yvy=$yv-$oschr{$_}*$osratio/2+2;
    my $yvx=$xstart-10;
    my $ynote =$svg->text(
           x=>$yvx, y=>$yvy,
           style=>{'text-anchor'=>'end','font-size'=>'xx-small','font-weight'=>100,'stroke-width'=>0.1}     
    )->cdata($notey); ###write note for each chromosome except for the last
    my $yy=$yheight/2+$ystart;
    my $yx=$xstart-35;
    my $yano=$svg->text(
       x=>$yx,y=>$yy,
       style=>{'font-size'=>'small','font-style'=>'italic','font-weight'=>100,'stroke-width'=>0.1},
       transform=>"rotate(-90,$yx,$yy)"
    ) ->cdata($yname);## write y anotation for y axis
    my $yline=$svg->line(
           x1=>$xstart, y1=>$yv,
           x2=>$xend,   y2=>$yv,  
           style=>{stroke=>'red','stroke-width'=>0.4}  
    );##draw internal line to seperate chromosome
}

###draw line for x and write note, similar with y
my @sbstart;
push (@sbstart,$xstart);
my $xv=$xstart;
my @array2=sort keys %sbchr;
my $narray2=@array2;
my $lastkeyx=pop @array2; 
my $lastxvy=$yend+20;
my $lastxvx=$xend-$sbchr{$lastkeyx}*$sbratio/2-2;
my $lastnotex="C$narray2";
my $xlastnote=$svg->text(
           x=>$lastxvx, y=>$lastxvy,
           style=>{'font-size'=>'xx-small','font-weight'=>100,'stroke-width'=>0.1}
)->cdata($lastnotex);
my $counterx;
foreach (@array2){
    $counterx++;
    my $notex="C$counterx";
    #print "$_\n"; 
    $xv+=$sbchr{$_}*$sbratio;  
    push (@sbstart,$xv); 
    #print "Sb\t$sbchr{$_}\n";
    #print "P\t$xv\n";
    my $xvx=$xv-$sbchr{$_}*$sbratio/2-2;
    my $xvy=$yend+20;
    my $xnote =$svg->text(
           x=>$xvx, y=> $xvy,
           style=>{'font-size'=>'xx-small','font-weight'=>100,'stroke-width'=>0.1}
    )->cdata($notex);
    my $xx=$xlength/2+$xstart-5;
    my $xy=$yend+40;
    my $xano=$svg->text(
         x=>$xx, y=>$xy,
         style=>{fill=>'black','font-size'=>'small','font-style'=>'italic','font-weight'=>100,'stroke-width'=>0.1},
    )->cdata($xname);
    my $xline=$svg->line(
           x1=>$xv,  y1=>$ystart,
           x2=>$xv,  y2=>$yend,
           style=>{stroke=>'red','stroke-width'=>0.4}
    );
}
#####

#######read in os_sb4dotplot and draw dot in each cell###
$/="##";
my $dotfile="$ARGV[0]"."4dotplot";
open DOT, "$dotfile" or die "can not open my dotfile";
while(<DOT>){
     my @unit=split("\n",$_);
     shift @unit;
     $unit[@unit]=~s/#//g;
     my @word=split("\t",$unit[0]);
     my $sub1=substr($word[1],0,2);
     my $sub2=substr($word[5],0,2);
     my $chromosome1=$word[0];
     my $chromosome2=$word[4];
     #print "$sub1\t$sub2\n"; 
     unless($sub1 eq $sub2){ ##if the blast result is from a self comparition,us if replace unless and delete the elsif steps 
           foreach(@unit){
                 unless($_=~/\w+/){next};
                 my @para=split("\t",$_);
      
                 my ($color,$x1, $x2, $y1, $y2);
                 if ($para[2] < $para[3]){
                     $color = '#DAA520';
                 }else{
                     $color = 'blue';
                 }
                 $y1 = $para[2] * $osratio + $osstart[$para[0]-1];
                 $y2 = $para[3] * $osratio + $osstart[$para[0]-1];               
                 $x1 = $para[6] * $sbratio + $sbstart[$para[4]-1];
                 $x2 = $para[7] * $sbratio + $sbstart[$para[4]-1];
                 my $line =$svg->line(
                         x1=>$x1,  y1=>$y1,
                         x2=>$x2,  y2=>$y2,
                         style=>{stroke=>$color,'stroke-width'=>0.4}
                 );
=pod
                 my $dotx;
                 if ($para[6]>$para[7]){
                     $dotx=$sbratio*(($para[6]-$para[7])/2+$para[7]);
                 }else{
                     $dotx=$sbratio*(($para[7]-$para[6])/2+$para[6]);
                 } 
                 my $doty;
                 if ($para[2]>$para[3]){
                     $doty=$osratio*(($para[2]-$para[3])/2+$para[3]);
                 }else{
                     $doty=$osratio*(($para[3]-$para[2])/2+$para[2]);
                 }
                 $dotx+=$sbstart[$para[4]-1];
                 $doty+=$osstart[$para[0]-1];
                 my $dot=$svg->circle(cx=>$dotx,cy=>$doty,r=>0.2,'fill'=>'#DAA520');
=cut
           }
=pod
     }elsif($sub1=~/$spec1/i){
          
          foreach(@unit){
                 unless($_=~/\w+/){next};
                 my @para=split("\t",$_);
                 my $dotx;
                 if ($para[6]>$para[7]){
                     $dotx=$sbratio*(($para[6]-$para[7])/2+$para[7]);
                 }else{
                     $dotx=$sbratio*(($para[7]-$para[6])/2+$para[6]);
                 }
                 my $doty;
                 if ($para[2]>$para[3]){
                     $doty=$osratio*(($para[2]-$para[3])/2+$para[3]);
                 }else{
                     $doty=$osratio*(($para[3]-$para[2])/2+$para[2]);
                 }
				 #print "$_\n";
                 #print "PARA4\t$para[4]\t$oschr{$para[4]}\n";
                 
                 my $rx=$sbchr{$os2sb{$para[4]}}/$oschr{$para[4]}; 
                 $dotx=$dotx*$rx;
                 $dotx+=$sbstart[$os2sb{$para[4]}-1];
                 
                 $doty+=$osstart[$para[0]-1];
                 my $dot=$svg->circle(cx=>$dotx,cy=>$doty,r=>0.2,'fill'=>'black');
          } 
     }elsif($sub1=~/$spec2/i){
          foreach(@unit){
	         unless($_=~/\w+/){next};
                 my @para=split("\t",$_);
                 my $dotx;
                 if ($para[6]>$para[7]){
                     $dotx=$sbratio*(($para[6]-$para[7])/2+$para[7]);
                 }else{
                     $dotx=$sbratio*(($para[7]-$para[6])/2+$para[6]);
                 }
                 my $doty;
                 if ($para[2]>$para[3]){
                     $doty=$osratio*(($para[2]-$para[3])/2+$para[3]);
                 }else{
                     $doty=$osratio*(($para[3]-$para[2])/2+$para[2]);
                 }
                 
                 $dotx+=$sbstart[$para[4]-1];
                 
				 
                 $doty=$doty*$oschr{$sb2os{$para[0]}}/$sbchr{$para[0]};
                 $doty+=$osstart[$sb2os{$para[0]}-1];
                 
                 my $dot=$svg->circle(cx=>$dotx,cy=>$doty,r=>0.2,'fill'=>'red');
          }
=cut
     } 




}

close DOT;
$/="\n";
######
my $outfile="$ARGV[0].svg";
writefile($outfile);



####sub writefile#####
sub writefile{
my ($file)=@_;
open OUT, ">$file";
print OUT $svg->xmlify;
close OUT;
system "/rhome/cjinfeng/software/tools/draw/svg2xxx_release/svg2xxx $file -t pdf";
#system "/home/jfchen/bgitraining/draw/svg2xxx_release/svg2xxx $file -t png";
}

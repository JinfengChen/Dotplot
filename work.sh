mummer -mum -b -c ../input/seq1.fa ../input/seq2.fa > mummer.mums
mummerplot -x "[0,275287]" -y "[0,265111]" -postscript -p mummer mummer.mums

echo "dotplot gene-based"
perl /rhome/cjinfeng/BigData/00.RD/ChromosomeAlign/bin/dotplot.pl ob_os

echo "sequence based"
mummer -mum -l 200 -b -c ../input/O.glaberrima.v1.0.fa ../input/MSU_r7.fa > og_os.mummer.mums 2> log2
python Mummer2dotplot.py --input og_os.mummer.mums --single 0 > og_os4dotplot
perl /rhome/cjinfeng/BigData/00.RD/ChromosomeAlign/bin/dotplot.pl og_os

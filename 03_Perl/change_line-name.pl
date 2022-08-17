#!/usr/bin/perl -w
use strict;
unless (@ARGV==2) {
    die"Usage: perl $0 <input.fa>  <output.fa> \n";
    }

   my ($infile,$outfile) = @ARGV;
open IN,$infile || die"error: can't open infile: $infile";
open OUT,">$outfile" || die$!;

    while(<IN>){
         chomp;
         if($_=~/^>KYQ/){
     
       s/(KYQ[0-9]+\.[0-9]).*$/$1/; #用括号先把整个匹配出来后，$1替换时，保留第一个括号里面的内容，如有2个括号则用$2表示,此时保留第二个括号内容。

       }
    print OUT " $_\n" ;
}
close IN;
close OUT; 








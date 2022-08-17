#!/usr/bin/perl
use strict;
use Getopt::Long;

my $usage = <<USAGE;
Usage:
    $0 [options] input.phy > input.nex

    --data-type <string>    defaut: PROTEIN
    设置数据类型。PROTEIN、DNA等。

USAGE
if (@ARGV==0){die $usage}

my $dataType;
GetOptions(
    "data-type:s" => \$dataType,
);
$dataType ||= "PROTEIN";

open IN, $ARGV[0] or die "Can not open file $ARGV[0], $!\n";
$_ = <IN>;
chomp;
@_ = split /\s+/;
print "#NEXUS\nBEGIN DATA;\n\tDIMENSIONS NTAX=$_[0] NCHAR=$_[1];\n\tFORMAT DATATYPE=PROTEIN MISSING=? GAP=-;\n\tMATRIX\n";
while (<IN>) {
    print;
}
print ";\nEND;\n";

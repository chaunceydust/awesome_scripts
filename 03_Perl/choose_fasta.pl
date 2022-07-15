# perl choose_fasta.pl gene_id.list gene_fasta 
open FH, $ARGV[0];
while(<FH>){
    chomp;
    $_=">".$_ if $_ !~ /^>/;
    $h{$_}=1;
}
close FH;

open FH, $ARGV[1];
$sign=0;
while(<FH>){
        chomp;
        if(/^>/){
            @F=split;
            exists $h{$F[0]} ? ($sign=1):($sign=0);
        }
        print "$_\n" if $sign==1;
}
close FH;
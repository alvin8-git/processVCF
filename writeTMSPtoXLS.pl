#!/usr/bin/perl -w

use strict;
use Excel::Writer::XLSX;
use File::Basename;
 
my ($file)=@ARGV;
my $sample = basename("$file",".vcf");

open(FILTER, "$sample.Filter.txt" ) or die "$sample.Filter.txt: $!";
open(ANNOTATION, "$sample.annotation.txt" ) or die "$sample.annotation.txt: $!";
open(CHECK, "$sample.annoCheck.txt" ) or die "$sample.annoCheck.txt: $!";
open(COMPARE, "$sample.compare.txt" ) or die "$sample.compare.txt: $!";
open(COMPAREVAF, "$sample.compareVAF.txt" ) or die "$sample.compareVAF.txt: $!";

#rename just in case
my $sampleName=$sample=~s/^(.{20}).*/$1/gr;

my $workbook  = Excel::Writer::XLSX->new( "$sample.xlsx" );
my $filter = $workbook->add_worksheet("$sampleName.filter");
my $annotation = $workbook->add_worksheet("$sampleName");
my $check = $workbook->add_worksheet("$sampleName.check");
my $compare = $workbook->add_worksheet("$sampleName.comp");
my $compareVAF = $workbook->add_worksheet("$sampleName.compVF");
 
# Row and column are zero indexed
my $a = 0;
 
while ( <FILTER> ) {
    chomp;
 
    # Split on single tab
    my @fields = split( '\t', $_ );
 
    my $b = 0;
    for my $c ( @fields ) {
        $filter->write( $a, $b, $c );
        $b++;
    }
    $a++;
}

my $d =0;
while ( <ANNOTATION> ) {
    chomp;
 
    # Split on single tab
    my @fields = split( '\t', $_ );
 
    my $e = 0;
    for my $f ( @fields ) {
        $annotation->write( $d, $e, $f );
        $e++;
    }
    $d++;
}

my $g =0;
while ( <CHECK> ) {
    chomp;

    # Split on single tab
    my @fields = split( '\t', $_ );

    my $h = 0;
    for my $i ( @fields ) {
        $check->write( $g, $h, $i );
        $h++;
    }
    $g++;
}

my $j =0;
while ( <COMPARE> ) {
    chomp;

    # Split on single tab
    my @fields = split( '\t', $_ );

    my $k = 0;
    for my $l ( @fields ) {
        $compare->write( $j, $k, $l );
        $k++;
    }
    $j++;
}

my $m =0;
while ( <COMPAREVAF> ) {
    chomp;

    # Split on single tab
    my @fields = split( '\t', $_ );

    my $n = 0;
    for my $o ( @fields ) {
        $compareVAF->write( $m, $n, $o );
        $n++;
    }
    $m++;
}

 
$workbook->close();

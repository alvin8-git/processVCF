#!/usr/bin/perl -w

#Usage Write1WStoXLS.pl <XLSname> <WS1 data>

use strict;
use Excel::Writer::XLSX;
use File::Basename;
 
my ($XLS,$file)=@ARGV;
my $data = basename("$file",".txt");

open(WS1, "$data.txt" ) or die "$data.txt: $!";
 
my $workbook  = Excel::Writer::XLSX->new( "$XLS.xlsx" );
my $worksheet1 = $workbook->add_worksheet("$data");
 
# Row and column are zero indexed
my $a = 0; 
while ( <WS1> ) {
    chomp;
 
    # Split on single tab
    my @fields = split( '\t', $_ );
 
    my $b = 0;
    for my $c ( @fields ) {
        $worksheet1->write( $a, $b, $c );
        $b++;
    }
    $a++;
}
 
$workbook->close();

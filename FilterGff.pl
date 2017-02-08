#!/usr/bin/perl
#Takes a GFF file, removes entries based on filter
#
# Usage: perl filterGFF.pl [options] my_gff_all.gff3
#
# options:
#    -f filter criteria with regex capture
#       e.g. 'alnlen=(\d+\.\d)'
#
#    -v filter value
#
#Travis Saari 1-20-2017

use warnings;
use strict;
use Getopt::Std;

my %options = ();
getopts('f:v:', \%options);

my $gffFile = $ARGV[0];

my $entry_count = 0; #Initialize count of GFF entries
my @entrylines; #Initialize array for entry lines

my $startCapture = '\tgene\t';
my $endCapture = '^###$';

open (GFF, "<$gffFile") or die "Can't open $gffFile!: $!\n";

while(<GFF>) {
    #Use flip-flop operator to push lines between & including the regular expressions
    if (my $rc = /$startCapture/../$endCapture/) {
        
        push @entrylines, $_;
        
        if ($rc =~ /E\d+/) { 
        #flip-flop operator trick. Do the following after last entryline
            my $entry_slurp = join('', @entrylines);
			$entry_slurp =~ /$options{'f'}/;
			my $captured_val = $1;
			if ($captured_val >= $options{'v'}) {
				print $entry_slurp;
			}
            @entrylines=(); #Clear array for next entry
        }
    }
}

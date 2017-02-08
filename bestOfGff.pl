#!/usr/bin/perl
#Takes a GFF file, outputs a GFF file with one entry (best score) for each ID
#
# usage: perl bestOfGff.pl [options] my_gff_all.gff3
#
# options:
#    -s regex found at start of entry 
#        Defaults to '\tgene\t'
#    -e regex found at end of entry
#        Defaults to '^###$'
#    -c regex to capture criterion (must capture numeric expression for comparison to work)
#        e.g. 'score=(\d+)'
#        Defaults to contents of 6th column of parent feature GFF line
#
#Travis Saari 4-22-2016

use warnings;
use strict;
use Getopt::Std;

my $gffFile = $ARGV[0];

my %options = ();
getopts('s:e:c:', \%options); # For customizing start, end, and criterion reg. expns.

my %scores;
my %lines;
my @entrylines;
my $entry_count = 0;

#Defaults for regex capture
my $startCapture = '\tgene\t';
my $endCapture = '^###$';
my $criterion = '';

if ($options{'s'}) {
    $startCapture = $options{'s'};
}
if ($options{'e'}) {
    $endCapture = $options{'e'};
}
if ($options{'c'}) {
    $criterion = $options{'c'};
}

open (GFF, "<$gffFile") or die "Can't open $gffFile!: $!\n";

#Step through the GFF file, process one entry at a time
while(<GFF>) {
    #Use flip-flop operator to push lines between & including the regular expression
    if (my $rc = /$startCapture/../$endCapture/) {

        push @entrylines, $_;
        if ($rc =~ /E\d+/) {
        #flip-flop operator trick. Do the following after last entryline
            $entry_count++;
            process_entry(@entrylines);
            @entrylines=(); #Clear array for next entry
        }
    }
}


foreach my $key (sort keys %lines) {
    print $lines{$key};
}


#sub process_entry($firstline, $lastline, @lines) {


# Description: For each entry, take score and ID, and update hash if score
# is greater than the current value.
# Inputs: Array of GFF lines 
# Outputs: None
# Side-effects: Update %scores and %lines
sub process_entry {
    #Exctract score and ID from lines
    my $parent = $_[0];
    my $score = 0;
    my $ID = "";
    if ($criterion) {
        if ($parent =~ /$criterion/){
            $score = $1;
        }
        else {
            warn "Criterion provided \'$criterion\' not found. Skipping\n";
        }
    }
    else {
        my @parentAttr = split("\t", $parent);
        $score = $parentAttr[5];
    }
    if ($parent =~ /ID=(.*?);/){
        $ID=$1;
        $ID =~ s/(\.\d+)\.\d+/$1/;
    }

    #Join the lines, correct the ID's if q
    my $str = join ("", @_);
    $str =~ s/$ID\.\d+/$ID/g;
    
    #If ID in hash, overwrite old entries if new score is higher.
    if($scores{$ID}){
        if($score > $scores{$ID}){
            #print STDERR "Old score: $scores{$ID}\nNew score: $score\n"; #TWS DEBUG
            $scores{$ID} = $score;
            $lines{$ID} = $str;
        }
    }
    #Otherwise, add the entry to the hash
    else {
        $scores{$ID} = $score;
        $lines{$ID} = $str;
    }

}

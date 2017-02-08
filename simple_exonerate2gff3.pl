#!/usr/bin/perl
#
# Convert exonerate output (with gff enabled) to gff3 format
#
# usage: simple_exonerate2gff3.pl [options] file1 file2 file3 > outfile.gff3
#
# options:
#    -p Parent feature 
#        Defaults to 'gene'
#        e.g. 'similarity'
#    -r Regex for ID capture
#        Defaults to 'Query (.*?) ;'
#        e.g. '\tgene\t.*sequence (.*?) ;'
#    -P Process RYO entry? boolean (Have to manually edit this program to process it how you wish)
#        no error handling on this part - user discretion is advised
#    -R RYO string for capture - instead of ending after GFF dump, pickup until specified string
#        e.g. '^---END RYO---'
#        no error handling on this part - user discretion is advised
#
# 6-30-2016 Travis Saari

use warnings;
use strict;
use Getopt::Std;

my %options = ();
getopts('p:r:R:P', \%options);


my $entry_count = 0; #Initialize count of GFF entries
my @entrylines; #Initialize array for entry lines

my $startCapture = '^# --- START OF GFF DUMP ---';
my $endCapture = '^# --- END OF GFF DUMP ---';

if ($options{'R'}) {
	$endCapture = $options{'R'};
}

while(<>) {
    #Use flip-flop operator to push lines between & including the regular expressions
    if (my $rc = /$startCapture/../$endCapture/) {
        
        push @entrylines, $_;
        
        if ($rc =~ /E\d+/) { 
        #flip-flop operator trick. Do the following after last entryline
            $entry_count++;
            my $gff3_entry = process_entry(@entrylines);
            print "$gff3_entry\n###\n" if $gff3_entry;
            @entrylines=(); #Clear array for next entry
        }
    }
}

################################
# Subroutines
################################

sub process_entry {
    my @lines = @_;
	my @other_lines; #Initialize array for other lines (e.g. RYO)
    my @lines_new; #Initialize array for processed lines
    
    #Set up for extracting ID by regex
    my $joined_lines = join('', @lines);
    my $Regex_ID = 'Query (.*?) ;' ; # Default regex
    if ( $options{'r'} ) {
        $Regex_ID = $options{'r'}; # User-def regex
    }

    #Extract ID from joined lines using regex
    my $ID;
    $joined_lines =~ m/$Regex_ID/ ;
    if (not (defined $1)) {
        print STDERR "Regex string failed to match an ID in entry $entry_count. Defaulting to \'entry_$entry_count\'\n";
        $ID = "entry_$entry_count";
    } else {
        $ID = $1; #Captured ID
    }

    #Set up for determining parent-child relationship
    my $parent = 'gene'; # Default parent
    if ( $options{'p'} ) {
        $parent = $options{'p'}; # User-def parent
    }

    #Convert GFF lines
    my $parent_found = 0; #Keep track of parent features found in gff lines

	#Keep track of GFF line number where parent feature is found
	my $current_gff_line = 0; 
	my $parent_line_no;
    foreach my $line (@lines) {

        my @attr = split("\t", $line);
        if ($line !~ /^#/ && scalar(@attr) == 9) { #Only GFF lines
            my $prepend; # Initialize variable for if/else result
            if ($attr[2] eq $parent) {
                $parent_found++; #Parent found
				$parent_line_no = $current_gff_line;
                $prepend = "ID=" . $ID . ";";
            } else {
                $prepend = "Parent=" . $ID . ";";
            }

            #Modify 9th column, then push result to @lines_new
            $attr[8] = process_gff_tags($prepend, $attr[8]);
            push(@lines_new, join("\t", @attr));

			$current_gff_line++;
        }
		else {
			push(@other_lines, $line);
		}
    }

	#If boolean flag to process RYO
	if ($options{'P'}) {
		my $newTags = process_RYO_to_tags(@other_lines);
		#append new tags to parent line
		$lines_new[$parent_line_no] .= ";$newTags";
	}

    #If exactly one parent feature found in entry, return processed entry
    if ($parent_found == 1) {
        my $processed_entry = join("\n", @lines_new);
        return $processed_entry;
    } else { #Return empty otherwise
        print STDERR "Parent feature \'$parent\' found $parent_found times in entry: $ID. Skipping...\n";
        return;
    }


}

# Trim whitespace, use equal sign for tag-value pairing,
# and prepend something onto the tag string (e.g. ID=Foo;)
sub process_gff_tags {
    my ($prepend, $tag_str) = @_;
    my @tags_orig = split(";", $tag_str);
    my @tags_new = ();
    foreach my $tag (@tags_orig) {
        $tag = trim($tag);
        $tag =~ s/ /=/;
        push(@tags_new, $tag);
    }
    $tag_str = join(';', @tags_new);
    $tag_str = $prepend . $tag_str;
    return $tag_str;
}

# This is specific to my RYO string and my desired use case
# It is kind of hack-y
# Edit this subroutine to change the output
sub process_RYO_to_tags {
	my @allLines = @_;
	my @RYOlines; #Initialize variable for RYOlines
	foreach(@allLines) {
		push (@RYOlines, $_) if $_ =~ /^>/;
	}
	my $subjectLine = $RYOlines[0];
	my @attr = split(' ', $subjectLine);
	my @keepThese = splice(@attr,1,2);
	my $len = $keepThese[0];
	$len =~ s/.*?=//;
	my $alnlen = $keepThese[1];
	$alnlen =~ s/.*?=//;
	my $pct_aln = $alnlen/$len*100;
	$pct_aln = sprintf("%.1f", $pct_aln);
	substr($pct_aln, 0, 0) = 'aln_pct='; #prepend 'aln-pct' to value
	$pct_aln .= ';'; #append semicolon
	push(@keepThese, $pct_aln);
	my $newTags = join(';', @keepThese);
	return $newTags;
}

# Simple regex trim whitespace
sub trim {
    my $s = shift; 
    $s =~ s/^\s+|\s+$//g; 
    return $s;
}

#!/usr/bin/perl
#
# Runs BLAST to find hits, then uses exonerate to perform spliced alignment
#
# usage: exonerate_afterblast.pl -c config.txt -q /path/to/queries/ -t /path/to/targets/ 
# 9-1-2016 Travis Saari

use warnings;
use strict;
use Getopt::Std;

use IPC::Run qw (run timeout);
use File::Basename;
use Data::Dumper;

my $startTime = time();

my %options = ();
getopts('c:q:t:', \%options);

#Load in configuration options from a file
my ($blastopts_ref, $exonopts_ref); #scope vars
if ($options{'c'}) {
    ($blastopts_ref, $exonopts_ref) = load_config($options{'c'});
} else { die "You must specify a configuration file. Exiting...\n"; }

#Get array of filenames in query dir
my @queries; #Scope var
if ($options{'q'}) {
    @queries = load_dir($options{'q'});
} else { die "You must specify a directory of query files. Exiting...\n"; }

#Get hash of filenames in target dir
my %targets; #Scope var
if ($options{'t'}) {
    my @target_arr = load_dir($options{'t'});
    foreach (@target_arr) {$targets{$_} = 1;}
} else { die "You must specify a directory of target files. Exiting...\n"; }

foreach (@queries) {
    my $queryFileName = $options{'q'} . $_;
    my @blast_results =  eval { blastJob($queryFileName, $blastopts_ref) };
    if ($@) {
        warn "Blast job failed: $@";
    } else {
        foreach (@blast_results) {
            my $targetFileName = check_for_target($_);
            if ($targetFileName) {
                $targetFileName = $options{'t'} . $targetFileName; #prepend directory

                my $result = exonerateJob($queryFileName, $targetFileName, $exonopts_ref) ;

                my ($Qname, $Qpath, $Qsuffix) = fileparse($queryFileName);
                my ($Tname, $Tpath, $Tsuffix) = fileparse($targetFileName);

                my $outFN = $Qname . ".exonerate_out";
                if ($result){
                    if (open(OUTFH, ">>$outFN")) {
                        print OUTFH $result;
                    } else { 
                        warn "Can't open $outFN: $!";
                    }
                }

            }
        }
    }
}

my $endTime = time();
my $timeElapsed = $endTime - $startTime;
print "exonerate_afterblast.pl took $timeElapsed seconds to complete\n";

#########################################
# Subroutines
#########################################
# Load in config file. Each line after a header (e.g. "### Blast Opts ###")
# contains an option for the corresponding program in the header (blast or exonerate)
sub load_config {
    my $config_file = shift;

    open( CONFIG, "<", $config_file ) or die "Can't open $config_file: $!\n";

    my $which_program = "";
    my @blast_opts;
    my @exonerate_opts;

    #Load config options into appropriate array based on header
    while (<CONFIG>) {
        if ($_ =~ /^###(.*?)###/i) {
            my $head_text = $1;
            if ($head_text =~ m/Blast/i) {
                $which_program = "blast";
            } elsif ($head_text =~ m/Exonerate/i) {
                $which_program = "exonerate";
            } else {
                die "Error parsing config file.\nCheck the following line for the words \'Blast\' or \'Exonerate\':\n$_\nExiting...\n";
            }
        } else {
            next if $_ =~ m/^\s+$/; # Skip empty lines

            if ($which_program eq "blast") {
                push(@blast_opts, $_);
            } elsif ($which_program eq "exonerate") {
                push(@exonerate_opts, $_);
            } else {
                die "Error parsing config file.\nIncorrectly formatted header.\nExiting...\n";
            }
        }
    }

    #Clean comments and leading/trailing whitespace from options
    my @cleaned_blastOpts = clean_opts(@blast_opts);
    my @cleaned_exonOpts = clean_opts(@exonerate_opts);

    return (\@cleaned_blastOpts, \@cleaned_exonOpts);
}

#Get a list of filenames in a directory
sub load_dir {
    my $dir_path = shift;
    opendir( DIR, $dir_path ) or die "can't open directory $dir_path";
    my @files = grep { !/^\./ } readdir(DIR);
    return @files;    
}

sub blastJob {
    my ($query, $blastopts_ref) = @_;
    #De-reference blast options
    my @blastCmd = @$blastopts_ref;
    my ($blastOut, $err); #for blast output
    #Add query information after first blast option
    splice(@blastCmd, 1, 0, ("-query", $query));
    my ($Qname, $Qpath, $Qsuffix) = fileparse($query);
    print "Running $blastCmd[0]: query $Qname...\n";
    run \@blastCmd, \undef, \$blastOut, \$err;
    if ($err) {
        warn "Error in BLAST query $Qname. $err\n";
    }
    my @results = split("\n", $blastOut);
    return uniq(@results);
}

sub exonerateJob {
    my ($query, $target, $exonopts_ref) = @_;
    #De-reference exonerate options
    my @exonCmd = @$exonopts_ref;
    my ($exonOut, $err); #for exonerate output
    #Add program, query, and target information to exonerate options
    unshift (@exonCmd, ("exonerate", "-q", $query, "-t", $target));
    my ($Qname, $Qpath, $Qsuffix) = fileparse($query);
    my ($Tname, $Tpath, $Tsuffix) = fileparse($target);
    eval {
        print "Running exonerate: query $Qname, target $Tname...\n";
        run \@exonCmd, \undef, \$exonOut, \$err, timeout(240);
        if ($err) {
            warn "Error in exonerate query $Qname, target $Tname. $err\n";
        }
    };
    if ($@ =~ /timeout/) {
        warn "Error: timeout in exonerate query $Qname, target $Tname\n";
    }

    return $exonOut;
}

sub check_for_target {
    my $targetName = $_;
    my $regex = qr/\Q$targetName\E/;
    my ($file_match) = grep { $_ =~ $regex } keys %targets;
    if ($file_match) {
        return $file_match;
    }
    else {
        warn "Error: Target file not found: $targetName. Can't run exonerate";
        return;
    }
}


# Clean options by removing comments and trimming whitespace
# Afterwards, split on first space
sub clean_opts {
    my @opts = @_;
    my @cleaned_opts;
    foreach my $opt (@opts) {
        $opt =~ s/(.*?)#.*/$1/;
        $opt = trim($opt);
        my @opt_pair = split(" ", $opt, 2);
        push(@cleaned_opts, @opt_pair);
    }
    return @cleaned_opts;
}

# Simple regex trim leading/trailing whitespace
sub trim {
    my $s = shift; 
    $s =~ s/^\s+|\s+$//g; 
    return $s;
}

#Simple test for unique elements in array
sub uniq {
    my %seen;
    return grep { !$seen{$_}++ } @_;
}

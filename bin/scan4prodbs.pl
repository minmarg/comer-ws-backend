#!/usr/bin/env perl
BEGIN {$^W=1}

## (C) 2020 Mindaugas Margelevicius, Institute of Biotechnology, Vilnius University
## Scanner for profile databases; the name of the directory scanned is read 
## from the configuration; appropriate links to profile dbs found are created.

use strict;
use FindBin;
use lib "$FindBin::Bin";
use readcfg;
use File::Spec;
use File::Basename;
use Getopt::Long;
use POSIX qw(strftime);

##constants:
my  $MYPROGNAME = basename($0);
my  $CFGFILE = File::Spec->catfile($FindBin::Bin,File::Spec->updir(),"var","comer-ws-backend.conf");
my  $devnull = File::Spec->devnull();
my  $dbsfx = '__comer2_db';
my  $cotherdbsfx = '__cother_db';

my  $usage = <<EOIN;

Scanner for profile databases. The name of the directory scanned is read 
from the configuration; appropriate links to profile dbs found are created.
(C) 2020 Mindaugas Margelevicius, Institute of Biotechnology, Vilnius University

Usage:
$0 [<Options>]

Options:

--help           This help text.

EOIN


my  $Fail = 0;

my  $result = GetOptions(
               'help|h'     => sub {print $usage; exit(0);}
);

print(STDERR "\n\n".GetDatetime()."\n\n");

unless(-f $CFGFILE) {
    print(STDERR "ERROR: $MYPROGNAME: Config file not found: '$CFGFILE'\n");
    exit(1);
}

## =============================================================================

my  $cfgvar = readcfg->new($CFGFILE);
my  %optionvalues;

unless(-d $cfgvar->PathComerDb_PDB()) {
    print(STDERR "ERROR: $MYPROGNAME: Database directory not found: '".$cfgvar->PathComerDb_PDB()."'\n");
    exit(1);
}
unless(-d $cfgvar->PathComerDb_SCOP()) {
    print(STDERR "ERROR: $MYPROGNAME: Database directory not found: '".$cfgvar->PathComerDb_SCOP()."'\n");
    exit(1);
}
unless(-d $cfgvar->PathComerDb_PFAM()) {
    print(STDERR "ERROR: $MYPROGNAME: Database directory not found: '".$cfgvar->PathComerDb_PFAM()."'\n");
    exit(1);
}

unless(-d $cfgvar->PathCotherDb_PDB()) {
    print(STDERR "ERROR: $MYPROGNAME: Database directory not found: '".$cfgvar->PathCotherDb_PDB()."'\n");
    exit(1);
}

## =============================================================================
## MAIN
##

print(STDERR "Processing pdb dbs...\n");
my  $pdbdbspat = '^pdb.+'.$dbsfx.'$';
exit(1) unless MakeLinks($cfgvar->PathComerDb_PDB(), $pdbdbspat);

print(STDERR "\nProcessing pfam dbs...\n");
my  $pfamdbspat = '^pfam.+'.$dbsfx.'$';
exit(1) unless MakeLinks($cfgvar->PathComerDb_PFAM(), $pfamdbspat);

print(STDERR "\nProcessing scop dbs...\n");
my  $scopdbspat = '^scop.+'.$dbsfx.'$';
exit(1) unless MakeLinks($cfgvar->PathComerDb_SCOP(), $scopdbspat);


print(STDERR "\n\nProcessing COTHER pdb dbs...\n");
my  $cpdbdbspat = '^pdb.+'.$cotherdbsfx.'$';
exit(1) unless MakeLinks($cfgvar->PathCotherDb_PDB(), $cpdbdbspat);

print(STDERR "\nFinished.\n");

exit(0);

## =============================================================================
##
##

sub MakeLinks
{
    my  $dirname = shift;
    my  $pattern = shift;
    my  $mysubname = (caller(0))[3];
    my  $pardir = dirname($dirname);
    my  $anyfilepat = '[^\.].+';
    my  $symlink_exists = eval { symlink("",""); 1 };
    my (@entries, @files);

    if(-l $dirname) {
        unless(unlink($dirname)) {
            print(STDERR "ERROR: $mysubname: Failed to unlink $dirname\n");
            return 0;
        }
    }
    unless(-d $dirname || mkdir($dirname)) {
        print(STDERR "ERROR: $mysubname: Failed to create directory: $dirname\n");
        return 0;
    }
    ##flush the directory
    if(ReadEntries($dirname, $anyfilepat, \@files, 1) < 0) {
        print(STDERR "ERROR: $mysubname: Failed to read directory: $dirname\n");
        return 0;
    }
    foreach my $f(@files) {
        my $fullf = File::Spec->catfile($dirname,$f);
        next if -d $fullf;
        unless(unlink($fullf)) {
            print(STDERR "ERROR: $mysubname: Failed to unlink $fullf\n");
            return 0;
        }
    }
    ##read subdirs of the parent directory
    if(ReadEntries($pardir, $pattern, \@entries, 0) < 0) {
        print(STDERR "ERROR: $mysubname: Failed to read directory: $pardir\n");
        return 0;
    }
    foreach my $sdir(@entries) {
        my @sdfiles;
        my $fullsdir = File::Spec->catfile($pardir,$sdir);
        ##read files in each subdirectory
        if(ReadEntries($fullsdir, $anyfilepat, \@sdfiles, 1) < 0) {
            print(STDERR "ERROR: $mysubname: Failed to read subdirectory: $fullsdir\n");
            return 0;
        }
        print(STDERR "   $sdir\n");
        foreach my $f(@sdfiles) {
            my $fullf = File::Spec->catfile($fullsdir,$f);
            my $source = File::Spec->catfile('..',$sdir,$f);##$fullf
            my $destin = File::Spec->catfile($dirname,$f);
            print(STDERR "      $f\n");
            unless($symlink_exists) {
                print(STDERR "WARNING: $mysubname: Symlink cannot be created for: $fullf\n");
                next;
            }
            unless(symlink($source,$destin)) {
                print(STDERR "ERROR: $mysubname: Failed to create symlink for: $fullf\n");
                return 0;
            }
        }
    }
    return 1;
}

## =============================================================================
## -----------------------------------------------------------------------------
## GetTime: get time string
##
sub GetTime
{
    return strftime("%H:%M:%S ",localtime());
}

sub GetDatetime
{
    return strftime("%a %b %d %H:%M:%S %Z %Y",localtime());
}

## =============================================================================
## -------------------------------------------------------------------
## read files or subdirectories in the given directory and return 1 
## if at least 1 entry has been read, 0 otherwise, -1 on error
##
sub ReadEntries
{
    my  $dirname = shift;
    my  $pattern = shift;  ## pattern of files to look for as locks
    my  $refents = shift; ## reference to vector of entries
    my  $readfiles = shift; ## flag of whether files have to be read
    my  @files;
    my  $locref = defined( $refents )? $refents: \@files;

    unless( opendir(DIR, $dirname)) {
        print(STDERR "ERROR: Cannot open directory $dirname.\n");
        return -1;
    }
    if($readfiles) {
        @{$locref} = grep { /$pattern/ && -f File::Spec->catfile($dirname,$_) } readdir(DIR);
    } else {
        @{$locref} = grep { /$pattern/ && -d File::Spec->catfile($dirname,$_) } readdir(DIR);
    }
    closedir(DIR);
    return 1 if 0 <= $#{$locref};
    return 0;
}

## -------------------------------------------------------------------
## ExecCommand: execute system command
##
sub CheckStatus
{
    ExecCommand();
}

sub ExecCommand
{
    my  $cmdline = shift;
    my  $returnstatus = shift;##if defined, return exit status of the command

    system($cmdline) if $cmdline;

    if($? == -1) {
        print(STDERR GetTime()."ERROR: Failed to execute command: $!\n");
        return($returnstatus? 1: 0);
    }
    if($? & 127) {
        printf(STDERR GetTime()."ERROR: Command terminated with signal %d, %s coredump.\n",
            ($? & 127), ($? & 128)? 'with': 'without');
        return($returnstatus? ($? & 127): 0);
    }
    else {
        if(($? >> 8) != 0) {
            printf( STDERR GetTime()."ERROR: Command failed and exited with status %d\n",
                    $? >> 8);
            return($returnstatus? ($? >> 8): 0);
        }
    }
    return($returnstatus? 0: 1);
}

## <<>>


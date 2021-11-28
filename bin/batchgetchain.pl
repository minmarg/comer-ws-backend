#!/usr/bin/perl -w

##
## (C)2021 Mindaugas Margelevicius
## Institute of Biotechnology, Vilnius University
##

use strict;
use Config;
use threads;
use threads::shared;
use FindBin;
use lib "$FindBin::Bin";
use File::Spec;
use File::Basename;
use Getopt::Long;

my  $MYPROGNAME = basename($0);
my  $LOCPDBDIR = glob("/data/databases/pdb");
my  $nCPUs = 1;
my  $maxncpus = 0;

my  $usage = <<EOIN;

Download pdb structures given by a list and extract their respective chains.
(C)2019-2020 Mindaugas Margelevicius, Vilnius University

Usage:
$MYPROGNAME <Parameters>

Parameters:

-i <list>          Input list of pdb ids with chain identifiers 
                   (e.g., 1ZEE_A), which can be given as a filename 
                   (with one item per line) or in-line with 
                   comma-separated entries.

-o <directory>     Output directory of resulting files.

--pdb <directory>  Local directory of pdb structure files.
           default=$LOCPDBDIR

-c <#cpus>         Number of CPUs to use for parallelizing the process.
           Default=$nCPUs

-h                 This text.

EOIN

my  $INLIST;
my  $OUTDIR;
my  $string;
my  $Fail = 0;

my  $result = GetOptions(
               'i=s'      => \$INLIST,
               'o=s'      => \$OUTDIR,
               'pdb=s'    => \$LOCPDBDIR,
               'c=i'      => \$nCPUs,
               'help|h'   => sub { print $usage; exit( 0 ); }
);

do { print $usage; $Fail = 1; }  unless $result;
do { print STDERR "ERROR: Input missing.\n$usage"; $Fail = 1; } unless($Fail || $INLIST);
do { print STDERR "ERROR: Directory missing.\n$usage"; $Fail = 1; } unless($Fail || $LOCPDBDIR);
do { print STDERR "ERROR: Directory of structure files not found: $LOCPDBDIR\n"; $Fail = 1; } unless($Fail || -d $LOCPDBDIR);
do { print STDERR "ERROR: Invalid number of CPUs.\n$usage"; $Fail = 1; } unless($Fail || $nCPUs >= 1);

if(GetNCPUs(\$maxncpus)) {
    if($maxncpus < $nCPUs) {
        print("ERROR: Specified #CPUs > #CPUs on system: $nCPUs > $maxncpus.\n");
        $Fail = 1;
    }
}

##switch to python3
#scl enable rh-python36 bash

if( $Fail ) {
    exit(1);
}

##programs:
my  $GETCHAIN = "$FindBin::Bin/getchain.py";
my  $UNZIP = "gunzip";

##remote addresses:
my  %FTPPDBDIR = (
    CIF => "ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/mmCIF",
    CIFOBS => "ftp://ftp.wwpdb.org/pub/pdb/data/structures/obsolete/mmCIF",
    PDB => "ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb",
    PDBOBS => "ftp://ftp.wwpdb.org/pub/pdb/data/structures/obsolete/pdb"
);

##check
do { print( STDERR "ERROR: Program not found: $GETCHAIN\n"); exit(1); } unless($GETCHAIN && -f $GETCHAIN);
do { print( STDERR "ERROR: Program not found: $UNZIP\n"); exit(1); } unless( RunCommand("which $UNZIP",0,\$string));


##go
my  $querydir = File::Spec->rel2abs( dirname($INLIST));
my  $querybasename = basename($INLIST);
my  $queryname = $querybasename; $queryname =~ s/\.[^\.]*$//;
my  $curdir = File::Spec->rel2abs(File::Spec->curdir());

do { print( STDERR "ERROR: Output dirname not given.\n"); exit(1); }unless($OUTDIR);

$OUTDIR = File::Spec->rel2abs($OUTDIR);

unless( -d $OUTDIR || mkdir($OUTDIR)) {
    print( STDERR "ERROR: Failed to create directory: $OUTDIR\n");
    exit(1);
}

my  @tids;

if(-f $INLIST) {
    ##the list given as a filename
    unless(open(I, $INLIST)) {
        print( STDERR "ERROR: Failed to open input file: $INLIST\n");
        exit(1);
    }
    while(<I>) {
        next if /^\s*#/;
        next if /^\s*$/;
        my @a = split(/\s+/);
        push @tids, $a[0];
    }
    close(I);
} else {
    ##the list provided directly
    @tids = split(',',$INLIST);
}



if($Config{useithreads}) {
    unless(LaunchThreads()) {
        print( STDERR "ERROR: Launching threads failed.\n");
        exit(1);
    }
} else {
    $nCPUs = 1;
    print(STDERR "WARNING: Perl compiled WITHOUT thread support: ".
        "Using only one cpu: #cpus=$nCPUs\n");
    unless(ProcessPartofTids(0)) {
        print( STDERR "ERROR: Extracting chains failed.\n");
        exit(1);
    }
}



##unless( chdir($curdir)) {
##    print( STDERR "ERROR: Failed to change directory to: $curdir\n");
##    exit(1);
##}

print("Finished.\n");
exit(0);

## ===================================================================
sub LaunchThreads
{
    my @workers;
    my @joinable = ();
    for(my $t = 0; $t <= $nCPUs; $t++) {
        $workers[$t] = threads->create(
                {'context' => 'list'},
                \&ProcessPartofTids, 
                $t
        );
        unless($workers[$t]) {
            print(STDERR "ERROR: Thread $t initialization failed.\n");
        }
    }
    while(threads->list(threads::all)) {
        @joinable = threads->list(threads::joinable);##...or use @workers
        foreach my $thr(@joinable) {
            my $thretcode = $thr->join();
            unless($thretcode) {
                printf(STDERR "ERROR: Thread %d returned the fail code.\n", $thr->tid());
                next;
            }
            printf(STDERR "\nI: Thread %d joined.\n", $thr->tid());
       }
    }
    return 1;
}

## ===================================================================
## Entry point for threads: process part of given PDB chains
##
sub ProcessPartofTids
{
    my $mythid = shift;
    for(my $i = $mythid; $i <= $#tids; $i += $nCPUs) {
        my $tmpl = $tids[$i];
        my $outfile;
        next if -f "$OUTDIR/$tmpl";
        unless( GetStructure($UNZIP, $GETCHAIN, $LOCPDBDIR, $OUTDIR, \%FTPPDBDIR, $tmpl, \$outfile)) {
            print( STDERR "ERROR: Failed to obtain structure for: $tmpl\n");
            next;##return(0);
        }
        ##if($OUTDIR cmp $LOCPDBDIR) {
        ##    unless( RunCommandV("mv $LOCPDBDIR/$outfile $OUTDIR/")) {
        ##       print( STDERR "ERROR: Failed to mv $outfile to $OUTDIR\n");
        ##       return(0);
        ##    }
        ##}
        print "Obtained: $outfile\n\n\n";
    }
    return(1);
}

## -------------------------------------------------------------------
sub GetNCPUs
{
    my $rncpus = shift;##ref to #cpus
    return 0 unless(open(F, "/proc/cpuinfo"));
    $$rncpus = scalar(map /^processor/, <F>); 
    close(F);
    return 1;
}

## -------------------------------------------------------------------
## Get the template name for output files
##
sub GetOutputTemplateName
{
    my $ltmplname = shift;##template name
    my $ltmplstrbasename = $ltmplname;##output template basename
    ##$ltmplstrbasename =~ s/\./_/g;
    return $ltmplstrbasename;
}

## -------------------------------------------------------------------
## download the structure of the given template from PDB if it is not 
## present yet
##
sub GetStructure
{
    my $unzipprog = shift;
    my $getchainprog = shift;
    my $lpdbdir = shift;##local directory of pdb structures
    my $outdir = shift;##output directory for chains
    my $rrmtdirs = shift;##ref to remote directories of structures (hash)
    my $ltmplname = shift;##template name
    my $rtmplstrfilename = shift;##ref to the filename of template structure (to return)
    my $lcurdir = File::Spec->rel2abs(File::Spec->curdir());
    my ($ltmplstruct, $ltmplchain) = ($ltmplname);
    my ($middle, $ciffilename, $pdbfilename, $fname);
    my $chre = qr/_([\da-zA-Z\.]+)$/;
    my $fail;

    ##unless( chdir($lpdbdir)) {
    ##    print( STDERR "ERROR: Failed to change directory to: $lpdbdir\n");
    ##    return(0);
    ##}

    $ltmplstruct =~ s/$chre//;
    $ltmplchain = $1 if $ltmplname =~ /$chre/;
    $$rtmplstrfilename = GetOutputTemplateName($ltmplname).".ent";

    if(-f File::Spec->catfile($outdir, $$rtmplstrfilename)) {
        ##unless( chdir($lcurdir)) {
        ##    print( STDERR "ERROR: Failed to change directory to: $lcurdir\n");
        ##    return(0);
        ##}
        return 1;
    }

    $middle = lc(substr($ltmplstruct,1,2));
    $ciffilename = lc(${ltmplstruct}).".cif";
    $pdbfilename = "pdb".lc(${ltmplstruct}).".ent";

    unless( -f File::Spec->catfile($lpdbdir,$ciffilename) || 
            -f File::Spec->catfile($lpdbdir,$pdbfilename)) {
        print("MSG: Downloading the structure for $ltmplname ...\n");
        ##first, try .cif file
        $fname = $ciffilename;
        unless( RunCommandV("wget -P $lpdbdir $$rrmtdirs{CIF}/$middle/${ciffilename}.gz")) {
            ##try obsolete .cif file
            unless( RunCommandV("wget -P $lpdbdir $$rrmtdirs{CIFOBS}/$middle/${ciffilename}.gz")) {
                ##then, .pdb file
                $fname = $pdbfilename;
                unless( RunCommandV("wget -P $lpdbdir $$rrmtdirs{PDB}/$middle/${pdbfilename}.gz")) {
                    ##obsolete .pdb file
                    unless( RunCommandV("wget -P $lpdbdir $$rrmtdirs{PDBOBS}/$middle/${pdbfilename}.gz")) {
                        print( STDERR "ERROR: Failed to download the structure for: $ltmplname\n");
                        $fail = 1;
                    }
                }
            }
        }
        print("\n") unless($fail);
    }

    if( !$fail && $fname ) {
        print("MSG: Unzipping...\n");
        $fail = 1 unless RunCommandV("$unzipprog ".File::Spec->catfile($lpdbdir,"${fname}.gz"));
        print("\n") unless($fail);
    }

    unless( $fail ) {
        $fname = (-f File::Spec->catfile($lpdbdir,$ciffilename))? 
                     File::Spec->catfile($lpdbdir,$ciffilename): 
                     File::Spec->catfile($lpdbdir,$pdbfilename);
        my $outfilename = File::Spec->catfile($outdir,$$rtmplstrfilename);
        print("MSG: Extracting chain of $fname ...\n");
        my $cmd = "python3 $getchainprog -i $fname -o $outfilename";
        $cmd .= " -c $ltmplchain" if $ltmplchain !~ /\./;
        $fail =1 unless RunCommandV($cmd);
        print("\n") unless($fail);
    }

    ##unless( chdir($lcurdir)) {
    ##    print( STDERR "ERROR: Failed to change directory to: $lcurdir\n");
    ##    return(0);
    ##}

    return !$fail;
}

## -------------------------------------------------------------------
## run system command
##
sub CheckStatus
{
    return RunCommand();
}

sub RunCommandV
{
    my  $cmdline = shift;
    my  $retstatus = shift;
    my  $routput = shift;##ref
    print( STDERR "CMD: $cmdline\n") if $cmdline;
    return RunCommand($cmdline, $retstatus, $routput);
}

sub RunCommand
{
    my  $cmdline = shift;
    my  $retstatus = shift;
    my  $routput = shift;##ref

    if($cmdline) {
        $$routput = `$cmdline 2>&1` if $routput;
        system( "$cmdline" ) unless $routput;
    }

    if( $? == -1 ) {
        printf( STDERR "ERROR: Failed to execute command: $!\n" );
        return 0;
    }
    if( $? & 127 ) {
        printf( STDERR "ERROR: Command terminated with signal %d (%s coredump).\n",
            ($? & 127), ($? & 128)? 'with' : 'without' );
        return 0;
    }
    else {
        if(( $? >> 8 ) != 0 ) {
            unless( $retstatus ) {
                printf( STDERR "ERROR: Command failed and exited with status %d\n", $? >> 8 );
                return 0
            }
            return( $? >> 8 );
        }
    }
    return 1;
}

##<<>>

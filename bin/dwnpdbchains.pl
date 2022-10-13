#!/usr/bin/perl -w

##
## (C)2019-2020 Mindaugas Margelevicius
## Institute of Biotechnology, Vilnius University
##

use strict;
use FindBin;
use lib "$FindBin::Bin";
use File::Spec;
use File::Basename;
use Getopt::Long;

my  $MYPROGNAME = basename($0);
my  $LOCPDBDIR = glob("/data/databases/pdb");

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

-h                 This text.

EOIN

my  $INLIST;
my  $OUTDIR;
my  $string;
my  $Fail = 0;

my  $result = GetOptions(
               'i=s'      => \$INLIST,
               'o=s'      => \$OUTDIR,
               'pdb=s'      => \$LOCPDBDIR,
               'help|h'   => sub { print $usage; exit( 0 ); }
);

do { print $usage; $Fail = 1; }  unless $result;
do { print STDERR "ERROR: Input missing.\n$usage"; $Fail = 1; } unless($Fail || $INLIST);
do { print STDERR "ERROR: Directory missing.\n$usage"; $Fail = 1; } unless($Fail || $LOCPDBDIR);
do { print STDERR "ERROR: Directory of structure files not found: $LOCPDBDIR\n"; $Fail = 1; } unless($Fail || -d $LOCPDBDIR);

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

my  $command;
my  @tids;

if(-f $INLIST) {
    ##the list given as a filename
    unless( open(I, $INLIST)) {
        print( STDERR "ERROR: Failed to open input file: $INLIST\n");
        exit(1);
    }
    while(<I>) {
        next if /^\s*#/;
        my @a = split(/\s+/);
        push @tids, $a[0];
    }
    close(I);
} else {
    ##the list provided directly
    @tids = split(',',$INLIST);
}

foreach my $tmpl(@tids) {
    my $outfile;
    next if -f "$OUTDIR/$tmpl";
    unless( GetStructure($UNZIP, $GETCHAIN, $LOCPDBDIR, \%FTPPDBDIR, $tmpl, \$outfile)) {
        print( STDERR "ERROR: Failed to obtain structure for: $tmpl\n");
        exit(1);
    }
    if(($OUTDIR cmp $LOCPDBDIR) && (! -f "$OUTDIR/$outfile")) {
        unless( RunCommandV("cp $LOCPDBDIR/$outfile $OUTDIR/")) {
           print( STDERR "ERROR: Failed to cp $outfile to $OUTDIR\n");
            exit(1);
        }
    }
    print "Obtained: $outfile\n\n\n";
}

unless( chdir($curdir)) {
    print( STDERR "ERROR: Failed to change directory to: $curdir\n");
    exit(1);
}

print("Finished.\n");
exit(0);

## ===================================================================
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
    my $rrmtdirs = shift;##ref to remote directories of structures (hash)
    my $ltmplname = shift;##template name
    my $rtmplstrfilename = shift;##ref to the filename of template structure (to return)
    my $lcurdir = File::Spec->rel2abs(File::Spec->curdir());
    my ($ltmplstruct, $ltmplchain) = ($ltmplname);
    my ($middle, $ciffilename, $pdbfilename, $fname);
    my $chre = qr/_([\da-zA-Z\.]+)$/;
    my $fail;

    unless( chdir($lpdbdir)) {
        print( STDERR "ERROR: Failed to change directory to: $lpdbdir\n");
        return(0);
    }

    $ltmplstruct =~ s/$chre//;
    $ltmplchain = $1 if $ltmplname =~ /$chre/;
    $$rtmplstrfilename = GetOutputTemplateName($ltmplname).".ent";

    if( -f $$rtmplstrfilename ) {
        unless( chdir($lcurdir)) {
            print( STDERR "ERROR: Failed to change directory to: $lcurdir\n");
            return(0);
        }
        return 1;
    }

    $middle = lc(substr($ltmplstruct,1,2));
    $ciffilename = lc(${ltmplstruct}).".cif";
    $pdbfilename = "pdb".lc(${ltmplstruct}).".ent";

    unless( -f $ciffilename || -f $pdbfilename) {
        print("MSG: Downloading the structure for $ltmplname ...\n");
        ##first, try .cif file
        $fname = $ciffilename;
        unless( RunCommandV("wget $$rrmtdirs{CIF}/$middle/${ciffilename}.gz")) {
            ##try obsolete .cif file
            unless( RunCommandV("wget $$rrmtdirs{CIFOBS}/$middle/${ciffilename}.gz")) {
                ##then, .pdb file
                $fname = $pdbfilename;
                unless( RunCommandV("wget $$rrmtdirs{PDB}/$middle/${pdbfilename}.gz")) {
                    ##obsolete .pdb file
                    unless( RunCommandV("wget $$rrmtdirs{PDBOBS}/$middle/${pdbfilename}.gz")) {
                        print( STDERR "ERROR: Failed to download the structure for: $ltmplname\n");
                        $fail = 1;
                    }
                }
            }
        }
        print("\n") unless($fail);
    }

    sleep(2+int(rand(5))) if(!$fail);

    if( !$fail && $fname && (-f "${fname}.gz") && (! -f ${fname})) {
        print("MSG: Unzipping...\n");
        $fail = 1 unless RunCommandV("$unzipprog ${fname}.gz");
        print("\n") unless($fail);
        $fail = 0 if(-f ${fname});
    }

    unless( $fail ) {
        $fname = (-f $ciffilename)? $ciffilename: $pdbfilename;
        print("MSG: Extracting chain of $fname ...\n");
        my $cmd = "python3 $getchainprog -i $fname -o $$rtmplstrfilename";
        $cmd .= " -c $ltmplchain" if $ltmplchain !~ /\./;
        $fail =1 unless RunCommandV($cmd);
        print("\n") unless($fail);
    }

    unless( chdir($lcurdir)) {
        print( STDERR "ERROR: Failed to change directory to: $lcurdir\n");
        return(0);
    }

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

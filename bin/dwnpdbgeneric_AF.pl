#!/usr/bin/perl -w

##
## (C)2019-2022 Mindaugas Margelevicius
## Institute of Biotechnology, Vilnius University
##

use strict;
use FindBin;
use lib "$FindBin::Bin";
use File::Spec;
use File::Basename;
use Getopt::Long;

my  $MYPROGNAME = basename($0);
my  $LOCPDBDIR = glob("/data/databases/pdbstyle");
my  $WEBADDR = 'https://alphafold.ebi.ac.uk/files';

my  $usage = <<EOIN;

Download pdb-style structures/models from given address.
(C)2019-2022 Mindaugas Margelevicius, Vilnius University

Usage:
$MYPROGNAME <Parameters>

Parameters:

--id <ID>          Filename with extension (e.g., AF-Q9PIM5-F1-model_v2.pdb).

-o <directory>     Output directory to locate the resulting file.

--pdb <directory>  Local directory of pdb-style structure files.
           default=$LOCPDBDIR

--add <web-address> Web address to the page of pdb-style structures.
           default=$WEBADDR

-h                 This text.

EOIN

my  $ID;
my  $OUTDIR;
my  $string;
my  $Fail = 0;

my  $result = GetOptions(
               'id=s'     => \$ID,
               'o=s'      => \$OUTDIR,
               'pdb=s'    => \$LOCPDBDIR,
               'add=s'    => \$WEBADDR,
               'help|h'   => sub { print $usage; exit( 0 ); }
);

do { print $usage; $Fail = 1; }  unless $result;
do { print STDERR "ERROR: Input ID not given.\n$usage"; $Fail = 1; } unless($Fail || $ID);
do { print STDERR "ERROR: Directory not given.\n$usage"; $Fail = 1; } unless($Fail || $LOCPDBDIR);
do { print STDERR "ERROR: Directory of structure files not found: $LOCPDBDIR\n"; $Fail = 1; } unless($Fail || -d $LOCPDBDIR);

exit(1) if $Fail;

##programs:
#my  $GETCHAIN = "$FindBin::Bin/getchain.py";
#my  $UNZIP = "gunzip";

##check
#do { print( STDERR "ERROR: Program not found: $GETCHAIN\n"); exit(1); } unless($GETCHAIN && -f $GETCHAIN);
#do { print( STDERR "ERROR: Program not found: $UNZIP\n"); exit(1); } unless( RunCommand("which $UNZIP",0,\$string));


##go
my  $querybasename = $ID;
my  $queryname = $querybasename;
my  $curdir = File::Spec->rel2abs(File::Spec->curdir());

do { print( STDERR "ERROR: Output dirname not given.\n"); exit(1); }unless($OUTDIR);

$OUTDIR = File::Spec->rel2abs($OUTDIR);

unless( -d $OUTDIR || mkdir($OUTDIR)) {
    print( STDERR "ERROR: Failed to create directory: $OUTDIR\n");
    exit(1);
}

my  $command;
my  @tids = ($ID);

foreach my $tmpl(@tids) {
    my $outfile = $tmpl;
    my $outfullfilename = File::Spec->catfile($OUTDIR, $outfile);
    next if -f $outfile;
    unless( GetStructure($LOCPDBDIR, $WEBADDR, $tmpl, $outfile, $outfullfilename)) {
        print( STDERR "ERROR: Failed to obtain structure for: $tmpl\n");
        exit(1);
    }
    printf("Obtained: %s\n\n",basename($outfile));
}

print("Finished.\n");
exit(0);

## ===================================================================
## -------------------------------------------------------------------
## download structure with the given filename from given site if 
## it is not present yet
##
sub GetStructure
{
    my $lpdbdir = shift;##local directory of pdb structures
    my $rmtaddr = shift;##remote (web) address to structures
    my $ltmplname = shift;##template name
    my $ltmplfilename = shift;##filename of output template structure
    my $ltmplfullfilename = shift;##full filename of output template structure
    my $pdbfilename;
    my @pdbinfo;
    my $fail;

    return 1 if(-f $ltmplfullfilename);

    $pdbfilename = File::Spec->catfile($lpdbdir, $ltmplfilename);

    unless(-f $pdbfilename) {
        print("MSG: Downloading the structure $ltmplname ...\n");
        $fail = 1;
        foreach(0..1) {
            ##try several times
            if(RunCommandV("wget -O ${ltmplfullfilename} ${rmtaddr}/${ltmplname}")) {
                $fail = 0;
                last;
            }
            sleep(1);
        }
        print(STDERR "ERROR: Failed to download the structure: $ltmplname\n") if($fail);
        print("\n") unless($fail);
    }

    return 0 if $fail;

#    unless(${pdbfilename} eq ${ltmplfullfilename}) {
#        unless(RunCommandV("mv ${pdbfilename} ${ltmplfullfilename}")) {
#            print(STDERR "ERROR: Failed to move and rename: ${pdbfilename}\n");
#            return 0;
#        }
#    }

#    unless(ProcessStructure($ltmplname, $pdbfilename, $ltmplfullfilename)) {
#        unless( RunCommandV("mv ${pdbfilename} ${ltmplfullfilename}.corrupted")) {
#            print(STDERR "ERROR: Failed to move and rename: ${pdbfilename}\n");
#        }
#        return 0;
#    }

    return !$fail;
}

## -------------------------------------------------------------------
## ProcessStructure: verify the structure and ensure chain consistency
## chain; write the final structure to file
##
sub ProcessStructure
{
    my $ltmplname = shift;##template name
    my $inputfile = shift;##input file of downloaded structure
    my $outputfile = shift;##name of the output file to be written to
    my @pdbinfo;
    my $fail;

    unless(open(F,$inputfile)) {
        print(STDERR "ERROR: Failed to open file: $inputfile\n");
        return 0;
    }
    while(<F>) {
        if(/^(?:HEADER|REMARK|MODEL|END)/) {
            push @pdbinfo, $_;
            next;
        }
        if(/^(?:ATOM|HETATM|TER)/) {
            if(length($_) < 22) {
                print(STDERR "ERROR: Invalid file format: line: $. file: $inputfile\n");
                $fail = 1;
                last;
            }
            substr($_,21,1) = ' ';#change chain ID to space
            push @pdbinfo, $_;
            next;
        }
        if(/^ENDMDL/) {
            push @pdbinfo, $_;
            push @pdbinfo,"END\n";
            last;
        }
    }
    close(F);
    return 0 if $fail;
    if($#pdbinfo<0) {
        print(STDERR "ERROR: No structure in file: $inputfile\n");
        return 0;
    }
    unless(open(F,'>',$outputfile)) {
        print(STDERR "ERROR: Failed to open file for writing: $outputfile\n");
        return 0;
    }
    print(F) foreach(@pdbinfo);
    close(F);
    return 1;
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

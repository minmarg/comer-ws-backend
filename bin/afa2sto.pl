#!/usr/bin/env perl
BEGIN {$^W=1}

## (C) 2020 Mindaugas Margelevicius, Institute of Biotechnology, Vilnius University

use strict;
use FindBin;
use lib "$FindBin::Bin";
use File::Basename;
use Getopt::Long;


my  $MYPROGNAME = basename($0);
my  $WIDTH = -1;

my  $usage = <<EOIN;

Convert the format of an MSA from FASTA to STOCKHOLM.
(C) 2020 Mindaugas Margelevicius, Institute of Biotechnology, Vilnius University

Usage:
$MYPROGNAME -i <input> -o <output> [OPTIONS]

Options:

-i <input>         Input MSA file in FASTA.
-o <output>        Output MSA file in STOCKHOLM.
-h                 This text

EOIN

my  $INFILE;
my  $OUTFILE;

my  $result = GetOptions(
               'i=s'      => \$INFILE,
               'o=s'      => \$OUTFILE,
               'help|h'   => sub { print $usage; exit( 0 ); }
);


do { print $usage; exit(1); }  unless $result;

die "ERROR: No input filename given." unless $INFILE;
die "ERROR: No output filename given." unless $OUTFILE;
die "ERROR: Input alingment file not found: $INFILE" unless($INFILE && -f $INFILE);

$INFILE =~ s/^\s*//; $INFILE =~ s/\s*$//; $INFILE =~ s/"//;
$OUTFILE =~ s/^\s*//; $OUTFILE =~ s/\s*$//; $OUTFILE =~ s/"//;

my  $myinput = $INFILE;
my  $myoutput = $OUTFILE;
my  $inbase = basename($INFILE);
my  $outdir = dirname($OUTFILE);
my  %alns;

unless(ReadFasta($myinput, \%alns)) {
    exit(1);
}

unless(WriteStockholm1($myoutput, \%alns, $WIDTH)) {
    exit(1);
}

exit(0);

## -------------------------------------------------------------------
## ReadFasta: Read an MSA in FASTA
##
sub ReadFasta
{
    my  $input = shift;
    my  $ralns = shift;#ref to alignment hash
    my ($desc,$seqn);
    my ($n,$len) = (0,-1);
    my  $ret = 1;

    unless( open( FF, $input )) {
        printf(STDERR "ERROR: Failed to open input file: %s\n", $input);
        return 0;
    }
    while(<FF>) {
        chomp;
        if( /^\s*>(.*)\s*$/ ) { 
            $$ralns{$1}{N} = $n++;
            $$ralns{$1}{S} = '';
            $seqn = \$$ralns{$1}{S};
            next;
        }
        s/\s//g;
        unless($seqn) {
            printf(STDERR "ERROR: Input file is not in FASTA: %s\n", $input);
            $ret = 0;
        }
        $$seqn .= $_;
    }
    close(FF);
    return $ret unless $ret;

    ##verify length of sequences
    foreach $desc(keys %$ralns) {
        $seqn = \$$ralns{$desc}{S};
        $len = length($$seqn) if $len < 0;
        if(length($$seqn) != $len || $len <= 0) {
            printf(STDERR "ERROR: Invalid length of the sequence: %s\n", $desc);
            $ret = 0;
            return $ret;
        }
    }

    return $ret;
}

## -------------------------------------------------------------------
## Write the alignment to file in STOCKHOLM
##
sub WriteStockholm1
{
    my  $output = shift;
    my  $ralns = shift;#ref to alignment hash
    my  $width = shift;#output width for aligned sequences
    my  @srtkeys;
    my ($hdn,$unkdesc,$desc,$seqn,$key,$chr) = ('','','','','','');
    my ($len,$counter) = (0,0);
    my (@matchndxs);##, @insndxs);##match,insert state indices
    my  $ret = 1;

    ##write sequences to file
    unless( open( OF, ">$output" )) {
        printf( STDERR "ERROR: Failed to open file for writing: %s\n", $output );
        $ret = 0;
        return $ret;
    }

    @srtkeys = sort {$$ralns{$a}{N} <=> $$ralns{$b}{N}} keys %$ralns;

    if(0 <= $#srtkeys) {
        $key = $srtkeys[0];
        if($key =~ /^(\S+)\s*(.*)$/) {
            $hdn = $1;
            $desc = $2;
        }
        $hdn = "Unknown" unless $hdn;
        printf(OF "# STOCKHOLM 1.0\n");
        printf(OF "#=GF ID %s\n", $hdn);
        printf(OF "##=GF DE %s\n", $desc);
        printf(OF "##=GF AU The author (%s)\n", $MYPROGNAME);
        printf(OF "\n");
        $seqn = \$$ralns{$key}{S};
        $$seqn = uc($$seqn);
        $$seqn =~ s/\-/\./g;
        push @matchndxs, $-[0] while $$seqn =~ /[a-zA-Z]/g;
        ##push @insndxs while $$seqn =~ /[\.]/g;
    }

    foreach $key( @srtkeys ) {
        $hdn = $desc = '';
        if($key =~ /^(\S+)\s*(.*)$/) {
            $hdn = "$1_$$ralns{$key}{N}";
            $desc = $2;
        }
        $hdn = sprintf("Unknown_seq_%d", $counter++) unless $hdn;
        printf(OF "#=GS %-26s DE %s\n", $hdn, $desc);
        $$ralns{$key}{H} = $hdn;
    }

    print(OF "\n");

    foreach $key( @srtkeys ) {
        $hdn = $$ralns{$key}{H};
        $seqn = \$$ralns{$key}{S};
        $$seqn = lc($$seqn);
        $$seqn =~ s/\-/\./g;
        foreach(@matchndxs) {
            $chr = substr($$seqn, $_, 1);
            if($chr eq '.') {substr($$seqn, $_, 1) = '-'}
            else {substr($$seqn, $_, 1) = uc($chr)}
        }
        printf(OF "%-34s ", $hdn);
        unless( WrapSequence(\*OF, $seqn, $width)) {
            $ret = 0;
            last;
        }
    }

    ##last state-annotated record
    printf(OF "%-34s %s\n", "#=GC RF", $$ralns{$srtkeys[0]}{S});
    print(OF "//\n");

    close(OF);
    return $ret;
}

## -------------------------------------------------------------------
## wrap sequence in fragments of equal length
##
sub WrapSequence {
    my  $reffile = shift;   ## reference to file descriptor
    my  $reffasta = shift;  ## reference to sequence
    my  $width = shift;     ## width of fragment per line
    my  $padding = 0;       ## padding at the beginning of each line
    my  $line;

    $width = 999999 if $width <= 0;

    if( ref( $reffile ) ne 'GLOB' && ref( $reffile ) ne 'SCALAR' ) {
        printf( STDERR "ERROR: WrapFasta: Wrong reference.\n" );
        return 0;
    }

##    $$reffile = '' if( ref( $reffile ) eq 'SCALAR' );

    for( my $n = 0; $n < length( $$reffasta ); $n += $width ) {
        if( $n && $padding ) {
            $line = sprintf( "%${padding}s%s\n", ' ', substr( $$reffasta, $n, $width ));
        } else {
            $line = sprintf( "%s\n", substr( $$reffasta, $n, $width ));
        }
        if( ref( $reffile ) eq 'SCALAR' ) {
                 $$reffile .= $line;
        } else { printf( $reffile $line );
        }
    }
    return 1;
}


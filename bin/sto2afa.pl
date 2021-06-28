#!/usr/bin/env perl
BEGIN {$^W=1}

## (C) 2020 Mindaugas Margelevicius, Institute of Biotechnology, Vilnius University

use strict;
use FindBin;
use lib "$FindBin::Bin";
use File::Basename;
use Getopt::Long;


my  $MYPROGNAME = basename($0);
my  $CONSENSUS = 0;
my  $WIDTH = -1;

my  $usage = <<EOIN;

Convert the format of an MSA from STOCKHOLM to FASTA.
(C) 2020 Mindaugas Margelevicius, Institute of Biotechnology, Vilnius University

Usage:
$MYPROGNAME -i <input> -o <output> [OPTIONS]

Options:

-i <input>         Input MSA file in STOCKHOLM.
-o <output>        Output MSA file in FASTA.
-c                 Calculate a consensus sequence to be the first
                   sequence in the resulting MSA.
-w <width>         Number of characters to wrap aligned sequences at.
                   No wrapping is used by default.
-h                 This text

EOIN

my  $INFILE;
my  $OUTFILE;

my  $result = GetOptions(
               'i=s'      => \$INFILE,
               'o=s'      => \$OUTFILE,
               'c'        => sub { $CONSENSUS = 1; },
               'w=i'      => \$WIDTH,
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

unless(ReadStockholm1($myinput, \%alns)) {
    exit(1);
}

if($CONSENSUS) {
    exit(1) unless MakeConsensusSeq(\%alns);
}

unless(WriteFasta($myoutput, \%alns, $WIDTH)) {
    exit(1);
}

exit(0);

## -------------------------------------------------------------------
## ReadStockholm1: Read an MSA in STOCKHOLM
##
sub ReadStockholm1
{
    my  $input = shift;
    my  $ralns = shift;#ref to alignment hash
    my  @srtkeys;
    my ($hdn,$desc,$seqn,$fstseqn) = ('','','','');
    my ($supseqn,$suptitle) = ('','');
    my ($n,$len) = (0,0);
    my  $ret = 1;

    unless( open( FF, $input )) {
        printf( STDERR "ERROR: Failed to open input file: %s\n", $input );
        return 0;
    }
    $_ = <FF>;
    unless(/^#\s+STOCKHOLM/) {
        close(FF);
        printf(STDERR "ERROR: Input not in STOCKHOLM format.\n");
        $ret = 0;
        return $ret;
    }
    while(<FF>) {
        chomp;
        $suptitle .= $1 if(/^#=GF\s+(?:AC|DE)\s+(.+)$/);
        if(/^#=GS\s+(\S+)\s+\S+\s+(.+)$/) {
            $hdn = $1;
            $desc = $2;
            $$ralns{$hdn}{D} .= "${hdn} ${desc}";
        }
        next if /^#/;
        last if /^\/\//;
        if(/^(\S+)\s+(\S+)$/) {
            $hdn = $1;
            $seqn = $2;
            $n++;
            $$ralns{$hdn}{N} = $n unless exists $$ralns{$hdn}{N};
            $$ralns{$hdn}{S} .= $seqn;
        }
    }
    close(FF);
    return $ret unless $ret;

    @srtkeys = sort {$$ralns{$a}{N} <=> $$ralns{$b}{N}} keys %$ralns;
    if(0 <= $#srtkeys) {
        $fstseqn = $$ralns{$srtkeys[0]}{S};
        ##change the description of the first sequence:
        $$ralns{$srtkeys[0]}{D} = $suptitle if($suptitle);
    }
    $len = length($fstseqn);

    ##verify length of sequences
    foreach $hdn( @srtkeys ) {
        $seqn = \$$ralns{$hdn}{S};
        unless( $hdn && $$seqn ) {
            printf( STDERR "ERROR: ReadStockholm1: Invalid sequence %s\n", $hdn );
            $ret = 0;
            last;
        }
        if( length($$seqn) != $len ) {
            printf(STDERR "ERROR: ReadStockholm1: Invalid length of the sequence: %s\n", $hdn);
            $ret = 0;
            return $ret;
        }
    }

    return $ret;
}

## -------------------------------------------------------------------
## MakeConsensusSeq: Make the consensus sequence from the MSA 
## encoded in STOCKHOLM
##
sub MakeConsensusSeq
{
    my  $ralns = shift;#ref to alignment hash
    my (@srtkeys,$fstseqn);
    my (%RES,@sorted,$slcted);
    my ($hdn,$seqn,$supseqn) = ('','','');
    my ($ch,$rr) = ('','');
    my ($n,$s,$len,$restot) = (0,0,0,0);
    my  $ret = 1;

    ##{{constants:
    ##minimum fraction of delete state symbols per position to consider 
    ##  position to be in delete state
    my  $frcDEL = 0.7;
    ##calculate delete states using $frcDEL as a criterion for the minimum 
    ##  fraction of gaps per column; 
    ##  see InputMultipleAlignment::TranslateSTOStates for discussion
    my  $cbCALCDELSTATES = 0;
    ##}}

    @srtkeys = sort {$$ralns{$a}{N} <=> $$ralns{$b}{N}} keys %$ralns;

    if($#srtkeys < 0) {
        printf(STDERR "ERROR: MakeConsensusSeq: No aligned sequences.\n");
        $ret = 0;
        return $ret;
    }

    $fstseqn = \$$ralns{$srtkeys[0]}{S};
    $len = length($$fstseqn);

    ##{{find out state at each position and make a consensus
    for( $n = 0; $n < $len; $n++ ) {
        $ch = substr($$fstseqn, $n, 1);
        unless( $ch =~ /[A-Z\-]/ ) {
            $supseqn .= '-';
            next;
        }
        ##find out residue distribution at the position
        undef %RES;
        $restot = 0;
        foreach $hdn( @srtkeys ) {
            $seqn = \$$ralns{$hdn}{S};
            $rr = substr($$seqn, $n, 1);
            $RES{$rr}++;
            $restot++;
        }
        ##sort by count
        @sorted = sort { $RES{$b} <=> $RES{$a} } keys %RES;
        ##if most observed is GAP
        $rr = $sorted[0];
        unless( '-' cmp $rr || !$restot || !$cbCALCDELSTATES ) {
            if( $frcDEL < $RES{$rr}/$restot ) {
                $supseqn .= '-';
                next;
            }
        }
        ##select a residue observed max number of times
        $slcted = 0;
        for( $s = 0; $s <= $#sorted; $s++ ) {
            $rr = $sorted[$s];
            unless( '-' cmp $rr ) {
                next;
            }
            if( $RES{$rr} < 1 ) {
                last;
            }
            if( $rr =~ /[XBZJO]/ && $s+1 <= $#sorted &&
                $sorted[$s+1] !~ /[XBZJO]/ && $RES{$rr} == $RES{$sorted[$s+1]}) {
                ##count of actual amino acid is the same
                next;
            }
            $supseqn .= $rr;
            $slcted = 1;
            last;
        }
        ##if no residue selected
        unless( $slcted ) {
            ##delete state columns containing only GAP may occur
            $supseqn .= '-';
            next;
            printf( STDERR "ERROR: MakeConsensusSeq: Unmatched state at position %d.", $n );
            $ret = 0;
            return $ret;
        }
    }
    ##}}

    if( length($supseqn) != $len ) {
        printf( STDERR "ERROR: MakeConsensusSeq: Invalid length of the resulting consensus sequence.\n");
        $ret = 0;
        return $ret;
    }

    ##save the information for the consensus sequence
    $hdn = "THISISCONSENSUSSEQUENCE_by_${MYPROGNAME}";
    $$ralns{$hdn}{N} = 0;
    ##NOTE: when "_consensus" appended to $srtkeys[0], hhblist crashes!
    $$ralns{$hdn}{D} = "$srtkeys[0]_cons. ";
    $$ralns{$hdn}{D} .= $$ralns{$srtkeys[0]}{D} if $$ralns{$srtkeys[0]}{D};
    $$ralns{$hdn}{D} =~ s/consensus/cons./g;
    $$ralns{$hdn}{D} .= " [Reference consensus sequence]";
    $$ralns{$hdn}{S} = $supseqn;

    return $ret;
}

## -------------------------------------------------------------------
## Write the alignment to file in FASTA
##
sub WriteFasta
{
    my  $output = shift;
    my  $ralns = shift;#ref to alignment hash
    my  $width = shift;#output width for aligned sequences
    my  @srtkeys;
    my ($hdn,$unkdesc,$desc,$seqn) = ('','','','');
    my ($len,$counter) = (0,0);
    my  $ret = 1;

    ##write sequences to file
    unless( open( OF, ">$output" )) {
        printf( STDERR "ERROR: Failed to open file for writing: %s\n", $output );
        $ret = 0;
        return $ret;
    }

    @srtkeys = sort {$$ralns{$a}{N} <=> $$ralns{$b}{N}} keys %$ralns;

    foreach $hdn( @srtkeys ) {
        $unkdesc = sprintf("Unknown sequence %d", $counter++);
        $desc = \$unkdesc;
        if(exists $$ralns{$hdn}{D}) { $desc = \$$ralns{$hdn}{D} }
        elsif($hdn) { $desc = \$hdn }
        $seqn = \$$ralns{$hdn}{S};
        unless( $hdn && $$seqn ) {
            printf( STDERR "ERROR: WriteFasta: Invalid sequence %s\n", $hdn );
            $ret = 0;
            last;
        }
        $$seqn =~ s/[\._~]/-/g;
        $$seqn = uc($$seqn);
        printf(OF ">%s\n", $$desc);
        unless( WrapSequence(\*OF, $seqn, $width)) {
            $ret = 0;
            last;
        }
    }

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


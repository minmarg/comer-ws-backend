#!/usr/bin/env perl
BEGIN {$^W=1}

## (C) 2020 Mindaugas Margelevicius, Institute of Biotechnology, Vilnius University

use strict;
use File::Basename;
use Getopt::Long;

my  $MYPROGNAME = basename( $0 );
my  $X = 'X';
my  $EVALUE = 1e-4;
my  $ROUND = -1;
my  $WIDTH = -1;
my  $ALTALN = 1;
my  $QUERYBEGSHIFT = 0;
my  $usage = <<EOIN;

Convert HHblits/HHsearch output of pairwise alignments to pairwise sequence 
alignments in FASTA-like format.
(C)2020 Mindaugas Margelevicius, Institute of Biotechnology, Vilnius University


Usage:
$MYPROGNAME <Parameters>

Parameters:

-i <filename>   HHblits/HHsearch output file of pairwise alignments.

-o <filename>   Name of output file of multiple alignment in FASTA.

-e <e-value>    Read alignments with expectation value less or equal
                than the value given.
        default=$EVALUE

-b <pos_shift>  Number of positions to add to the query alignment boundaries.
        default=$QUERYBEGSHIFT

-r <round>      Number of iteration to extract alignments from.
        By default, alignments from the last round will be processed.

-a              Do not include alternative alignments for same
                target sequences.

-w <width>      Column width to wrap sequence data within.
        By default, no wrapping is used.

-h              Short description.

EOIN


my  $INPUT;
my  $OUTPUT;
my  $Fail = 0;
my  $processed = 0;

my  $result = GetOptions(
               'i=s'      => \$INPUT,
               'o=s'      => \$OUTPUT,
               'e=f'      => \$EVALUE,
               'b=i'      => \$QUERYBEGSHIFT,
               'r=i'      => \$ROUND,
               'w=i'      => \$WIDTH,
               'a'        => sub { $ALTALN = 0; },
               'help|h'   => sub { print $usage; exit( 0 ); }
);


do { print $usage; $Fail = 1; }  unless $result;
do { print STDERR "ERROR: Missing arguments.\n$usage"; $Fail = 1; } unless( $Fail || ( $INPUT && $OUTPUT ));
do { print STDERR "ERROR: File $INPUT does not exist.\n"; $Fail = 1; } unless( $Fail || -f $INPUT );


## ===================================================================

if( $Fail ) {
    exit( 1 );
}

unless( FormatTopHits( $INPUT, $OUTPUT, $EVALUE, $QUERYBEGSHIFT, $ROUND, $WIDTH, \$processed )) {
    printf( STDERR "Failed.\n" );
    exit( 1 );
}

printf( STDERR "\nTotal sequences, %d\nDone.\n", $processed );
exit( 0 );

## ===================================================================
## align top hits from blast output
##

sub FormatTopHits
{
    my  $input = shift;
    my  $output = shift;
    my  $evalue = shift;
    my  $querybegshift = shift;
    my  $round = shift;
    my  $width = shift;
    my  $rproc = shift;
    my  %hithash;
    my  %descrip;
    my  $alnmask;
    my  $query = 'theonly';
    my ($rec, $len ) = ( 1, 0 ); ## the first record is reserved

    if( $round <= 0 ) {
        return 0 unless DetermineNoRounds( $input, \$round );
    }

    return 0 unless ExtractAlignments( $input, $evalue, $querybegshift, $round, $query, \%hithash, $rec, \%descrip );
    return 0 unless PrintAligned( \$descrip{$query}, \%hithash, $output, $width, $rproc );
    return 1;
}

## -------------------------------------------------------------------
## determine number of rounds in HHblits/HHsearch output
##

sub DetermineNoRounds
{
    my  $input = shift;
    my  $refround = shift;

    unless( open( IN, $input )) {
        printf( STDERR "ERROR: Failed to open $input: $!\n" );
        return 0;
    }
    while( <IN> ) {
        $$refround = $1 if /^Results\s+from\s+round\s+(\d+)/;
    }
    close( IN );
    return 1;
}

## -------------------------------------------------------------------
## extract top alignments from HHblits/HHsearch output
##

sub ExtractAlignments
{
    my  $input = shift;
    my  $evalue = shift;
    my  $querybegshift = shift;
    my  $round = shift;
    my  $query = shift;
    my  $refhash = shift; ## reference to hash of hits
    my  $stindex = shift; ## start index
    my  $refdesc = shift; ## reference to hash of query descriptions

    my  $hitnum = 0;
    my  $rec = $stindex;
    my  $last;
    my  $e_val;
    my  $score;

    my ($sbjct, $titletmp, $skip );

    my  $queryfasta;
    my  $sbjctfasta;

    my  $querystart;
    my  $sbjctstart;
    my  $queryend;
    my  $sbjctend;
    my ($sbjctlen)=(-1);

    my ($curound, $desc ) = ( 0, 0 );
    $$refdesc{$query} = '';

    unless( open( IN, $input )) {
        printf( STDERR "ERROR: Failed to open $input: $!\n" );
        return 0;
    }
    while( <IN> ) {
##        next if /^$/; ## comment to make eof( IN ) to work
        chomp;
        $last = $_;

        if( $last =~ /^Query\s+(.*)$/ ) {
            $$refdesc{$query} = $1;
            $desc = 1;
            next;
        }
        if( $last =~ /^Match_columns\s+\d+/ && $desc ) {
            $desc = 0;
            next;
        }
        if( $desc ) {
            $last =~ s/^\s*(.+)\s*$/$1/;
            $$refdesc{$query} .= " $last";
            next;
        }


        if( $last =~ /^Results\s+from\s+round\s+(\d+)/ ) {
            $curound = $1;
        }

        next if 0 < $round && $curound < $round;

        if( $last =~ /^>(.*)\s*$/ || 
           ($last =~ /^Probab=/ && $queryfasta ) ||
            eof( IN ))
        {
            $titletmp = $1;
            $titletmp = '' if $last =~ /^Probab=/;
            $hitnum++ if $last =~ /^>/;
            $skip = 0;

            if( defined( $sbjct ) && defined( $e_val )) {
                if( !$sbjct ) {
                    printf( STDERR "WARNING: No subject description (Hit no. %d). Skipped.\n", $hitnum - 1 ) if $ALTALN;
                    $skip = 1;
                }
                elsif( $e_val =~ /^\-1/ ) {
                    printf( STDERR "WARNING: No expect value (Hit no. %d: %s...). Skipped.\n",
                            $hitnum - 1, substr( $sbjct, 0, 20 ));
                    $skip = 1;
                }
                elsif( length( $queryfasta ) != length( $sbjctfasta ) || length( $sbjctfasta ) == 0 ) {
                    printf( STDERR "WARNING: Invalid alignment (sequence lengths, Hit no. %d: %s...). Skipped.\n",
                            $hitnum - 1, substr( $sbjct, 0, 20 ));
                    $skip = 1;
                }
                elsif( $querystart == -1 || $sbjctstart == -1 ) {
                    printf( STDERR "WARNING: No alignment start positions (Hit no. %d: %s...). Skipped.\n",
                            $hitnum - 1, substr( $sbjct, 0, 20 ));
                    $skip = 1;
                }
                elsif( $queryend == -1 || $sbjctend == -1 ) {
                    printf( STDERR "WARNING: No alignment end positions (Hit no. %d: %s...). Skipped.\n",
                            $hitnum - 1, substr( $sbjct, 0, 20 ));
                    $skip = 1;
                }

                $skip = 1 if !defined($e_val) || $evalue < $e_val;

                unless( $skip ) {
                    $$refhash{$query}[$rec][0] = $sbjct;
                    $$refhash{$query}[$rec][1] = $e_val;

                    $$refhash{$query}[$rec][2] = ''; ## reserved
                    $$refhash{$query}[$rec][3] = ''; ## reserved
                    $$refhash{$query}[$rec][4] = $queryfasta;
                    $$refhash{$query}[$rec][5] = $sbjctfasta;

                    $$refhash{$query}[$rec][6] = $querystart + $querybegshift;
                    $$refhash{$query}[$rec][7] = $queryend + $querybegshift;
                    $$refhash{$query}[$rec][8] = $sbjctstart;
                    $$refhash{$query}[$rec][9] = $sbjctend;
                    $$refhash{$query}[$rec][10]= $sbjctlen;

                    if($$refhash{$query}[$rec][6] < 1 || $$refhash{$query}[$rec][7] < 1) {
                        printf( STDERR "ERROR: Invalid positional shift parameter: Obtained for Hit no. %d: %s...\n",
                            $hitnum-1, substr($sbjct, 0, 20));
                        close( IN );
                        return 0;
                    }

                    $rec++;
                }
                last if 0 < $round && $round < $curound;
            }

            $sbjct = $titletmp if $last =~ /^>/ || !$ALTALN;
            $sbjctlen = -1 if $last =~ /^>/;
            $queryfasta = '';
            $sbjctfasta = '';
            $querystart = -1;
            $sbjctstart = -1;
            $queryend = -1;
            $sbjctend = -1;
            $e_val = -1;
            next if $last =~ /^>/;
        }

        if( $last =~ /^>/ ) {
            printf( STDERR "ERROR: Missed hit (Hit no. %d: %s...). Probably wrong format. Terminating.\n",
                    $hitnum, substr( $sbjct, 0, 20 ));
            close( IN );
            return 0;
        }

        if( $last =~ /^Probab=\s*([\d\.eE\-\+]+)\s+E-value=\s*([\d\.eE\-\+]+)\s+
                       Score=\s*([\d\.eE\-\+]+)/x )
        {
            $score = $3;
            $e_val = $2;
            $e_val = "1$e_val" if $e_val =~ /^e/i;
            next;
        }

        if( $last =~ /^Q\s+(?!Consensus)\S+\s+(\d+)\s+([a-zA-Z\-]+)\s+(\d+)/ ) {
            $querystart = $1 if $querystart == -1;
            $queryend = $3;
            $queryfasta .= $2; #uc("$2");
            next;
        }

        if( $last =~ /^T\s+(?!Consensus)\S+\s+(\d+)\s+([a-zA-Z\-]+)\s+(\d+)\s+\((\d+)\)/ ) {
            $sbjctstart = $1 if $sbjctstart == -1;
            $sbjctend = $3;
            $sbjctlen = $4 if !defined($sbjctlen) || $sbjctlen < 0;
            $sbjctfasta .= $2; #uc("$2");
            next;
        }
    }

    close( IN );
    return 1;
}

## -------------------------------------------------------------------
## print aligned top hits
##

sub PrintAligned
{
    my  $rquerydescription = shift; ##ref
    my  $rhithash = shift; ## reference
    my  $filename = shift;
    my  $width = shift;
    my  $rproc = shift;

    my  $rec;
    my  $e_val;
    my  $query;
    my  $sbjct;
    my  $qufasta;
    my  $sbfasta;

    my  $querydesc;
    my  $querystart;
    my  $sbjctstart;
    my  $queryend;
    my  $sbjctend;

    unless( open( OUT, ">$filename" )) {
        printf( STDERR "ERROR: Failed to open $filename for writing.\n" );
        return 0;
    }

    foreach $query( keys %{$rhithash} ) {
        for( $rec = 1; $rec <= $#{$$rhithash{$query}}; $rec++ ) {
            $sbjct      = $$rhithash{$query}[$rec][0];
            $e_val      = $$rhithash{$query}[$rec][1];
            $qufasta   = \$$rhithash{$query}[$rec][4];
            $sbfasta   = \$$rhithash{$query}[$rec][5];

            $querystart = $$rhithash{$query}[$rec][6];
            $queryend   = $$rhithash{$query}[$rec][7];
            $sbjctstart = $$rhithash{$query}[$rec][8];
            $sbjctend   = $$rhithash{$query}[$rec][9];

            $querydesc = $$rquerydescription;
            $querydesc =~ s/^(\S+)\s*(.*)$/$1 (ALN: $querystart-$queryend) $2  Expect=$e_val/;
            printf( OUT ">%s\n", $querydesc );
            WrapFasta( \*OUT, $qufasta, $width );
            $sbjct =~ s/^(\S+)\s*(.*)$/$1 (ALN: $sbjctstart-$sbjctend) $2  Expect=$e_val/;
            printf( OUT ">%s\n", $sbjct );
            WrapFasta( \*OUT, $sbfasta, $width );
            printf( OUT "//\n");
            $$rproc++;
        }
    }

    close( OUT );
    return 1;
}



## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## -------------------------------------------------------------------
## wrap sequence to fragments of equal length
##

sub WrapFasta {
    my  $reffile = shift;   ## reference to file descriptor
    my  $reffasta = shift;  ## reference to sequence
    my  $width = shift;     ## width of fragment per line
    my  $padding = 0;       ## padding at the beginning of each line
    my  $line;

    $width = 99999 if $width <= 0;

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

## <<>>


package readcfg;

## (C) 2020 Mindaugas Margelevicius, Institute of Biotechnology, Vilnius University
## Perl package for reading backend configuration

use strict;
use FindBin;
use lib "$FindBin::Bin";

## -------------------------------------------------------------------

BEGIN {
    my  $MAILADDRESSEE = '';
    my  $MAILSENDER = '';
    my  $MAILSERVER = '';

    sub GetMailAddressee { return $MAILADDRESSEE; }
    sub GetMailSender  { return $MAILSENDER; }
    sub GetMailServer  { return $MAILSERVER; }
}

## ===================================================================

sub new {
    my $that = shift;
    my $confile = shift;
    my $class = ref( $that ) || $that;
    my $self;

    $self->{LOCPDBDIR} = '';
    $self->{LOCSCOPeDIR} = '';
    $self->{WEBSCOPePDB} = '';
    $self->{MODELLER} = '';

    $self->{JOB_NUM_CPUS} = 1;
    $self->{JOB_MAX_NQUERIES} = 1;
    $self->{JOB_MAX_NQUERIES_COTHER} = 1;
    $self->{JOB_MAX_NMODELS} = 1;
    $self->{JOB_MAX_N3DTEMPLATES} = 1;

    $self->{PATHSEQDB_UNIREF} = '';
    $self->{PATHSEQDB_UNIREFHHS} = '';

    $self->{PATHCPRODB_PDB} = '';
    $self->{PATHCPRODB_SCOP} = '';
    $self->{PATHCPRODB_PFAM} = '';
    $self->{PATHCPRODB_SwissProt} = '';
    $self->{PATHCOTHERPRODB_PDB} = '';

    $self->{INSTALL_COMER} = '';
    $self->{INSTALL_COTHER} = '';
    $self->{INSTALL_ROPIUS0} = '';
    $self->{INSTALL_HHsuite} = '';
    $self->{INSTALL_HMMER} = '';
    $self->{INSTALL_BLAST} = '';
    $self->{INSTALL_PSIPRED} = '';
    $self->{INSTALL_MODPLUS} = '';

    $self->{MAILADDRESSEE} = GetMailAddressee();
    $self->{MAILSENDER}  = GetMailSender();
    $self->{MAILSERVER}  = GetMailServer();

    while( scalar( @_ )) {
        $self->{uc( $_[0] )} = $_[1];
        shift, shift;
    }

    bless( $self, $class );
    $self->Read( $confile );
    return $self;
}

## -------------------------------------------------------------------
## read/write member methods
##

sub Path3dDb_PDB { my $self = shift; if (@_) { $self->{LOCPDBDIR} = shift } return $self->{LOCPDBDIR}; }
sub Path3dDb_SCOPe { my $self = shift; if (@_) { $self->{LOCSCOPeDIR} = shift } return $self->{LOCSCOPeDIR}; }
sub WebAddr3dDb_SCOPe { my $self = shift; if (@_) { $self->{WEBSCOPePDB} = shift } return $self->{WEBSCOPePDB}; }
sub Executable_MODELLER { my $self = shift; if (@_) { $self->{MODELLER} = shift } return $self->{MODELLER}; }

sub JobMaxNoQueries { my $self = shift; if (@_) { $self->{JOB_MAX_NQUERIES} = shift } return $self->{JOB_MAX_NQUERIES}; }
sub JobMaxNoQueriesCother { my $self = shift; if (@_) { $self->{JOB_MAX_NQUERIES_COTHER} = shift } return $self->{JOB_MAX_NQUERIES_COTHER}; }
sub JobMaxNo3DModels { my $self = shift; if (@_) { $self->{JOB_MAX_NMODELS} = shift } return $self->{JOB_MAX_NMODELS}; }
sub JobMaxNo3DTemplates { my $self = shift; if (@_) { $self->{JOB_MAX_N3DTEMPLATES} = shift } return $self->{JOB_MAX_N3DTEMPLATES}; }

sub PathSeqDb_UniRef { my $self = shift; if (@_) { $self->{PATHSEQDB_UNIREF} = shift } return $self->{PATHSEQDB_UNIREF}; }
sub PathSeqDb_UniRefHHS { my $self = shift; if (@_) { $self->{PATHSEQDB_UNIREFHHS} = shift } return $self->{PATHSEQDB_UNIREFHHS}; }

sub PathComerDb_PDB { my $self = shift; if (@_) { $self->{PATHCPRODB_PDB} = shift } return $self->{PATHCPRODB_PDB}; }
sub PathComerDb_SCOP { my $self = shift; if (@_) { $self->{PATHCPRODB_SCOP} = shift } return $self->{PATHCPRODB_SCOP}; }
sub PathComerDb_PFAM { my $self = shift; if (@_) { $self->{PATHCPRODB_PFAM} = shift } return $self->{PATHCPRODB_PFAM}; }
sub PathComerDb_SwissProt { my $self = shift; if (@_) { $self->{PATHCPRODB_SwissProt} = shift } return $self->{PATHCPRODB_SwissProt}; }
sub PathCotherDb_PDB { my $self = shift; if (@_) { $self->{PATHCOTHERPRODB_PDB} = shift } return $self->{PATHCOTHERPRODB_PDB}; }

sub InstallDir_COMER { my $self = shift; if (@_) { $self->{INSTALL_COMER} = shift } return $self->{INSTALL_COMER}; }
sub InstallDir_COTHER { my $self = shift; if (@_) { $self->{INSTALL_COTHER} = shift } return $self->{INSTALL_COTHER}; }
sub InstallDir_ROPIUS0 { my $self = shift; if (@_) { $self->{INSTALL_ROPIUS0} = shift } return $self->{INSTALL_ROPIUS0}; }
sub InstallDir_HHsuite { my $self = shift; if (@_) { $self->{INSTALL_HHsuite} = shift } return $self->{INSTALL_HHsuite}; }
sub InstallDir_HMMER { my $self = shift; if (@_) { $self->{INSTALL_HMMER} = shift } return $self->{INSTALL_HMMER}; }
sub InstallDir_BLAST { my $self = shift; if (@_) { $self->{INSTALL_BLAST} = shift } return $self->{INSTALL_BLAST}; }
sub InstallDir_PSIPRED { my $self = shift; if (@_) { $self->{INSTALL_PSIPRED} = shift } return $self->{INSTALL_PSIPRED}; }
sub InstallDir_MODPLUS { my $self = shift; if (@_) { $self->{INSTALL_MODPLUS} = shift } return $self->{INSTALL_MODPLUS}; }

sub MailAddressee { my $self = shift; if (@_) { $self->{MAILADDRESSEE} = shift } return $self->{MAILADDRESSEE}; }
sub MailSender  { my $self = shift; if (@_) { $self->{MAILSENDER} = shift }  return $self->{MAILSENDER}; }
sub MailServer  { my $self = shift; if (@_) { $self->{MAILSERVER} = shift }  return $self->{MAILSERVER}; }


sub Exists     { my $self = shift; if( @_ ) { return exists $self->{$_[0]} }return 0; }
sub GetValue   { my $self = shift; if( @_ && $self->Exists( $_[0] )) { return $self->{$_[0]} }return 0; }
sub SetValue   { my $self = shift; if( @_ && defined( $_[0] )&& defined( $_[1] )) { $self->{$_[0]} = $_[1]; } }

sub ExistsKeyByValue  { my $self = shift; if( @_ ) { my $value = shift; foreach( keys %$self ){ return 1 if $self->{$_} eq $value; } }return 0; }
sub GetKeyByValue     { my $self = shift; if( @_ ) { my $value = shift; foreach( keys %$self ){ return $_ if $self->{$_} eq $value; } }return 0; }

## -------------------------------------------------------------------
## read configuration variables
##

sub Read
{
    my  $self = shift;
    my  $confile = shift;
    my  $class = ref( $self ) || die( "ERROR: Read: Should be called by object." );
    my ($key, $value );

    unless( open( F, $confile )) {
        printf( STDERR "ERROR: Read: Failed to open %s\n", $confile );
        return( 0 );
    }

    while( <F> ) {
        chomp;
        next if /^$/ || /^\s*#/;
        s/\s*(.*)$/$1/;
        next unless /^\s*([\w\.]+)\s*=\s*(.+)$/;
        $key = $1;
        $value = $2;

        ##unless( grep {/^$key$/} keys %$self ) {
        ##    printf( STDERR "WARNING: Read: Unrecognized variable %s; skipped.\n", $key );
        ##    next;
        ##}

        if( $value =~ /^\s*['"]([^'"]*)['"]/ ) { $value = $1; }
        elsif( $value =~ /^\s*['"]?([^\s#]*)/ ) { $value = $1; }
        else {
            printf( STDERR "WARNING: Read: Unable to read value of variable.\n" );
            next;
        }
        $self->{$key} = $value;
    }

    close( F );

    return 1;
}

## -------------------------------------------------------------------

1;


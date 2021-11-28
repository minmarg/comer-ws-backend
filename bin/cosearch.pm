package cosearch;

## (C) 2020-2021 Mindaugas Margelevicius, Institute of Biotechnology, Vilnius University
## Package for conducting COMER engine-powered search given input sequences or MSAs: 
## from sequences to profile-profile alignments

use strict;
use Config;
use threads;
use threads::shared;
use FindBin;
use lib "$FindBin::Bin";
use readcfg;
use readopt;
use File::Spec;
use File::Copy;
use File::Basename;
use Scalar::Util qw(looks_like_number);
use POSIX qw(strftime);

## =============================================================================

sub new {
    my $that = shift;
    my $class = ref($that) || $that;
    my $self;

    ##constants:
    $self->{MYPROGNAME} = __PACKAGE__;
    $self->{HHBLITShelper} = File::Spec->catfile($FindBin::Bin,"hhblits_helper.sh");
    $self->{HMMERhelper} = File::Spec->catfile($FindBin::Bin,"hmmsearch_iterated.sh");
    $self->{PWFA2MSAprog} = File::Spec->catfile($FindBin::Bin,"pwfa2msa.pl");
    $self->{SENDMAILprog} = File::Spec->catfile($FindBin::Bin,"sendemail.pl");
    $self->{CFGFILE} = File::Spec->catfile($FindBin::Bin,File::Spec->updir(),"var","comer-ws-backend.conf");
    $self->{devnull} = File::Spec->devnull();
    $self->{TARPROG} = `which tar 2>$self->{devnull}`; chomp($self->{TARPROG});
    $self->{GZPROG} = `which gzip 2>$self->{devnull}`; chomp($self->{GZPROG});
    $self->{MAXNCPUs} = 6;##maximum number of CPU cores assigned for a job
    ##threshold for the number of queries which, if exceeded, implies creating threads, each assigned 1 cpu:
    ##(use a large number for serialization with multiple cores)
    $self->{MAXNQUERIES_multicore} = 6;
    ##maximum number of sequences per search program permitted to be included in the results
    $self->{MAXNSEQS_perprog} = 20000;
    $self->{MAXNQUERIES} = 100;##maximum number of queries allowed in the input
    $self->{MAXNQUERIESCOTHER} = 10;##maximum number of queries allowed in the input for COTHER

    $self->{MAXSEQLENCOTHER} = 1000;##maximum sequence length COTHER queries are allowed of

    ($self->{QRYEXT}, $self->{FASEXT}, $self->{AFAEXT}, $self->{STOEXT}, $self->{A3MEXT}, 
     $self->{PWFEXT}, $self->{COVEXT}, $self->{PROEXT}, $self->{TPROEXT}, $self->{DBEXT}, 
     $self->{CTHDBEXT}, $self->{OUTEXT}) = 
       ('qry','fa','afa','sto','a3m', 'pwfa','cov','pro','tpro','bin','tbin','json');

    $self->{NEFFEXT} = 'neff';

    $self->{INPFILENAME} = '';
    $self->{OPTFILENAME} = '';
    $self->{STAFILENAME} = '';
    $self->{COMLOGFILENAME} = '';
    $self->{RESLSTFILENAME} = '';
    $self->{RESULTSFILENAME} = '';
    $self->{ERRFILENAME} = '';
    $self->{NOPRO} = 0;
    $self->{NORUN} = 0;
    $self->{USINGSSSCORING} = 1;##using SS scoring

    while( scalar(@_)) {
        $self->{uc( $_[0] )} = $_[1];
        shift, shift;
    }

    bless( $self, $class );
    return $self;
}

## =============================================================================

sub SetOptions
{
    my  $self = shift;
    my  $mysubname = (caller(0))[3];
    my  $class = ref($self) || die("ERROR: $mysubname: Should be called by object.");
    $self->{INPFILENAME} = shift;
    $self->{OPTFILENAME} = shift;
    $self->{STAFILENAME} = shift;
    $self->{COMLOGFILENAME} = shift;
    $self->{RESLSTFILENAME} = shift;
    $self->{RESULTSFILENAME} = shift;
    $self->{ERRFILENAME} = shift;
    $self->{NOPRO} = shift;
    $self->{NORUN} = shift;
    $self->{METHOD} = shift;
}

## =============================================================================
##
sub Initialize
{
    my  $self = shift;
    my  $mysubname = (caller(0))[3];
    my  $class = ref($self) || die("ERROR: $mysubname: Should be called by object.");
    $self->Preinitialize();
    $self->ValidateOptionsFile();
    $self->CheckConfig();
    $self->InitializeOptions();
}

## =============================================================================
## preinitilalize and check paths and files
##
sub Preinitialize
{
    my  $self = shift;
    my  $mysubname = (caller(0))[3];
    my  $class = ref($self) || die("ERROR: $mysubname: Should be called by object.");

    print(STDERR "\n\n".GetDatetime()."\n\n");

    my $suffix;
    ($self->{inpbasename},$self->{inpdirname},$suffix) = fileparse($self->{INPFILENAME}, qr/\.[^.]*/);
    $self->{inpdirname} = File::Spec->rel2abs($self->{inpdirname});

    unless($self->{STAFILENAME}) {
        $self->{STAFILENAME} = File::Spec->catfile($self->{inpdirname},"$self->{inpbasename}.status");
        print(STDERR "WARNING: $self->{MYPROGNAME}: Filename for computation progress mesages not given; ".
              "it has been set to: '$self->{STAFILENAME}'\n");
    }

    unless($self->{NORUN}) {
        unless($self->{COMLOGFILENAME}) {
            $self->{COMLOGFILENAME} = File::Spec->catfile($self->{inpdirname},"$self->{inpbasename}__comer_out.log");
            print(STDERR "WARNING: $self->{MYPROGNAME}: Filename for COMER log mesages not given; ".
                  "it has been set to: '$self->{COMLOGFILENAME}'\n");
        }
        unless($self->{RESLSTFILENAME}) {
            $self->{RESLSTFILENAME} = File::Spec->catfile($self->{inpdirname},"$self->{inpbasename}__comer_out.lst");
            print(STDERR "WARNING: $self->{MYPROGNAME}: Filename for a results list not given; ".
                  "it has been set to: '$self->{RESLSTFILENAME}'\n");
        }
        unless($self->{RESULTSFILENAME}) {
            $self->{RESULTSFILENAME} = File::Spec->catfile($self->{inpdirname},"$self->{inpbasename}__comer_out.tar");
            print(STDERR "WARNING: $self->{MYPROGNAME}: Filename for compressed results not given; ".
                  "it has been set to: '$self->{RESULTSFILENAME}'\n");
        }
    }
    unless($self->{ERRFILENAME}) {
        $self->{ERRFILENAME} = File::Spec->catfile($self->{inpdirname},"$self->{inpbasename}.err");
        print(STDERR "WARNING: $self->{MYPROGNAME}: Filename for high-level error mesages not given; ".
              "it has been set to: '$self->{ERRFILENAME}'\n");
    }
    ## truncate error file before using Error!
    if(open(F,'>',$self->{ERRFILENAME})){close(F);}
    ##unless($result) {
    ##    $self->Error("ERROR: $self->{MYPROGNAME}: Error in command-line arguments.\n",
    ##            "Command-line arguments error.\n");## h-l error message
    ##    $self->MyExit(1);
    ##}
    unless($self->{INPFILENAME} && -f $self->{INPFILENAME}) {
        $self->Error("ERROR: $self->{MYPROGNAME}: Input file not found: '$self->{INPFILENAME}'\n",
                "Input file not found.\n");## h-l error message
        $self->MyExit(1);
    }
    unless($self->{OPTFILENAME} && -f $self->{OPTFILENAME}) {
        $self->Error("ERROR: $self->{MYPROGNAME}: Input job options file not found: '$self->{OPTFILENAME}'\n",
                "Job options file not found.\n");## h-l error message
        $self->MyExit(1);
    }
    unless(-f $self->{CFGFILE}) {
        $self->Error("ERROR: $self->{MYPROGNAME}: Config file not found: '$self->{CFGFILE}'\n",
                "Configuration file not found.\n");## h-l error message
        $self->MyExit(1);
    }
    unless(-f $self->{HHBLITShelper}) {
        $self->Error("ERROR: $self->{MYPROGNAME}: Program file not found: '$self->{HHBLITShelper}'\n",
                "Some of the required program files not found.\n");## h-l error message
        $self->MyExit(1);
    }
    unless(-f $self->{HMMERhelper}) {
        $self->Error("ERROR: $self->{MYPROGNAME}: Program file not found: '$self->{HMMERhelper}'\n",
                "Some of the required program files not found.\n");## h-l error message
        $self->MyExit(1);
    }
    unless(-f $self->{PWFA2MSAprog}) {
        $self->Error("ERROR: $self->{MYPROGNAME}: Program file not found: '$self->{PWFA2MSAprog}'\n",
                "Some of the required program files not found.\n");## h-l error message
        $self->MyExit(1);
    }
    unless($self->{TARPROG} && $self->{GZPROG}) {
        $self->Error("ERROR: $self->{MYPROGNAME}: System programs 'tar' and/or 'gzip' not found.\n",
                "Some of the system programs not found.\n");## h-l error message
        $self->MyExit(1);
    }
}

## =============================================================================
## check configuration variables and initilalize corresponding paths
##
sub CheckConfig
{
    my  $self = shift;
    my  $mysubname = (caller(0))[3];
    my  $class = ref($self) || die("ERROR: $mysubname: Should be called by object.");

    $self->{cfgvar} = readcfg->new($self->{CFGFILE});##backend configuration object

    if( $self->{METHOD} eq 'cother') {
        if($self->{cfgvar}->JobMaxNoQueriesCother() < 1) {
            print(STDERR "WARNING: $self->{MYPROGNAME}: Maximum number of queries for COTHER not given; ".
                  "it has been set to: '$self->{MAXNQUERIESCOTHER}'\n");
            $self->{MAXNQUERIES} = $self->{MAXNQUERIESCOTHER};
        } else {
            $self->{MAXNQUERIES} = $self->{cfgvar}->JobMaxNoQueriesCother();
        }
    } else {
        if($self->{cfgvar}->JobMaxNoQueries() < 1) {
            print(STDERR "WARNING: $self->{MYPROGNAME}: Maximum number of queries not given; ".
                  "it has been set to: '$self->{MAXNQUERIES}'\n");
        } else {
            $self->{MAXNQUERIES} = $self->{cfgvar}->JobMaxNoQueries();
        }
    }

    unless(-d $self->{cfgvar}->PathSeqDb_UniRef()) {
        $self->Error("ERROR: $self->{MYPROGNAME}: Database directory not found: '".$self->{cfgvar}->PathSeqDb_UniRef()."'\n",
                "Some of the database directories not found.\n");## h-l error message
        $self->MyExit(1);
    }
    unless(-d $self->{cfgvar}->PathSeqDb_UniRefHHS()) {
        $self->Error("ERROR: $self->{MYPROGNAME}: Database directory not found: '".$self->{cfgvar}->PathSeqDb_UniRefHHS()."'\n",
                "Some of the database directories not found.\n");## h-l error message
        $self->MyExit(1);
    }

    unless(-d $self->{cfgvar}->PathComerDb_PDB()) {
        $self->Error("ERROR: $self->{MYPROGNAME}: Database directory not found: '".$self->{cfgvar}->PathComerDb_PDB()."'\n",
                "Some of the database directories not found.\n");## h-l error message
        $self->MyExit(1);
    }
    unless(-d $self->{cfgvar}->PathComerDb_SCOP()) {
        $self->Error("ERROR: $self->{MYPROGNAME}: Database directory not found: '".$self->{cfgvar}->PathComerDb_SCOP()."'\n",
                "Some of the database directories not found.\n");## h-l error message
        $self->MyExit(1);
    }
    unless(-d $self->{cfgvar}->PathComerDb_PFAM()) {
        $self->Error("ERROR: $self->{MYPROGNAME}: Database directory not found: '".$self->{cfgvar}->PathComerDb_PFAM()."'\n",
                "Some of the database directories not found.\n");## h-l error message
        $self->MyExit(1);
    }
    unless(-d $self->{cfgvar}->PathComerDb_SwissProt()) {
        $self->Error("ERROR: $self->{MYPROGNAME}: Database directory not found: '".$self->{cfgvar}->PathComerDb_SwissProt()."'\n",
                "Some of the database directories not found.\n");## h-l error message
        $self->MyExit(1);
    }


    if( $self->{METHOD} eq 'cother') {
        unless(-d $self->{cfgvar}->PathCotherDb_PDB()) {
            $self->Error("ERROR: $self->{MYPROGNAME}: Database directory not found: '".$self->{cfgvar}->PathCotherDb_PDB()."'\n",
                "Some of the database directories not found.\n");## h-l error message
            $self->MyExit(1);
        }
    }


    unless($self->{NOPRO} && $self->{NORUN}) {
        unless(-d $self->{cfgvar}->InstallDir_COMER()) {
            $self->Error("ERROR: $self->{MYPROGNAME}: COMER installation directory not found: '".$self->{cfgvar}->InstallDir_COMER()."'\n",
                "Some of the installation directories not found.\n");## h-l error message
            $self->MyExit(1);
        }
        unless(-d $self->{cfgvar}->InstallDir_BLAST()) {
            $self->Error("ERROR: $self->{MYPROGNAME}: BLAST installation directory not found: '".$self->{cfgvar}->InstallDir_BLAST()."'\n",
                "Some of the installation directories not found.\n");## h-l error message
            $self->MyExit(1);
        }
        unless(-d $self->{cfgvar}->InstallDir_PSIPRED()) {
            $self->Error("ERROR: $self->{MYPROGNAME}: PSIPRED installation directory not found: '".$self->{cfgvar}->InstallDir_PSIPRED()."'\n",
                "Some of the installation directories not found.\n");## h-l error message
            $self->MyExit(1);
        }

        if( $self->{METHOD} eq 'cother') {
            unless(-d $self->{cfgvar}->InstallDir_COTHER()) {
                $self->Error("ERROR: $self->{MYPROGNAME}: COTHER installation directory not found: '".$self->{cfgvar}->InstallDir_COTHER()."'\n",
                    "Some of the installation directories not found.\n");## h-l error message
                $self->MyExit(1);
            }
            unless(-d $self->{cfgvar}->InstallDir_ROPIUS0()) {
                $self->Error("ERROR: $self->{MYPROGNAME}: ROPIUS0 installation directory not found: '".$self->{cfgvar}->InstallDir_ROPIUS0()."'\n",
                    "Some of the installation directories not found.\n");## h-l error message
                $self->MyExit(1);
            }
        }
    }
    unless(-d $self->{cfgvar}->InstallDir_HHsuite()) {
        $self->Error("ERROR: $self->{MYPROGNAME}: HHsuite installation directory not found: '".$self->{cfgvar}->InstallDir_HHsuite()."'\n",
                "Some of the installation directories not found.\n");## h-l error message
        $self->MyExit(1);
    }
    unless(-d $self->{cfgvar}->InstallDir_HMMER()) {
        $self->Error("ERROR: $self->{MYPROGNAME}: HMMER installation directory not found: '".$self->{cfgvar}->InstallDir_HMMER()."'\n",
                "Some of the installation directories not found.\n");## h-l error message
        $self->MyExit(1);
    }
    #unless(-d $self->{cfgvar}->InstallDir_MODPLUS()) {
    #    $self->Error("ERROR: $self->{MYPROGNAME}: ModPlus installation directory not found: '".$self->{cfgvar}->InstallDir_MODPLUS()."'\n",
    #            "Some of the installation directories not found.\n");## h-l error message
    #    $self->MyExit(1);
    #}


    ##check particular programs
    $self->{optionvalues}->{prog_comer_comer} = '';
    $self->{optionvalues}->{prog_comer_makepro} = '';
    $self->{optionvalues}->{prog_comer_makecov} = '';
    $self->{optionvalues}->{prog_comer_neff} = '';

    $self->{optionvalues}->{prog_cother_cother} = '';
    $self->{optionvalues}->{prog_cother_adddist} = '';
    $self->{optionvalues}->{prog_cother_batchadddistpy} = '';

    $self->{optionvalues}->{prog_ropius0_promage4cother_519py} = '';
    $self->{optionvalues}->{prog_ropius0_batchpred4segm_519py} = '';
    $self->{optionvalues}->{prog_ropius0_combinepredictionspy} = '';
    $self->{optionvalues}->{prog_ropius0_distopred} = '';

    unless($self->{NOPRO} && $self->{NORUN}) {
        $self->{optionvalues}->{prog_comer_comer} = File::Spec->catfile($self->{cfgvar}->InstallDir_COMER(),'bin','comer');
        $self->{optionvalues}->{prog_comer_makepro} = File::Spec->catfile($self->{cfgvar}->InstallDir_COMER(),'bin',($self->{USINGSSSCORING})? 'makepro.sh': 'makepro');
        $self->{optionvalues}->{prog_comer_makecov} = File::Spec->catfile($self->{cfgvar}->InstallDir_COMER(),'bin','makecov');
        $self->{optionvalues}->{prog_comer_neff} = File::Spec->catfile($self->{cfgvar}->InstallDir_COMER(),'bin','neff');
        unless(-f $self->{optionvalues}->{prog_comer_comer}) {
            $self->Error("ERROR: $self->{MYPROGNAME}: COMER executable not found: '".$self->{optionvalues}->{prog_comer_comer}."'\n",
                "Incomplete software installed on the system.\n");## h-l error message
            $self->MyExit(1);
        }
        unless(-f $self->{optionvalues}->{prog_comer_makepro}) {
            $self->Error("ERROR: $self->{MYPROGNAME}: COMER executable not found: '".$self->{optionvalues}->{prog_comer_makepro}."'\n",
                "Incomplete software installed on the system.\n");## h-l error message
            $self->MyExit(1);
        }
        unless(-f $self->{optionvalues}->{prog_comer_makecov}) {
            $self->Error("ERROR: $self->{MYPROGNAME}: COMER executable not found: '".$self->{optionvalues}->{prog_comer_makecov}."'\n",
                "Incomplete software installed on the system.\n");## h-l error message
            $self->MyExit(1);
        }
        unless(-f $self->{optionvalues}->{prog_comer_neff}) {
            $self->Error("ERROR: $self->{MYPROGNAME}: COMER executable not found: '".$self->{optionvalues}->{prog_comer_neff}."'\n",
                "Incomplete software installed on the system.\n");## h-l error message
            $self->MyExit(1);
        }

        if( $self->{METHOD} eq 'cother') {
            $self->{optionvalues}->{prog_cother_cother} = File::Spec->catfile($self->{cfgvar}->InstallDir_COTHER(),'bin','cother');
            $self->{optionvalues}->{prog_cother_adddist} = File::Spec->catfile($self->{cfgvar}->InstallDir_COTHER(),'bin','adddist');
            $self->{optionvalues}->{prog_cother_batchadddistpy} = File::Spec->catfile($self->{cfgvar}->InstallDir_COTHER(),'bin','batchadddist.py');
            unless(-f $self->{optionvalues}->{prog_cother_cother}) {
                $self->Error("ERROR: $self->{MYPROGNAME}: COTHER executable not found: '".$self->{optionvalues}->{prog_cother_cother}."'\n",
                    "Incomplete software installed on the system.\n");## h-l error message
                $self->MyExit(1);
            }
            unless(-f $self->{optionvalues}->{prog_cother_adddist}) {
                $self->Error("ERROR: $self->{MYPROGNAME}: COTHER executable not found: '".$self->{optionvalues}->{prog_cother_adddist}."'\n",
                    "Incomplete software installed on the system.\n");## h-l error message
                $self->MyExit(1);
            }
            unless(-f $self->{optionvalues}->{prog_cother_batchadddistpy}) {
                $self->Error("ERROR: $self->{MYPROGNAME}: COTHER executable not found: '".$self->{optionvalues}->{prog_cother_batchadddistpy}."'\n",
                    "Incomplete software installed on the system.\n");## h-l error message
                $self->MyExit(1);
            }

            $self->{optionvalues}->{prog_ropius0_promage4cother_519py} = File::Spec->catfile($self->{cfgvar}->InstallDir_ROPIUS0(),'infer','promage4cother_519.py');
            $self->{optionvalues}->{prog_ropius0_batchpred4segm_519py} = File::Spec->catfile($self->{cfgvar}->InstallDir_ROPIUS0(),'infer','batchpred4segm_519.py');
            $self->{optionvalues}->{prog_ropius0_combinepredictionspy} = File::Spec->catfile($self->{cfgvar}->InstallDir_ROPIUS0(),'bin','combinepredictions.py');
            $self->{optionvalues}->{prog_ropius0_distopred} = File::Spec->catfile($self->{cfgvar}->InstallDir_ROPIUS0(),'srvs','distopred.sh');
            unless(-f $self->{optionvalues}->{prog_ropius0_promage4cother_519py}) {
                $self->Error("ERROR: $self->{MYPROGNAME}: ROPIUS0 executable not found: '".$self->{optionvalues}->{prog_ropius0_promage4cother_519py}."'\n",
                    "Incomplete software installed on the system.\n");## h-l error message
                $self->MyExit(1);
            }
            unless(-f $self->{optionvalues}->{prog_ropius0_batchpred4segm_519py}) {
                $self->Error("ERROR: $self->{MYPROGNAME}: ROPIUS0 executable not found: '".$self->{optionvalues}->{prog_ropius0_batchpred4segm_519py}."'\n",
                    "Incomplete software installed on the system.\n");## h-l error message
                $self->MyExit(1);
            }
            unless(-f $self->{optionvalues}->{prog_ropius0_combinepredictionspy}) {
                $self->Error("ERROR: $self->{MYPROGNAME}: ROPIUS0 executable not found: '".$self->{optionvalues}->{prog_ropius0_combinepredictionspy}."'\n",
                    "Incomplete software installed on the system.\n");## h-l error message
                $self->MyExit(1);
            }
            unless(-f $self->{optionvalues}->{prog_ropius0_distopred}) {
                $self->Error("ERROR: $self->{MYPROGNAME}: ROPIUS0 executable not found: '".$self->{optionvalues}->{prog_ropius0_distopred}."'\n",
                    "Incomplete software installed on the system.\n");## h-l error message
                $self->MyExit(1);
            }
        }
    }

    $self->{optionvalues}->{prog_hhsuite_hhblits} = File::Spec->catfile($self->{cfgvar}->InstallDir_HHsuite(),'bin','hhblits');
    $self->{optionvalues}->{prog_hhsuite_reformat} = File::Spec->catfile($self->{cfgvar}->InstallDir_HHsuite(),'scripts','reformat.pl');
    unless(-f $self->{optionvalues}->{prog_hhsuite_hhblits}) {
        $self->Error("ERROR: $self->{MYPROGNAME}: HHsuite executable not found: '".$self->{optionvalues}->{prog_hhsuite_hhblits}."'\n",
            "Incomplete software installed on the system.\n");## h-l error message
        $self->MyExit(1);
    }
    unless(-f $self->{optionvalues}->{prog_hhsuite_reformat}) {
        $self->Error("ERROR: $self->{MYPROGNAME}: HHsuite executable not found: '".$self->{optionvalues}->{prog_hhsuite_reformat}."'\n",
            "Incomplete software installed on the system.\n");## h-l error message
        $self->MyExit(1);
    }

    $self->{optionvalues}->{prog_hmmer_jackhmmer} = File::Spec->catfile($self->{cfgvar}->InstallDir_HMMER(),'bin','jackhmmer');
    $self->{optionvalues}->{prog_hmmer_hmmsearch} = File::Spec->catfile($self->{cfgvar}->InstallDir_HMMER(),'bin','hmmsearch');
    $self->{optionvalues}->{prog_hmmer_hmmbuild} = File::Spec->catfile($self->{cfgvar}->InstallDir_HMMER(),'bin','hmmbuild');
    unless(-f $self->{optionvalues}->{prog_hmmer_jackhmmer}) {
        $self->Error("ERROR: $self->{MYPROGNAME}: HMMER executable not found: '".$self->{optionvalues}->{prog_hmmer_jackhmmer}."'\n",
            "Incomplete software installed on the system.\n");## h-l error message
        $self->MyExit(1);
    }
    unless(-f $self->{optionvalues}->{prog_hmmer_hmmsearch}) {
        $self->Error("ERROR: $self->{MYPROGNAME}: HMMER executable not found: '".$self->{optionvalues}->{prog_hmmer_hmmsearch}."'\n",
            "Incomplete software installed on the system.\n");## h-l error message
        $self->MyExit(1);
    }
    unless(-f $self->{optionvalues}->{prog_hmmer_hmmbuild}) {
        $self->Error("ERROR: $self->{MYPROGNAME}: HMMER executable not found: '".$self->{optionvalues}->{prog_hmmer_hmmbuild}."'\n",
            "Incomplete software installed on the system.\n");## h-l error message
        $self->MyExit(1);
    }
}

## =============================================================================
#$ helper function for validating values in the options file
sub ValidateHelper
{
    my  $self = shift;
    my  $optname = shift;
    my  $optcurvalue = shift;
    my  $optminvalue = shift;
    my  $optmaxvalue = shift;
    my  $optdefvalue = shift;
    my  $mysubname = (caller(0))[3];
    my  $class = ref($self) || die("ERROR: $mysubname: Should be called by object.");
    my  $warn = '';
    if(!looks_like_number($optcurvalue) || $optcurvalue < $optminvalue || $optcurvalue > $optmaxvalue) {
        $warn = "\nWARNING: Disallowed values: Option $optname changed: $optcurvalue -> $optdefvalue\n";
        $self->Warning($warn, $warn, 1);##no e-mail
        $_ = "$optname = $optdefvalue\n";
    }
}
## validate options file and modify illegal values
##
sub ValidateOptionsFile
{
    my  $self = shift;
    my  $mysubname = (caller(0))[3];
    my  $class = ref($self) || die("ERROR: $mysubname: Should be called by object.");
    my ($contents, $warn, $value) = ('','',0);

    return unless(open(F, $self->{OPTFILENAME}));
    while(<F>) {
        if(/\s*EVAL\s*=\s*(\S+)/) { $self->ValidateHelper("EVAL", $1, 0, 100, 10); }
        elsif(/\s*NOHITS\s*=\s*(\S+)/) { $self->ValidateHelper("NOHITS", $1, 1, 2000, 700); }
        elsif(/\s*NOALNS\s*=\s*(\S+)/) { $self->ValidateHelper("NOALNS", $1, 1, 2000, 700); }
        elsif(/\s*ADJWGT\s*=\s*(\S+)/) { $self->ValidateHelper("ADJWGT", $1, 0.001, 0.999, 0.33); }
        elsif(/\s*CVSWGT\s*=\s*(\S+)/) { $self->ValidateHelper("CVSWGT", $1, 0, 0.999, 0.15); }
        elsif(/\s*SSSWGT\s*=\s*(\S+)/) { 
            if($1=~/^0*\.?0+$/) { $self->{USINGSSSCORING} = 0; }
            else { $self->ValidateHelper("SSSWGT", $1, 0, 0.999, 0.12); }
        }
        elsif(/\s*DDMSWGT\s*=\s*(\S+)/) { $self->ValidateHelper("DDMSWGT", $1, 0, 0.999, 0.2); }
        elsif(/\s*LCFILTEREACH\s*=\s*(\S+)/) { $self->ValidateHelper("LCFILTEREACH", $1, 0, 1, 1); }
        elsif(/\s*MINPP\s*=\s*(\S+)/) { $self->ValidateHelper("MINPP", $1, 0, 0.999, 0.28); }
        ##
        elsif(/\s*(hhsuite_opt_niterations)\s*=\s*(\S+)/) { $self->ValidateHelper($1, $2, 1, 4, 2); }
        elsif(/\s*(hhsuite_opt_evalue)\s*=\s*(\S+)/) { $self->ValidateHelper($1, $2, 0, 1, 0.001); }
        elsif(/\s*(hmmer_opt_niterations)\s*=\s*(\S+)/) { $self->ValidateHelper($1, $2, 1, 4, 2); }
        elsif(/\s*(hmmer_opt_evalue)\s*=\s*(\S+)/) { $self->ValidateHelper($1, $2, 0, 1, 0.001); }
        $contents .= $_
    }
    close(F);

    ## no error check
    move($self->{OPTFILENAME},"$self->{OPTFILENAME}.org");

    unless(open(F, ">", $self->{OPTFILENAME})) {
        $self->Error("ERROR: $self->{MYPROGNAME}: $mysubname: ".
                "Failed to open job options file: '$self->{OPTFILENAME}'\n",
                "Faled to open job options file.\n");## h-l error message
        $self->MyExit(1);
    }
    print(F $contents);
    close(F);
}

## =============================================================================
## initialize options for the job
##
sub InitializeOptions
{
    my  $self = shift;
    my  $mysubname = (caller(0))[3];
    my  $class = ref($self) || die("ERROR: $mysubname: Should be called by object.");

    $self->{options} = readopt->new($self->{OPTFILENAME});##options object
    $self->{optionkeys} = {
        job_num_cpus => 'JOB_NUM_CPUS',
        comer_db => 'comer_db',
        cother_db => 'cother_db',
        sequence_db => 'sequence_db',
        hhsuite_db => 'hhsuite_db',

        hhsuite_in_use => 'hhsuite_in_use',
        hhsuite_opt_niterations => 'hhsuite_opt_niterations',
        hhsuite_opt_evalue => 'hhsuite_opt_evalue',

        hmmer_in_use => 'hmmer_in_use',
        hmmer_opt_niterations => 'hmmer_opt_niterations',
        hmmer_opt_evalue => 'hmmer_opt_evalue'
    };
    my  $nsyscpus = $self->GetNCPUs();
    $self->{ncpus} = $self->{MAXNCPUs};

    if($Config{useithreads}) {
        if($self->{cfgvar}->Exists($self->{optionkeys}->{job_num_cpus})) {
            $self->{ncpus} = $self->{cfgvar}->GetValue($self->{optionkeys}->{job_num_cpus});
            $self->{ncpus} = $nsyscpus if($nsyscpus > 0 && $nsyscpus < $self->{ncpus});
        }
        else {
            print(STDERR "WARNING: $self->{MYPROGNAME}: Job option $self->{optionkeys}->{job_num_cpus} not specified: ".
                "Using default #cpus=$self->{ncpus}\n");
        }
    }
    else {
        $self->{ncpus} = 1;
        print(STDERR "WARNING: $self->{MYPROGNAME}: Perl compiled WITHOUT thread support: ".
            "Using only one cpu: #cpus=$self->{ncpus}\n");
    }
}

## =============================================================================
##







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

## -----------------------------------------------------------------------------
## GetNCPUs: get the number of CPU cores in the system; return <0 on error;
##
sub GetNCPUs
{
    my  $self = shift;
    my  $mysubname = (caller(0))[3];
    my  $class = ref($self) || die("ERROR: $mysubname: Should be called by object.");

    unless(open(H, "/proc/cpuinfo")) {
        print(STDERR "WARNING: $self->{MYPROGNAME}: Unable to determine #cpu cores.\n");
        return -1;##failed to open cpuinfo: $!
    }
    my $ncpus = scalar(grep(/^processor/,<H>)); 
    close(H);
    return $ncpus;
}

## -----------------------------------------------------------------------------
## ProgressMsg: Print progress message to file
##
sub ProgressMsg
{
    my  $self = shift;
    my  $msg = shift;
    my  $mysubname = (caller(0))[3];
    my  $class = ref($self) || die("ERROR: $mysubname: Should be called by object.");
    return 0 unless(open(F,'>>',$self->{STAFILENAME}));
    print(F $msg);
    close(F);
    return 1;
}

## -----------------------------------------------------------------------------
## GetWarnings: extract warnings from a log file
##
sub GetWarnings
{
    my  $self = shift;
    my  $comerlogfile = shift;##comer search log file
    my  $rlist = shift;##ref to the output list of warnings
    my  $mysubname = (caller(0))[3];
    my  $class = ref($self) || die("ERROR: $mysubname: Should be called by object.");
    unless(open(F, $comerlogfile)) {
        print(STDERR "WARNING: $self->{MYPROGNAME}: Unable to open COMER based tool's log file: '$comerlogfile'.\n");
        return 0;
    }
    @$rlist = grep(/\sWARNING/,<F>);
    close(F);
    return 1;
}

## -----------------------------------------------------------------------------
## GetProfileLength: extract profile length from profile file
##
sub GetProfileLength
{
    my  $self = shift;
    my  $profile = shift;##profile pathname
    my  $rlength = shift;##ref to the profile length
    my  $mysubname = (caller(0))[3];
    my  $class = ref($self) || die("ERROR: $mysubname: Should be called by object.");
    $$rlength = 0;
    unless(open(F, $profile)) {
        print(STDERR "WARNING: $self->{MYPROGNAME}: Unable to open profile: '$profile'.\n");
        return 0;
    }
    while(<F>) {
        do {$$rlength = $1; last} if(/^LEN:\s+(\d+)/);
    }
    close(F);
    return 1;
}

## -----------------------------------------------------------------------------
## ChangeProfileFileField: change the FILE field of a COMER profile
##
sub ChangeProfileFileField
{
    my  $self = shift;
    my  $profile = shift;##comer profile
    my  $outprofile = shift;##output comer profile
    my  $fieldvalue = shift;##new value of the FILE field
    my  $mysubname = (caller(0))[3];
    my  $class = ref($self) || die("ERROR: $mysubname: Should be called by object.");
    my  $procontents = '';
    unless(open(F, $profile)) {
        print(STDERR "ERROR: $self->{MYPROGNAME}: ChangeProfileFileField: Failed to open file: '$profile'.\n\n");
        return 0;
    }
    while(<F>) {
        s/^(FILE:).*$/$1 $fieldvalue/;
        $procontents .= $_;
    }
    close(F);
    unless(open(F, '>', $outprofile)) {
        print(STDERR "ERROR: $self->{MYPROGNAME}: ChangeProfileFileField: Failed to open file for writing: '$profile'.\n\n");
        return 0;
    }
    print F $procontents;
    close(F);
    return 1;
}

## -----------------------------------------------------------------------------
## AddFileToArchive: add a file to a given archive
##
sub AddFileToArchive
{
    my  $self = shift;
    my  $filepathname = shift;##full pathname to the file
    my  $dirname = shift;##name of directory where $filepathname is; may be empty
    my  $rerrmsg = shift;##ref to the error message string to be put in logs
    my  $rhlerrmsg = shift;##ref to the h-l error message string
    my  $create = shift;##flag of whether the archive is to be created

    my  $mysubname = (caller(0))[3];
    my  $class = ref($self) || die("ERROR: $mysubname: Should be called by object.");
    my  $preamb = "[ ${mysubname} ] ";
    my  $archive = $self->{RESULTSFILENAME};
    my  $command = '';
    my  $opt = $create?'c':'r';
    my  $ret = 1;

    if($dirname) {
        $command = "$self->{TARPROG} -${opt}f \"${archive}\" -C \"${dirname}\" \"${filepathname}\"";
    }
    else {
        my $filedirname = dirname($filepathname);
        my $filename = basename($filepathname);
        $command = "$self->{TARPROG} -${opt}f \"${archive}\" -C \"${filedirname}\" \"${filename}\"";
    }

    print(STDERR GetTime()."${preamb} ${command}\n");

    unless($self->ExecCommand($command)) {
        $$rerrmsg = "ERROR: $self->{MYPROGNAME}: $mysubname: Failed to add file to archive: '${filepathname}'\n";
        $$rhlerrmsg = "Failed to archive some of the results files.\n";
        return 0;
    }

    return $ret;
}

## -----------------------------------------------------------------------------
## VerifyOptionValues: verify directory and filename information given in the 
## job options file
##
sub VerifyOptionValues 
{
    my  $self = shift;
    my  $rerrmsg = shift;##ref to the error message string
    my  $rhlerrmsg = shift;##ref to the h-l error message string

    my  $mysubname = (caller(0))[3];
    my  $class = ref($self) || die("ERROR: $mysubname: Should be called by object.");
    my ($filename, $fullname, $fullnameext);
    my  $ret = 1;


    unless($self->{NORUN}) {
        if( $self->{METHOD} eq 'cother')
        {
            unless($self->{options}->Exists($self->{optionkeys}->{cother_db})) {
                $$rerrmsg = "ERROR: $self->{MYPROGNAME}: $mysubname: ".
                    "Option $self->{optionkeys}->{cother_db} not specified in the job options file.\n";
                $$rhlerrmsg = "COTHER profile database not specified.\n";
                return 0;
            }

            $self->{optionvalues}->{$self->{optionkeys}->{cother_db}} = '';
            $filename = $self->{options}->GetValue($self->{optionkeys}->{cother_db});
            my @fnames = split(',', $filename);
            my @fullnamelist;
            foreach $filename(@fnames) {
                $fullname = File::Spec->catfile($self->{cfgvar}->PathCotherDb_PDB(),${filename});
                $fullnameext = "${fullname}.$self->{CTHDBEXT}";

                if(-f $fullnameext) {
                    push @fullnamelist, $fullname;
                } else {
                    my $text = "WARNING: COTHER profile database not found: $filename\n";
                    $self->Warning($text, $text, 1);##no e-mail
                }
            }

            if($#fullnamelist < 0) {
                $$rerrmsg = "ERROR: $self->{MYPROGNAME}: $mysubname: ".
                  "All COTHER db(s) '${filename}' not found in the COTHER db directories.\n";
                $$rhlerrmsg = "COTHER profile database(s) not found.\n";
                return 0;
            }

            $self->{optionvalues}->{$self->{optionkeys}->{cother_db}} = join(',',@fullnamelist);
        }
        else {## method comer2
            unless($self->{options}->Exists($self->{optionkeys}->{comer_db})) {
                $$rerrmsg = "ERROR: $self->{MYPROGNAME}: $mysubname: ".
                    "Option $self->{optionkeys}->{comer_db} not specified in the job options file.\n";
                $$rhlerrmsg = "COMER profile database not specified.\n";
                return 0;
            }

            $self->{optionvalues}->{$self->{optionkeys}->{comer_db}} = '';
            $filename = $self->{options}->GetValue($self->{optionkeys}->{comer_db});
            my @fnames = split(',', $filename);
            my @fullnamelist;
            foreach $filename(@fnames) {
                $fullname = File::Spec->catfile($self->{cfgvar}->PathComerDb_PDB(),${filename});
                $fullnameext = "${fullname}.$self->{DBEXT}";

                if(-f $fullnameext) {
                    push @fullnamelist, $fullname;
                } else {
                    $fullname = File::Spec->catfile($self->{cfgvar}->PathComerDb_SCOP(),${filename});
                    $fullnameext = "${fullname}.$self->{DBEXT}";
                    if(-f $fullnameext) {
                        push @fullnamelist, $fullname;
                    } else {
                        $fullname = File::Spec->catfile($self->{cfgvar}->PathComerDb_PFAM(),${filename});
                        $fullnameext = "${fullname}.$self->{DBEXT}";
                        if(-f $fullnameext) {
                            push @fullnamelist, $fullname;
                        } else {
                            $fullname = File::Spec->catfile($self->{cfgvar}->PathComerDb_SwissProt(),${filename});
                            $fullnameext = "${fullname}.$self->{DBEXT}";
                            if(-f $fullnameext) {
                                push @fullnamelist, $fullname;
                            } else {
                                my $text = "WARNING: COMER profile database not found: $filename\n";
                                $self->Warning($text, $text, 1);##no e-mail
                            }
                        }
                    }
                }
            }

            if($#fullnamelist < 0) {
                $$rerrmsg = "ERROR: $self->{MYPROGNAME}: $mysubname: ".
                  "All COMER db(s) '${filename}' not found in the COMER db directories.\n";
                $$rhlerrmsg = "COMER profile database(s) not found.\n";
                return 0;
            }

            $self->{optionvalues}->{$self->{optionkeys}->{comer_db}} = join(',',@fullnamelist);
        }
    }


    $self->{optionvalues}->{$self->{optionkeys}->{hhsuite_in_use}} = 
        ($self->{options}->Exists($self->{optionkeys}->{hhsuite_in_use})? 
            $self->{options}->GetValue($self->{optionkeys}->{hhsuite_in_use}): 0);

    if($self->{optionvalues}->{$self->{optionkeys}->{hhsuite_in_use}}) {
        if(!$self->{options}->Exists($self->{optionkeys}->{hhsuite_opt_niterations}) || 
            $self->{options}->GetValue($self->{optionkeys}->{hhsuite_opt_niterations}) !~ /^0*[1-9]\d*\s*/)
        {
            $$rerrmsg = "ERROR: $self->{MYPROGNAME}: $mysubname: ".
                "Option $self->{optionkeys}->{hhsuite_opt_niterations} either undefined or invalid.\n";
            $$rhlerrmsg = "A job option, number of iterations, is invalid.\n";
            return 0;
        }

        $self->{optionvalues}->{$self->{optionkeys}->{hhsuite_opt_niterations}} = 
            $self->{options}->GetValue($self->{optionkeys}->{hhsuite_opt_niterations});

        if(!$self->{options}->Exists($self->{optionkeys}->{hhsuite_opt_evalue}) || 
            $self->{options}->GetValue($self->{optionkeys}->{hhsuite_opt_evalue}) !~ /^\d*\.?\d*(?:e[\+\-])?\d+\s*/)
        {
            $$rerrmsg = "ERROR: $self->{MYPROGNAME}: $mysubname: ".
                "Option $self->{optionkeys}->{hhsuite_opt_evalue} either undefined or invalid.\n";
            $$rhlerrmsg = "A job option, E-value, is invalid.\n";
            return 0;
        }

        $self->{optionvalues}->{$self->{optionkeys}->{hhsuite_opt_evalue}} =
            $self->{options}->GetValue($self->{optionkeys}->{hhsuite_opt_evalue});

        unless($self->{options}->Exists($self->{optionkeys}->{hhsuite_db})) {
            $$rerrmsg = "ERROR: $self->{MYPROGNAME}: $mysubname: ".
                "Option $self->{optionkeys}->{hhsuite_db} not specified in the job options file.\n";
            $$rhlerrmsg = "HHsuite database not specified.\n";
            return 0;
        }
        
        $self->{optionvalues}->{$self->{optionkeys}->{hhsuite_db}} = '';
        $filename = $self->{options}->GetValue($self->{optionkeys}->{hhsuite_db});
        $fullname = File::Spec->catfile($self->{cfgvar}->PathSeqDb_UniRefHHS(),$filename);
        
        if(scalar(glob("${fullname}*"))) {
            $self->{optionvalues}->{$self->{optionkeys}->{hhsuite_db}} = $fullname;
        }
        
        unless($self->{optionvalues}->{$self->{optionkeys}->{hhsuite_db}}) {
            $$rerrmsg = "ERROR: $self->{MYPROGNAME}: $mysubname: ".
                "HHsuite db '${filename}*' not found in all HHsuite db directories.\n";
            $$rhlerrmsg = "HHsuite database not found.\n";
            return 0;
        }
    }


    $self->{optionvalues}->{$self->{optionkeys}->{hmmer_in_use}} = 
        ($self->{options}->Exists($self->{optionkeys}->{hmmer_in_use})? 
            $self->{options}->GetValue($self->{optionkeys}->{hmmer_in_use}): 0);

    if($self->{optionvalues}->{$self->{optionkeys}->{hmmer_in_use}}) {
        if(!$self->{options}->Exists($self->{optionkeys}->{hmmer_opt_niterations}) || 
            $self->{options}->GetValue($self->{optionkeys}->{hmmer_opt_niterations}) !~ /^0*[1-9]\d*\s*/)
        {
            $$rerrmsg = "ERROR: $self->{MYPROGNAME}: $mysubname: ".
                "Option $self->{optionkeys}->{hmmer_opt_niterations} either undefined or invalid.\n";
            $$rhlerrmsg = "A job option, number of iterations, is invalid.\n";
            return 0;
        }

        $self->{optionvalues}->{$self->{optionkeys}->{hmmer_opt_niterations}} = 
            $self->{options}->GetValue($self->{optionkeys}->{hmmer_opt_niterations});

        if(!$self->{options}->Exists($self->{optionkeys}->{hmmer_opt_evalue}) || 
            $self->{options}->GetValue($self->{optionkeys}->{hmmer_opt_evalue}) !~ /^\d*\.?\d*(?:e[\+\-])?\d+\s*/)
        {
            $$rerrmsg = "ERROR: $self->{MYPROGNAME}: $mysubname: ".
                "Option $self->{optionkeys}->{hmmer_opt_evalue} either undefined or invalid.\n";
            $$rhlerrmsg = "A job option, E-value, is invalid.\n";
            return 0;
        }

        $self->{optionvalues}->{$self->{optionkeys}->{hmmer_opt_evalue}} =
            $self->{options}->GetValue($self->{optionkeys}->{hmmer_opt_evalue});

        unless($self->{options}->Exists($self->{optionkeys}->{sequence_db})) {
            $$rerrmsg = "ERROR: $self->{MYPROGNAME}: $mysubname: ".
                "Option $self->{optionkeys}->{sequence_db} not specified in the job options file.\n";
            $$rhlerrmsg = "A sequence database is not specified.\n";
            return 0;
        }
        
        $self->{optionvalues}->{$self->{optionkeys}->{sequence_db}} = '';
        $filename = $self->{options}->GetValue($self->{optionkeys}->{sequence_db});
        $fullname = File::Spec->catfile($self->{cfgvar}->PathSeqDb_UniRef(),$filename);
        
        if(-f $fullname) {
            $self->{optionvalues}->{$self->{optionkeys}->{sequence_db}} = $fullname;
        }
        
        unless($self->{optionvalues}->{$self->{optionkeys}->{sequence_db}}) {
            $$rerrmsg = "ERROR: $self->{MYPROGNAME}: $mysubname: ".
                "Sequence db '${filename}' not found in all sequence db directories.\n";
            $$rhlerrmsg = "Sequence database not found.\n";
            return 0;
        }
    }

    return $ret;
}

## -----------------------------------------------------------------------------
## DistrInputToSubdirs: parse the input and distribute individual queries to 
## subdirectories;
##
sub DistrInputToSubdirs
{
    my  $self = shift;
    my  $rnqueries = shift;##ref to the number of queries
    my  $rinputs = shift;##ref to the hash of individual inputs
    my  $rinfmsg = shift;##ref to an information message
    my  $rerrmsg = shift;##ref to the error message string to be put in logs
    my  $rhlerrmsg = shift;##ref to the h-l error message string

    my  $mysubname = (caller(0))[3];
    my  $class = ref($self) || die("ERROR: $mysubname: Should be called by object.");
    my  $inputsep = qr/^\/\//;
    my ($ndescs) = (0);
    my ($qrycontents, $qrysubdir, $qrybasename, $qryfullname, $ext) = ('','','','','');
    my ($lastext, $seqn, $prevseqn) = ('','','');
    my  $ret = 1;

    unless(open(F, $self->{INPFILENAME})) {
        $$rerrmsg = "ERROR: $self->{MYPROGNAME}: $mysubname: Failed to open the input file: '$self->{INPFILENAME}'\n";
        $$rhlerrmsg = "Input file not found.\n";
        return 0;
    }
    $$rnqueries = 0;
    while(<F>) {
        next if(!eof(F) && /^\s*$/);
        if(eof(F) || !/$inputsep/) {
            $ext = $self->{STOEXT} if(!$qrycontents && /^#\s+STOCKHOLM/);
            $ext = $self->{PROEXT} if(!$qrycontents && /^COMER\s+profile/);
            $ext = $self->{TPROEXT} if(!$qrycontents && /^COTHER\s+profile/);
            if(/^>(?:ss_dssp|ss_pred|ss_conf)/) {
                $ext = $self->{A3MEXT};
            }
            $qrycontents .= $_ unless(/$inputsep/);
            $ndescs++ if(/^\s*>/);
            if(eof(F) || !/^\s*>/) {
                ##NOTE: $_ changes here
                chomp;
                s/^\s+(.+)\s+$/$1/;
                $seqn .= $_ unless(/$inputsep/);
            }
        }
        if(eof(F) || /$inputsep/ || /^\s*>/) {
            $prevseqn = $seqn unless $prevseqn;
            if(!$ext && length($prevseqn) != length($seqn)) {
                $ext = $self->{A3MEXT};
            }
            $seqn = '';
        }
        if(eof(F) || /$inputsep/) {
            next unless $qrycontents;
            if($self->{MAXNQUERIES} <= $$rnqueries) {
                my $text = "\nWARNING: Number of queries ".
                    "reduced to the maximum allowed: $self->{MAXNQUERIES}\n";
                $self->Warning($text, $text, 1);##no e-mail
                last;
            }
            if($ext && $ext eq $self->{STOEXT}) {##STOCKHOLM format
                $qrycontents .= $_;
                $$rinfmsg = 'MSA STOCKHOLM format' unless $lastext;
            } elsif($ext && $ext eq $self->{PROEXT}) {##COMER profile
                $$rinfmsg = 'profile COMER format' unless $lastext;
            } elsif($ext && $ext eq $self->{TPROEXT}) {##COTHER profile
                $$rinfmsg = 'profile COTHER format' unless $lastext;
            } elsif($ext && $ext eq $self->{A3MEXT}) {##A3M profile
                $qrycontents =~ s/^\s*>\s*\n/>Query_${$rnqueries} (unnamed)\n/;##added to avoid empty ids
                $$rinfmsg = 'MSA A3M format' unless $lastext;
            } elsif($ndescs < 1) {##assume a plain sequence: make FASTA
                $qrycontents = ">Query_${$rnqueries} (unnamed)\n".$qrycontents;
                $ext = $self->{FASEXT};
                $$rinfmsg = 'sequence FASTA/plain format' unless $lastext;
            } elsif($ndescs == 1) {##sequence in FASTA
                $qrycontents =~ s/^\s*>\s*\n/>Query_${$rnqueries} (unnamed)\n/;##added to avoid empty ids
                $ext = $self->{FASEXT};
                $$rinfmsg = 'sequence FASTA format' unless $lastext;
            } else {##MSA in FASTA
                $ext = $self->{AFAEXT};
                $qrycontents =~ s/^\s*>\s*\n/>Query_${$rnqueries} (unnamed)\n/;##added to avoid empty ids
                $$rinfmsg = 'MSA aligned FASTA format' unless $lastext;
            }
            $lastext = $ext unless $lastext;
            $$rinfmsg = '  The queries in the input have different formats' if $lastext cmp $ext;
            $qrysubdir = File::Spec->catfile($self->{inpdirname},"$self->{inpbasename}__${$rnqueries}");
            $qrybasename = File::Spec->catfile($qrysubdir,"$self->{inpbasename}__${$rnqueries}");
            $qryfullname = "${qrybasename}.${ext}";
            ;;
            unless( -d $qrysubdir || mkdir($qrysubdir)) {
                $$rerrmsg = "ERROR: $self->{MYPROGNAME}: $mysubname: Failed to create directory: '$qrysubdir'\n";
                $$rhlerrmsg = "Creating a directory failed.\n";
                $ret = 0;
                last;
            }
            unless(open(FIN, ">", $qryfullname)) {
                $$rerrmsg = "ERROR: $self->{MYPROGNAME}: $mysubname: Failed to open file for writing: '$qryfullname'\n";
                $$rhlerrmsg = "Failed to create a file.\n";
                $ret = 0;
                last;
            }
            unless(print(FIN $qrycontents)) {
                $$rerrmsg = "ERROR: $self->{MYPROGNAME}: $mysubname: Failed to write to file: '$qryfullname'\n";
                $$rhlerrmsg = "Write to a file failed.\n";
                close(FIN);
                $ret = 0;
                last;
            }
            close(FIN);
            $$rinputs{"${$rnqueries}_rcode"} = 0;
            $$rinputs{"${$rnqueries}_error"} = '';##error message
            $$rinputs{"${$rnqueries}_errhl"} = '';##high-level error message
            $$rinputs{"${$rnqueries}_isize"} = length($qrycontents);
            $$rinputs{"${$rnqueries}_input"} = $qryfullname;
            $$rinputs{"${$rnqueries}_bname"} = $qrybasename;
            $$rinputs{"${$rnqueries}_logfl"} = "${qrybasename}.log";
            $$rinputs{"${$rnqueries}_msafl"} = '';
            $$rinputs{"${$rnqueries}_profl"} = '';
            $$rinputs{"${$rnqueries}_nefffl"} = '';##neff file
            $$rinputs{"${$rnqueries}_tprof"} = '';##COTHER profile when it is in use
            $$rinputs{"${$rnqueries}_outfl"} = '';
            $$rinputs{"${$rnqueries}_hhsrun"} = 0;##hhsuite has not been run
            $$rinputs{"${$rnqueries}_hmmrun"} = 0;##hmmer has not been run
            ;;
            $qrycontents = $qrysubdir = $qrybasename = $qryfullname = $ext = '';
            $seqn = $prevseqn = '';
            $ndescs = 0;
            $$rnqueries++;
            next;
        }
    }
    close(F);

    if($$rnqueries < 1) {
        $$rerrmsg = $$rhlerrmsg = '' if $ret;
        $$rerrmsg .= "ERROR: $self->{MYPROGNAME}: $mysubname: No queries found in the input.\n";
        $$rhlerrmsg .= "Invalid input format: No queries.\n";
        return 0;
    }

    if($$rnqueries == 1) {
        $$rinfmsg =~ s/^(MSA|sequence|profile)(.+)$/  The $1 is in$2/;
    } else {
        $$rinfmsg =~ s/^(MSA|sequence|profile)(.+)$/  The $1s in the input are all in$2/;
    }

    return $ret;
}

## -----------------------------------------------------------------------------
## CheckForErrorsandWarnings: check for error codes and any messages that could 
## have been returned from threads
##
sub CheckForErrorsandWarnings
{
    my  $self = shift;
    my  $nqueries = shift;##number of queries
    my  $rinputs = shift;##ref to the hash of individual inputs

    my  @qrynums = 0..$nqueries-1;

    for(my $q = 0; $q <= $#qrynums; $q++) {
        my $qnum = $qrynums[$q];
        ##print errors issued by threads if any
        if($$rinputs{"${qnum}_rcode"}) { ##don't send an e-mail
            $self->Error($$rinputs{"${qnum}_error"}, $$rinputs{"${qnum}_errhl"}, 1);
            next;
        }
        elsif($$rinputs{"${qnum}_errhl"}) {
            ##no error, but messages can be present; print them
            $self->Warning($$rinputs{"${qnum}_error"}, $$rinputs{"${qnum}_errhl"}, 1);##no e-mail
        }
    }

}

## -----------------------------------------------------------------------------
## RunCOMER: put all profiles together in one file and run COMER search
##
sub RunCOMER
{
    my  $self = shift;
    my  $nqueries = shift;##number of queries
    my  $rinputs = shift;##ref to the hash of individual inputs
    my  $rerrmsg = shift;##ref to the error message string to be put in logs
    my  $rhlerrmsg = shift;##ref to the h-l error message string

    my  $mysubname = (caller(0))[3];
    my  $class = ref($self) || die("ERROR: $mysubname: Should be called by object.");
    my  $preamb = "[ ${mysubname} ] ";
    my  @qrynums = 0..$nqueries-1;
    my  $fileofprofiles = File::Spec->catfile($self->{inpdirname},"$self->{inpbasename}__nqries${nqueries}.$self->{PROEXT}");
    my  $outputpfx = File::Spec->catfile($self->{inpdirname},"$self->{inpbasename}__comer_out");
    my  $outputsubdir = $outputpfx;
    my  $comeropt_pass2memp = 50;
    ##no. trials to launch comer; NOTE: may fail when two or more processes request resources at the same time:
    my  $ntrials = 3;
    my  $command = '';
    my  @warnings;
    my  $ret = 1;

    for(my $q = 0; $q <= $#qrynums; $q++) {
        my $qnum = $qrynums[$q];
        ##print errors issued by threads if any
        if($$rinputs{"${qnum}_rcode"}) { ##don't send an e-mail
            ##$self->Error($$rinputs{"${qnum}_error"}, $$rinputs{"${qnum}_errhl"}, 1);
            next;
        }
        $command .= " \"".$$rinputs{"${qnum}_profl"}."\"";
    }

    unless($command) {
        $$rerrmsg = "ERROR: $self->{MYPROGNAME}: $mysubname: Making profiles for all queries failed.\n";
        $$rhlerrmsg = "Making profiles for all queries failed. Please check your input.\n";
        return 0;
    }

    ##put all profiles into one file
    $command = "cat ${command} >\"${fileofprofiles}\"";

    print(STDERR GetTime()."${preamb} ${command}\n\n");

    unless($self->ExecCommand($command)) {
        $$rerrmsg = "ERROR: $self->{MYPROGNAME}: $mysubname: Failed to combine profiles into one file\n";
        $$rhlerrmsg = "Combining profiles into one file failed.\n";
        return 0;
    }

    ##run comer
    my $gpumemopt = '';
    if($self->{cfgvar}->Exists('JOB_GPU_MEM')) {
        my $valuemb = $self->{cfgvar}->GetValue('JOB_GPU_MEM')*1024;
        $gpumemopt = "--dev-mem=${valuemb}" if $valuemb > 256;
    }
    $command = "$self->{optionvalues}->{prog_comer_comer} ".
            "-v -i \"${fileofprofiles}\" -d \"".$self->{optionvalues}->{$self->{optionkeys}->{comer_db}}.
            "\" -o \"${outputsubdir}\" -p \"$self->{OPTFILENAME}\" -f 1 ".
            "--dev-pass2memp=${comeropt_pass2memp} ${gpumemopt} >\"$self->{COMLOGFILENAME}\" 2>&1";

    for(my $i=1;; $i++) {
        print(STDERR GetTime()."${preamb} ${command}\n\n");
        unless($self->ExecCommand($command)) {
            do{ sleep(2); next} if $i < $ntrials;
            $$rerrmsg = "ERROR: $self->{MYPROGNAME}: $mysubname: Unable to conduct a COMER search.\n";
            $$rhlerrmsg = "Unable to conduct a COMER search. If you submitted profile(s) as query(-ies), please make sure they have correct format.\n";
            return 0;
        }
        last;
    }

    ##extract warnings if any; no check of return code
    $self->GetWarnings($self->{COMLOGFILENAME}, \@warnings);

    if(0 <= $#warnings) {
        my $lst = join("", @warnings);
        my $text = "\nWarnings from the COMER search:\n\n$lst\n\n";
        $self->Warning($text, $text, 1);##no e-mail
    }

    ##associate comer results with particular queries
    $command = '';
    for(my ($q,$c) = (0,0); $q <= $#qrynums; $q++) {
        my $qnum = $qrynums[$q];
        ##skip queries with pending errors and which have not been searched for
        next if($$rinputs{"${qnum}_rcode"});

        my $namepattern = "${outputsubdir}/".basename($$rinputs{"${qnum}_bname"})."*__${c}.$self->{OUTEXT}";
        my @outfiles = glob($namepattern);
        if($#outfiles < 0) {
            $$rinputs{"${qnum}_rcode"} = 1;
            $$rinputs{"${qnum}_error"} = "ERROR: $self->{MYPROGNAME}: $mysubname: COMER results for query No.${qnum} not found.\n";
            $$rinputs{"${qnum}_errhl"} = "COMER results for query No.${qnum} not found.\n";
            $self->Error($$rinputs{"${qnum}_error"}, $$rinputs{"${qnum}_errhl"}, 1);##no e-mail
        }
        else {
            $$rinputs{"${qnum}_outfl"} = $command = $outfiles[0];
        }
        $c++;
    }

    unless($command) {
        $$rerrmsg = "ERROR: $self->{MYPROGNAME}: $mysubname: No COMER results for all queries.\n";
        $$rhlerrmsg = "No COMER results for all queries.\n";
        return 0;
    }

    return $ret;
}

## -----------------------------------------------------------------------------
## MakeResultsList: make and write to a file the list of the results files over 
## all queries
##
sub MakeResultsList
{
    my  $self = shift;
    my  $rfilelist = shift;##ref to the list of files listed in the results list file
    my  $nqueries = shift;##number of queries
    my  $rinputs = shift;##ref to the hash of individual inputs
    my  $rerrmsg = shift;##ref to the error message string to be put in logs
    my  $rhlerrmsg = shift;##ref to the h-l error message string

    my  $mysubname = (caller(0))[3];
    my  $class = ref($self) || die("ERROR: $mysubname: Should be called by object.");
    my  $preamb = "[ ${mysubname} ] ";
    my  @qrynums = 0..$nqueries-1;
    my  $ret = 1;

    unless(open(F,'>',$self->{RESLSTFILENAME})){
        $$rerrmsg = "ERROR: $self->{MYPROGNAME}: $mysubname: Failed to open file for writing: '$self->{RESLSTFILENAME}'.\n";
        $$rhlerrmsg = "Opening file for writing failed.\n";
        return 0;
    }

    print(F "# Search_results Profile MSA Query neff_file logfile\n");

    for(my $q = 0; $q <= $#qrynums; $q++) {
        my $qnum = $qrynums[$q];
        ##skip queries with pending errors
        next if($$rinputs{"${qnum}_rcode"});

        my $profile = $$rinputs{"${qnum}_profl"};
        $profile = $$rinputs{"${qnum}_tprof"} if $$rinputs{"${qnum}_tprof"};

        my @files = ( $$rinputs{"${qnum}_outfl"},
                      $profile,
                      $$rinputs{"${qnum}_msafl"},
                      $$rinputs{"${qnum}_input"},
                      $$rinputs{"${qnum}_nefffl"},
                      $$rinputs{"${qnum}_logfl"}
        );

        do{$_ = '' unless -f $_} foreach @files;

        $_ =~ s/^$self->{inpdirname}\/*(.+)$/$1/ foreach @files;

        print(F "\"$files[0]\"\t\"$files[1]\"\t\"$files[2]\"\t\"$files[3]\"\t\"$files[4]\"\t\"$files[5]\"\n");

        push @$rfilelist, @files;
    }

    close(F);
    return $ret;
}

## -----------------------------------------------------------------------------
## CompressResults: compress all required results files to a single archive
##
sub CompressResults
{
    my  $self = shift;
    my  $rfilelist = shift;##ref to the list of ($dirname-prefix-removed) files to include in the archive
    my  $rerrmsg = shift;##ref to the error message string to be put in logs
    my  $rhlerrmsg = shift;##ref to the h-l error message string

    my  $mysubname = (caller(0))[3];
    my  $class = ref($self) || die("ERROR: $mysubname: Should be called by object.");
    my  $preamb = "[ ${mysubname} ] ";
    my  $command = '';
    my  $ret = 1;

    return 0 unless($self->AddFileToArchive($self->{COMLOGFILENAME}, '', $rerrmsg, $rhlerrmsg, 1)); #no dirname; 1==create

    return 0 unless($self->AddFileToArchive($self->{RESLSTFILENAME}, '', $rerrmsg, $rhlerrmsg)); #no dirname

    for(my $i = 0; $i <= $#{$rfilelist}; $i++) {
        next unless $$rfilelist[$i];
        return 0 unless($self->AddFileToArchive($$rfilelist[$i], $self->{inpdirname}, $rerrmsg, $rhlerrmsg)); #dirname given
    }

    unless($self->AddFileToArchive($self->{STAFILENAME}, '', $rerrmsg, $rhlerrmsg)){ #no dirname
        ##record an error, don't stop and  don't send an e-mail
        $self->Error($$rerrmsg, $$rhlerrmsg, 1);
    }

    unless($self->AddFileToArchive($self->{ERRFILENAME}, '', $rerrmsg, $rhlerrmsg)){ #no dirname
        ##record an error, don't stop and  don't send an e-mail
        $self->Error($$rerrmsg, $$rhlerrmsg, 1);
    }

    $command = "$self->{GZPROG} -f \"$self->{RESULTSFILENAME}\"";

    print(STDERR GetTime()."${preamb} ${command}\n\n");

    unless($self->ExecCommand($command)) {
        $$rerrmsg = "ERROR: $self->{MYPROGNAME}: $mysubname: Failed to gzip the archive: '$self->{RESULTSFILENAME}'\n";
        $$rhlerrmsg = "Failed to compress the archive of results files.\n";
        return 0;
    }

    return $ret;
}



## -----------------------------------------------------------------------------
## CreateLink: create a symbolic link to file
##
sub CreateLink
{
    my  $self = shift;
    my  $oldfile = shift;##source file
    my  $newfile = shift;##destination/new file (link)
    my  $rerrmsg = shift;##ref to the error message string to be put in logs
    my  $rhlerrmsg = shift;##ref to the h-l error message string

    my  $mysubname = (caller(0))[3];
    my  $class = ref($self) || die("ERROR: $mysubname: Should be called by object.");
    my  $preamb = "[ ${mysubname} ] ";
    my  $command = '';

    my $symlink_exists = eval {symlink("",""); 1};

    ##create a link to the msa file
    unless(-f $newfile) {
        ##if($symlink_exists) {
            ##unless(symlink($oldfile, $newfile)) {
            $command = "ln -s $oldfile $newfile";
            print(STDERR GetTime()."${preamb} ${command}\n\n");
            unless($self->ExecCommand($command)) {
                $$rerrmsg = "ERROR: $self->{MYPROGNAME}: $mysubname: Failed to create a link to a file: '$oldfile' -> '$newfile'\n";
                $$rhlerrmsg = "Creating a link failed.\n";
                return 0;
            }
        ##}
    }
    return 1;
}

## -----------------------------------------------------------------------------
## PredictDistancesAndMakeCOTHERProfiles: predict inter-residue distances using 
## ROPIUS0 and integrate them into the COTHER profile
##
sub PredictDistancesAndMakeCOTHERProfiles
{
    my  $self = shift;
    my  $nqueries = shift;##number of queries
    my  $rinputs = shift;##ref to the hash of individual inputs
    my  $rerrmsg = shift;##ref to the error message string to be put in logs
    my  $rhlerrmsg = shift;##ref to the h-l error message string

    my  $mysubname = (caller(0))[3];
    my  $class = ref($self) || die("ERROR: $mysubname: Should be called by object.");
    my  $preamb = "[ ${mysubname} ] ";
    my  @qrynums = 0..$nqueries-1;
    my  $command = '';
    my  @warnings;
    my  $ret = 1;

    for(my $q = 0; $q <= $#qrynums; $q++) {
        my $qnum = $qrynums[$q];
        ##print errors issued by threads if any
        next if $$rinputs{"${qnum}_rcode"};

        my $logfile = $$rinputs{"${qnum}_logfl"};
        my $redirecttext = (-f $logfile)? ">>\"${logfile}\" 2>&1": '';

        my $msafile = $$rinputs{"${qnum}_msafl"};
        my $profile = $$rinputs{"${qnum}_profl"};
        ##my $thisdir = dirname($profile);
        ##my $proname = basename($profile,$self->{PROEXT});
        my ($proname, $thisdir, $prosfx) = fileparse($profile, qr/\.[^.]*/);
        my $covfile = File::Spec->catfile($thisdir, "$proname.$self->{COVEXT}");
        my $msafilecopy = File::Spec->catfile($thisdir, $proname);

        ##profile is in COTHER format; no need for distance prediction and profile construction
        if($prosfx eq ".$self->{TPROEXT}") {
           $$rinputs{"${qnum}_tprof"} = $profile;
           next;
        }

        unless(-f $profile) {
            $$rinputs{"${qnum}_rcode"} = 1;
            $$rinputs{"${qnum}_error"} = "ERROR: $self->{MYPROGNAME}: $preamb Profile No.${qnum} not found: '${profile}'\n";
            $$rinputs{"${qnum}_errhl"} = "Profile No.${qnum} not found.\n";
            next;
        }

        unless(-f $covfile) {
            $$rinputs{"${qnum}_rcode"} = 1;
            $$rinputs{"${qnum}_error"} = "ERROR: $self->{MYPROGNAME}: $preamb Xcov file No.${qnum} not found: '${covfile}'\n";
            $$rinputs{"${qnum}_errhl"} = "Xcov file No.${qnum} not found (Query in COMER rather than COTHER format?).\n";
            next;
        }

        my $ropsubdir = File::Spec->catfile($thisdir,"ropius0_${qnum}");
        my $profilecopy = File::Spec->catfile($ropsubdir, $proname . $prosfx);
        my $covfilecopy = File::Spec->catfile($ropsubdir,"$proname.$self->{COVEXT}");

        unless( -d $ropsubdir || mkdir($ropsubdir)) {
            $$rerrmsg = "ERROR: $self->{MYPROGNAME}: $mysubname: Failed to create directory: '$ropsubdir'\n";
            $$rhlerrmsg = "Creating a directory failed.\n";
            $ret = 0;
            last;
        }

        ##create a link to the msa file
        unless($self->CreateLink($msafile, $msafilecopy, $rerrmsg, $rhlerrmsg)) {
            $ret = 0;
            last;
        }

        ##create a link to profile
        unless($self->CreateLink($profile, $profilecopy, $rerrmsg, $rhlerrmsg)) {
            $ret = 0;
            last;
        }
        ##create a link to xcov file
        unless($self->CreateLink($covfile, $covfilecopy, $rerrmsg, $rhlerrmsg)) {
            $ret = 0;
            last;
        }

        ##predict distances
        $command = $self->{optionvalues}->{prog_ropius0_distopred} .
            " -i \"${msafilecopy}\" -o \"${ropsubdir}\"    ${redirecttext}";

        if(open(F,'>>',$logfile)){print(F GetTime()."${preamb} ${command}\n\n");close(F);}

        unless($self->ExecCommand($command)) {
            $$rinputs{"${qnum}_rcode"} = 1;
            $$rinputs{"${qnum}_error"} = "ERROR: $self->{MYPROGNAME}: $preamb Failed to predict distances for query No.${qnum}: '$msafilecopy'\n";
            $$rinputs{"${qnum}_errhl"} = "Distance prediction for query No.${qnum} failed.\n";
            next;
        }

        ##incorporate predictions into a COTHER profile
        my $outtprofile = File::Spec->catfile($thisdir, "$proname.$self->{TPROEXT}");
        my $dstprbfile = File::Spec->catfile($ropsubdir,"${proname}__pred__nonavg.prb");

        unless(-f $dstprbfile) {
            $$rinputs{"${qnum}_rcode"} = 1;
            $$rinputs{"${qnum}_error"} = "ERROR: $self->{MYPROGNAME}: $preamb Distance predictions for query No.${qnum} not found: '$dstprbfile'\n";
            $$rinputs{"${qnum}_errhl"} = "Distance predictions for query No.${qnum} not found.\n";
            next;
        }

        $command = $self->{optionvalues}->{prog_cother_adddist} .
            " -v -i \"${dstprbfile}\" -j \"${profile}\" -o \"${outtprofile}\" -p \"$self->{OPTFILENAME}\" --dst=3,5,7 --prb=0.05    ${redirecttext}";

        if(open(F,'>>',$logfile)){print(F GetTime()."${preamb} ${command}\n\n");close(F);}

        unless($self->ExecCommand($command)) {
            $$rinputs{"${qnum}_rcode"} = 1;
            $$rinputs{"${qnum}_error"} = "ERROR: $self->{MYPROGNAME}: $preamb Failed to incorporate distance predictions for query No.${qnum}: '$dstprbfile'\n";
            $$rinputs{"${qnum}_errhl"} = "Incorporation of distance predictions for query No.${qnum} failed.\n";
            next;
        }

        $$rinputs{"${qnum}_tprof"} = $outtprofile;
    }

    return $ret unless $ret;

    return $ret;
}

## -----------------------------------------------------------------------------
## RunCOTHER: put all profiles together in one file and run COTHER search
## TODO: Make one general subroutine for RunCOMER and RunCOTHER
##
sub RunCOTHER
{
    my  $self = shift;
    my  $nqueries = shift;##number of queries
    my  $rinputs = shift;##ref to the hash of individual inputs
    my  $rerrmsg = shift;##ref to the error message string to be put in logs
    my  $rhlerrmsg = shift;##ref to the h-l error message string

    my  $mysubname = (caller(0))[3];
    my  $class = ref($self) || die("ERROR: $mysubname: Should be called by object.");
    my  $preamb = "[ ${mysubname} ] ";
    my  @qrynums = 0..$nqueries-1;
    my  $fileofprofiles = File::Spec->catfile($self->{inpdirname},"$self->{inpbasename}__nqries${nqueries}.$self->{TPROEXT}");
    my  $outputpfx = File::Spec->catfile($self->{inpdirname},"$self->{inpbasename}__cother_out");
    my  $outputsubdir = $outputpfx;
    my  $cotheropt_pass2memp = 100;
    ##no. trials to launch the program; NOTE: may fail when two or more processes request resources at the same time:
    my  $ntrials = 1;
    my  $command = '';
    my  @warnings;
    my  $ret = 1;

    for(my $q = 0; $q <= $#qrynums; $q++) {
        my $qnum = $qrynums[$q];
        ##print errors issued by threads if any
        if($$rinputs{"${qnum}_rcode"}) { ##don't send an e-mail
            ##$self->Error($$rinputs{"${qnum}_error"}, $$rinputs{"${qnum}_errhl"}, 1);
            next;
        }
        $command .= " \"".$$rinputs{"${qnum}_tprof"}."\"";
    }

    unless($command) {
        $$rerrmsg = "ERROR: $self->{MYPROGNAME}: $mysubname: Making profiles for all queries failed.\n";
        $$rhlerrmsg = "Making profiles for all queries failed. Please check your input.\n";
        return 0;
    }

    ##put all profiles into one file
    $command = "cat ${command} >\"${fileofprofiles}\"";

    print(STDERR GetTime()."${preamb} ${command}\n\n");

    unless($self->ExecCommand($command)) {
        $$rerrmsg = "ERROR: $self->{MYPROGNAME}: $mysubname: Failed to combine profiles into one file\n";
        $$rhlerrmsg = "Combining profiles into one file failed.\n";
        return 0;
    }

    ##run cother
    my $gpumemopt = '';
    if($self->{cfgvar}->Exists('JOB_GPU_MEM')) {
        my $valuemb = $self->{cfgvar}->GetValue('JOB_GPU_MEM')*1024;
        $gpumemopt = "--dev-mem=${valuemb}" if $valuemb > 256;
    }
    $command = "$self->{optionvalues}->{prog_cother_cother} ".
            "-v -i \"${fileofprofiles}\" -d \"".$self->{optionvalues}->{$self->{optionkeys}->{cother_db}}.
            "\" -o \"${outputsubdir}\" -p \"$self->{OPTFILENAME}\" -f 1 ".
            "--dev-pass2memp=${cotheropt_pass2memp} ${gpumemopt} >\"$self->{COMLOGFILENAME}\" 2>&1";

    for(my $i=1;; $i++) {
        print(STDERR GetTime()."${preamb} ${command}\n\n");
        unless($self->ExecCommand($command)) {
            do{ sleep(2); next} if $i < $ntrials;
            $$rerrmsg = "ERROR: $self->{MYPROGNAME}: $mysubname: Unable to conduct a COTHER search.\n";
            $$rhlerrmsg = "Unable to conduct a COTHER search. If you submitted profile(s) as query(-ies), please make sure they have correct format.\n";
            return 0;
        }
        last;
    }

    ##extract warnings if any; no check of return code
    $self->GetWarnings($self->{COMLOGFILENAME}, \@warnings);

    if(0 <= $#warnings) {
        my $lst = join("", @warnings);
        my $text = "\nWarnings from the COTHER search:\n\n$lst\n\n";
        $self->Warning($text, $text, 1);##no e-mail
    }

    ##associate cother results with particular queries
    $command = '';
    for(my ($q,$c) = (0,0); $q <= $#qrynums; $q++) {
        my $qnum = $qrynums[$q];
        ##skip queries with pending errors and which have not been searched for
        next if($$rinputs{"${qnum}_rcode"});

        my $namepattern = "${outputsubdir}/".basename($$rinputs{"${qnum}_bname"})."*__${c}.$self->{OUTEXT}";
        my @outfiles = glob($namepattern);
        if($#outfiles < 0) {
            $$rinputs{"${qnum}_rcode"} = 1;
            $$rinputs{"${qnum}_error"} = "ERROR: $self->{MYPROGNAME}: $mysubname: COTHER results for query No.${qnum} not found.\n";
            $$rinputs{"${qnum}_errhl"} = "COTHER results for query No.${qnum} not found.\n";
            $self->Error($$rinputs{"${qnum}_error"}, $$rinputs{"${qnum}_errhl"}, 1);##no e-mail
        }
        else {
            $$rinputs{"${qnum}_outfl"} = $command = $outfiles[0];
        }
        $c++;
    }

    unless($command) {
        $$rerrmsg = "ERROR: $self->{MYPROGNAME}: $mysubname: No COTHER results for all queries.\n";
        $$rhlerrmsg = "No COTHER results for all queries.\n";
        return 0;
    }

    return $ret;
}







## THREADS =====================================================================
##
## RunThreadsAndWait: launch threads for producing MSAs and making profiles, and 
## wait for them to finish
##
sub RunThreadsAndWait
{
    my  $self = shift;
    my  $subroutine_t = shift;##thread subroutine address
    my  $oneworker = shift;##flag instructing to use only one worker with spec. #cpus
    my  $nqueries = shift;##number of queries
    my  $rinputs = shift;##ref to the hash of individual inputs
    my  $rerrmsg = shift;##ref to the error message string to be put in logs
    my  $rhlerrmsg = shift;##ref to the h-l error message string

    my  $mysubname = (caller(0))[3];
    my  $class = ref($self) || die("ERROR: $mysubname: Should be called by object.");

    my  $ncpus = $self->{ncpus};
    my  $nworkers = 1;##number of threads
    my ($tmMSA, $tmPro, $procd) = (0,0,0);##timespans and fraction of processed queries
    ##query serial numbers sorted by query size:
    my  @qrynums = sort {$$rinputs{"${b}_isize"} <=> $$rinputs{"${a}_isize"}} 0..$nqueries-1;
    my  @workers;
    my  $ret = 1;

    if(!$oneworker && 
        $#qrynums+1 > $self->{MAXNQUERIES_multicore})
    {
        $nworkers = $ncpus;
        $ncpus = 1;
    }

    printf(STDERR "\nINFO: #queries= %d  dedicated #workers= %d #cpus= %d/wrk.\n\n",
        $#qrynums+1,$nworkers,$ncpus);

    my @joinable = ();
    my ($nsuccess, $nqries, $nvalidqries) = (0,$#qrynums+1,$#qrynums+1);

    for(my $q = 0; $q <= $#qrynums || threads->list(threads::all); ) {
        my $w = threads->list(threads::running);
        for( ; $w < $nworkers && $q <= $#qrynums; $w++, $q++) {
            my $qnum = $qrynums[$q];
            ##skip if errors from the previous run are pending
            next if $$rinputs{"${qnum}_rcode"};
            $workers[$qnum] = threads->create(
                    {'context' => 'list'},
                    sub{$self->$subroutine_t(@_)}, 
                    $self->{NOPRO}, $self->{NORUN}, $ncpus, $qnum, $self->{optionkeys}, $self->{optionvalues}, $rinputs
            );
            unless($workers[$qnum]) {
                $$rinputs{"${qnum}_rcode"} = 1;##code of fail
                $$rinputs{"${qnum}_error"} = "Creation of thread $qnum (input $qnum) failed.\n";##error message
                $$rinputs{"${qnum}_errhl"} = "Creation of thread $qnum (input $qnum) failed.\n";##high-level msg
            }
            ##print(STDERR "[test] --> w= $w q= $q\n");
        }

        @joinable = threads->list(threads::joinable);##...or use @workers

        foreach my $thr(@joinable) {
            my @thretlst = $thr->join();
            my $qnum = $thretlst[0];
            if($#qrynums < $qnum) {
                $$rerrmsg = "ERROR: $self->{MYPROGNAME}: 
                    Invalid query/input number returned by thread ".$thr->tid().".\n";
                $$rhlerrmsg = '';##not to be shown
                $ret = 0;
                next;
            }
            $$rinputs{"${qnum}_rcode"} = $thretlst[1];##code of fail
            $$rinputs{"${qnum}_error"} = $thretlst[2];##error message
            $$rinputs{"${qnum}_errhl"} = $thretlst[3];##high-level msg
            $$rinputs{"${qnum}_msafl"} = $thretlst[4] if $thretlst[4];##MSA filename
            $$rinputs{"${qnum}_profl"} = $thretlst[5] if $thretlst[5];##profile filename
            $$rinputs{"${qnum}_nefffl"} = $thretlst[6] if $thretlst[6];##.neff filename
            $tmMSA += ($$rinputs{"${qnum}_tmMSA"} = $thretlst[7]);##MSA building timespan
            $tmPro += ($$rinputs{"${qnum}_tmPro"} = $thretlst[8]);##profile construction timespan
            $$rinputs{"${qnum}_hhsrun"} = $thretlst[9] if $thretlst[9];##whether hhsuite has been run
            $$rinputs{"${qnum}_hmmrun"} = $thretlst[10] if $thretlst[10];##whether hmmer has been run

            ##next if $$rinputs{"${qnum}_rcode"};
            do{$nvalidqries--; next} if $$rinputs{"${qnum}_rcode"};

            $nsuccess++;

            next if $procd == int($nsuccess*10/$nvalidqries);

            $procd = int($nsuccess*10/$nvalidqries);

            my $tmdmsg = '';
            if($tmMSA && $tmPro) {
                $tmdmsg = sprintf("(time distr.: %.0f%% MSA, %.0f%% profile construction)",
                        $tmMSA*100./($tmMSA+$tmPro),$tmPro*100./($tmMSA+$tmPro));
            } elsif($tmMSA) {
                $tmdmsg = sprintf("(time distr.: %.0f%% MSA)", $tmMSA*100./($tmMSA+$tmPro));
            } elsif($tmPro) {
                $tmdmsg = sprintf("(time distr.: %.0f%% profile construction)", $tmPro*100./($tmMSA+$tmPro));
            }

            $self->ProgressMsg("  $nsuccess/$nqries queries done $tmdmsg\n");
        }

        sleep(2) if $#joinable < 0;
    }


    return $ret;
}

## -----------------------------------------------------------------------------
##
##
sub ProcessQuery_hhsuite_t
{
    my  $self = shift;
    my  $nopro = shift;##unused flag
    my  $norun = shift;##unsed flag
    my  $ncpus = shift;##number of cpus
    my  $qrynum = shift;##query serial number
    my  $roptionkeys_t = shift;##ref to the keys of options
    my  $roptionvalues_t = shift;## ref to the values of options
    my  $rinputs_t = shift;##ref to the hash of individual inputs

    my  $mysubname = (caller(0))[3];
    my  $class = ref($self) || die("ERROR: $mysubname: Should be called by object.");

    my  $preamb = "[ ${mysubname} ".threads->tid()." ] ";
    my  $rcode = 0;
    ## These are general return fields that should be returned even if not computed!
    ## error strings, output MSA and profile filenames, time for building the MSA, 
    ##   time for profile construction, and flags of hhsuite and hmmer having run:
    my ($error,$errhl,$msafile,$profile,$timeMSA,$timePro,$hhsrun,$hmmrun) = ('','','','',0,0,0,0);
    my  $nefffile = '';

    my  $combine = $$roptionvalues_t{$$roptionkeys_t{hhsuite_in_use}} &&
                   $$roptionvalues_t{$$roptionkeys_t{hmmer_in_use}};

    ## use these contants over all thread subprograms:
    my ($resfxhhs) = ('_resulthhs');
    my  $inputfile = $$rinputs_t{"${qrynum}_input"};
    my  $outputreformatfile = $$rinputs_t{"${qrynum}_bname"}."_rfm.$self->{AFAEXT}";
    my  $outputhhsfile = $$rinputs_t{"${qrynum}_bname"}.${resfxhhs};
    my  $logfile = $$rinputs_t{"${qrynum}_logfl"};
    my  $inputisprofile = ($inputfile =~ /\.(?:$self->{PROEXT}|$self->{TPROEXT})$/)? 1: 0;
    my  $inputisa3m = ($inputfile =~ /$self->{A3MEXT}$/)? 1: 0;
    my  $opta = $combine? '': '-a';
    my  $command = '';

    $hhsrun = $$rinputs_t{"${qrynum}_hhsrun"};
    $hmmrun = $$rinputs_t{"${qrynum}_hmmrun"};

    if($$roptionvalues_t{$$roptionkeys_t{hhsuite_in_use}} && !$inputisprofile) {
        $timeMSA = time();

        my $hhinputfile = $inputfile;

        if($inputisa3m) {
            my $tmpfile = "${outputreformatfile}.f";
            $command = "$$roptionvalues_t{prog_hhsuite_reformat} a3m fas \"${inputfile}\" \"${tmpfile}\" && ".
                "(perl -e 'while(<>){\$d=1 if /^>/; \$d=0 if /^>ss_/;print if \$d}' ${tmpfile} >${outputreformatfile})";
            ##$command = "$$roptionvalues_t{prog_hhsuite_reformat} a3m fas \"${inputfile}\" \"${outputreformatfile}\"";

            unless(-f $outputreformatfile) {
                if(open(F,'>>',$logfile)){print(F GetTime()."${preamb} ${command}\n\n");close(F);}

                unless($self->ExecCommand($command)) {
                    $error = "ERROR: $self->{MYPROGNAME}: $preamb Failed to reformat MSA No.${qrynum}: '${inputfile}'\n";
                    $errhl = "Reformatting the MSA (query No.${qrynum}) failed. Please check your input.\n";
                    return ($qrynum, 1, $error, $errhl, $msafile, $profile, $nefffile, time()-$timeMSA, $timePro, $hhsrun, $hmmrun);
                }
            }

            $hhinputfile = $outputreformatfile;
        }

        $command = "$self->{HHBLITShelper} -i \"${hhinputfile}\" -o \"${outputhhsfile}\" -d \"".
                $$roptionvalues_t{$$roptionkeys_t{hhsuite_db}}."\" -n ".
                $$roptionvalues_t{$$roptionkeys_t{hhsuite_opt_niterations}}." -e ".
                $$roptionvalues_t{$$roptionkeys_t{hhsuite_opt_evalue}}." -p ${ncpus} ${opta} -b \"".
                $self->{cfgvar}->InstallDir_HHsuite()."\" -N $self->{MAXNSEQS_perprog} >>\"${logfile}\" 2>&1";

        my $condfile = $combine? "${outputhhsfile}.$self->{PWFEXT}": "${outputhhsfile}.$self->{AFAEXT}";

        unless($hhsrun || -f $condfile) {
            if(open(F,'>>',$logfile)){print(F GetTime()."${preamb} ${command}\n\n");close(F);}

            unless($self->ExecCommand($command)) {
                $error = "ERROR: $self->{MYPROGNAME}: $preamb Failed to make MSA No.${qrynum} by hhsuite for '${inputfile}'\n";
                $errhl = "Building an MSA by searching for query No.${qrynum} failed. Please check your input.\n";
                return ($qrynum, 1, $error, $errhl, $msafile, $profile, $nefffile, time()-$timeMSA, $timePro, $hhsrun, $hmmrun);
            }
            $hhsrun = 1;
        }

        $msafile = $condfile;
        $timeMSA = time()-$timeMSA;
    }

    return ($qrynum, $rcode, $error, $errhl, $msafile, $profile, $nefffile, $timeMSA, $timePro, $hhsrun, $hmmrun);
}

## -----------------------------------------------------------------------------
## TODO: Modularize!
##
sub ProcessQuery_t
{
    my  $self = shift;
    my  $nopro = shift;##flag of no construction of COMER profiles
    my  $norun = shift;##flag of no COMER run
    my  $ncpus = shift;##number of cpus
    my  $qrynum = shift;##query serial number
    my  $roptionkeys_t = shift;##ref to the keys of options
    my  $roptionvalues_t = shift;## ref to the values of options
    my  $rinputs_t = shift;##ref to the hash of individual inputs

    my  $mysubname = (caller(0))[3];
    my  $class = ref($self) || die("ERROR: $mysubname: Should be called by object.");

    my  $preamb = "[ ${mysubname} ".threads->tid()." ] ";
    my  $rcode = 0;
    ## These are general return fields that should be returned even if not computed!
    ## error strings, output MSA and profile filenames, time for building the MSA, 
    ##   time for profile construction, and flags of hhsuite and hmmer having run:
    my ($error,$errhl,$msafile,$profile,$timeMSA,$timePro,$hhsrun,$hmmrun) = ('','','','',0,0,0,0);
    my  $nefffile = '';

    my  $search = $$roptionvalues_t{$$roptionkeys_t{hhsuite_in_use}} ||
                  $$roptionvalues_t{$$roptionkeys_t{hmmer_in_use}};
    my  $combine = $$roptionvalues_t{$$roptionkeys_t{hhsuite_in_use}} &&
                   $$roptionvalues_t{$$roptionkeys_t{hmmer_in_use}};

    my ($resfxhhs,$resfxhmm,$resfx) = ('_resulthhs','_resulthmmer','_result');
    my  $inputfile = $$rinputs_t{"${qrynum}_input"};
    my  $outputreformatfile = $$rinputs_t{"${qrynum}_bname"}."_rfm.$self->{AFAEXT}";
    my  $rfminputfile = $inputfile;
    my  $outputhhsfile = $$rinputs_t{"${qrynum}_bname"}.${resfxhhs};
    my  $outputhmmfile = $$rinputs_t{"${qrynum}_bname"}.${resfxhmm};
    my  $reshhspwfafile = "${outputhhsfile}.$self->{PWFEXT}";
    my  $reshmmpwfafile = "${outputhmmfile}.$self->{PWFEXT}";
    my  $logfile = $$rinputs_t{"${qrynum}_logfl"};
    my  $opta = $combine? '': '-a';
    my  $respwfafile = $$rinputs_t{"${qrynum}_bname"}."${resfx}.$self->{PWFEXT}";
    my  $resfilepat = $$rinputs_t{"${qrynum}_bname"}.${resfx};##pattern of result files
    my  $resafafile = "${resfilepat}.$self->{AFAEXT}";##resulting MSA file
    my  $resprofile = "${resfilepat}.$self->{PROEXT}";##profile name
    my  $rescovfile = "${resfilepat}.$self->{COVEXT}";##xcov file name
    my  $resqryfile = '';##query filename
    my  $inputiscomer2profile = ($inputfile =~ /\.$self->{PROEXT}$/)? 1: 0;
    my  $inputiscotherprofile = ($inputfile =~ /\.$self->{TPROEXT}$/)? 1: 0;
    my  $inputisprofile = ($inputiscomer2profile | $inputiscotherprofile);
    ##name of profile as given by the query (comer2 or cother):
    my  $resorgprofile = "${resfilepat}.".($inputiscotherprofile? $self->{TPROEXT}: $self->{PROEXT});
    my  $inputisa3m = ($inputfile =~ /$self->{A3MEXT}$/)? 1: 0;
    my ($command, $cmbcmd) = ('','');

    $hhsrun = $$rinputs_t{"${qrynum}_hhsrun"};
    $hmmrun = $$rinputs_t{"${qrynum}_hmmrun"};


    if($inputisa3m) {
        my $tmpfile = "${outputreformatfile}.f";
        $command = "$$roptionvalues_t{prog_hhsuite_reformat} a3m fas \"${inputfile}\" \"${tmpfile}\" &&".
            "(perl -e 'while(<>){\$d=1 if /^>/; \$d=0 if /^>ss_/;print if \$d}' ${tmpfile} >${outputreformatfile})";

        unless(-f $outputreformatfile) {
            if(open(F,'>>',$logfile)){print(F GetTime()."${preamb} ${command}\n\n");close(F);}

            unless($self->ExecCommand($command)) {
                $error = "ERROR: $self->{MYPROGNAME}: $preamb Failed to reformat MSA No.${qrynum}: '${inputfile}'\n";
                $errhl = "Reformatting the MSA (query No.${qrynum}) failed. Please check your input.\n";
                return ($qrynum, 1, $error, $errhl, $msafile, $profile, $nefffile, time()-$timeMSA, $timePro, $hhsrun, $hmmrun);
            }
        }

        $rfminputfile = $outputreformatfile;
    }


    if($search && !$inputisprofile)
    {
        $timeMSA = time();

        if($$roptionvalues_t{$$roptionkeys_t{hhsuite_in_use}}) {
            $command = "$self->{HHBLITShelper} -i \"${rfminputfile}\" -o \"${outputhhsfile}\" -d \"".
                    $$roptionvalues_t{$$roptionkeys_t{hhsuite_db}}."\" -n ".
                    $$roptionvalues_t{$$roptionkeys_t{hhsuite_opt_niterations}}." -e ".
                    $$roptionvalues_t{$$roptionkeys_t{hhsuite_opt_evalue}}." -p ${ncpus} ${opta} -b \"".
                    $self->{cfgvar}->InstallDir_HHsuite()."\" -N $self->{MAXNSEQS_perprog} >>\"${logfile}\" 2>&1";

            my $condfile = $combine? "${outputhhsfile}.$self->{PWFEXT}": "${outputhhsfile}.$self->{AFAEXT}";

            unless($hhsrun || -f $condfile) {
                if(open(F,'>>',$logfile)){print(F GetTime()."${preamb} ${command}\n\n");close(F);}

                unless($self->ExecCommand($command)) {
                    $error = "ERROR: $self->{MYPROGNAME}: $preamb Failed to make MSA No.${qrynum} by hhsuite for '${inputfile}'\n";
                    $errhl = "Building an MSA by searching for query No.${qrynum} failed. Please check your input.\n";
                    return ($qrynum, 1, $error, $errhl, $msafile, $profile, $nefffile, time()-$timeMSA, $timePro, $hhsrun, $hmmrun);
                }
                $hhsrun = 1;
            }

            $resqryfile = "${outputhhsfile}.$self->{QRYEXT}";
            $resfilepat =  ${outputhhsfile};
            $resafafile = "${outputhhsfile}.$self->{AFAEXT}";
            $cmbcmd .= " \"${outputhhsfile}.$self->{PWFEXT}\"";
        }

        if($$roptionvalues_t{$$roptionkeys_t{hmmer_in_use}}) {
            $command = "$self->{HMMERhelper} -i \"${rfminputfile}\" -o \"${outputhmmfile}\" -d \"".
                    $$roptionvalues_t{$$roptionkeys_t{sequence_db}}."\" -n ".
                    $$roptionvalues_t{$$roptionkeys_t{hmmer_opt_niterations}}." -e ".
                    $$roptionvalues_t{$$roptionkeys_t{hmmer_opt_evalue}}." -p ${ncpus} ${opta} -b \"".
                    $self->{cfgvar}->InstallDir_HMMER()."\" -N $self->{MAXNSEQS_perprog} >>\"${logfile}\" 2>&1";

            my $condfile = $combine? "${outputhmmfile}.$self->{PWFEXT}": "${outputhmmfile}.$self->{AFAEXT}";

            unless($hmmrun || -f $condfile) {
                if(open(F,'>>',$logfile)){print(F GetTime()."${preamb} ${command}\n\n");close(F);}

                unless($self->ExecCommand($command)) {
                    $error = "ERROR: $self->{MYPROGNAME}: $preamb Failed to make MSA No.${qrynum} by HMMER for '${inputfile}'\n";
                    $errhl = "Building an MSA by searching for query No.${qrynum} failed. Please check your input.\n";
                    return ($qrynum, 1, $error, $errhl, $msafile, $profile, $nefffile, time()-$timeMSA, $timePro, $hhsrun, $hmmrun);
                }
                $hmmrun = 1;
            }

            $resqryfile = "${outputhmmfile}.$self->{QRYEXT}";
            $resfilepat =  ${outputhmmfile};
            $resafafile = "${outputhmmfile}.$self->{AFAEXT}";
            $cmbcmd .= " \"${outputhmmfile}.$self->{PWFEXT}\"";
        }

        if($combine) {
            $resfilepat = $$rinputs_t{"${qrynum}_bname"}.${resfx};
            $resafafile = "${resfilepat}.$self->{AFAEXT}";
        }

        $msafile = $resafafile;

        if($combine && !-f "${resfilepat}.$self->{AFAEXT}") {
            $resfilepat = $$rinputs_t{"${qrynum}_bname"}.${resfx};
            $resafafile = "${resfilepat}.$self->{AFAEXT}";

            my $evalue = $$roptionvalues_t{$$roptionkeys_t{hhsuite_opt_evalue}};
            $evalue = $$roptionvalues_t{$$roptionkeys_t{hmmer_opt_evalue}}
                    if($evalue < $$roptionvalues_t{$$roptionkeys_t{hmmer_opt_evalue}});

            $command = "cat ${cmbcmd} >\"${respwfafile}\" 2>>\"${logfile}\"";

            if(open(F,'>>',$logfile)){print(F GetTime()."${preamb} ${command}\n\n");close(F);}

            unless($self->ExecCommand($command)) {
                $error = "ERROR: $self->{MYPROGNAME}: $preamb Failed to combine MSAs for query No.${qrynum}: '${inputfile}'\n";
                $errhl = "Combining MSAs obtained for query No.${qrynum} failed.\n";
                return ($qrynum, 1, $error, $errhl, $msafile, $profile, $nefffile, time()-$timeMSA, $timePro, $hhsrun, $hmmrun);
            }

            ##no matter which $resqryfile (by hhs or hmmer)
            $command = "$self->{PWFA2MSAprog} -i \"${respwfafile}\" -o \"${resafafile}\" -f 0 ".
                    "-q \"${resqryfile}\" -e ${evalue} ".
                    "-N ".($self->{MAXNSEQS_perprog}*2)." >>\"${logfile}\" 2>&1";##x2 to account for two search programs
 
            if(open(F,'>>',$logfile)){print(F GetTime()."${preamb} ${command}\n\n");close(F);}

            unless($self->ExecCommand($command)) {
                $error = "ERROR: $self->{MYPROGNAME}: $preamb Failed to convert combined MSAs for query No.${qrynum}: '${inputfile}'\n";
                $errhl = "Combining MSAs obtained for query No.${qrynum} failed.\n";
                return ($qrynum, 1, $error, $errhl, $msafile, $profile, $nefffile, time()-$timeMSA, $timePro, $hhsrun, $hmmrun);
            }
        }

        $timeMSA = time()-$timeMSA;
    }


    if($inputisprofile) {
        $profile = $resorgprofile;
        unless($self->ChangeProfileFileField($inputfile, $profile, basename($inputfile))) {
            $error = "ERROR: $self->{MYPROGNAME}: $preamb Failed to preprocess profile No.${qrynum}: '${inputfile}'\n";
            $errhl = "Profile preprocessing for query No.${qrynum} failed. Please check your input.\n";
            return ($qrynum, 1, $error, $errhl, $msafile, $profile, $nefffile, $timeMSA, $timePro, $hhsrun, $hmmrun);
        }
    }
    elsif(!$nopro)
    {
        $msafile = $search? $resafafile: $rfminputfile;

        $nefffile = "$msafile.".$self->{NEFFEXT};

        unless(-f $nefffile)
        {
            $command = "$$roptionvalues_t{prog_comer_neff} -v -i \"${msafile}\" >\"${nefffile}\" 2>>\"${logfile}\"";
            if(open(F,'>>',$logfile)){print(F GetTime()."${preamb} ${command}\n\n");close(F);}
            unless($self->ExecCommand($command)) {
                $error = "ERROR: $self->{MYPROGNAME}: $preamb Failed to calculate Neff for query No.${qrynum} MSA '${msafile}'\n";
                $errhl = "WARNING: Calculating Neff for query No.${qrynum} failed.\n";
                $nefffile = '';
            }
        }

        $profile = $resprofile;

        $timePro = time();

        if($self->{USINGSSSCORING}) {
            $command = "$$roptionvalues_t{prog_comer_makepro} ".
                "-v -i \"${msafile}\" -o \"${profile}\" -p \"$self->{OPTFILENAME}\" -P \"".
                $self->{cfgvar}->InstallDir_PSIPRED()."\" -B \"".
                $self->{cfgvar}->InstallDir_BLAST()."\" >>\"${logfile}\" 2>&1";
        } else {
            $command = "$$roptionvalues_t{prog_comer_makepro} ".
                "-v -i \"${msafile}\" -o \"${profile}\" -p \"$self->{OPTFILENAME}\" >>\"${logfile}\" 2>&1";
        }

        unless(-f $profile)
        {
            unless(-f $msafile) {
                $error = "ERROR: $self->{MYPROGNAME}: $preamb MSA file for profile construction not found '${msafile}'\n";
                $errhl = "An MSA file for profile construction not found.\n";
                return ($qrynum, 1, $error, $errhl, $msafile, $profile, $nefffile, $timeMSA, time()-$timePro, $hhsrun, $hmmrun);
            }

            if(open(F,'>>',$logfile)){print(F GetTime()."${preamb} ${command}\n\n");close(F);}

            unless($self->ExecCommand($command)) {
                $error = "ERROR: $self->{MYPROGNAME}: $preamb Failed to construct profile No.${qrynum} by makepro for '${msafile}'\n";
                $errhl = "Constructing a profile for query No.${qrynum} failed. Please check your input. ".
                         "You might consider turning off secondary structure scoring (by setting \"Weight of SS scores\" to 0) or ".
                         "low-complexity filtering (by unchecking the \"Invoke low-complexity filtering for each ".
                           "sequence in alignment\" checkbox).\n";
                return ($qrynum, 1, $error, $errhl, $msafile, $profile, $nefffile, $timeMSA, time()-$timePro, $hhsrun, $hmmrun);
            }
        }

        if( $self->{METHOD} eq 'cother' && !-f $rescovfile) {
            my $prolen = 0;
            unless($self->GetProfileLength($profile, \$prolen)) {
            }
            if($self->{MAXSEQLENCOTHER} < $prolen) {
                $error = "ERROR: $self->{MYPROGNAME}: $preamb Query No.${qrynum} length exceeds max allowed ($self->{MAXSEQLENCOTHER}).\n";
                $errhl = "Length of query No.${qrynum} exceeds max allowed: $prolen > $self->{MAXSEQLENCOTHER}. Please check your input.\n";
                return ($qrynum, 1, $error, $errhl, $msafile, $profile, $nefffile, $timeMSA, time()-$timePro, $hhsrun, $hmmrun);
            }

            $command = "$$roptionvalues_t{prog_comer_makecov} ".
                    "-v -i \"${msafile}\" -o \"${rescovfile}\" -p \"$self->{OPTFILENAME}\" --scale --mi ".
                    " >>\"${logfile}\" 2>&1";

            ##NOTE: assume that the existence of msafile has been checked above

            if(open(F,'>>',$logfile)){print(F GetTime()."${preamb} ${command}\n\n");close(F);}

            unless($self->ExecCommand($command)) {
                $error = "ERROR: $self->{MYPROGNAME}: $preamb Failed to make xcov file No.${qrynum} by makecov for '${msafile}'\n";
                $errhl = "Constructing an xcov file for query No.${qrynum} failed. Please check your input.\n";
                return ($qrynum, 1, $error, $errhl, $msafile, $profile, $nefffile, $timeMSA, time()-$timePro, $hhsrun, $hmmrun);
            }
        }

        $timePro = time()-$timePro;
    }

    return ($qrynum, $rcode, $error, $errhl, $msafile, $profile, $nefffile, $timeMSA, $timePro, $hhsrun, $hmmrun);
}

## General =====================================================================
##
##







## =============================================================================
## MyExit: print a code to file and exit
##
sub MyExit
{
    my $self = shift;
    my $ecode = shift;##exit code
    my $mysubname = (caller(0))[3];
    my $class = ref($self) || die("ERROR: $mysubname: Should be called by object.");
    if(open(F,">>",$self->{ERRFILENAME})){print(F "\n${ecode}\n");close(F);}
    exit($ecode);
}

## -------------------------------------------------------------------
## Warning: output a warning message and optionally, send an email
##
sub Warning
{
    my ($self, $msg, $hlerrmsg, $nosend) = @_;
    my $mysubname = (caller(0))[3];
    my $class = ref($self) || die("ERROR: $mysubname: Should be called by object.");
    return $self->Issue($msg, $hlerrmsg, $nosend,
        '',
        "Warning from comer-ws ($self->{MYPROGNAME})");
}

## Error: output an error message and optionally, send an email
##
sub Error
{
    my ($self, $msg, $hlerrmsg, $nosend) = @_;
    my $mysubname = (caller(0))[3];
    my $class = ref($self) || die("ERROR: $mysubname: Should be called by object.");
    return $self->Issue($msg, $hlerrmsg, $nosend,
        "ERROR: Issue on the server's backend side: ",
        "Error message from comer-ws ($self->{MYPROGNAME})");
}

## Issue: output an issue message and send an email
##
sub Issue
{
    my ($self, $msg, $hlerrmsg, $nosend,  $leadtxt1, $sbjct) = @_;
    my $mysubname = (caller(0))[3];
    my $class = ref($self) || die("ERROR: $mysubname: Should be called by object.");
    print(STDERR GetTime().$msg);
    if($self->{ERRFILENAME} && $hlerrmsg) {
        if(open(EF,">>",$self->{ERRFILENAME})) {
            print(EF $leadtxt1.$hlerrmsg);
            close(EF);
        } else {
            print(STDERR "ERROR: $self->{MYPROGNAME}: Write to error file skipped ".
                "due to the fail to open the file: '$self->{ERRFILENAME}'\n");
        }
    }
    return 1 if $nosend;
    unless(-f $self->{SENDMAILprog}) {
        print(STDERR "ERROR: $self->{MYPROGNAME}: Send-mail program not found: $self->{SENDMAILprog}\n");
        return 0;
    }
    my  $command = "$self->{SENDMAILprog} --sub \"$sbjct\" --body \"$msg\" 2>&1";
    print(STDERR "$command\n");
    return $self->ExecCommand($command);
}

## -------------------------------------------------------------------
## try read lock files in a directory and return 1 if locks exist,
## 0 otherwise, -1 on error
##

sub ReadLocks
{
    my  $self = shift;
    my  $mysubname = (caller(0))[3];
    my  $class = ref($self) || die("ERROR: $mysubname: Should be called by object.");
    return $self->ReadFiles(shift, shift);
}

sub ReadFiles
{
    my  $self = shift;
    my  $dirname = shift;
    my  $pattern = shift;  ## pattern of files to look for as locks
    my  $reffiles = shift; ## reference to vector of files
    my  $mysubname = (caller(0))[3];
    my  $class = ref($self) || die("ERROR: $mysubname: Should be called by object.");
    my  @files;
    my  $locref = defined($reffiles)? $reffiles: \@files;

    unless( opendir(DIR, $dirname)) {
        printf( STDERR "ERROR: Cannot open directory $dirname.\n" );
        return -1;
    }
    @{$locref} = grep { /$pattern$/ && -f File::Spec->catfile($dirname,$_) } readdir(DIR);
    closedir(DIR);
    return 1 if 0 <= $#{$locref};
    return 0;
}

## -------------------------------------------------------------------
## ExecCommand: execute system command
##
sub CheckStatus
{
    my  $self = shift;
    my  $mysubname = (caller(0))[3];
    my  $class = ref($self) || die("ERROR: $mysubname: Should be called by object.");
    $self->ExecCommand();
}

sub ExecCommand
{
    my  $self = shift;
    my  $cmdline = shift;
    my  $returnstatus = shift;##if defined, return exit status of the command
    my  $mysubname = (caller(0))[3];
    my  $class = ref($self) || die("ERROR: $mysubname: Should be called by object.");

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
## -----------------------------------------------------------------------------

1;


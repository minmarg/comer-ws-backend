#!/usr/bin/env perl
BEGIN {$^W=1}

## (C) 2020-2021 Mindaugas Margelevicius, Institute of Biotechnology, Vilnius University
## Helper script for modeling the 3D structure of the query using pairwise 
## alignments given for it

use strict;
use Config;
use threads;
use threads::shared;
use FindBin;
use lib "$FindBin::Bin";
use readcfg;
use readopt;
use File::Spec;
use File::Basename;
use Getopt::Long;
use POSIX qw(strftime);

##constants:
my  $MYPROGNAME = basename($0);
my  $RENUMhelper = File::Spec->catfile($FindBin::Bin,"renumchain.py");
my  $DWNPDBSTRhelper = File::Spec->catfile($FindBin::Bin,"dwnpdbchains.pl");
my  $DWNSCPSTRhelper = File::Spec->catfile($FindBin::Bin,"dwnscopedomain.pl");
my  $PWFA2MSAprog = File::Spec->catfile($FindBin::Bin,"pwfa2msa.pl");
my  $SENDMAILprog = File::Spec->catfile($FindBin::Bin,"sendemail.pl");
my  $CFGFILE = File::Spec->catfile($FindBin::Bin,File::Spec->updir(),"var","comer-ws-backend.conf");
my  $devnull = File::Spec->devnull();
my  $TARPROG = `which tar 2>$devnull`; chomp($TARPROG);
my  $GZPROG = `which gzip 2>$devnull`; chomp($GZPROG);
my  $MAXNCPUs = 6;##maximum number of CPU cores assigned for a job
##threshold for the number of queries which, if exceeded, implies creating threads, each assigned 1 cpu:
##(use a large number for serialization with multiple cores)
my  $MAXNQUERIES_multicore = 0;
##maximum number of sequences permitted to be included in the results when each pair represents a model
my  $MAXNSEQS_perprog = 100;
##maximum number of sequences permitted to be included in the results when producing a single model using all alignments
my  $MAXNSEQS_multim = 7;

my ($FASEXT, $AFAEXT, $PWFEXT, $PIREXT, $ENTEXT, $PDBEXT) = 
   ('fa','afa','pwfa','pir','ent','pdb');

my  $usage = <<EOIN;

Helper script for modeling the 3D structure of the query using pairwise 
alignments given for it.
(C) 2021 Mindaugas Margelevicius, Institute of Biotechnology, Vilnius University

Usage:
$0 <Options>

Options:

--in <input>     Input file of pairwise alignments. Format is provided below.

--opt <options>  Filename of options for this job.

--status <file>  Name of the output status file of computation progress
                 messages.

--log <logfile>  Name of the output log file generated during the process.

--reslst <file>  Name of the output file listing the results files.

--results <file> Name of the output compressed file containing the results
                 files. (Extension .gz will be added.)

--err <errfile>  Name of the file of high-level error messages to write to on
                 error.

--help           This help text.


Format of the input file of pairwise alignments:

><query_description> (ALN:<query_start>-<query_end>)
<aligned_query_sequence>
><db_sequence_description> (ALN:<dbseq_start>-<dbseq_end>)
<aligned_database_sequence>
//
...

EOIN


my  $INPFILENAME = '';
my  $OPTFILENAME = '';
my  $STAFILENAME = '';
my  $LOGFILENAME = '';
my  $RESLSTFILENAME = '';
my  $RESULTSFILENAME = '';
my  $ERRFILENAME = '';
my  $Fail = 0;

my  $result = GetOptions(
               'in=s'       => \$INPFILENAME,
               'opt=s'      => \$OPTFILENAME,
               'status=s'   => \$STAFILENAME,
               'log=s'      => \$LOGFILENAME,
               'reslst=s'   => \$RESLSTFILENAME,
               'results=s'  => \$RESULTSFILENAME,
               'err=s'      => \$ERRFILENAME,
               'help|h'     => sub {print $usage; exit(0);}
);

print(STDERR "\n\n".GetDatetime()."\n\n");

my ($inpbasename,$inpdirname,$suffix) = fileparse($INPFILENAME, qr/\.[^.]*/);
$inpdirname = File::Spec->rel2abs($inpdirname);

unless($STAFILENAME) {
    $STAFILENAME = File::Spec->catfile(${inpdirname},"${inpbasename}.status");
    print(STDERR "WARNING: $MYPROGNAME: Filename for computation progress mesages not given; ".
          "it has been set to: '$STAFILENAME'\n");
}
unless($LOGFILENAME) {
    $LOGFILENAME = File::Spec->catfile(${inpdirname},"${inpbasename}__3d_out.log");
    print(STDERR "WARNING: $MYPROGNAME: Filename for COMER log mesages not given; ".
          "it has been set to: '$LOGFILENAME'\n");
}
unless($RESLSTFILENAME) {
    $RESLSTFILENAME = File::Spec->catfile(${inpdirname},"${inpbasename}__3d_out.lst");
    print(STDERR "WARNING: $MYPROGNAME: Filename for a results list not given; ".
          "it has been set to: '$RESLSTFILENAME'\n");
}
unless($RESULTSFILENAME) {
    $RESULTSFILENAME = File::Spec->catfile(${inpdirname},"${inpbasename}__3d_out.tar");
    print(STDERR "WARNING: $MYPROGNAME: Filename for compressed results not given; ".
          "it has been set to: '$RESULTSFILENAME'\n");
}
unless($ERRFILENAME) {
    $ERRFILENAME = File::Spec->catfile(${inpdirname},"${inpbasename}.err");
    print(STDERR "WARNING: $MYPROGNAME: Filename for high-level error mesages not given; ".
          "it has been set to: '$ERRFILENAME'\n");
}
## truncate error file before using Error!
if(open(F,'>',$ERRFILENAME)){close(F);}
unless($result) {
    Error("ERROR: $MYPROGNAME: Error in command-line arguments.\n",
            "Command-line arguments error.\n");## h-l error message
    MyExit(1);
}
unless($INPFILENAME && -f $INPFILENAME) {
    Error("ERROR: $MYPROGNAME: Input file not found: '$INPFILENAME'\n",
            "Input file not found.\n");## h-l error message
    MyExit(1);
}
unless($OPTFILENAME && -f $OPTFILENAME) {
    Error("ERROR: $MYPROGNAME: Input job options file not found: '$OPTFILENAME'\n",
            "Job options file not found.\n");## h-l error message
    MyExit(1);
}
unless(-f $CFGFILE) {
    Error("ERROR: $MYPROGNAME: Config file not found: '$CFGFILE'\n",
            "Configuration file not found.\n");## h-l error message
    MyExit(1);
}
$RENUMhelper = '' unless(-f $RENUMhelper);
unless(-f $DWNPDBSTRhelper) {
    Error("ERROR: $MYPROGNAME: Program file not found: '$DWNPDBSTRhelper'\n",
            "Some of the required program files not found.\n");## h-l error message
    MyExit(1);
}
unless(-f $DWNSCPSTRhelper) {
    Error("ERROR: $MYPROGNAME: Program file not found: '$DWNSCPSTRhelper'\n",
            "Some of the required program files not found.\n");## h-l error message
    MyExit(1);
}
unless(-f $PWFA2MSAprog) {
    Error("ERROR: $MYPROGNAME: Program file not found: '$PWFA2MSAprog'\n",
            "Some of the required program files not found.\n");## h-l error message
    MyExit(1);
}
unless($TARPROG && $GZPROG) {
    Error("ERROR: $MYPROGNAME: System programs 'tar' and/or 'gzip' not found.\n",
            "Some of the system programs not found.\n");## h-l error message
    MyExit(1);
}


## =============================================================================

my  $cfgvar = readcfg->new($CFGFILE);
my  %optionvalues;

if($cfgvar->JobMaxNo3DModels() < 1) {
    print(STDERR "WARNING: $MYPROGNAME: Maximum number of 3D models not given; ".
          "it has been set to: '$MAXNSEQS_perprog'\n");
} else {
    $MAXNSEQS_perprog = $cfgvar->JobMaxNo3DModels();
}
if($cfgvar->JobMaxNo3DTemplates() < 1) {
    print(STDERR "WARNING: $MYPROGNAME: Maximum number of 3D templates not given; ".
          "it has been set to: '$MAXNSEQS_multim'\n");
} else {
    $MAXNSEQS_multim = $cfgvar->JobMaxNo3DTemplates();
}


unless(-d $cfgvar->Path3dDb_PDB()) {
    Error("ERROR: $MYPROGNAME: Database directory not found: '".$cfgvar->Path3dDb_PDB()."'\n",
            "Some of the database directories not found.\n");## h-l error message
    MyExit(1);
}
unless(-d $cfgvar->Path3dDb_SCOPe()) {
    Error("ERROR: $MYPROGNAME: Database directory not found: '".$cfgvar->Path3dDb_SCOPe()."'\n",
            "Some of the database directories not found.\n");## h-l error message
    MyExit(1);
}

unless($cfgvar->WebAddr3dDb_SCOPe()=~/^http/) {
    Error("ERROR: $MYPROGNAME: Invalid web address for accessing SCOPe entries: '".$cfgvar->WebAddr3dDb_SCOPe()."'\n",
            "Some of the web addresses specified incorrectly.\n");## h-l error message
    MyExit(1);
}


unless(-f $cfgvar->Executable_MODELLER()) {
    Error("ERROR: $MYPROGNAME: Modeller executable not found: '".$cfgvar->Executable_MODELLER()."'\n",
            "Some of the executables not found.\n");## h-l error message
    MyExit(1);
}
unless(-d $cfgvar->InstallDir_MODPLUS()) {
    Error("ERROR: $MYPROGNAME: ModPlus installation directory not found: '".$cfgvar->InstallDir_MODPLUS()."'\n",
            "Some of the installation directories not found.\n");## h-l error message
    MyExit(1);
}


##check particular programs
##
$optionvalues{pdbdb_dirname} = $cfgvar->Path3dDb_PDB();
$optionvalues{scopedb_dirname} = $cfgvar->Path3dDb_SCOPe();
$optionvalues{web_address_scope} = $cfgvar->WebAddr3dDb_SCOPe();
$optionvalues{prog_modeller_mod} = $cfgvar->Executable_MODELLER();
$optionvalues{prog_modplus_modplus} = File::Spec->catfile($cfgvar->InstallDir_MODPLUS(),'modplus.pl');
unless(-f $optionvalues{prog_modplus_modplus}) {
    Error("ERROR: $MYPROGNAME: ModPlus executable not found: '".$optionvalues{prog_modplus_modplus}."'\n",
        "Incomplete software installed on the system.\n");## h-l error message
    MyExit(1);
}


## =============================================================================
## MAIN
##

my  $options = readopt->new($OPTFILENAME);
my  %optionkeys = (
    job_num_cpus => 'JOB_NUM_CPUS',
    comer_db => 'comer_db',
    modeller_key => 'modeller_key',
    model_all_pairs => 'model_all_pairs'
);
my  $nsyscpus = GetNCPUs();
my  $ncpus = $MAXNCPUs;

my  $errmsg = "ERROR: $MYPROGNAME: Invalid job options.\n";
my  $hlerrmsg = "Invalid job options.\n";

unless(VerifyOptionValues($cfgvar, $options,
        \%optionkeys, \%optionvalues, \$errmsg, \$hlerrmsg)) {
    Error($errmsg, $hlerrmsg);## err msg and h-l error message
    MyExit(1);
}

if($Config{useithreads}) {
    if($cfgvar->Exists($optionkeys{job_num_cpus})) {
        $ncpus = $cfgvar->GetValue($optionkeys{job_num_cpus});
        $ncpus = $nsyscpus if($nsyscpus > 0 && $nsyscpus < $ncpus);
    }
    else {
        print(STDERR "WARNING: Job option $optionkeys{job_num_cpus} not specified: ".
            "Using default #cpus=${ncpus}\n");
    }
}
else {
    $ncpus = 1;
    print(STDERR "WARNING: Perl compiled WITHOUT thread support: ".
        "Using only one cpu: #cpus=${ncpus}\n");
}

my  $modelallpairs = !($optionvalues{$optionkeys{model_all_pairs}} cmp '1');
my  $nqueries = 0;## :shared = 0;##number of individual queries/inputs
my  %inputs;## :shared;##all inputs divided

$errmsg = "ERROR: $MYPROGNAME: Distribution of the input to subdirs failed.\n";
$hlerrmsg = "Parsing or distribution of the input failed.\n";

ProgressMsg("Preparing input...\n");

unless(PreprocessInputFile($INPFILENAME, \$errmsg, \$hlerrmsg)) {
    Error($errmsg, $hlerrmsg);## err msg and h-l error message
    MyExit(1);
}

unless(DistrInputToSubdirs($INPFILENAME,
        $MAXNSEQS_perprog, $MAXNSEQS_multim, $modelallpairs,
        $inpdirname, $inpbasename, \$nqueries, \%inputs, \$errmsg, \$hlerrmsg)) {
    Error($errmsg, $hlerrmsg);## err msg and h-l error message
    MyExit(1);
}

$errmsg = '';
$hlerrmsg = '';
my $oneworker;

if($modelallpairs) {
    if(1<$nqueries) {
        ProgressMsg("Preparing alignments and modeling multiple query fragments...\n");
    } else {
        ProgressMsg("Preparing the alignment and modeling the query...\n");
    }
} else {
    ProgressMsg("Preparing the MSA and modeling the query...\n");
}

unless(RunThreadsAndWait(\&ProcessQuery_t,
        $oneworker=0, $ncpus, $nqueries, \%inputs,
        \%optionkeys, \%optionvalues, \$errmsg, \$hlerrmsg)) {
    Error($errmsg);## err msg and h-l error message
    ##continue on...
}

ProgressMsg("Verifying modeling output...\n");

unless(VerifyOutput($inpdirname, $inpbasename, $LOGFILENAME, $nqueries,
    \%inputs, \%optionkeys, \%optionvalues, \$errmsg, \$hlerrmsg)) {
    Error($errmsg, $hlerrmsg);
    MyExit(1);
}

ProgressMsg("Finalizing results...\n");

my  @resfilelist;

unless(MakeResultsList(\@resfilelist, $inpdirname, $RESLSTFILENAME, $nqueries,
    \%inputs, \$errmsg, \$hlerrmsg)) {
    Error($errmsg, $hlerrmsg);
    MyExit(1);
}

unless(CompressResults($TARPROG, $GZPROG, $inpdirname, $RESULTSFILENAME, 
    \@resfilelist, 
    $LOGFILENAME, $RESLSTFILENAME, $STAFILENAME, $ERRFILENAME,
    \$errmsg, \$hlerrmsg)) {
    Error($errmsg, $hlerrmsg);
    MyExit(1);
}

ProgressMsg("Finished.\n");

MyExit(0);

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
    unless(open(H, "/proc/cpuinfo")) {
        print(STDERR "WARNING: Unable to determine #cpu cores.\n");
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
    my  $msg = shift;
    return 0 unless(open(F,'>>',$STAFILENAME));
    print(F $msg);
    close(F);
    return 1;
}

## -----------------------------------------------------------------------------
## GetWarnings: extract warnings (if any) from file
##
sub GetWarnings
{
    my  $logfile = shift;##comer search log file
    my  $rlist = shift;##ref to the output list of warnings
    unless(open(F, $logfile)) {
        print(STDERR "WARNING: Unable to open log file: '$logfile'.\n");
        return 0;
    }
    @$rlist = grep(/WARNING/i,<F>);
    close(F);
    return 1;
}

## -----------------------------------------------------------------------------
## ProcessStructureInline: filter out unuseful records from the structure file
##
sub ProcessStructureInline
{
    my  $filename = shift;##full filename of the structure
    my  $rmsg = shift;##ref to a message formed on error
    my @pdbinfo;
    unless(open(F,$filename)) {
        $$rmsg = "Failed to open structure file: $filename";
        return 0;
    }
    while(<F>) {
        next if(/^(?:HETATM)/);
        push @pdbinfo, $_;
    }
    close(F);
    unless(open(F,'>',$filename)) {
        $$rmsg = "Failed to open structure file for writing: $filename";
        return 0;
    }
    print(F) foreach(@pdbinfo);
    close(F);
    return 1;
}

## -----------------------------------------------------------------------------
## AddFileToArchive: add a file to a given archive
##
sub AddFileToArchive
{
    my  $tarprog = shift;##full pathname to the tar program
    my  $archive = shift;##fullname of the resulting archive file
    my  $filepathname = shift;##full pathname to the file
    my  $dirname = shift;##name of directory where $filepathname is; may be empty
    my  $rerrmsg = shift;##ref to the error message string to be put in logs
    my  $rhlerrmsg = shift;##ref to the h-l error message string
    my  $create = shift;##flag of whether the archive is to be created

    my  $mysubname = (caller(0))[3];
    my  $preamb = "[ ${mysubname} ] ";
    my  $command = '';
    my  $opt = $create?'c':'r';
    my  $ret = 1;

    if($dirname) {
        $command = "${tarprog} -${opt}f \"${archive}\" -C \"$dirname\" \"${filepathname}\"";
    }
    else {
        my $filedirname = dirname($filepathname);
        my $filename = basename($filepathname);
        $command = "${tarprog} -${opt}f \"${archive}\" -C \"${filedirname}\" \"${filename}\"";
    }

    print(STDERR GetTime()."${preamb} ${command}\n");

    unless(ExecCommand($command)) {
        $$rerrmsg = "ERROR: $MYPROGNAME: $mysubname: Failed to add file to archive: '${filepathname}'\n";
        $$rhlerrmsg = "Failed to archive some of the results files.\n";
        return 0;
    }

    return $ret;
}

## -----------------------------------------------------------------------------
## VerifyOptionValues: verify filename and related information given in the 
## job options file
##
sub VerifyOptionValues 
{
    my  $cfgobj = shift;##backend configuration object
    my  $optobj = shift;##options object
    my  $roptionkeys = shift;##ref to the keys of options
    my  $roptionvalues = shift;## ref to the values of options
    my  $rerrmsg = shift;##ref to the error message string
    my  $rhlerrmsg = shift;##ref to the h-l error message string

    my  $mysubname = (caller(0))[3];
    my  $ret = 1;


    unless($optobj->Exists($$roptionkeys{comer_db})) {
        $$rerrmsg = "ERROR: $MYPROGNAME: $mysubname: ".
            "Option $$roptionkeys{comer_db} not specified in the job options file.\n";
        $$rhlerrmsg = "COMER profile database not specified.\n";
        return 0;
    }

    $$roptionvalues{$$roptionkeys{comer_db}} = '';
    $$roptionvalues{$$roptionkeys{comer_db}} = 
        $optobj->GetValue($$roptionkeys{comer_db});

    $$roptionvalues{$$roptionkeys{modeller_key}} = '';

    if($optobj->Exists($$roptionkeys{modeller_key})) {
        $$roptionvalues{$$roptionkeys{modeller_key}} = 
            $optobj->GetValue($$roptionkeys{modeller_key});
    }

    $$roptionvalues{$$roptionkeys{model_all_pairs}} = '0';

    if($optobj->Exists($$roptionkeys{model_all_pairs})) {
        $$roptionvalues{$$roptionkeys{model_all_pairs}} =
            $optobj->GetValue($$roptionkeys{model_all_pairs});
    }

    return $ret;
}

## -----------------------------------------------------------------------------
## PreprocessInputFile: parse the input file and change the name of each query 
## inline so that duplicate names are avoided (Modeller will crash if the 
## sequence name and a template name are the same)
##
sub PreprocessInputFile
{
    my  $inpfilename = shift;##the input to be parsed for multiple individual inputs
    my  $rerrmsg = shift;##ref to the error message string to be put in logs
    my  $rhlerrmsg = shift;##ref to the h-l error message string

    my  $mysubname = (caller(0))[3];
    my  $inputsep = qr/^\/\//;
    my  $suffix = '_to_model';
    my ($filecontents, $nseq) = ('',0);
    my  $ret = 1;

    unless(open(F, $inpfilename)) {
        $$rerrmsg = "ERROR: $MYPROGNAME: $mysubname: Failed to open the input file: '$inpfilename'\n";
        $$rhlerrmsg = "Input file not found.\n";
        return 0;
    }
    while(<F>) {
        next if(!eof(F) && /^\s*$/);
        $nseq = 0 if(/$inputsep/);
        if(/^>(\S*)(.*)$/ && $nseq == 0) {
            my ($id, $rem) = ($1, $2);
            $id =~ s/[^0-9a-zA-Z\.\-\+\^~=_]/_/g;
            $id .= $suffix unless $id =~ /$suffix$/;
            $_ = ">$id$rem\n";
        }
        $nseq++ if /^>/;
        $filecontents .= $_;
    }
    close(F);
    unless(open(F, '>', $inpfilename)) {
        $$rerrmsg = "ERROR: $MYPROGNAME: $mysubname: Failed to open the input file for writing: '$inpfilename'\n";
        $$rhlerrmsg = "Preprocessing of the input file failed.\n";
        return 0;
    }
    print(F $filecontents);
    close(F);
}

## -----------------------------------------------------------------------------
## DistrInputToSubdirs: parse the input and distribute individual queries to 
## subdirectories;
##
sub DistrInputToSubdirs
{
    my  $inpfilename = shift;##the input to be parsed for multiple individual inputs
    my  $maxnseqs_perprog = shift;##max #pairs when they are all to be modeled
    my  $maxnseqs_multim = shift;##max #sequences when one model is to be produced
    my  $modelallpairs = shift;##flag of modeling all pairs listed in the input file
    my  $dirname = shift;##name of directory where to create subdirectories
    my  $basename = shift;##name pattern for to-be-created files in subdirectories
    my  $rnqueries = shift;##ref to the number of queries
    my  $rinputs = shift;##ref to the hash of individual inputs
    my  $rerrmsg = shift;##ref to the error message string to be put in logs
    my  $rhlerrmsg = shift;##ref to the h-l error message string

    my  $mysubname = (caller(0))[3];
    my  $inputsep = qr/^\/\//;
    my ($qrycontents, $qrysubdir, $qrybasename, $qryfullname) = ('','','','');
    my  $nalnpairs = 0;
    my  $ret = 1;

    unless(open(F, $inpfilename)) {
        $$rerrmsg = "ERROR: $MYPROGNAME: $mysubname: Failed to open the input file: '$inpfilename'\n";
        $$rhlerrmsg = "Input file not found.\n";
        return 0;
    }
    $$rnqueries = 0;
    while(<F>) {
        next if(!eof(F) && /^\s*$/);
        $qrycontents .= $_;
        if(/$inputsep/) {
            $nalnpairs++;
            if($modelallpairs) {
                next unless $qrycontents;
                if($maxnseqs_perprog <= $$rnqueries) {
                    my $text = "\nWARNING: Number of alignment pairs (models) ".
                        "reduced to the maximum allowed: $maxnseqs_perprog\n";
                    Warning($text, $text, 1);##no e-mail
                    last;
                }
            } else {
                last if $$rnqueries;
                if($maxnseqs_multim <= $nalnpairs) {
                    my $text = "\nWARNING: Number of alignment templates ".
                        "reduced to the maximum allowed: $maxnseqs_multim\n";
                    Warning($text, $text, 1);##no e-mail
                } else {
                    next unless eof(F);
                }
            }
            $qrysubdir = File::Spec->catfile($dirname,"${basename}__${$rnqueries}");
            $qrybasename = File::Spec->catfile($qrysubdir,"${basename}__${$rnqueries}");
            $qryfullname = "${qrybasename}.${PWFEXT}";
            ;;
            unless( -d $qrysubdir || mkdir($qrysubdir)) {
                $$rerrmsg = "ERROR: $MYPROGNAME: $mysubname: Failed to create directory: '$qrysubdir'\n";
                $$rhlerrmsg = "Creating a directory failed.\n";
                $ret = 0;
                last;
            }
            unless(open(FIN, ">", $qryfullname)) {
                $$rerrmsg = "ERROR: $MYPROGNAME: $mysubname: Failed to open file for writing: '$qryfullname'\n";
                $$rhlerrmsg = "Failed to create a file.\n";
                $ret = 0;
                last;
            }
            unless(print(FIN $qrycontents)) {
                $$rerrmsg = "ERROR: $MYPROGNAME: $mysubname: Failed to write to file: '$qryfullname'\n";
                $$rhlerrmsg = "Write to a file failed.\n";
                close(FIN);
                $ret = 0;
                last;
            }
            close(FIN);
            $$rinputs{"${$rnqueries}_rcode"} = 0;
            $$rinputs{"${$rnqueries}_error"} = '';##error message
            $$rinputs{"${$rnqueries}_errhl"} = '';##high-level error message
            $$rinputs{"${$rnqueries}_warng"} = '';##warning message
            $$rinputs{"${$rnqueries}_wrnhl"} = '';##high-level warning message
            $$rinputs{"${$rnqueries}_isize"} = length($qrycontents);
            $$rinputs{"${$rnqueries}_input"} = $qryfullname;
            $$rinputs{"${$rnqueries}_bname"} = $qrybasename;
            $$rinputs{"${$rnqueries}_logfl"} = "${qrybasename}.log";
            $$rinputs{"${$rnqueries}_pirfl"} = '';
            $$rinputs{"${$rnqueries}_outfl"} = '';
            $$rinputs{"${$rnqueries}_tmpls"} = '';##templates used
            ;;
            $qrycontents = $qrysubdir = $qrybasename = $qryfullname = '';
            $$rnqueries++;
            next;
        }
    }
    close(F);

    if($$rnqueries < 1) {
        $$rerrmsg = "ERROR: $MYPROGNAME: $mysubname: No queries found in input.\n";
        $$rhlerrmsg = "Invalid input format: No records.\n";
        return 0;
    }

    return $ret;
}

## -----------------------------------------------------------------------------
## VerifyOutput: check for warning and error messages generated during the 
## execution of the theads
##
sub VerifyOutput
{
    my  $dirname = shift;##name of the job directory
    my  $basename = shift;##basename of the (whole) input transfered to the backend
    my  $logfile = shift;##fullname of the log file
    my  $nqueries = shift;##number of queries
    my  $rinputs = shift;##ref to the hash of individual inputs
    my  $roptionkeys = shift;##ref to the keys of options
    my  $roptionvalues = shift;## ref to the values of options
    my  $rerrmsg = shift;##ref to the error message string to be put in logs
    my  $rhlerrmsg = shift;##ref to the h-l error message string

    my  $mysubname = (caller(0))[3];
    my  $preamb = "[ ${mysubname} ] ";
    my  @qrynums = 0..$nqueries-1;
    my  $strhlp = ($nqueries < 2)? "the query": "all queries";
    my  $hlstrhlp = ($nqueries < 2)? "the input": "all inputs";
    my  $command = '';
    my  $ret = 1;

    for(my $q = 0; $q <= $#qrynums; $q++) {
        my $qnum = $qrynums[$q];
        ##print warnings and errors issued by threads if any
        if($$rinputs{"${qnum}_warng"} || $$rinputs{"${qnum}_wrnhl"}) {
            Warning($$rinputs{"${qnum}_warng"}, $$rinputs{"${qnum}_wrnhl"}, 1);##no e-mail
        }
        if($$rinputs{"${qnum}_rcode"}) { ##don't send an e-mail
            Error($$rinputs{"${qnum}_error"}, $$rinputs{"${qnum}_errhl"}, 1);
            next;
        }
        $command .= "1";
    }

    unless($command) {
        $$rerrmsg = "ERROR: $MYPROGNAME: $mysubname: Generating models for ${strhlp} failed.\n";
        $$rhlerrmsg = "Model generation for ${hlstrhlp} failed.\n";
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
    my  $rfilelist = shift;##ref to the list of files listed in the results list file
    my  $dirname = shift;##name of the job directory
    my  $reslstfile = shift;##name of the output results listing file to be generated
    my  $nqueries = shift;##number of queries
    my  $rinputs = shift;##ref to the hash of individual inputs
    my  $rerrmsg = shift;##ref to the error message string to be put in logs
    my  $rhlerrmsg = shift;##ref to the h-l error message string

    my  $mysubname = (caller(0))[3];
    my  $preamb = "[ ${mysubname} ] ";
    my  @qrynums = 0..$nqueries-1;
    my  $ret = 1;

    unless(open(F,'>',$reslstfile)){
        $$rerrmsg = "ERROR: $MYPROGNAME: $mysubname: Failed to open file for writing: '${reslstfile}'.\n";
        $$rhlerrmsg = "Opening file for writing failed.\n";
        return 0;
    }

    print(F "# Model_file PIR_file Template_list\n");

    for(my $q = 0; $q <= $#qrynums; $q++) {
        my $qnum = $qrynums[$q];
        ##skip queries with pending errors
        next if($$rinputs{"${qnum}_rcode"});

        my @files = ( $$rinputs{"${qnum}_outfl"},
                      $$rinputs{"${qnum}_pirfl"}
        );
        my $tmplids = $$rinputs{"${qnum}_tmpls"};

        $_ =~ s/^${dirname}\/*(.+)$/$1/ foreach @files;

        print(F "\"$files[0]\"\t\"$files[1]\"\t${tmplids}\n");

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
    my  $tarprog = shift;##full pathname to the tar program
    my  $gzprog = shift;##full pathname to the gzip program
    my  $dirname = shift;##name of the job directory
    my  $resultsfile = shift;##name of the output results archive file to be created
    my  $rfilelist = shift;##ref to the list of ($dirname-prefix-removed) files to include in the archive
    my  $comerlogfile = shift;##fullname of the log file to include in the archive
    my  $reslstfile = shift;##fullname of the results listing file to include in the archive
    my  $statusfile = shift;##fullname of the status messages file to include in the archive
    my  $errorfile = shift;##fullname of the high-level error messages file to include in the archive
    my  $rerrmsg = shift;##ref to the error message string to be put in logs
    my  $rhlerrmsg = shift;##ref to the h-l error message string

    my  $mysubname = (caller(0))[3];
    my  $preamb = "[ ${mysubname} ] ";
    my  $command = '';
    my  $ret = 1;

    return 0 unless(AddFileToArchive($tarprog, 
        $resultsfile, $reslstfile, '', $rerrmsg, $rhlerrmsg, 1)); #no dirname; 1==create

    for(my $i = 0; $i <= $#{$rfilelist}; $i++) {
        return 0 unless(AddFileToArchive($tarprog, 
            $resultsfile, $$rfilelist[$i], $dirname, $rerrmsg, $rhlerrmsg)); #dirname given
    }

    unless(AddFileToArchive($tarprog, 
        $resultsfile, $statusfile, '', $rerrmsg, $rhlerrmsg)){ #no dirname
        ##record an error, don't stop and  don't send an e-mail
        Error($$rerrmsg, $$rhlerrmsg, 1);
    }

    unless(AddFileToArchive($tarprog, 
        $resultsfile, $errorfile, '', $rerrmsg, $rhlerrmsg)){ #no dirname
        ##record an error, don't stop and  don't send an e-mail
        Error($$rerrmsg, $$rhlerrmsg, 1);
    }

    $command = "${gzprog} -f \"${resultsfile}\"";

    print(STDERR GetTime()."${preamb} ${command}\n\n");

    unless(ExecCommand($command)) {
        $$rerrmsg = "ERROR: $MYPROGNAME: $mysubname: Failed to gzip the archive: '${resultsfile}'\n";
        $$rhlerrmsg = "Failed to compress the archive of results files.\n";
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
    my  $subroutine_t = shift;##thread subroutine address
    my  $oneworker = shift;##flag instructing to use only one worker with spec. #cpus
    my  $ncpus = shift;##number of cpus
    my  $nqueries = shift;##number of queries
    my  $rinputs = shift;##ref to the hash of individual inputs
    my  $roptionkeys = shift;##ref to the keys of options
    my  $roptionvalues = shift;## ref to the values of options
    my  $rerrmsg = shift;##ref to the error message string to be put in logs
    my  $rhlerrmsg = shift;##ref to the h-l error message string

    my  $nworkers = 1;##number of threads
    my ($tmPIR, $tmOut, $procd) = (0,0,0);##timespans and fraction of processed queries
    ##query serial numbers sorted by query size:
    my  @qrynums = sort {$$rinputs{"${b}_isize"} <=> $$rinputs{"${a}_isize"}} 0..$nqueries-1;
    my  @workers;
    my  $ret = 1;

    if(!$oneworker && 
        $#qrynums+1 > $MAXNQUERIES_multicore)
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
                    $subroutine_t, 
                    $ncpus, $qnum, $roptionkeys, $roptionvalues, $rinputs
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
                $$rerrmsg = "ERROR: $MYPROGNAME: 
                    Invalid query/input number returned by thread ".$thr->tid().".\n";
                $$rhlerrmsg = '';##not to be shown
                $ret = 0;
                next;
            }
            $$rinputs{"${qnum}_rcode"} = $thretlst[1];##code of fail
            $$rinputs{"${qnum}_error"} = $thretlst[2];##error message
            $$rinputs{"${qnum}_errhl"} = $thretlst[3];##high-level msg
            $$rinputs{"${qnum}_warng"} = $thretlst[4];##warning message
            $$rinputs{"${qnum}_wrnhl"} = $thretlst[5];##high-level warning message
            $$rinputs{"${qnum}_pirfl"} = $thretlst[6] if $thretlst[6];##PIR filename
            $$rinputs{"${qnum}_outfl"} = $thretlst[7] if $thretlst[7];##model filename
            $$rinputs{"${qnum}_tmpls"} = $thretlst[8] if $thretlst[8];##comma-sep. template ids
            $tmPIR += ($$rinputs{"${qnum}_tmPIR"} = $thretlst[9]);##PIR file building timespan
            $tmOut += ($$rinputs{"${qnum}_tmOut"} = $thretlst[10]);##model generation timespan

            do{$nvalidqries--; next} if $$rinputs{"${qnum}_rcode"};

            $nsuccess++;

            next if $procd == int($nsuccess*10/$nvalidqries);

            $procd = int($nsuccess*10/$nvalidqries);

            my $tmdmsg = '';
            if($tmPIR && $tmOut) {
                $tmdmsg = sprintf("(time distr.: %.0f%% Data, %.0f%% model generation)",
                        $tmPIR*100./($tmPIR+$tmOut),$tmOut*100./($tmPIR+$tmOut));
            } elsif($tmPIR) {
                $tmdmsg = sprintf("(time distr.: %.0f%% Data)", $tmPIR*100./($tmPIR+$tmOut));
            } elsif($tmOut) {
                $tmdmsg = sprintf("(time distr.: %.0f%% model generation)", $tmOut*100./($tmPIR+$tmOut));
            }

            ProgressMsg("  $nsuccess/$nqries queries done $tmdmsg\n");
        }

        sleep(2) if $#joinable < 0;
    }


    return $ret;
}

## -----------------------------------------------------------------------------
## ProcessQuery_t: process one query: from preparing the alignment to model 
## generation
##
sub ProcessQuery_t
{
    my  $ncpus = shift;##number of cpus
    my  $qrynum = shift;##query serial number
    my  $roptionkeys_t = shift;##ref to the keys of options
    my  $roptionvalues_t = shift;## ref to the values of options
    my  $rinputs_t = shift;##ref to the hash of individual inputs

    my  $mysubname = (caller(0))[3];
    my  $preamb = "[ ${mysubname} ".threads->tid()." ] ";
    my  $rcode = 0;
    ## These are general return fields that should be returned even if not computed!
    ## error strings, warning messages, output MSA and model filenames, template id list, 
    ##   time for building the MSA, and time for model generation:
    my ($error,$errhl,$warng,$wrnhl,$pirfile,$outfile,$tmpls,$timePIR,$timeOut) = ('','','','','','','',0,0);

    my  $pdbdbdir = $$roptionvalues_t{pdbdb_dirname};
    my  $scopedbdir = $$roptionvalues_t{scopedb_dirname};
    my  $webaddrscope = $$roptionvalues_t{web_address_scope};
    my  $modeller = $$roptionvalues_t{prog_modeller_mod};
    my  $modplus = $$roptionvalues_t{prog_modplus_modplus};
    my  $comerdb = $$roptionvalues_t{$$roptionkeys_t{comer_db}};
    my  $modkey = $$roptionvalues_t{$$roptionkeys_t{modeller_key}};

    my  $inputfile = $$rinputs_t{"${qrynum}_input"};
    my  $fullnameonly = $$rinputs_t{"${qrynum}_bname"};
    my  $logfile = $$rinputs_t{"${qrynum}_logfl"};

    my ($resmsasfx,$resmsampsfx) = ('_resmsa','_resmsamodplus');

    my  $outputmsafile = "${fullnameonly}${resmsasfx}.${AFAEXT}";
    my  $outputmsampfile = "${fullnameonly}${resmsampsfx}.${AFAEXT}";
    my  $templatelstfile = "${fullnameonly}_templates.lst";
    my  $templatesdir = "${fullnameonly}_templates";
    my  $modelsdir = "${fullnameonly}_models";
    my  $target = '';

    my ($command, $cmbcmd) = ('','');
    my  $qrystartpos = -1;##query model start position


    $timePIR = time();

    unless(MakeAfaFromPwfa($qrynum, $inputfile, $outputmsafile, $logfile, \$qrystartpos,
            \$error, \$errhl, \$warng, \$wrnhl)) {
        return ($qrynum, 1, $error, $errhl, $warng, $wrnhl, $pirfile, $outfile, $tmpls, time()-$timePIR, $timeOut);
    }

    unless(PrepareAlnDataForModPlus($qrynum, $outputmsafile, $outputmsampfile, \$target, \$tmpls, $logfile, 
            $templatelstfile, $templatesdir, $comerdb, $pdbdbdir, $scopedbdir, $webaddrscope, 
            \$error, \$errhl, \$warng, \$wrnhl)) {
        return ($qrynum, 1, $error, $errhl, $warng, $wrnhl, $pirfile, $outfile, $tmpls, time()-$timePIR, $timeOut);
    }

    $timePIR = time() - $timePIR;
    $timeOut = time();

    unless(RunModPlus($qrynum, $outputmsampfile, \$pirfile, \$outfile, $target, $tmpls, $logfile, $qrystartpos,
            $templatesdir, $modelsdir, $modplus, $modeller, $modkey,
            \$error, \$errhl, \$warng, \$wrnhl)) {
        return ($qrynum, 1, $error, $errhl, $warng, $wrnhl, $pirfile, $outfile, $tmpls, $timePIR, time()-$timeOut);
    }

    $timeOut = time() - $timeOut;

    return ($qrynum, $rcode, $error, $errhl, $warng, $wrnhl, $pirfile, $outfile, $tmpls, $timePIR, $timeOut);
}

## -----------------------------------------------------------------------------
## MakeAfaFromPwfa: build the MSA in the afa format from pairwise alignments
##
sub MakeAfaFromPwfa
{
    my  $qrynum = shift;##query serial number
    my  $inpfilename = shift;##full name of input
    my  $ouputFILENAME = shift;##name of output file
    my  $logfilename = shift;##name of log file
    my  $rqrystartpos = shift;##query model start position
    my  $rerrmsg = shift;##ref to the error message string to be put in logs
    my  $rhlerrmsg = shift;##ref to the h-l error message string
    my  $rwrnmsg = shift;##ref to the warning message string to be put in logs
    my  $rhlwrnmsg = shift;##ref to the h-l warning message string

    my  $mysubname = (caller(0))[3];
    my  $preamb = "[ ${mysubname} ".threads->tid()." ] ";
    my  $command = '';
    my  @warnings;
    my  $ret = 1;
        
    ##run builder
    $command = "${PWFA2MSAprog} ".
            "-i \"${inpfilename}\" -o \"${ouputFILENAME}\" -a -f 0 -e 1.e19 ".
            "-N ${MAXNSEQS_perprog}  >>\"${logfilename}\" 2>&1";
    
    if(open(F,'>',$logfilename)){print(F GetTime()."${preamb} ${command}\n\n");close(F);}

    unless(ExecCommand($command)) {
        $$rerrmsg = "ERROR: $MYPROGNAME: $mysubname: Failed to consruct the MSA for model No.${qrynum}.\n";
        $$rhlerrmsg = "Consruction of the MSA for model No.${qrynum} failed.\n";
        return 0;
    }
    
    ##extract warnings if any; no check of return code
    GetWarnings($logfilename, \@warnings);

    if(0 <= $#warnings) {
        my $lst = join("", @warnings);
        my $text = "\nWarnings from building the MSA for model No.${qrynum}:\n\n$lst\n\n";
        $rwrnmsg .= $text;
        $rhlwrnmsg .= $text;
    }

    if(open(F,$ouputFILENAME)) {
        $_ = <F>;
        $$rqrystartpos = $1 if /^\S+\s+\((\d+)\-\d+\)/;
        close(F);
    }

    return $ret;
}

## -----------------------------------------------------------------------------
## PrepareAlnDataForModPlus: prepare the input MSA in the format suitable to 
## modplus and download required template structures for modeller
##
sub PrepareAlnDataForModPlus
{
    my  $qrynum = shift;##query serial number
    my  $inpfilename = shift;##full name of input
    my  $ouputFILENAME = shift;##name of output file
    my  $rtarget = shift;##ref to the target name (to return)
    my  $rtemplateids = shift;##ref to the string of template ids (to return)
    my  $logfilename = shift;##name of log file
    my  $templatelstfile = shift;##template listing filename
    my  $templatesdir = shift;##name of the directory of templates
    my  $comerdb = shift;##name of a COMER profile database
    my  $pdbdbdir = shift;##directory of pdb structures database
    my  $scopedbdir = shift;##directory of SCOPe structures database
    my  $webaddrscope = shift;##web address for downloading SCOPe domain structures
    my  $rerrmsg = shift;##ref to the error message string to be put in logs
    my  $rhlerrmsg = shift;##ref to the h-l error message string
    my  $rwrnmsg = shift;##ref to the warning message string to be put in logs
    my  $rhlwrnmsg = shift;##ref to the h-l warning message string

    my  $mysubname = (caller(0))[3];
    my  $preamb = "[ ${mysubname} ".threads->tid()." ] ";
    my ($command, $cmdbase, $cmdbasepdb, $cmdbasescop, $dirname) = ('','','','','');
    my ($outcontents, $text) = ('','');
    my  $ret = 1;

    unless(open(F,$inpfilename)) {
        $$rerrmsg = "ERROR: $MYPROGNAME: $mysubname: Failed to open the MSA file for model No.${qrynum}.\n";
        $$rhlerrmsg = "Consructed MSA file for model No.${qrynum} not found.\n";
        return 0;
    }
    my ($tmpl, $seqn)=('','');
    my (@templids, @records);##template ids
    while(<F>) {
        $seqn .= $_ unless /^\s*>/;
        if(eof(F) || /^\s*>(\S+)/) {
            if($tmpl && $seqn) {
                $tmpl =~ s/\|/_/g;
                push @templids, $tmpl;
                push @records, ">$tmpl ;\n$seqn";
            }
            $tmpl = $1;
            $seqn = '';
            next;
        }
    }
    close(F);

    if($#templids < 1 || $#records < 1) {
        $$rerrmsg = "ERROR: $MYPROGNAME: $mysubname: No valid pair of aligned sequences ".
            "found in the MSA file for model No.${qrynum}.\n";
        $$rhlerrmsg = "No valid pair of aligned sequences found in the MSA file for ".
            "model No.${qrynum} not found.\n";
        return 0;
    }

    $dirname = $pdbdbdir;
    $cmdbasepdb = "${DWNPDBSTRhelper} -o \"${templatesdir}\" --pdb \"${pdbdbdir}\" -i ";

    unless( -d $dirname || mkdir($dirname)) {
        $$rerrmsg = "ERROR: $MYPROGNAME: $mysubname: Failed to create a directory of structures: '$dirname' ".
            "(model No.${qrynum}).\n";
        $$rhlerrmsg = "Creating a directory failed (model No.${qrynum}).\n";
        return 0;
    }

    $dirname = $scopedbdir;
    $cmdbasescop = "${DWNSCPSTRhelper} -o \"${templatesdir}\" --pdb \"${scopedbdir}\" --add \"${webaddrscope}\" --id ";

##    } else {
##        $$rerrmsg = "ERROR: $MYPROGNAME: $mysubname: No or invalid COMER database specified: '$comerdb' ".
##            "(model No.${qrynum}).\n";
##        $$rhlerrmsg = "No or invalid database of protein structures specified: '$comerdb' ".
##            "(model No.${qrynum}).\n";
##        return 0;
##    }

    unless( -d $dirname || mkdir($dirname)) {
        $$rerrmsg = "ERROR: $MYPROGNAME: $mysubname: Failed to create a directory of structures: '$dirname' ".
            "(model No.${qrynum}).\n";
        $$rhlerrmsg = "Creating a directory failed (model No.${qrynum}).\n";
        return 0;
    }

    $$rtarget = $templids[0];

    ##start with template 1 as 0th is the target
    for(my $i = 1; $i <= $#templids; $i++) {
        if($templids[$i] =~ /^[\da-zA-Z]{4}_[^\s_]+$/) {
            $cmdbase = $cmdbasepdb;
        } elsif($templids[$i] =~ /^[deg][\da-z]{4}[\da-zA-Z_\.][\da-zA-Z_]$/) {
            $cmdbase = $cmdbasescop;
        } else {
            $text = "WARNING: ID $templids[$i] is not recognized to have a structure (model No.${qrynum}). ".
                "Alignment skipped.\n";
            $$rwrnmsg .= $text;
            $$rhlwrnmsg .= $text;
            next;
        }
        $command = "$cmdbase $templids[$i]  >>\"${logfilename}\" 2>&1";
        if(open(F,'>>',$logfilename)){print(F GetTime()."${preamb} ${command}\n\n");close(F);}
        unless(ExecCommand($command)) {
            $text = "WARNING: Failed to download template $templids[$i] for model No.${qrynum}. ".
                "Alignment skipped.\n";
            $$rwrnmsg .= $text;
            $$rhlwrnmsg .= $text;
            next;
        }
        my $templatefile = File::Spec->catfile($templatesdir,$templids[$i].".${ENTEXT}");
        unless(ProcessStructureInline($templatefile, \$text)) {
            $$rerrmsg = "ERROR: $MYPROGNAME: $mysubname: $text (model No.${qrynum}).\n";
            $$rhlerrmsg = "Processing of template $templids[$i] failed (model No.${qrynum}).\n";
            return 0;
        }
        $outcontents .= $records[$i];
        $$rtemplateids .= ',' if $$rtemplateids;
        $$rtemplateids .= $templids[$i];
    }

    unless($outcontents) {
        $$rerrmsg = "ERROR: $MYPROGNAME: $mysubname: No templates obtained for model No.${qrynum}.\n";
        $$rhlerrmsg = "No templates obtained for model No.${qrynum}.\n";
        return 0;
    }

    unless(open(F,'>',$ouputFILENAME)) {
        $$rerrmsg = "ERROR: $MYPROGNAME: $mysubname: Failed to open ModPlus MSA file for writing: model No.${qrynum}.\n";
        $$rhlerrmsg = "File open for writing an MSA for model No.${qrynum} failed.\n";
        return 0;
    }
    print(F $records[0]);
    print(F $outcontents);
    close(F);

    return $ret;
}

## -----------------------------------------------------------------------------
## RunModPlus: run the ModPlus wrapper for modeller and arrange results files
##
sub RunModPlus
{
    my  $qrynum = shift;##query serial number
    my  $inpfilename = shift;##full name of input
    my  $rouputPIRFILE = shift;##ref to the name of output PIR file
    my  $rouputPDBFILE = shift;##ref to the name of output model file
    my  $target = shift;##target name
    my  $templateids = shift;##string of template ids
    my  $logfilename = shift;##name of log file
    my  $qrystartpos = shift;##query model start position
    my  $templatesdir = shift;##name of the directory of templates
    my  $modelsdir = shift;##directory of locate output models
    my  $modplus = shift;##full pathname to modplus executable
    my  $modeller = shift;##full pathname to modeller executable
    my  $modkey = shift;##modeller key
    my  $rerrmsg = shift;##ref to the error message string to be put in logs
    my  $rhlerrmsg = shift;##ref to the h-l error message string
    my  $rwrnmsg = shift;##ref to the warning message string to be put in logs
    my  $rhlwrnmsg = shift;##ref to the h-l warning message string

    my  $mysubname = (caller(0))[3];
    my  $preamb = "[ ${mysubname} ".threads->tid()." ] ";
    my  $tmpl1 = (split(',',$templateids))[0];
    my  $command = '';
    my  $ret = 1;

    ##run ModPlus
    $command = "${modplus} ".
            "--in \"${inpfilename}\" --pdb \"${templatesdir}\" --dir \"${modelsdir}\" ".
            "--modeller \"${modeller}\" >>\"${logfilename}\" 2>&1";

    if(open(F,'>>',$logfilename)){print(F GetTime()."${preamb} ${command}\n\n");close(F);}

    unless(ExecCommand($command)) {
        $$rerrmsg = "ERROR: $MYPROGNAME: $mysubname: Failed to generate model No.${qrynum}.\n";
        $$rhlerrmsg = "Generation of model No.${qrynum} failed.\n";
        return 0;
    }

    $$rouputPIRFILE = File::Spec->catfile($modelsdir,"${target}-${tmpl1}","${target}.${PIREXT}");
    $$rouputPDBFILE = File::Spec->catfile($modelsdir,"${target}-${tmpl1}","${target}.B99990001.${PDBEXT}");

    unless(-f $$rouputPDBFILE) {
        $$rerrmsg = "ERROR: $MYPROGNAME: $mysubname: Final model No.${qrynum} not found: '$$rouputPDBFILE'.\n";
        $$rhlerrmsg = "No final model No.${qrynum}.\n";
        return 0;
    }

    if(-f $RENUMhelper && 1 < $qrystartpos) {
        ##renumber the residues in the model
        my $renumoutputfile = File::Spec->catfile($modelsdir,"${target}-${tmpl1}","${target}.renum.${PDBEXT}");
        $command = "python3 ${RENUMhelper} ".
            "-i \"${$rouputPDBFILE}\" -r ${qrystartpos} -o \"${renumoutputfile}\" >>\"${logfilename}\" 2>&1";

        if(open(F,'>>',$logfilename)){print(F "\n".GetTime()."${preamb} ${command}\n\n");close(F);}

        unless(ExecCommand($command)) {
            $$rwrnmsg = "WARNING: $MYPROGNAME: $mysubname: Residue renumbering failed for model No.${qrynum}.\n";
            $$rhlwrnmsg = "WARNING: Residue renumbering failed for model No.${qrynum}.\n";
        } else {
            my $outcontents = '';
            if(open(F,$$rouputPDBFILE)) {
                while(<F>) {last if /^(?:ATOM|HETATM)/; $outcontents .= $_}
                close(F);
            }
            if(open(F,$renumoutputfile)) {
                $outcontents .= $_ while(<F>);
                close(F);
                if(open(F,'>',$renumoutputfile)) {
                    print F $outcontents;
                    close(F);
                    $$rouputPDBFILE = $renumoutputfile;
                }
            }
        }
    }

    return $ret;
}

## General =====================================================================
##
##







## =============================================================================
## MyExit: print a code to file and exit
##
sub MyExit
{
    my $ecode = shift;##exit code
    if(open(F,">>",$ERRFILENAME)){print(F "\n${ecode}\n");close(F);}
    exit($ecode);
}

## -------------------------------------------------------------------
## Warning: output a warning message and optionally, send an email
##
sub Warning
{
    my ($msg, $hlerrmsg, $nosend) = @_;
    return Issue($msg, $hlerrmsg, $nosend,
        '',
        "Warning from comer-ws ($MYPROGNAME)");
}

## Error: output an error message and optionally, send an email
##
sub Error
{
    my ($msg, $hlerrmsg, $nosend) = @_;
    return Issue($msg, $hlerrmsg, $nosend,
        "ERROR: Issue on the server's backend side: ",
        "Error message from comer-ws ($MYPROGNAME)");
}

## Issue: output an issue message and send an email
##
sub Issue
{
    my ($msg, $hlerrmsg, $nosend,  $leadtxt1, $sbjct) = @_;
    print(STDERR GetTime().$msg);
    if($ERRFILENAME && $hlerrmsg) {
        if(open(EF,">>",$ERRFILENAME)) {
            print(EF $leadtxt1.$hlerrmsg);
            close(EF);
        } else {
            print(STDERR "ERROR: $MYPROGNAME: Write to error file skipped ".
                "due to the fail to open the file: '$ERRFILENAME'\n");
        }
    }
    return 1 if $nosend;
    unless(-f $SENDMAILprog) {
        print(STDERR "ERROR: $MYPROGNAME: Send-mail program not found: $SENDMAILprog\n");
        return 0;
    }
    my  $command = "$SENDMAILprog --sub \"$sbjct\" --body \"$msg\" 2>&1";
    print(STDERR "$command\n");
    return ExecCommand($command);
}

## -------------------------------------------------------------------
## try read lock files in a directory and return 1 if locks exist,
## 0 otherwise, -1 on error
##

sub ReadLocks
{
    return ReadFiles( shift, shift );
}

sub ReadFiles
{
    my  $dirname = shift;
    my  $pattern = shift;  ## pattern of files to look for as locks
    my  $reffiles = shift; ## reference to vector of files
    my  @files;
    my  $locref = defined( $reffiles )? $reffiles: \@files;

    unless( opendir( DIR, $dirname )) {
        printf( STDERR "ERROR: Cannot open directory $dirname.\n" );
        return -1;
    }
    @{$locref} = grep { /$pattern$/ && -f File::Spec->catfile($dirname,$_) } readdir( DIR );
    closedir( DIR );
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


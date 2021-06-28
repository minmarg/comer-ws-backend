#!/bin/bash

## (C) 2020 Mindaugas Margelevicius, Institute of Biotechnology, Vilnius University
## Helper script for HHblits search

dirname="$(dirname $0)"
[[ "${dirname:0:1}" != "/" ]] && dirname="$(pwd)/$dirname"
basename="$(basename $0)"

ROUNDS=2
EVALUE=1e-3
CPUs=2

## make .afa (aligned fasta) on normal exit
MAKEAFA=0
## maximum number of iterations allowed
MAXNROUNDS=8

## predefined and constant values used for the corresponding options:
PAR_Z=35000
PAR_B=35000
PAR_maxfilt=20000

usage="
Initiate HHblits search and write results in aligned FASTA format.
(C)2020 Mindaugas Margelevicius, Institute of Biotechnology, Vilnius University

Usage:
$basename <Options>

Options:

-i <input>     Input sequence (plain or in FASTA) or MSA in FASTA or 
               STOCKHOLM format, used as the query to search a database.
-o <output>    Filename pattern for output files of pairwise alignments 
               and MSA in FASTA (option -a).
-d <database>  HHsuite database prefix.
-n <num_its>   Number of iterations to search a database.
           Default=${ROUNDS}
-N <count>     Maximum number of sequences permitted to be included for
               the next iteration.
           Default=${PAR_maxfilt}
-e <e_value>   E-value threshold for output and inclusion of sequences 
               for the next iteration.
           Default=${EVALUE}
-p <num_cores> Number of CPU cores to use for search.
           Default=${CPUs}
-a             Produce aligned-fasta results file from the alignments 
               obtained in the last iteration.
-b <path>      HHsuite installation directory.
-h             short description.
"

echo

while getopts "i:o:d:n:N:e:p:ab:h" Option
do
    case $Option in
        i ) INPUT=${OPTARG} ;;
        o ) OUTPUT=${OPTARG} ;;
        d ) DB=${OPTARG} ;;
        n ) ROUNDS=${OPTARG} ;;
        N ) PAR_maxfilt=${OPTARG}
            PAR_Z=${PAR_maxfilt}; PAR_B=${PAR_maxfilt}
            ;;
        e ) EVALUE=${OPTARG} ;;
        p ) CPUs=${OPTARG} ;;
        a ) MAKEAFA=1 ;;
        b ) HHSDIR=${OPTARG} ;;
        h ) echo "$usage"; exit 0 ;;
        * ) echo ERROR: Unrecognized argument. >&2; exit 1 ;;
    esac
done
shift $(($OPTIND - 1))

fail=""
DBdir="$(dirname ${DB})"

if [[ -z "$fail" && ( -z "${INPUT}" || ! -f "${INPUT}" ) ]]; then 
  fail+="ERROR: Input file not found: \"${INPUT}\"\n"
fi
if [[ -z "$fail" && -z "${OUTPUT}" ]]; then
  fail+="ERROR: Output filename not provided.\n"
fi
if [[ -z "$fail" && ( -z "${DBdir}" || ! -d "${DBdir}" ) ]]; then
  fail+="ERROR: Database directory not found: \"${DBdir}\".\n"
fi
if [[ -z "$fail" && !( "${ROUNDS}" =~ ^[1-9][0-9]*$ ) ]]; then
  fail+="ERROR: Invalid number of iterations specified: ${ROUNDS}\n"
fi
if [[ -z "$fail" && !( "${PAR_maxfilt}" =~ ^[1-9][0-9]*$ ) ]]; then
  fail+="ERROR: Invalid number of sequences specified: ${PAR_maxfilt}\n"
fi
if [[ -z "$fail" && !( "${EVALUE}" =~ ^[0-9\.]+[eE]?[\-\+]?[0-9]+$ ) ]]; then
  fail+="ERROR: Invalid E-value threshold specified: ${EVALUE}\n"
fi
if [[ -z "$fail" && !( "${CPUs}" =~ ^[1-9][0-9]*$ ) ]]; then
  fail+="ERROR: Invalid number of CPU cores specified: ${CPUs}\n"
fi
if [[ -z "$fail" && ( -z "${HHSDIR}" || ! -d "${HHSDIR}" ) ]]; then 
  fail+="ERROR: HHsuite installation directory not found: \"${HHSDIR}\".\n"
fi

hhblits="${HHSDIR}/bin/hhblits"

if [[ -z "$fail" && ! -f "${hhblits}" ]]; then fail+="ERROR: HHblits program not found: \"${hhblits}\"\n"; fi

hhr2pwfa="${dirname}/hhr2pwfa.pl"
pwfa2msa="${dirname}/pwfa2msa.pl"
sto2afa="${dirname}/sto2afa.pl"

datecmd="date +%H:%M:%S"

if [[ -z "$fail" && ! -f "${hhr2pwfa}" ]]; then fail+="ERROR: Program not found: \"${hhr2pwfa}\"\n"; fi
if [[ -z "$fail" && ! -f "${pwfa2msa}" ]]; then fail+="ERROR: Program not found: \"${pwfa2msa}\"\n"; fi
if [[ -z "$fail" && ! -f "${sto2afa}" ]]; then fail+="ERROR: Program not found: \"${sto2afa}\"\n"; fi

if [ -n "$fail" ]; then
  echo -e "$fail" >&2
  exit 1
fi

if [ ${ROUNDS} -gt ${MAXNROUNDS} ]; then
    echo -e "WARNING: Decreasing #iterations to a maximum allowed: ${ROUNDS} -> ${MAXNROUNDS}\n" >&2
    ROUNDS=${MAXNROUNDS}
fi

## -----------------------------------------------------------------------------
## Functions
## =========
##
## MakeQueryFile: make the query file given an MSA in FASTA or sequence file
##
function MakeQueryFile()
{
    local locmqfinput="$1"
    local locmqfoutputfas="$2"
    cmd="perl -e '"
    cmd+='
      while(<>){if(/^\s*>(.+)$/){last if $desc;$desc=$1;next} chomp; s/[\-\.]//g; $seqn.=$_}
      $desc="Undefined" unless $desc;
      if($seqn=~/[^\-\.a-zA-Z]/){print STDERR "ERROR: Illegal symbols in query sequence.\n";exit(1)}
      print ">$desc\n$seqn\n"'
    cmd+="' \"${locmqfinput}\" >\"${locmqfoutputfas}\""
    echo -e "$(${datecmd}) ${cmd}\n"; eval "${cmd}"
    if [ $? -ne 0 ]; then echo -e "ERROR: MakeQueryFile: Command failed.\n" >&2; return 1; fi
    return 0
}
## MakeQueryAndInput: make the query file and prepare input for hmmsearch if 
## needed
##
function MakeQueryAndInput()
{
    local locnamepat="$1"
    local locext=$2
    local __locqueryvarname=$3 ##name of variable to contain the query pathname
    local __locinputvarname=$4 ##name of variable to contain the processed input filename

    local locfmt=$(perl -e ' ##<<possible formats: STO, AFA, FAS, and SEQ>>
      while(<>){if($.<2 && /^#\s+STOCKHOLM/){$f="STO";last}elsif(/^>/){$c++;if($c>1){$f="AFA";last}}}
      if($f){print "$f";exit(0)} if($c){print "FAS";exit(0)} print "SEQ";exit(0)
    ' "${locnamepat}.${locext}")

    local locqueryvalue="${locnamepat}_${locext}.qry"
    ##eval locqueryvalue="\"\${${__locqueryvarname}}\""
    local locinputvalue="${locnamepat}.${locext}"
    ##make input equal to the query file in FASTA
    locinputvalue="${locqueryvalue}"

    ##values for the case of AFA
    local locafaname="${locnamepat}.${locext}"

    case "${locfmt}" in
      STO ) ##convert STOCKHOLM1 MSA to FASTA with the consensus sequence written first (reference)
            locafaname="${locnamepat}_${locext}.afa"
            cmd="${sto2afa} -i \"${locnamepat}.${locext}\" -o \"${locafaname}\" -c"
            echo -e "$(${datecmd}) ${cmd}\n"; eval "${cmd}"
            if [ $? -ne 0 ]; then echo -e "ERROR: MakeQueryAndInput: Command failed.\n" >&2; return 1; fi
            ;&

      AFA ) ##input will be an aligned-fasta MSA file; assign appropriately
            locinputvalue="${locafaname}"
            ;&

      FAS ) ;&

      SEQ ) ##make the query file; <<continued from the above>>
            ##${locafaname} is the input and contains (presumably) only a plain sequence or that in FASTA
            cmd="MakeQueryFile \"${locafaname}\" \"${locqueryvalue}\""
            echo -e "$(${datecmd}) ${cmd}\n"; eval "${cmd}"
            if [ $? -ne 0 ]; then return 1; fi
            ;;

      * ) echo ERROR: MakeQueryAndInput: Unrecognized format of the input file. >&2
          return 1 ;;
    esac

    eval ${__locqueryvarname}="'${locqueryvalue}'"
    eval ${__locinputvarname}="'${locinputvalue}'"
    return 0
}
## -----------------------------------------------------------------------------

DNAME="$(dirname ${INPUT})"
BNAME="$(basename ${INPUT})"
NAMEPAT="${DNAME}/${BNAME%.*}" #longest name pattern w/o extension
EXT="${BNAME##*.}" #shortest matching extension

queryseq="${OUTPUT}.qry"
input="${INPUT}"

cmd="MakeQueryAndInput \"${NAMEPAT}\" ${EXT} queryseq input"
echo -e "$(${datecmd}) ${cmd}\n"; eval "${cmd}"

if [ $? -ne 0 ]; then echo -e "ERROR: Command failed.\n" >&2; exit 1; fi

DNAME="$(dirname ${input})"
BNAME="$(basename ${input})"
NAMEPAT="${DNAME}/${BNAME%.*}"
EXT="${BNAME##*.}"

sfx="hhb_n${ROUNDS}"
hhbnhhr="${NAMEPAT}.${sfx}.hhr"
##NOTE: written in psi, a3m, and hmm format are the results of one before the last iteration:
##hhbnpsi="${NAMEPAT}.${sfx}.psi"
hhbnpwfa="${OUTPUT}.pwfa"
hhbnafa="${OUTPUT}.afa"

if [ "${queryseq}" != "${OUTPUT}.qry" ]; then
    cmd="cp \"${queryseq}\" \"${OUTPUT}.qry\""
    echo -e "$(${datecmd}) ${cmd}\n"; eval "${cmd}"
    if [ $? -ne 0 ]; then echo -e "ERROR: Failed to copy the query file.\n" >&2; exit 1; fi
fi

cmd="largerE=\$(perl -e \"print (${EVALUE}<1)? ${EVALUE}*10.: ${EVALUE}*2.\")"
echo -e "$(${datecmd}) ${cmd}\n"; eval "${cmd}"

if [ $? -ne 0 ]; then echo -e "ERROR: Increasing e-value threshold failed.\n" >&2; exit 1; fi

cmd="${hhblits} -v 1 -i \"${input}\" -d \"${DB}\" -o \"${hhbnhhr}\" " #-opsi \"${hhbnpsi}\" "
cmd+="-n ${ROUNDS} -e ${EVALUE} -E ${largerE} -cpu ${CPUs} "
cmd+="-M first -Z ${PAR_Z} -B ${PAR_B} -maxfilt ${PAR_maxfilt}"
echo -e "$(${datecmd}) ${cmd}\n"; eval "${cmd}"

if [ $? -ne 0 ]; then echo -e "ERROR: Command failed.\n" >&2; exit 1; fi

## extract pairwise alignments from the output; 
cmd="${hhr2pwfa} -i \"${hhbnhhr}\" -o \"${hhbnpwfa}\" -e ${EVALUE}"
echo -e "$(${datecmd}) ${cmd}\n"; eval "${cmd}"

if [ $? -ne 0 ]; then echo -e "ERROR: Command failed.\n" >&2; exit 1; fi

if [ ${MAKEAFA} -eq 1 ]; then
    ## convert pairwise fasta-like format to aligned fasta
    cmd="${pwfa2msa} -f 0 -i \"${hhbnpwfa}\" -o \"${hhbnafa}\" -q \"${queryseq}\" -e ${EVALUE}"
    echo -e "$(${datecmd}) ${cmd}\n"; eval "${cmd}"

    if [ $? -ne 0 ]; then echo -e "ERROR: Command failed.\n" >&2; exit 1; fi
fi

exit 0


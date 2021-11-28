#!/bin/bash

## (C) 2020 Mindaugas Margelevicius, Institute of Biotechnology, Vilnius University
## Perform iterative HMM search for the same input sequence (i.e., HMMs constructed 
## at every iteration describe the same input) using HMMER

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
## maximum number of sequence to include for the next iteration
MAXNSEQS=20000

usage="
Perform iterative HMM search for the same input sequence/MSA (i.e., HMMs
constructed at every iteration describe the same input) using HMMER.
(C)2020 Mindaugas Margelevicius, Institute of Biotechnology, Vilnius University

Usage:
$basename <Options>

Options:

-i <input>     Input sequence (plain or in FASTA) or MSA in FASTA or 
               STOCKHOLM format, used as the query to search a database.
-o <output>    Filename pattern for output files of pairwise alignments 
               and MSA in FASTA (option -a).
-d <database>  Sequence database file.
-n <num_its>   Number of iterations using the result of an iteration 
               as the query for the next iteration of search.
           Default=${ROUNDS}
-N <count>     Maximum number of sequences permitted to be included for
               the next iteration.
           Default=${MAXNSEQS}
-e <e_value>   E-value threshold for output and inclusion of sequences 
               for the next iteration.
           Default=${EVALUE}
-p <num_cores> Number of CPU cores to use for search.
           Default=${CPUs}
-a             Produce aligned-fasta results file from the alignments 
               obtained in the last iteration.
-b <path>      HMMER installation directory.
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
        N ) MAXNSEQS=${OPTARG} ;;
        e ) EVALUE=${OPTARG} ;;
        p ) CPUs=${OPTARG} ;;
        a ) MAKEAFA=1 ;;
        b ) HMMERDIR=${OPTARG} ;;
        h ) echo "$usage"; exit 0 ;;
        * ) echo ERROR: Unrecognized argument. >&2; exit 1 ;;
    esac
done
shift $(($OPTIND - 1))

fail=""

if [[ -z "$fail" && ( -z "${INPUT}" || ! -f "${INPUT}" ) ]]; then 
  fail+="ERROR: Input file not found: \"${INPUT}\"\n"
fi
if [[ -z "$fail" && -z "${OUTPUT}" ]]; then
  fail+="ERROR: Output filename not provided.\n"
fi
if [[ -z "$fail" && ( -z "${DB}" || ! -f "${DB}" ) ]]; then
  fail+="ERROR: Database file not found: \"${DB}\".\n"
fi
if [[ -z "$fail" && !( "${ROUNDS}" =~ ^[1-9][0-9]*$ ) ]]; then
  fail+="ERROR: Invalid number of iterations specified: ${ROUNDS}\n"
fi
if [[ -z "$fail" && !( "${MAXNSEQS}" =~ ^[1-9][0-9]*$ ) ]]; then
  fail+="ERROR: Invalid number of sequences specified: ${MAXNSEQS}\n"
fi
if [[ -z "$fail" && !( "${EVALUE}" =~ ^[0-9\.]+[eE]?[\-\+]?[0-9]+$ ) ]]; then
  fail+="ERROR: Invalid E-value threshold specified: ${EVALUE}\n"
fi
if [[ -z "$fail" && !( "${CPUs}" =~ ^[1-9][0-9]*$ ) ]]; then
  fail+="ERROR: Invalid number of CPU cores specified: ${CPUs}\n"
fi
if [[ -z "$fail" && ( -z "${HMMERDIR}" || ! -d "${HMMERDIR}" ) ]]; then 
  fail+="ERROR: HMMER installation directory not found: \"${HMMERDIR}\".\n"
fi

jackhmmer="${HMMERDIR}/bin/jackhmmer"
hmmsearch="${HMMERDIR}/bin/hmmsearch"
hmmbuild="${HMMERDIR}/bin/hmmbuild"

if [[ -z "$fail" && ! -f "${jackhmmer}" ]]; then fail+="ERROR: HMMER program not found: \"${jackhmmer}\"\n"; fi
if [[ -z "$fail" && ! -f "${hmmsearch}" ]]; then fail+="ERROR: HMMER program not found: \"${hmmsearch}\"\n"; fi
if [[ -z "$fail" && ! -f "${hmmbuild}" ]]; then fail+="ERROR: HMMER program not found: \"${hmmbuild}\"\n"; fi

hmmer2pwfa="${dirname}/hmmer2pwfa.pl"
pwfa2msa="${dirname}/pwfa2msa.pl"
sto2afa="${dirname}/sto2afa.pl"
afa2sto="${dirname}/afa2sto.pl"

datecmd="date +%H:%M:%S"

if [[ -z "$fail" && ! -f "${hmmer2pwfa}" ]]; then fail+="ERROR: Program not found: \"${hmmer2pwfa}\"\n"; fi
if [[ -z "$fail" && ! -f "${pwfa2msa}" ]]; then fail+="ERROR: Program not found: \"${pwfa2msa}\"\n"; fi
if [[ -z "$fail" && ! -f "${sto2afa}" ]]; then fail+="ERROR: Program not found: \"${sto2afa}\"\n"; fi
if [[ -z "$fail" && ! -f "${afa2sto}" ]]; then fail+="ERROR: Program not found: \"${afa2sto}\"\n"; fi

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
## AdjustHMM: modify an HMM to make match state symbols correspond to the 
## reference states!
##
function AdjustHMM()
{
    local locahinputhmm="$1"
    local locahoutputhmm="$2"
    cmd="perl -e '"
    cmd+='while(<>){s/^(\s+\d+\s+.+\d+\s+)([\w\-]\s)([\w\-]\s)([\w\-]\s[\w\-])$/$1$3$3$4/;print}'
    cmd+="' \"${locahinputhmm}\" >\"${locahoutputhmm}\""
    echo -e "$(${datecmd}) ${cmd}\n"; eval "${cmd}"
    if [ $? -ne 0 ]; then echo -e "ERROR: AdjustHMM: Command failed.\n" >&2; return 1; fi
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
    local locstoname="${locnamepat}_${locext}.sto"

    case "${locfmt}" in
      STO ) ##convert STOCKHOLM1 MSA to FASTA with the consensus sequence written first
            locafaname="${locnamepat}_${locext}.afa"
            locstoname="${locnamepat}_${locext}_afa.sto"
            cmd="${sto2afa} -i \"${locnamepat}.${locext}\" -o \"${locafaname}\" -c"
            echo -e "$(${datecmd}) ${cmd}\n"; eval "${cmd}"
            if [ $? -ne 0 ]; then echo -e "ERROR: MakeQueryAndInput: Command failed.\n" >&2; return 1; fi
            ;&

      AFA ) ##convert aligned fasta back to STOCKHOLM1, which describes now the consensus
            cmd="${afa2sto} -i \"${locafaname}\" -o \"${locstoname}\""
            echo -e "$(${datecmd}) ${cmd}\n"; eval "${cmd}"
            if [ $? -ne 0 ]; then echo -e "ERROR: MakeQueryAndInput: Command failed.\n" >&2; return 1; fi

            ##make an HMM; use --hand to instruct it to use reference annotation for match states!
            local lochmmname="${locnamepat}_${locext}_afa_sto.hmm"
            cmd="${hmmbuild} --amino --hand \"${lochmmname}\" \"${locstoname}\""
            echo -e "$(${datecmd}) ${cmd}\n"; eval "${cmd}"
            if [ $? -ne 0 ]; then echo -e "ERROR: MakeQueryAndInput: Command failed.\n" >&2; return 1; fi

            ##adjust appropriately the HMM
            local lochmmCNname="${locnamepat}_${locext}_afa_sto_CN.hmm"
            cmd="AdjustHMM \"${lochmmname}\" \"${lochmmCNname}\""
            echo -e "$(${datecmd}) ${cmd}\n"; eval "${cmd}"
            if [ $? -ne 0 ]; then return 1; fi

            locinputvalue="${lochmmCNname}"
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

hmmeroutputpwfa="${OUTPUT}.pwfa"
hmmeroutputafa="${OUTPUT}.afa"

if [ "${queryseq}" != "${OUTPUT}.qry" ]; then
    cmd="cp \"${queryseq}\" \"${OUTPUT}.qry\""
    echo -e "$(${datecmd}) ${cmd}\n"; eval "${cmd}"
    if [ $? -ne 0 ]; then echo -e "ERROR: Failed to copy the query file.\n" >&2; exit 1; fi
fi

for ((i=0; i<ROUNDS; i++)); do
    sfx="hmm$((i+1))"
    outfile="${NAMEPAT}.${sfx}.out"
    pwfafile="${NAMEPAT}.${sfx}.pwfa"
    stofile="${NAMEPAT}.${sfx}.sto"
    hmmfile="${NAMEPAT}.${sfx}.hmm"
    hmmCNfile="${NAMEPAT}.${sfx}.CN.hmm"

    if [[ "${EXT}" == "fa" || "${EXT}" == "qry" ]]; then
        ## hmmer N1
        cmd="${jackhmmer} -o \"${outfile}\" -N 1 -E ${EVALUE} --domE ${EVALUE} --cpu ${CPUs} \"${queryseq}\" \"${DB}\""
        echo -e "$(${datecmd}) ${cmd}\n"; eval "${cmd}"
    elif [ "${EXT}" == "hmm" ]; then
        ## profile search
        cmd="${hmmsearch} -o \"${outfile}\" -E ${EVALUE} --domE ${EVALUE} --cpu ${CPUs} \"${input}\" \"${DB}\""
        echo -e "$(${datecmd}) ${cmd}\n"; eval "${cmd}"
    else
        echo -e "ERROR: Unsupported file extension!\n" >&2
        exit 1
    fi

    if [ $? -ne 0 ]; then echo -e "ERROR: Command failed.\n" >&2; exit 1; fi


    if [[ $((i+1)) -ge ${ROUNDS} ]]; then pwfafile="${hmmeroutputpwfa}"; fi  ##***

    ## extract pairwise alignments from the output; 
    ## consider using a larger value for E-value to include domains that individually score low
    cmd="${hmmer2pwfa} -i \"${outfile}\" -o \"${pwfafile}\" -r 1 -e ${EVALUE}" ##<<larger E here
    echo -e "$(${datecmd}) ${cmd}\n"; eval "${cmd}"

    if [ $? -ne 0 ]; then echo -e "ERROR: Command failed.\n" >&2; exit 1; fi

    if [[ ${MAKEAFA} != 1 && $((i+1)) -ge ${ROUNDS} ]]; then break; fi  ##***


    ## make STOCKHOLM1 MSA from the pairwise alignments, which will include the aligned query sequence
    cmd="${pwfa2msa} -f 2 -i \"${pwfafile}\" -o \"${stofile}\" -q \"${queryseq}\" -N ${MAXNSEQS} -e ${EVALUE}" ##<<larger E here
    echo -e "$(${datecmd}) ${cmd}\n"; eval "${cmd}"

    if [ $? -ne 0 ]; then echo -e "ERROR: Command failed.\n" >&2; exit 1; fi

    ## make an HMM; NOTE: use --hand to instruct it to use reference annotation for match states!
    cmd="${hmmbuild} --amino --hand \"${hmmfile}\" \"${stofile}\""
    echo -e "$(${datecmd}) ${cmd}\n"; eval "${cmd}"

    if [ $? -ne 0 ]; then echo -e "ERROR: Command failed.\n" >&2; exit 1; fi

    ## modify the profile to make match state symbols correspond to the reference states!
    cmd="AdjustHMM \"${hmmfile}\" \"${hmmCNfile}\""
    echo -e "$(${datecmd}) ${cmd}\n"; eval "${cmd}"

    if [ $? -ne 0 ]; then echo -e "ERROR: Command failed.\n" >&2; exit 1; fi

    input="${hmmCNfile}"
    EXT=hmm
done

if [ ${MAKEAFA} -eq 1 ]; then
    ## convert STOCKHOLM1 MSA obtained in the last round to FASTA
    cmd="${sto2afa} -i \"${stofile}\" -o \"${hmmeroutputafa}\""
    echo -e "$(${datecmd}) ${cmd}\n"; eval "${cmd}"

    if [ $? -ne 0 ]; then echo -e "ERROR: Command failed.\n" >&2; exit 1; fi
fi

exit 0


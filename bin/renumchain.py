#!/usr/bin/env python

##
## (C)2019 Mindaugas Margelevicius, Vilnius University
##

import sys, os
from Bio.PDB import PDBParser, PDBIO
from optparse import OptionParser

description = "Re-enumerate the residues of the specified chain of a PDB file."

parser = OptionParser(description=description)

parser.add_option("-i", "--infile", dest="input",
                  help="Input PDB file", metavar="FILE")
parser.add_option("-c", "--chain", dest="cid", default=" ",
                  help="PDB chain id to reenumerate", metavar="ID")
parser.add_option("-r", "--rnum", dest="resnum", default=1, type=int,
                  help="Residue number to start the reenumeration with", \
                    metavar="NUMBER")
parser.add_option("-o", "--outfile", dest="output",
                  help="Output PDB file", metavar="FILE")

(options, args) = parser.parse_args()

if not (options.input and options.output):
    sys.stderr.write("ERROR: Input and/or output file is not provided.\n")
    sys.exit(1)

if len(options.cid) != 1:
    sys.stderr.write("ERROR: Invalid chain id specified.\n")
    sys.exit(1)

if options.input and not os.path.isfile(options.input):
    sys.stderr.write("ERROR: Input file does not exist: "+options.input+"\n")
    sys.exit(1)


basename = os.path.basename(options.input)
name = os.path.splitext(basename)[0]

pdb = PDBParser().get_structure(name, options.input)

for model in pdb:
    for chain in model:
        if chain.get_id() != options.cid: continue
        resfst = chain.get_list()[0].get_id()
        #if options.resnum < resfst[1]:
        for residue in chain:
            resid = list(residue.id)
            resid[1] = options.resnum + (resid[1]-resfst[1])
            residue._id = tuple(resid)
        #elif resfst[1] < options.resnum:
        #    for residue in reversed(chain):
        #        resid = list(residue.id)
        #        resid[1] = options.resnum + (resid[1]-resfst[1])
        #        residue.id = tuple(resid)

io = PDBIO()
io.set_structure(pdb)
io.save(options.output)

#<<>>

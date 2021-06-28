#!/usr/bin/env python

import sys, os
from Bio.PDB import PDBParser, MMCIFParser, PDBIO
from Bio.PDB.PDBIO import Select
from optparse import OptionParser

description = "Extract a full chain from a PDB file."
chainNA = "n/a"

parser = OptionParser(description=description)

parser.add_option("-i", "--infile", dest="input",
                  help="input PDB file", metavar="FILE")
parser.add_option("-c", "--chain", dest="cid", default=chainNA,
                  help="PDB chain id to extract from file; " \
                      "by default, the first chain is selected", metavar="ID")
parser.add_option("-o", "--outfile", dest="output",
                  help="output PDB file", metavar="FILE")

(options, args) = parser.parse_args()

if not (options.input and options.output):
    sys.stderr.write("ERROR: Input and/or output file is not provided.\n")
    sys.exit()

if 4 < len(options.cid) or len(options.cid) < 1:
    sys.stderr.write("ERROR: Invalid chain id specified.\n")
    sys.exit()

if options.input and not os.path.isfile(options.input):
    sys.stderr.write("ERROR: Input file does not exist: "+options.input+"\n")
    sys.exit()


class ChainSelect(Select):
    def __init__(self, cid):
        self.cid = cid

    def accept_chain(self, chain):
        if chain.get_id() == self.cid:
            return 1
        else:
            return 0


basename = os.path.basename(options.input)
name, extension = os.path.splitext(basename)

if extension == '.cif':
    pdb = MMCIFParser().get_structure(name, options.input)
else:
    pdb = PDBParser().get_structure(name, options.input)

model = pdb.get_list()[0]

if options.cid == chainNA:
    ##get the first chain if not specified
    chain = model.get_list()[0]
    options.cid = chain.get_id()
else:
    chain = model[options.cid]

if 1 < len(options.cid):
    newid = options.cid[0]
    if model.has_id(newid):
        model.detach_child(newid)
    chain.id = newid

io = PDBIO()
io.set_structure(pdb) #chain)
io.save(options.output, ChainSelect(chain.get_id()))


COMER web server backend module for protein analysis by homology

(C)2020-2021 Mindaugas Margelevicius,
Institute of Biotechnology, Vilnius University

# Description

   [The COMER web server](https://bioinformatics.lt/comer) is a web 
   server for protein analysis based on the sequence alignment and 
   homology search methods [COMER2](https://github.com/minmarg/comer2) and 
   [COTHER](https://github.com/minmarg/cother). They exhibit sensitive, 
   accurate, and fast homology searches.

   The backend module represents computation logic without a graphical 
   interface. It can work independently as stand-alone software and be 
   useful for installing the core of the webserver locally, provided 
   the SLURM workload manager is installed in the system and variable 
   values and the paths to the external software are properly set in the 
   configuration file `var/comer-ws-backend.conf`. The programs in this 
   package do not require adjustments.

   The backend module provides these services as starting points:

  *  services/comersearch\_SRVC.sh

   &nbsp;&nbsp;&nbsp;&nbsp;
   The service of GPU-accelerated homology searches based on profile-profile 
   comparison using COMER2

  *  services/cothersearch\_SRVC.sh

   &nbsp;&nbsp;&nbsp;&nbsp;
   The service of GPU-accelerated homology (and analogy) searches by threading
   using COTHER

  *  services/comerMSA\_SRVC.sh

   &nbsp;&nbsp;&nbsp;&nbsp;
   The service of building a multiple sequence alignment (MSA) based on 
   (resulting accurate profile-profile) pairwise alignments

  *  services/comer3D\_SRVC.sh

   &nbsp;&nbsp;&nbsp;&nbsp;
   The service of generating 3D structural models using the structures of 
   identified proteins as templates and produced alignments as restraints

   The input parameters to all these services are (i) a filename pattern for 
   two input files, one of which (.in) includes input data and the other 
   (.options) is the options file (please find an example `var/job.options`),
   and (ii) the directory where the input files can be found. Input data is 
   sequence, MSA, and/or profile querie(s) separated by the line `//` for the 
   first two services, while it is a list of pairwise alignments for the last 
   two services. Please run the services with the `-h` option to see a more 
   detailed description of command-line options.

# License

   The COMER web server backend module is licensed under GNU General Public 
   License version 3. Please find the LICENSE and COPYING files 
   included in the software package.

# Funding

The work was supported by the European Regional Development Fund 
[grant number 01.2.2-LMT-K-718-01-0028]

---

Contact: <mindaugas.margelevicius@bti.vu.lt>


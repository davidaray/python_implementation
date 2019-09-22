import argparse
import shutil
import re
import os
import subprocess
from Bio import SeqIO
import pandas as pd
import numpy as np
from pybedtools import BedTool
from pyfaidx import Fasta
import logging
import time
import datetime
pd.options.mode.chained_assignment = None  # default='warn'

LOGGER = logging.getLogger(__name__)
'''The text below is from the Hubley, perl version of this script. It provides general information but some aspecs have changed as we converted to python.

   Given a .2bit file and either the batch count or batch size, create a 
   directory structure with subdirectories for each batch and a set of 
   scripts for use with David's Sun Grid Engine Cluster.  For details 
   on how this script works and how best to use David's cluster with
   RepeatMasker, see the "background" section below.

   An example run looks like this:

      % cd /lustre/work/daray/GENOMES/Hubley.genomes
      % mkdir mm10; cd mm10
      % ... scp mm10.unmasked.2bit from UCSC to mm10/
      % generateSGEClusterRun.pl -twoBit mm10.unmasked.2bit \
           -batch_count 24 \
           -species "Mus musculus" \
           -genomeDir /lustre/work/daray/GENOMES/Hubley.genomes/mm10
      % /bin/csh < qsub.sh
          # Submits 24 jobs to the queue
      % tail RMPart/*/run.log
          # view the end of the log file for each batch in the queue
          # The last few lines should look like this for a completed
          # batch:
          #   cycle 10 .......................................
          #   Generating output......................
          #   masking
          #   done
          # 
      % cat RMPart/*/run.log | egrep -v "identifying|Checking|refining|cycle|Generating|processing|masking|Search Engine:|Master RepeatMasker|RepeatMasker version"| egrep -v "^$" | more
          # When all batches are complete, check for errors.  This should 
          # display a header plus two lines per batch:
          #     analyzing file 000.fa
          #     done
          # If any messages appear in between these two lines take note 
          # and concat the RepeatMasker folks.
      % qsub doLift.sh
         
          # cluster to combine the results, and run some high-level
          # analsysis scripts on the results. When complete all final
          # files are gziped.  The RMPart directory may be deleted
          # once the results have been verifed.

   BACKGROUND
   ==========

   Sun Grid Engine (SGE) RepeatMasker Batch Template
   
   RepeatMasker is a shared memory multi-threaded ( forking ) 
   perl program.  Each invocation of RepeatMasker must run on a
   single node but can optionally use more than one processor on
   that node using the RepeatMasker "-pa #" parameter.  
  
   Site Specifics:
      Job Submission Node: hrothgar.hpcc.ttu.edu
      Cluster Configuration:
          Nodes:
              ray512cc@compute-17-1.local    40 slots ( large memory )
              raycc@compute-17-2.local       40 slots
              raycc@compute-17-3.local       40 slots
              raycc@compute-17-4.local       40 slots
              raycc@compute-17-5.local       20 slots
              raycc@compute-17-6.local       20 slots
              raycc@compute-17-7.local       20 slots
              raycc@compute-17-8.local       20 slots
                                            ---------
                                   Total    240 slots
  
         Queues:
          Two queues one for the large memory node and another for
          all other nodes
  
         Parallel Execution Environments
          On hrothgar SGE has two general parallel execution ( -pe )
          settings: fill and sm. 
   
          fill:  Allocate slots on the cluster on various nodes
                 in a round-robin fashion until the total number
                 requested is reached.  This option cannot be
                 used with shared-memory tasks like RepeatMasker
                 as slots are not gauranteed to be on the same
                 node.  
  
          sm:    Allocate all slots requested on a single node.
                 This is the general "shared memory" environment.
  
  
  
        RepeatMasker Run Recommendation:
           There is a tradeoff made when deciding how fine grained
           to make the job size of a particular RepeatMasker run.
           You could place the entire genome in one batch and use
           one invocation of RepeatMasker running with all the slots
           on one node. (ie. queue ray512cc, -pe sm 40, and 
           RepeatMasker -pa 40 ).  This would ensure that no one
           else got onto this machine while we were running and cause
           it to be overloaded.  It would also mean that we couldn't
           benefit from free slots on other nodes to complete this run
           even faster.  Also it means that if there are only 39
           slots available on ray512cc we will have to wait until all
           40 are available before we can begin running.  Lastly
           as a single computer running a genome could take quite
           some time, we have essentially hogged this node and made
           sure that no one else will be able to use this node for 
           perhaps many days.
  
           We could go to the other extreme and divide the genome into
           many small batches and allocate just a single slot per job.
           ( ie. queue raycc,ray512cc, -pe sm 1, and RepeatMasker -pa 1 ).
           This would mean that we have short running jobs allowing for
           others in the queue to jump up in priority and get some 
           slots.  Also this allows for us to use the entire cluster
           rather than just one node.  The downside is that results
           need to be combined after a run is complete to build a complete
           dataset, boundary ( edge cases ) can occur where repeats may
           be missed at boundaries of batches, and finally the overhead 
           many RepeatMasker invocations starts to reduce the benefit 
           of the parallelization. 
  
           A middle of the road approach is a good choice in these
           situations.  Here I decided to use a jobs size of 10 slots
           (-pe sm 10 ) per RepeatMasker invocation ( using 
           RepeatMasker -pa 9 ).  Why -pa 9 when I told SGE to use
           to slots ( 10 cpus in your case )?  This is because 
           RepeatMasker itself uses one slot and the other 9 
           are used by RepeatMasker to invoke 9 parallel runs of
           cross_match ( or whichever search engine you are using ).
           This reduces the ovearhead of too many RepeatMasker scripts
           running.  This also gives me between 2 and 4 batches 
           running per machine.  It also gives someone a chance
           to hold jobs and let another user jump ahead in a reasonable
           amount of runtime.
      
           I further decided to minimize the batch size to about
           the size that would fill the cluster ( should it be free ).
           240 slots / 10 = 24.

     A typical run works like this:
             - Break the genome up into 24 non-overlapping batches.
             - Submit 24 batches to the queue with the following
                 parameters:
                     o Use either ray512cc or raycc queue to
                       run ( -q raycc,ray512cc )
                     o Use a shared memory environment and request
                       10 slots per batch ( -pe sm 10 )
                     o Each batch runs RepeatMasker using -pa 9
               - When jobs are finished and no errors are 
                 found in the log files, the individual batch
                 runs are lifted up ( coordinates transformed back
                 from batch coordinates to global genomic coordinates )
                 and combined into a single result set.  Also I 
                 typically run some optional higher level analsysis
                 scripts on the results ( such as summarizing the 
                 *.out file, building repeat landscapes etc ).

      Batch Mechanism:

           The design of the batch file organisation is from
           UCSC with modifications.  
               - Gather information about the assembly ( size, 
                 sequences ) from the *.2bit file using the
                 UCSC utility "twoBitInfo".
               - Create subdirectories for each batch using 3
                 digit numbers ( zero padded ) starting from 
                 "000" and containing several files:
                            
	                  000.fa  : Sequence for the batch 
                                    pre-extracted from the .2bit
                                    file.
                          000.lft : A UCSC lift file which describes
                                    the genomic coordinates of the
                                    sequences in the batch file.
                          000.seqlst : A file which contains just the
                                       names of the sequences in the
                                       batch.
                          batch-000.sh : The SGE script to run
                                         RM on the batch


               - Once runs are complete in all the batch directories,
                 the results need to be adjusted so they can be combined.
                      * In each batch directory run the UCSC
                        liftUp program, specifying the file type ( *.out
                        or *.align ) and the *.lft file and the unlifted
                        RM output file.  
                      * One final run of liftUp is run to combine all the
                        lifted results from the previous step.
                - Finally the RepeatMasker buildSummary program is run
                  on the *.out file to create a summary of the results
                  ( similar to the *.tbl file ).
                        
                 
The options are:

=over 4

=item -version

Displays the version of the program

=back

=head1 DEPENDENCIES

RepeatMasker, UCSC: liftUp, twoBitInfo, twoBitToFa

=head1 SEE ALSO

Based initially on the simplePartition.pl UCSC script for breaking up
jobs for their cluster.  This version has a much improved batching 
mechanism and rolled in lots of script generation for the job run.
 
'''


## Set up input arguments
def get_args():
	parser = argparse.ArgumentParser(description="Given a fasta file and a batch count, create a directory structure with subdirectories for each batch and a set of scripts for use with TTU's Univa Grid Engine Cluster.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-g', '--genome_fasta', type=str, help='Name of the fasta formatted genome to be queried.', required=True)
	parser.add_argument('-b', '--batch_count', type=int, help="Integer number of batches into which to split the genome and run RepeatMasker.", required = True)
	parser.add_argument('-s', '--species', type=str, help='If using a Repbase derived library. Must use this if a custom library is not used.')
	parser.add_argument('-l', '--library', type=str, help='Custom repeat library in fasta format with #Class/Family designations. Must use this if a species library is not used.')
	parser.add_argument('-x', '--xsmall', type=str, help='Mask to lowercase letters. If not invoked, the default is to mask to Ns. Default = n', default='n')
	parser.add_argument('-n', '--nlow', type=str, help='Mask low-complexity regions? Default = no.', default = 'n')
	parser.add_argument('-i', '--inv', type=str, help='Return alignments in the orientation of the repeats. Default is y.', default = 'y')
	parser.add_argument('-p', '--path_to_genome', type=str, help='Path to the genome fasta? Default is the current working directory.', default = '.')
	parser.add_argument('-R', '--RMPart', type=str, help='If rerunning, delete the RMPart directory? Default = yes.', default = 'y')
	parser.add_argument('-q', '--qsub', type=str, help='Automatically submit all of the batch jobs? Default = yes.', default = 'y')
	parser.add_argument('-f', '--fast', type=str, help='Run in rapid mode? Default = no, run in sensitive mode.', default = 'n')
	parser.add_argument(-'qu', '--queue', type=str, help="Run on 'quanah' or 'hrothgar'? Default = quanah.", default = 'quanah')
	parser.add_argument("-log", "--log_level", default="INFO")

	args = parser.parse_args()
	GENOMEFA = args.genome_fasta
	BATCH = args.batch_count
	SPECIES = args.species
	LIBRARY = args.library
	XSMALL = args.xsmall
	NOLOW = args.nolow
	INV = args.inv
	GENOMEPATH = args.path_to_genome
	QSUB = args.qsub
	FAST = args.fast
	QUEUE = args.queue
	LOG = args.log_level

	return GENOMEFA, BATCH, SPECIES. LIBRARY, XSMALL, NOLOW, INV, GENOMEPATH, QSUB, FAST, QUEUE, LOG

##FUNCTIONS

##MAIN
def main():	
##Get input arguments
	GENOMEFA, BATCH, SPECIES. LIBRARY, XSMALL, NOLOW, INV, GENOMEPATH, QSUB, FAST, QUEUE, LOG = get_args()


#
# Wrap script functionality in main() to avoid automatic execution
# when imported ( e.g. when help is called on file )
#
if __name__ =="__main__":main()
		

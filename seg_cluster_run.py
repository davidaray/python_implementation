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



# Where RepeatMasker is stored
REPEATMASKER = "/lustre/work/daray/software/RepeatMasker"
# Where this script can find liftUp, twoBitInfo and twoBitToFa
BIN_DIR = "/lustre/work/daray/software"

# Define arguments
def get_args():
	#What this script does
	parser = argparse.ArgumentParser(description="Generate SGE cluster runs for RepeatMasker; built in RepeatMasker parameters are -xsmall [softmasks repetitive regions] -a [.align output file] -gff [generates a GFF format output] -pa [runs in parallel], please see RepeatMasker for details of these run options", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	required = parser.add_argument_group('required arguments')
	#Give input genome FASTA
	parser.add_argument('-i', '--input', type=str, help='genome file in FASTA format', required=True)
	#Argument of species name
	parser.add_argument('-sp', '--species', type=str, help='Source species of query DNA FASTA', required=False)
	# Desired batch number
	parser.add_argument('-b', '--batch_count', type=int, help='Batch count', default=50)
	# Input genome directory
	parser.add_argument('-dir', '--genome_dir', type=str, help='Path to genome FASTA', required=True)
	# Argument for output directory
	parser.add_argument('-od', '--outdir', type=str, help='Location of directory for the output subdirectory', default='.')
	# Which queue to use
	parser.add_argument('-q', '--queue', type=str, help='Select the queue to run RepeatMasker in [quanah|hrothgar] with the quanah option being the general quanah omni queue, and hrothgar being the communitycluster Chewie queue', choices=['quanah', 'hrothgar'], default='quanah')
	#Argument of RepeatMasker run parameter
	parser.add_argument('-lib', type=str, help='RepeatMasker run parameter custom library "-lib [filename]" option', required=False)
	#Argument of RepeatMasker run parameter
	parser.add_argument('-xsmall', type=str, help='Select a RepeatMasker masking option as lowercase bases [-xsmall], default is to mask as Ns', action='store_true')
	#Argument of RepeatMasker run parameter
	parser.add_argument('-engine', type=str, help='RepeatMasker run parameter "-engine <search_engine>" option; select a non-default search engine to use, otherwise RepeatMasker will used the default configured at install time; [crossmatch|abblast|rmblast|hmmer]', choices=['crossmatch', 'abblast', 'rmblast', 'hmmer'], required=False)
	#Argument of RepeatMasker run parameter
	parser.add_argument('-inv', type=str, help='RepeatMasker parameter flag "-inv" option; alignments are presented in the orientation of the repeat', action='store_true')
	#Argument of RepeatMasker run parameter
	parser.add_argument('-nolow', type=str, help='RepeatMasker parameter flag "-nolow" option; does not mask low complexity DNA or simple repeats', action='store_true')
	#Argument of RepeatMasker run parameter
	parser.add_argument('-s', '-speed', type=str, help='RepeatMasker run parameter "-q" or "-s" option; q=quick search; 5-10% less sensitive, 3-4 times faster than default; s=slow search; 0-5% more sensitive, 2.5 times slower than default', choices=['q', 's'], required=False)
	#Argument of RepeatMasker run parameter
	parser.add_argument('-div', type=int, help='RepeatMasker run parameter "-div [number]" option; masks only those repeats that are less than [number] percent diverged from the consensus sequence', required=False)
	
	args = parser.parse_args()
	GENOME = args.input
	SPECIES = args.species
	BATCH_COUNT = args.batch_count
	GENOME_DIR = args.genome_dir
	OUTDIR = args.outdir
	QUEUE = args.queue
	LIBRARY = args.lib
	XSMALL = args.xsmall
	ENGINE = args.engine
	INV = args.inv
	NOLOW = args.nolow
	SPEED = args.speed
	DIV = args.div
	
	return GENOME, SPECIES, BATCH_COUNT, GENOME_DIR, OUTDIR, QUEUE, LIBRARY, XSMALL, ENGINE, INV, NOLOW, SPEED, DIV
	
GENOME, SPECIES, BATCH_COUNT, GENOME_DIR, OUTDIR, QUEUE, LIBRARY, XSMALL, ENGINE, INV, NOLOW, SPEED, DIV = get_args()

# Sanity checks
print("The species is {}, the query genome is {}.\n").format(SPECIES, GENOME)
print("{} batches will be made.\n").format(str(BATCH_COUNT))
print("The genome FASTA is located in '{}'.\n").format(GENOME_DIR)
print("The output directory is '{}'.\n").format(OUTDIR)
print("The job queue is {}.\n").format(QUEUE)

if not SPECIES or LIBRARY:
	sys.exit("Must supply value for option 'species' or 'lib'!")
if SPECIES and LIBRARY:
	sys.exit("Only supply a value for one option: 'species' or 'lib'! Not both!")

FLAGS = [LIBRARY, XSMALL, ENGINE, INV, NOLOW, SPEED, DIV]
if not FLAGS:
	print("All default RepeatMasker parameters were used, no custom library.")
else:
	print("Custom parameters used:\n")
	if XSMALL:
		print("-xsmall flag used.\n")
	if INV:
		print("-inv flag used.\n")
	if NOLOW:
		print("-nolow flag used.\n")
	if LIBRARY:
		print("-lib flag used. Custom library is '{}'.\n").format(os.path.basename(LIBRARY))
	if ENGINE:
		print("-engine flag used. Changed search engine to {}.\n").format(ENGINE)
	if SPEED:
		print("-{} flag used. Search sensitivity has changed.\n").format(SPEED)
	if DIV:
		print("-div flag used. RepeatMasker will mask only repeats that are less than {}% diverged from the consensus sequence.\n").format(str(DIV))


if not os.path.isdir(GENOME_DIR):
	sys.exit("The given genome directory, '{}', does not exist.").format(GENOME_DIR)

GENOME_FASTA = os.path.join(GENOME_DIR, GENOME)

#if not os.path.isfile(GENOME_FASTA):
#	sys.exit("The given genome file '{}' does not exist.").format(GENOME_FASTA)
#if os.stat(GENOME).st_size==0:
#	sys.exit("The genome file, '{}', is empty.").format(GENOME_FASTA)
#if not os.path.isfile(LIBRARY):
#	sys.exit("The given library file '{}' does not exist.").format(LIBRARY)
#if os.stat(LIBRARY).st_size==0:
#	sys.exit("The library file, '{}', is empty.").format(LIBRARY)

try:
	if not os.path.getsize(GENOME_FASTA) > 0:
		sys.exit("The genome file, '{}', is empty.").format(GENOME_FASTA)
except OSError as e:
	sys.exit("The genome file '{}' does not exist or is inaccessible.").format(GENOME_FASTA)
	
try:
	if not os.path.getsize(LIBRARY) > 0:
		sys.exit("The library file, '{}', is empty.").format(LIBRARY)
except OSError as e:
	sys.exit("The library file '{}' does not exist or is inaccessible.").format(LIBRARY)
	
if not os.path.isdir(OUTDIR):
	sys.exit("The output directory '{}' does not exist.").format(OUTDIR)

PARTITION_DIR = os.path.join(GENOME_DIR, "RMPart")

SLOTS_PER_BATCH = 10
MAX_DIR_SIZE = 1000
NUM_BATCHES = BATCH_COUNT

PARTITION_DIR = os.path.abspath(PARTITION_DIR)
if not os.listdir(PARTITION_DIR):
	print("{} is empty. Continuing.").format(PARTITION_DIR)
else:
	print("{} is not empty. Removing contents and continuing.").format(PARTITION_DIR)
	os.remove(os.path.join(PARTITION_DIR, '*'))

def get_batches():
	# Return a 3-level list(ref): partitions -> chunks -> chunk properties (scaffold + coordinates)
	PARTS = []
	
	## Index the genome 
	LOGGER.info('Indexing the genome')
	GENOMEIDX = Fasta(GENOME)
	GENOMEPREFIX = os.path.splitext(GENOME)[0]
	FAIDX = pd.read_csv(GENOME	+ '.fai', sep='\t', names=['CONTIGHEADER', 'CONTIGSIZE', 'three', 'four', 'five'])
	FAIDX = FAIDX[['CONTIGHEADER', 'CONTIGSIZE']]
	FAIDX.to_csv(GENOMEPREFIX + '.fai', sep='\t', header=False, index=False)
	
	##Calculate total size of all contigs
	TOTALSIZE = FAIDX['CONTIGSIZE'].sum()
	## Get all of the headers into a list
	HEADERLIST = FAIDX['CONTIGHEADER'].tolist()
	## Remove index column and use all the contig names as the index
	NEWFAIDX = FAIDX.set_index('CONTIGHEADER')
	NEWFAIDX = NEWFAIDX.rename_axis(None)
		
	##Definitions --> BATCH = 000, 001, etc.
	##Definitions --> CHUNK = the peice of a contig to be added to a 000.fa, 001.fa, etc.
	## If batch_count is invoked...
	if BATCH_COUNT > 0:
		##Calculate the size of each chunk to be used in each batch
		TARGET_TOTAL_CHUNK_SIZE_PER_BATCH = int(TOTALSIZE / NUM_BATCHES) + 1
		## For each batch...
		for CURRENT_BATCH in range(0, NUM_BATCHES):
			## Set THIS_TOTAL_CHUNK_SIZE chunk size to zero
			THIS_TOTAL_CHUNK_SIZE = 0
			CURRENT_BATCH_HEADERLIST = []
			## Start putting together your list of headers, indicating how much of each contig you will use.
			## So, while THIS_TOTAL_CHUNK_SIZE (sum of all CHUNKS being added) is less than TARGET_TOTAL_CHUNK_SIZE_PER_BATCH
			CURRENT_CONTIG_SIZE = FAIDX.at[CONTIG_NAME, 'CONTIGSIZE']
			if CURRENT_CONTIG_SIZE + THIS_TOTAL_CHUNK_SIZE < TARGET_TOTAL_CHUNK_SIZE_PER_BATCH:
				CURRENT_HEADER_TO_ADD = HEADERLIST[0] + ':1 - ' + CURRENT_CONTIG_SIZE
				CURRENT_BATCH_HEADERLIST.append(CURRENT_HEADER_TO_ADD)
			
			while THIS_TOTAL_CHUNK_SIZE < TARGET_TOTAL_CHUNK_SIZE_PER_BATCH:
				AVAILABLE_SPACE = TARGET_TOTAL_CHUNK_SIZE_PER_BATCH - THIS_TOTAL_CHUNK_SIZE
				## Get the contig name for the first contig
				CONTIG_NAME = HEADERLIST[0]
				## Get the size of that contig 
				CONTIG_SIZE = FAIDX.at[CONTIG_NAME, 'CONTIGSIZE']
				## Check if this contig is larger than what remains available in THIS_TOTAL_CHUNK_SIZE
				## Add CONTIG_SIZE to THIS_TOTAL_CHUNK_SIZE
				THIS_TOTAL_CHUNK_SIZE += CONTIGSIZE
				##Remove the contig under consideration from HEADERLIST 
				del HEADERLIST[0]
				



	COMMAND = BIN + "/twoBitInfo " + GENOME_FASTA + " stdout |"
	#open(P, $COMMAND)
		#|| die "Couldn't open pipe ($COMMAND): $_\n"
	TOTAL_SIZE = 0
	SEQS = ()
	#while (<P>) {
	#	chomp;
	#	my ($seq, $seqsize) = split("\t");
	# ....
	
	if NUM_BATCHES > 0:
		CHUNK_SIZE = int(TOTAL_SIZE / NUM_BATCHES) + 1
		
	BATCHES = ()
	CURRENT_BATCH_SIZE = 0
	for SEQ in SEQS:
		SEQ_SIZE = SEQ_SIZES{SEQ}
		SEQ_IDX = 0
		while SEQ_SIZE > 0:
		
		
##MAIN
def main():	
##Get input arguments
	GENOMEFA, BATCH, SPECIES. LIBRARY, XSMALL, NOLOW, INV, GENOMEPATH, QSUB, FAST, QUEUE, LOG = get_args()


#
# Wrap script functionality in main() to avoid automatic execution
# when imported ( e.g. when help is called on file )
#
if __name__ =="__main__":main()
		
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


		

		

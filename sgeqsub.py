def generate_rmsubs(BATCH_COUNT, QUEUE, PREFIX, LIBRARY, SPECIES, XSMALL, NOLOW, SPEED, PROC):
	RMPATH = '/lustre/work/daray/software/RepeatMasker/'
	for BATCH in range(BATCH_COUNT):
		BATCHNUMBER = '{:03}'.format(BATCH)
		QSUBFILENAME = 'batch-'BATCH'.sh'
		if QUEUE = 'hrothgar':
			PROJECT = 'communitycluster'
			KUE = 'Chewie'
		elif QUEUE = 'quanah':
			PROJECT = 'quanah'
			KUE = 'omni'
		else:
			sys.exit('Please choose hrothgar or quanah as the queue to use.')
		with open("RMPart/" + BATCH + "/" + QSUB1FILENAME, 'w') as THISFILE:
			THISFILE.write('#!/bin/sh' + '\n')
			THISFILE.write('#$ -V' + '\n')
			THISFILE.write('#$ -cwd' + '\n')
			THISFILE.write('#$ -S /bin/bash' + '\n')
			THISFILE.write('#$ -N ' + PREFIX + '-' + BATCH + '\n')
			THISFILE.write('#$ -o $JOB_NAME.o$JOB_ID' + '\n')
			THISFILE.write('#$ -e $JOB_NAME.e$JOB_ID' + '\n')
			THISFILE.write('#$ -q ' + KUE + '\n')
			THISFILE.write('#$ -pe sm ' + str(PROC) + '\n')
			THISFILE.write('#$ -P ' + PROJECT + '\n')
			THISFILE.write('\n')
			THISFILE.write('cd RMPart/' + BATCH)
			FLAGS = [LIBRARY, SPECIES, XSMALL, NOLOW, SPEED]
			#Default RepeatMasker if no flags are raised.
			if not FLAGS:
				THISFILE.write(RMPATH + 'RepeatMasker -species ' + SPECIES + ' -a -inv -gff -s -pa ' + PROC + '. ' + BATCH + '.fa >& run.log\n'
			#All other possibilities. 
			else LIBRARY:
				THISFILE.write(RMPATH + 'RepeatMasker -lib ' + LIBRARY + ' -a -inv -gff -s -pa ' + PROC + '. ' + BATCH + '.fa >& run.log\n'
			else LIBRARY and XSMALL:
				THISFILE.write(RMPATH + 'RepeatMasker -lib ' + LIBRARY + ' -xsmall -a -inv -gff -s -pa ' + PROC + '. ' + BATCH + '.fa >& run.log\n'
			else SPECIES and XSMALL:
				THISFILE.write(RMPATH + 'RepeatMasker -species ' + SPECIES + ' -xsmall -a -inv -gff -s -pa ' + PROC + '. ' + BATCH + '.fa >& run.log\n'
			else LIBRARY and NOLOW:
				THISFILE.write(RMPATH + 'RepeatMasker -lib ' + LIBRARY + ' -nolow -a -inv -gff -s -pa ' + PROC + '. ' + BATCH + '.fa >& run.log\n'
			else SPECIES and NOLOW:
				THISFILE.write(RMPATH + 'RepeatMasker -species ' + SPECIES + ' -nolow -a -inv -gff -s -pa ' + PROC + '. ' + BATCH + '.fa >& run.log\n'
			else LIBRARY and NOLOW and XSMALL:
				THISFILE.write(RMPATH + 'RepeatMasker -lib ' + LIBRARY + ' -nolow -xsmall  -a -inv -gff -s -pa ' + PROC + '. ' + BATCH + '.fa >& run.log\n'
			else SPECIES and NOLOW and XSMALL:
				THISFILE.write(RMPATH + 'RepeatMasker -species ' + SPECIES + ' -nolow -xsmall -a -inv -gff -s -pa ' + PROC + '. ' + BATCH + '.fa >& run.log\n'
			else LIBRARY and SPEED:
				THISFILE.write(RMPATH + 'RepeatMasker -lib ' + LIBRARY + ' -s -a -inv -gff -s -pa ' + PROC + '. ' + BATCH + '.fa >& run.log\n'
			else SPECIES and SPEED:
				THISFILE.write(RMPATH + 'RepeatMasker -species ' + SPECIES + ' -s -a -inv -gff -s -pa ' + PROC + '. ' + BATCH + '.fa >& run.log\n'
			else LIBRARY and SPEED and XSMALL :
				THISFILE.write(RMPATH + 'RepeatMasker -lib ' + LIBRARY + ' -s -xsmall -a -inv -gff -s -pa ' + PROC + '. ' + BATCH + '.fa >& run.log\n'
			else SPECIES and SPEED and XSMALL:
				THISFILE.write(RMPATH + 'RepeatMasker -species ' + SPECIES + ' -s -xsmall -a -inv -gff -s -pa ' + PROC + '. ' + BATCH + '.fa >& run.log\n'
			else LIBRARY and SPEED and NOLOW :
				THISFILE.write(RMPATH + 'RepeatMasker -lib ' + LIBRARY + ' -s -nolow -a -inv -gff -s -pa ' + PROC + '. ' + BATCH + '.fa >& run.log\n'
			else SPECIES and SPEED and NOLOW:
				THISFILE.write(RMPATH + 'RepeatMasker -species ' + SPECIES + ' -s -nolow -a -inv -gff -s -pa ' + PROC + '. ' + BATCH + '.fa >& run.log\n'
			else LIBRARY and SPEED and NOLOW and XSMALL:
				THISFILE.write(RMPATH + 'RepeatMasker -lib ' + LIBRARY + ' -s -nolow -xsmall -a -inv -gff -s -pa ' + PROC + '. ' + BATCH + '.fa >& run.log\n'
			else SPECIES and SPEED and NOLOW and XSMALL:
				THISFILE.write(RMPATH + 'RepeatMasker -species ' + SPECIES + ' -s -nolow -xsmall -a -inv -gff -s -pa ' + PROC + '. ' + BATCH + '.fa >& run.log\n'



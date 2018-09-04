#!/usr/bin/env python

#merge all biological replicates into one bam file
#works from within the same directory as the bam files as is

import glob
import subprocess

g = glob.glob("*Rep1*.bam")

for sample in g:
	print("{}\n".format(sample))

	sampleName = sample.split("Rep1")[0]
	extension = sample.split("Rep1")[-1]
	reps = glob.glob(sampleName+"*.bam")

	if len(reps) == 1:
		print("There is only one replicate of {}\n".format(sample))
		cp_cmd = "cp {} {}".format(sample, sampleName+"mergedReps"+extension)
		print("{}\n".format(cp_cmd)) 
		subprocess.check_output(cp_cmd, shell=True)	

	else:
		samtools_cmd = "samtools merge -h {header_bam} {out_bam}".format(header_bam = sample,
		out_bam = sampleName+"mergedReps"+extension).split()

		for i in reps:
			samtools_cmd.append(i)

		cmd = " ".join(samtools_cmd)
		print("{}\n".format(cmd))
		subprocess.check_output(cmd, shell=True)		



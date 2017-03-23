#! /usr/bin/env python

import sys, subprocess
import utils

def main(args):


    utils.get_dnv_singleton(args.bam_file, args.output_file + ".tmp.dnv_1.txt")
    utils.enlarge_bed(args.output_file + ".tmp.dnv_1.txt", args.output_file + ".tmp.dnv_1.bed")

    hout = open(args.output_file + ".tmp.dnv_1.pileup", 'w')
    subprocess.check_call(["samtools", "mpileup", "-f", args.reference_genome, args.bam_file, "-l", args.output_file + ".tmp.dnv_1.bed"], stdout = hout)
    hout.close()

    utils.remove_lowdepth_snp(args.output_file + ".tmp.dnv_1.txt", args.output_file + ".tmp.dnv_1.pileup", args.output_file)


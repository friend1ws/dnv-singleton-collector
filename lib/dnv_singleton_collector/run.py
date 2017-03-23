#! /usr/bin/env python

import sys, subprocess
import utils

def main(args):

    utils.get_dnv_singleton(args.bam_file, args.output_prefix + ".tmp.dnv_inv_1.txt")
    utils.enlarge_bed(args.output_prefix + ".tmp.dnv_inv_1.txt", args.output_prefix + ".tmp.dnv_inv_1.bed")

    hout = open(args.output_prefix + ".tmp.dnv_inv_1.pileup", 'w')
    subprocess.check_call(["samtools", "mpileup", "-f", args.reference_genome, args.bam_file, "-l", args.output_prefix + ".tmp.dnv_inv_1.bed"], stdout = hout)
    hout.close()

    utils.remove_lowdepth_snp(args.output_prefix + ".tmp.dnv_inv_1.txt", args.output_prefix + ".tmp.dnv_inv_1.pileup", args.output_prefix + ".dnv_inv_list.txt")

    utils.get_dnv_inv_profile(args.output_prefix + ".dnv_inv_list.txt", args.output_prefix + ".dnv_profile.txt", args.output_prefix + ".inv_profile.txt")


    subprocess.check_call(["rm", "-rf", args.output_prefix + ".tmp.dnv_inv_1.txt"])
    subprocess.check_call(["rm", "-rf", args.output_prefix + ".tmp.dnv_inv_1.bed"])
    subprocess.check_call(["rm", "-rf", args.output_prefix + ".tmp.dnv_inv_1.pileup"])



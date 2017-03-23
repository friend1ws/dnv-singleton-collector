#! /usr/bin/env python


def get_dnv_singleton(bam_file, output_file, min_mapping_qual = 30, min_base_qual = 30, max_mismatch_num = 5, 
                      min_edge_mean_base_qual = 20, edge_size_for_qual_check = 5, ignore_edge_size = 12):

    import re, pysam

    bamfile = pysam.Samfile(bam_file, "rb")
    hout = open(output_file, "w")
 
    MD_tag_DNV_re = re.compile(r'[ACGT][01][ACGT]')
    MD_tag_split_re = re.compile(r'[\^ACGT]+|\d+')

    dnv2count = {}
    # maybe add the regional extraction of bam files
    for read in bamfile.fetch():

        if read.qname == "B00GFABXX110215:3:2106:20301:28994":
            pass

        # skip if there is no SA tags
        NM_tag = None 
        MD_tag = None
        for item in read.tags:
            if item[0] == "NM": NM_tag = item[1]
            if item[0] == "MD": MD_tag = item[1]

        if NM_tag <= 1 or NM_tag > max_mismatch_num: continue
        if MD_tag_DNV_re.search(MD_tag) is None: continue
        if len(MD_tag_DNV_re.findall(MD_tag)) > 1: continue

        qual_str = read.qual
        if sum([float(ord(x) - 33) for x in qual_str[:edge_size_for_qual_check]]) / edge_size_for_qual_check <= min_edge_mean_base_qual: continue
        if sum([float(ord(x) - 33) for x in qual_str[-edge_size_for_qual_check:]]) / edge_size_for_qual_check <= min_edge_mean_base_qual: continue

        # skip if there is indel or soft clippings
        if len(read.cigar) > 1: continue

        # skip if below the minimum mapping quality
        if read.mapq < min_mapping_qual: continue

        # get the flag information
        flags = format(int(read.flag), "#014b")[:1:-1]

        # skip unproper pair
        if flags[1] == "0": continue

        # skip if either of the read pair is unmapped
        if flags[2] == "1" or flags[3] == "1": continue

        # skip supplementary alignment
        if flags[8] == "1" or flags[11] == "1": continue

        # skip duplicated reads
        if flags[10] == "1": continue


        genomic_pos = int(read.pos + 1)
        sequence_pos = 0 
        pos1, pos2 = 0, 0
        qual1, qual2 = 0, 0
        ref1, ref2, alt1, alt2 = '', '', '', ''
        prev_base1, prev_base2, next_base1, next_base2 = '', '', '', ''
        dnv_flag = False
        finish_flag = False

        for MD_elm in MD_tag_split_re.findall(MD_tag):
            if MD_elm.isdigit(): 
                genomic_pos = genomic_pos + int(MD_elm)
                sequence_pos = sequence_pos + int(MD_elm)
                dnv_flag = True if int(MD_elm) <= 1 else False
            else:
                if ref1 != '' and ref2 == '' and dnv_flag == True and finish_flag == False and len(read.seq) - sequence_pos - 1 >= ignore_edge_size:
                    pos2 = genomic_pos
                    ref2 = MD_elm
                    alt2 = read.seq[sequence_pos]
                    next_base1 = read.seq[sequence_pos + 1]
                    next_base2 = read.seq[sequence_pos + 2]
                    qual2 = str(ord(read.qual[sequence_pos]) - 33) 
                    finish_flag = True
                elif finish_flag == False and sequence_pos >= ignore_edge_size:
                    pos1 = genomic_pos
                    ref1 = MD_elm
                    alt1 = read.seq[sequence_pos]
                    prev_base1 = read.seq[sequence_pos - 1]
                    prev_base2 = read.seq[sequence_pos - 2]
                    qual1 = str(ord(read.qual[sequence_pos]) - 33)
                genomic_pos = genomic_pos + 1
                sequence_pos = sequence_pos + 1


        if int(qual1) < min_base_qual or int(qual2) < min_base_qual: continue
        if flags[4] == "0" and prev_base2 == prev_base1 == alt1 == alt2: continue
        if flags[4] == "1" and alt1 == alt2 == next_base1 == next_base2: continue 

        dnv = bamfile.getrname(read.tid) +  '\t' + str(pos1) + '\t' + str(pos2) + '\t' + ref1 + '\t' + ref2 + '\t' + alt1 + '\t' + alt2
        if dnv not in dnv2count: dnv2count[dnv] = 0

        dnv2count[dnv] = dnv2count[dnv] + 1

        del_dnv_list = []
        for tdnv in dnv2count:
            tchr, tpos1, tpos2, tref1, tref2, talt1, talt2 = tdnv.split('\t')
            if bamfile.getrname(read.tid) != tchr or pos1 - int(tpos1) >= 10000:
                if dnv2count[tdnv] == 1: print >> hout, tdnv
                del_dnv_list.append(tdnv)

        for tdnv in del_dnv_list:
            del dnv2count[tdnv]
 

    for tdnv in dnv2count:
        tchr, tpos1, tpos2, tref1, tref2, talt1, talt2 = tdnv.split('\t')
        if bamfile.getrname(read.tid) != tchr or pos1 - int(tpos1) >= 10000:
            if dnv2count[tdnv] == 1: print >> hout, tdnv
            del dnv2count[tdnv]


    hout.close()


def enlarge_bed(input_file, output_file):
    hout = open(output_file, 'w')
    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            print >> hout, F[0] + '\t' + str(int(F[1]) - 2) + '\t' + str(int(F[2]) + 1) + '\t' + '\t'.join(F[3:])

    hout.close()



def remove_lowdepth_snp(input_file, pileup_file, output_file, min_depth = 10, max_mismatch_size = 0.1):

    remove_region = {}
    tmp_chr = ""
    tmp_start = 0
    tmp_pos = 0
    tmp_flag = False
    with open(pileup_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')

            if F[0] != tmp_chr or int(F[1]) - tmp_pos > 1:
                if tmp_flag == True:
                    for i in range(tmp_start, tmp_pos + 1):
                        remove_region[tmp_chr + '\t' + str(i)] = 1

                tmp_chr = F[0]
                tmp_start = int(F[1])
                tmp_flag = False

            tmp_pos = int(F[1])
            if int(F[3]) < min_depth: 
                tmp_flag = True
            else:
                num_ref = F[4].count(".") + F[4].count(",")
                if (float(F[3]) - float(num_ref)) / float(F[3]) > max_mismatch_size: tmp_flag = True

        if tmp_flag == True:
            for i in range(tmp_start, tmp_pos + 1):
                remove_region[tmp_chr + '\t' + str(i)] = 1
 

    hout = open(output_file, 'w')
    with open(input_file, 'r') as hin:
        for line in hin:
            F = line.rstrip('\n').split('\t')
            if F[0] + '\t' + F[1] in remove_region: continue
            print >> hout, '\t'.join(F)

    hout.close()



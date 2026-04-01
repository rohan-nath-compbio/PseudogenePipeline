################################################################################
#
# 10/14/10 By Shinhan
# This script does two things:
# 1) Based on the pseudogene file from step3, get the coordinates of pseudogenes
#    on the subject sequences.
# 2) Generate a pair list with prot_id, genome_seq_coord_as_id (based on 1).
#
#

import sys, ast

with open(sys.argv[1]) as inp, \
     open(sys.argv[1] + ".subj_coord", "w") as oup1, \
     open(sys.argv[1] + ".pairs", "w") as oup2:
    for inl in inp:
        L = inl.split("\t")		# genome_seq_id, prot_id, genome_region, prot_region
        gr = ast.literal_eval(L[2])
        gr = [gr[0][0], gr[-1][1]]  # Get the regions out

        oup1.write("%s\t%s\t%s\n" % (L[0], gr[0], gr[1]))
        oup2.write("%s\t%s|%s-%s\n" % (L[1], L[0], gr[0], gr[1]))

#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (https://opensource.org/licenses/MIT).
#

from hts import Bam
from hts.htsffi import libhts 

class COAlignment: pass

class COPileupRead:
    def __init__(self, alignment, pileup_pos):
        self.alignment = COAlignment()
        self.is_del = False
        self.qpos = self.get_qpos(alignment, pileup_pos)
        self.alignment.tags = self.get_tags(alignment)
        self.alignment.qual = ''.join(chr(qual+33) for qual in alignment.base_qualities)
        self.alignment.seq = alignment.seq
        self.alignment.mapq = alignment.mapping_quality
        self.alignment.is_reverse = alignment.flag & 16 == 16
        self.alignment.is_proper_pair = alignment.flag & 2 == 2
        self.alignment.is_read1 = alignment.flag & 64 == 64
        self.alignment.qname = alignment.qname

    def get_qpos(self, alignment, pileup_pos):
        qpos = pileup_pos - alignment.pos
        cigarcnt = 0
        for cnt, cigar in self.get_cigarlist(alignment):
            if qpos < cigarcnt: break
            if cigar == "S" or cigar == "I":
                qpos += cnt
            elif cigar == "D":
                if cigarcnt + cnt > qpos:
                    qpos = cigarcnt
                    self.is_del = True
                    break
                else:
                    qpos -= cnt
            if cigar in "MIS=X":
                cigarcnt += cnt
        return qpos

    def get_tags(self, alignment):
        tags = []
        for tag, typ, val in alignment.tags:
            if typ == "i":
                tags.append( (tag, int(val), ) )
            elif typ == "f":
                tags.append( (tag, float(val), ) )
            else:
                tags.append( (tag, val) )
        return tags

    def get_cigarlist(self, alignment):
        cigarlist = []
        cig = alignment.cigar._c
        for k in range(alignment.cigar.n_cigar):
            cigarlist.append( (libhts.bam_cigar_oplen(cig[k]), libhts.bam_cigar_opchr(cig[k]), ) )
        return cigarlist

class COPileupColumn:
    def __init__(self, pileup_pos):
        self.pileups = None
        self.pos = pileup_pos

class AlignmentFile:
    def __init__(self, fn, mode):
        self.alignmentfile = Bam(fn, mode)

    def pileup(self, chrom, pileup_start, pileup_end):
        BAM_FUNMAP, BAM_FSECONDARY, BAM_FQCFAIL, BAM_FDUP = 4, 256, 512, 1024

        for pileup_pos in range(pileup_start, pileup_end):
            mypileupcolumn = COPileupColumn(pileup_pos)
            mypileups = []
            for alignment in self.alignmentfile(chrom + ":" + str(pileup_pos+1) + "-" + str(pileup_pos+1)):
                if alignment.flag & BAM_FUNMAP or alignment.flag & BAM_FSECONDARY or alignment.flag & BAM_FQCFAIL or alignment.flag & BAM_FDUP:
                    continue
                mypileup = COPileupRead(alignment, pileup_pos)
                mypileups.append(mypileup)
            mypileupcolumn.pileups = tuple(mypileups)
            yield mypileupcolumn

    def close(self):
        self.alignmentfile.close()

    def __del__(self):
        self.close()

class Samfile(AlignmentFile): pass
#!/usr/bin/env python
"""
Simple python wrapper for SSW library
Please put the path of libssw.so into LD_LIBRARY_PATH or pass it explicitly as a parameter
By Yongan Zhao (March 2016)

Modify the original python wrapper to split PacBio subreads bam file based on a custom smrtbell sequence.
Basically it is to replace the primary analysis.
By Bo Han(2016/07)
"""

import sys
import os.path as op
import argparse as ap
import ctypes as ct
import timeit as ti
import gzip
import math

import ssw_lib
import pysam


def read(sFile):
    """
    read a sequence file
    @param  sFile   sequence file
    """

    def read_one_fasta(f):
        """
        read a fasta file
        @param  f   file handler
        """
        sId = ''
        sSeq = ''
        for l in f:
            if l.startswith('>'):
                if sSeq:
                    yield sId, sSeq, ''
                sId = l.strip()[1:].split()[0]
                sSeq = ''
            else:
                sSeq += l.strip()

        yield sId, sSeq, ''

    def read_one_fastq(f):
        """
        read a fastq file
        @param  f   file handler
        """
        sId = ''
        sSeq = ''
        s3 = ''
        sQual = ''
        for l in f:
            sId = l.strip()[1:].split()[0]
            sSeq = f.next()
            s3 = f.next()
            sQual = f.next()

            yield sId, sSeq, sQual

            # test if fasta or fastq

    bFasta = True
    ext = op.splitext(sFile)[1][1:].strip().lower()
    if ext == 'gz' or ext == 'gzip':
        with gzip.open(sFile, 'r') as f:
            l = f.next()
            if l.startswith('>'):
                bFasta = True
            elif l.startswith('@'):
                bFasta = False
            else:
                print >> sys.stderr, 'file format cannot be recognized'
                sys.exit()
    else:
        with open(sFile, 'r') as f:
            l = f.next()
            if l.startswith('>'):
                bFasta = True
            elif l.startswith('@'):
                bFasta = False
            else:
                print >> sys.stderr, 'file format cannot be recognized'
                sys.exit()

                # read
    if ext == 'gz' or ext == 'gzip':
        with gzip.open(sFile, 'r') as f:
            if bFasta == True:
                for sId, sSeq, sQual in read_one_fasta(f):
                    yield sId, sSeq, sQual
            else:
                for sId, sSeq, sQual in read_one_fastq(f):
                    yield sId, sSeq, sQual
    else:
        with open(sFile, 'r') as f:
            if bFasta == True:
                for sId, sSeq, sQual in read_one_fasta(f):
                    yield sId, sSeq, sQual
            else:
                for sId, sSeq, sQual in read_one_fastq(f):
                    yield sId, sSeq, sQual


def to_int(seq, lEle, dEle2Int):
    """
    translate a sequence into numbers
    @param  seq   a sequence
    """
    num_decl = len(seq) * ct.c_int8
    num = num_decl()
    for i, ele in enumerate(seq):
        try:
            n = dEle2Int[ele]
        except KeyError:
            n = dEle2Int[lEle[-1]]
        finally:
            num[i] = n

    return num


def align_one(ssw, qProfile, rNum, nRLen, nOpen, nExt, nFlag, nMaskLen):
    """
    align one pair of sequences
    @param  qProfile   query profile
    @param  rNum   number array for reference
    @param  nRLen   length of reference sequence
    @param  nFlag   alignment flag
    @param  nMaskLen   mask length
    """
    res = ssw.ssw_align(qProfile, rNum, ct.c_int32(nRLen), nOpen, nExt, nFlag, 0, 0, nMaskLen)

    nScore = res.contents.nScore
    nScore2 = res.contents.nScore2
    nRefBeg = res.contents.nRefBeg
    nRefEnd = res.contents.nRefEnd
    nQryBeg = res.contents.nQryBeg
    nQryEnd = res.contents.nQryEnd
    nRefEnd2 = res.contents.nRefEnd2
    lCigar = [res.contents.sCigar[idx] for idx in range(res.contents.nCigarLen)]
    nCigarLen = res.contents.nCigarLen
    ssw.align_destroy(res)

    return (nScore, nScore2, nRefBeg, nRefEnd, nQryBeg, nQryEnd, nRefEnd2, nCigarLen, lCigar)

def split_bam(sam_record, run_name, zmw, left_start, right_end, res, flag):
    new_bam = pysam.AlignedSegment()
    new_tags = []
    if flag == 0x2: # left
        new_bam.query_name = '/'.join([run_name, zmw, '_'.join([str(left_start), str(left_start + res[2])])])
        new_bam.query_sequence = sam_record.query_sequence[0:res[2]]
    elif flag == 0x1:  # right
        new_bam.query_name = '/'.join([run_name, zmw, '_'.join([str(left_start + res[3] + 1), str(right_end)])])
        new_bam.query_sequence = sam_record.query_sequence[res[3] + 1:]
    else:
        print >> sys.argv, "unknown tag {}".format(tag)
        sys.exit(1)
    new_bam.flag = sam_record.flag
    new_bam.reference_id  = -1
    new_bam.reference_start  = -1
    new_bam.mapping_quality  = 255
    new_bam.next_reference_id  = -1
    new_bam.next_reference_start  = -1
    for tag in sam_record.tags:
        if tag[0] in ("RG", "np", "rq", "sn", "zm"):
            new_tags.append(tag)
        elif tag[0] == 'qs':
            if flag == 0x2:
                new_tags.append((tag[0], left_start))
            else:
                new_tags.append((tag[0], left_start + res[3] + 1))
        elif tag[0] == "qe":
            if flag == 0x2:
                new_tags.append((tag[0], left_start + res[2]))
            else:
                new_tags.append((tag[0], right_end))
        elif tag[0] == "cx":
            new_tags.append((tag[0], tag[1] | flag))
        elif tag[0] in ("dq", "dt", "ip", "iq", "mq", "sq"):
            if flag == 0x2:
                new_tags.append((tag[0], tag[1][0:res[2]]))
            else:
                new_tags.append((tag[0], tag[1][res[3] + 1:]))
        else:
            print >> sys.argv, "unknown tag {}".format(tag[0])
    new_bam.tags = tuple(new_tags)
    return new_bam

def main(args):
    lEle = []
    dRc = {}
    dEle2Int = {}
    dInt2Ele = {}
    # init DNA score matrix
    if not args.sMatrix:
        lEle = ['A', 'C', 'G', 'T', 'N']
        dRc = {'A': 'C', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 'C', 'c': 'G', 'g': 'C', 't': 'A'}
        for i, ele in enumerate(lEle):
            dEle2Int[ele] = i
            dEle2Int[ele.lower()] = i
            dInt2Ele[i] = ele
        nEleNum = len(lEle)
        lScore = [0 for i in xrange(nEleNum ** 2)]
        for i in xrange(nEleNum - 1):
            for j in xrange(nEleNum - 1):
                if lEle[i] == lEle[j]:
                    lScore[i * nEleNum + j] = args.nMatch
                else:
                    lScore[i * nEleNum + j] = -args.nMismatch
    else:
        lEle, dEle2Int, dInt2Ele, lScore = ssw.read_matrix(args.sMatrix)

    # translate score matrix to ctypes
    mat = (len(lScore) * ct.c_int8)()
    mat[:] = lScore
    nFlag = 1  # always return CIGAR
    # (from high to low)
    # 00000000
    #     | bit 5: return the best alignment beginning position;
    #      | 6: if (ref_end1 - ref_begin1 <= filterd) && (read_end1 - read_begin1 <= filterd), return cigar;
    #       | 7: if max score >= filters, return cigar;
    #        | 8: always return cigar; if 6 & 7 are both setted, only return cigar when both filter fulfilled
    # if ((7&flag) == 0 || ((2&flag) != 0 && r->score1 < filters) || ((4&flag) != 0 && (r->ref_end1 - r->ref_begin1 > filterd || r->read_end1 - r->read_begin1 > filterd))) goto end;

    ssw = ssw_lib.CSsw(args.sLibPath)
    # iterate query sequence
    for sQId, sQSeq, sQQual in read(args.query):
        # build query profile
        qNum = to_int(sQSeq, lEle, dEle2Int)
        qProfile = ssw.ssw_init(qNum, ct.c_int32(len(sQSeq)), mat, len(lEle), 2)
        # set mask len
        if len(sQSeq) > 30:
            nMaskLen = len(sQSeq) / 2
        else:
            nMaskLen = 15

        # iter target sequence
        # for sRId, sRSeq, _ in read(args.target):
        in_fh = pysam.AlignmentFile(args.target, "rb")
        # out_fh = pysam.AlignmentFile("/home/hanb/PycharmProjects/untitled1/out.bam", "wb", header=in_fh.header)
        out_fh = pysam.AlignmentFile("-", "wb", header=in_fh.header)
        for sam_record in in_fh:
            # process sam info
            run_name, zmw, tmp = sam_record.query_name.split('/')
            left_start, right_end = map(int, tmp.split('_'))

            # prepare for alignment
            sRSeq = sam_record.query_sequence
            rNum = to_int(sRSeq, lEle, dEle2Int)
            res = align_one(ssw, qProfile, rNum, len(sRSeq), args.nOpen, args.nExt, nFlag, nMaskLen)

            # don't consider RC alignment
            # resRc = align_one(ssw, qRcProfile, rNum, len(sRSeq), args.nOpen, args.nExt, nFlag, nMaskLen)

            # optimal alignment score: res[0]
            # sub optimal alignment score: res[1]
            # target (query sequence) begin: res[2]
            # target end: res[3]
            # query (primer) begin: res[4]
            # query end: res[5]
            #                 0-based     0-based
            # left_start      res[2]      res[3]                    right_end
            # V               V           V                         V
            # ------------------------------------------------------
            #                 ------------

            # TODO: unaligned segment?
            if res[2] > 0:
                out_fh.write(split_bam(sam_record, run_name, zmw, left_start, right_end, res, 0x2))
            if left_start + res[3] + 1 < right_end:
                out_fh.write(split_bam(sam_record, run_name, zmw, left_start, right_end, res, 0x1))

        ssw.init_destroy(qProfile)


if __name__ == '__main__':
    parser = ap.ArgumentParser()
    parser.add_argument('-l', '--sLibPath', default='', help='path of libssw.so')
    parser.add_argument('-m', '--nMatch', type=int, default=2,
                        help='a positive integer as the score for a match in genome sequence alignment. [default: 2]')
    parser.add_argument('-x', '--nMismatch', type=int, default=2,
                        help='a positive integer as the score for a mismatch in genome sequence alignment. [default: 2]')
    parser.add_argument('-o', '--nOpen', type=int, default=3,
                        help='a positive integer as the penalty for the gap opening in genome sequence alignment. [default: 3]')
    parser.add_argument('-e', '--nExt', type=int, default=1,
                        help='a positive integer as the penalty for the gap extension in genome sequence alignment. [default: 1]')
    parser.add_argument('-a', '--sMatrix', default='',
                        help='a file for either Blosum or Pam weight matrix. [default: Blosum50]')
    parser.add_argument('-f', '--nThr', default=0,
                        help='a positive integer. Only output the alignments with the Smith-Waterman score >= N.')
    parser.add_argument('target', help='targe file')
    parser.add_argument('query', help='query file')
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    args = parser.parse_args()

    t1 = ti.default_timer()
    main(args)
    t2 = ti.default_timer()
    print >> sys.stderr, 'CPU time: {} seconds'.format(t2 - t1)

#!/usr/bin/env python
'''
this script takes two bam files
first bam is converted directly from bax.h5 files using bax2bam
second bam is converted from "re-farmed" bas.h5 file using bax2bam
because only one end of the smrtbells uses the canonical smrtbell sequence, so
sequences in the first bam contain the second smrtbell sequence

    ATCTCTCTCAATTTTTTTTTTTTTTTTTTTTTTTAAGAGAGAGAT

the second bam is generated from primary analysis that is aware of the
 heterogeneic smrtbells

this script then compare the sequences and see what was removed in the second bam,
providing information on the ligated sequences

'''

import sys, os
import pysam


class Interval(object):
    '''
    a simple class to represent a containuous segment in the genome
    chromosome name is not used in this specific application
    '''
    __slots__ = ["start", "end"]

    def __init__(self, s=0, e=0):
        '''
        :param s: start of the interval, int
        :param e: end of the interval, int
        '''
        if e < s:
            s, e = e, s
        self.start = s
        self.end = e

    @staticmethod
    def from_generator(gen):
        try:
            s = gen.next()
            e = gen.next()
        except StopIteration:
            s = 0
            e = 1
        return Interval(s, e)

    @staticmethod
    def from_list(lst):
        try:
            s = lst[0]
            e = lst[1]
        except OverflowError:
            s = 0
            e = 1
        return Interval(s, e)

    def __lt__(self, other):
        if self.start == other.start:
            return self.end < other.end
        else:
            return self.start < other.start

    def __gt__(self, other):
        if self.start == other.start:
            return self.end > other.end
        else:
            return self.start > other.start

    def __eq__(self, other):
        return self.start == other.start and self.end == other.end

    def is_valid(self):
        return self.start >= 0 and self.end > 0 and self.start < self.end

    def overlap(self, other):
        '''
        :param other: another Interval
        :return:  bool
        '''
        return min(self.end, other.end) > max(self.start, other.start)

    def substract(self, to_remove):
        '''
        :param to_remove: another Interval
        :return: generator of Intervals with to_remove deleted
        '''
        if not self.overlap(to_remove):
            yield Interval(self.start, self.end)
        else:
            if self.start < to_remove.start:
                yield Interval(self.start, to_remove.start)
            if self.end > to_remove.end:
                yield Interval(to_remove.end, self.end)


class Segment(object):
    '''
    Represent a polymerase read, which include multiple read of insert
    some of the ROI contains to adaptor to be removed
    '''
    __slots__ = ["run_name", "zmw", "intervals", "seqs"]

    def __init__(self, sam_query):
        self.run_name, zmw, tmp = sam_query.qname.split('/')
        self.zmw = int(zmw)
        self.intervals = []
        self.seqs = []
        self.intervals.append(Interval.from_list(map(int, tmp.split('_'))))
        self.seqs.append(sam_query.query_sequence)

    def insert(self, new_seg):
        new_interval = new_seg.intervals[0]
        i = 0
        while i < len(self.intervals):
            if new_interval < self.intervals[i]:
                break
            i += 1
        self.intervals.insert(i, new_interval)
        self.seqs.insert(i, new_seg.seqs[0])

    def substract(self, r):
        l = []
        s = []
        i = 0
        while i < len(self.intervals):
            interval = self.intervals[i]
            if interval.overlap(r):
                for frag in interval.substract(r):
                    l.append(frag)
                    s.append(self.seqs[i][frag.start - interval.start: frag.end - interval.start])
            else:
                l.append(interval)
                s.append(self.seqs[i])
            i += 1
        self.intervals = l  # shouldn't need to do sorting again
        self.seqs = s

    def substr(self):
        for seq in self.seqs:
            yield seq


def load_original_sam(sam_file):
    h = {}
    with pysam.AlignmentFile(sam_file, 'r') as f:
        for sam_entry in f:
            seg = Segment(sam_entry)
            if seg.zmw in h:
                h[seg.zmw].insert(seg)
            else:
                h[seg.zmw] = seg
    return h


def load_processed_sam(sam_file, hash):
    for sam_entry in pysam.AlignmentFile(sam_file, 'r'):
        seg = Segment(sam_entry)
        if seg.zmw in hash:
            hash[seg.zmw].substract(seg.intervals[0])
        else:
            print >> sys.stderr, "meet with new zmw {}".format(seg.zmw)


def output(data):
    for key, seg in data.items():
        for i in seg.substr():
            print i


def main(before_sam, after_sam):
    data = load_original_sam(before_sam)
    load_processed_sam(after_sam, data)
    output(data)


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print "usage: {}  default.subreads.bam  refarmed.subreads.bam".format(sys.argv[0])
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])

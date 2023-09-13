#!/usr/bin/env python

'''
Organizational structure for tracking individual bases from reads that cover the same genomic
position.
'''

class CompiledReads:
    '''
    Manager for CompiledPositions
    '''
    _strands = ('-', '+', '*')
    def __init__(self,
                 strand=0,
                 min_base_position=0,
                 max_base_position=float("inf"),
                 min_base_quality=0):
        self._nucleotides = {}
        if strand == 0:
            self.get_strand = lambda x: 2
        else:
            self.get_strand = self._get_strand
        self._strand_one = strand == 1
        self._ref = None
        self._ref_seq = self._get_ref_from_read

        self._qc = {
            "min_base_quality": min_base_quality,
            "min_base_position": min_base_position,
            "max_base_position": max_base_position
        }

    def add_reference(self, ref):
        '''
        Add a reference FASTA file to use.

        Parameters:
            ref (reditools.FastaFile)
        '''
        self._ref = ref
        self._ref_seq = self._get_ref_from_fasta

    def _get_ref_from_read(self, read):
        return list(read.get_reference_sequence().upper())

    def _get_ref_from_fasta(self, read):
        indices = [i for _, i in read.get_aligned_pairs(matches_only=True)]
        return self._ref.get_base(read.reference_name, *indices)

    def _prep_read(self, read):
        aligned_pairs = read.get_aligned_pairs(matches_only=True)
        seq = read.query_sequence.upper()
        qualities = read.query_qualities
        ref_seq = self._ref_seq(read)
        while aligned_pairs and aligned_pairs[0][0] < self._qc["min_base_position"]:
            aligned_pairs.pop(0)
            ref_seq.pop(0)
        if not aligned_pairs:
            return

        while aligned_pairs and \
                read.query_length - aligned_pairs[0][0] >= self._qc["max_base_position"]:
            i, j = aligned_pairs.pop(0)
            ref = ref_seq.pop(0)
            if qualities[i] < self._qc["min_base_quality"] or "N" in (ref, seq[i]):
                continue
            yield (j, seq[i], qualities[i], ref)

    def _get_strand(self, read):
        return read.is_read2 ^ self._strand_one ^ read.is_reverse

    def add_reads(self, reads):
        '''
        Adds an iterable of pysam reads to the object. The reads are broken down
        into individual nucleotides that are tracked by chromosomal location.

        Parameters:
            reads (iterable) - pysam reads
        '''

        for read in reads:
            strand = CompiledReads._strands[self.get_strand(read)]
            for pos, base, quality, ref in self._prep_read(read):
                try:
                    self._nucleotides[pos].add_base(quality, strand, base)
                except KeyError:
                    self._nucleotides[pos] = CompiledPosition(ref=ref)
                    self._nucleotides[pos].add_base(quality, strand, base)

    def pop(self, position):
        '''
        Removes and returns the CompiledPosition at position. Method returns None
        if the position is empty.

        Parameters:
            position (int) - The chromosomal location to pop

        Returns:
            A CompiledPosition or None if position is empty.
        '''
        return self._nucleotides.pop(position, None)

    def is_empty(self):
        '''
        Determines if there are any CompiledPositions still in the object.

        Returns:
            True if the object is empty, else False
        '''
        return len(self._nucleotides) == 0

class CompiledPosition:
    '''
    Keeps an ongoing list of details for bases and reads at a given position in the
    reference genome.
    '''
    _complement_map = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}

    def __init__(self, ref):
        self.qualities = []
        self.strands = []
        self.bases = []
        self.counter = False
        self.ref = ref

    def __len__(self):
        return len(self.bases)

    def add_base(self, quality, strand, base):
        '''
        Add details for a base at this position.

        Parameters:
            quality (int): The quality of the read
            strand (str): The strand the base is on (+, -, or *)
            base (str): The nucleotide at the position( A, C, G, or T)
        '''
        self.qualities.append(quality)
        self.strands.append(strand)
        self.bases.append(base)
        self.counter = False

    def filter_by_strand(self, vstrand):
        '''
        Filters the position content to only report bases on a specified strand.

        Parameters:
            vstrand (str): The strand to fitler by (+ or -)
        '''
        indices = [i for i in range(len(self.strands)) if self.strands[i] == vstrand]
        self._filter(indices)

    def filter_by_quality(self, min_quality):
        '''
        Filters the position content to only report bases with a minimum quality.

        Paramters:
            min_quality (int)
        '''
        indices = [i for i in range(len(self.qualities)) if self.qualities[i] >= min_quality]
        self._filter(indices)

    def _filter(self, indices):
        self.bases = [self.bases[i] for i in indices]
        self.strands = [self.strands[i] for i in indices]
        self.qualities = [self.qualities[i] for i in indices]
        self.counter = False

    def complement(self):
        '''
        Modifies all the summarized nucleotides to their complements.
        '''
        self.bases = [CompiledPosition._complement_map[b] for b in self.bases]
        self.ref = CompiledPosition._complement_map[self.ref]
        self.counter = False

    def mean_quality(self):
        '''
        Computers the mean quality across all reads.

        Returns:
            Mean read quality as a float
        '''
        return sum(self.qualities) / len(self.qualities) if self.qualities else 0

    def _count_bases(self):
        self.counter = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        for base in self.bases:
            self.counter[base] += 1

    def frequencies(self):
        '''
        Returns a list of all base frequencies in order A, C, G, T
        '''
        return [self.base_frequency(i) for i in ('A', 'C', 'G', 'T')]

    def base_frequency(self, base):
        '''
        Reports the frequency of a given nucleotide at this position.

        Parameters:
            base (str): The nucleotide (A, C, G, or T)

        Returns:
            The total number of reads with the given base as an int
        '''
        if not self.counter:
            self._count_bases()
        return self.counter[base]

    def ref_frequency(self):
        '''
        Determines the number of reads where the base is the same as a reference genome.

        Returns:
            Real number count of reference bases at the position.
        '''
        return self.base_frequency(self.ref)

    def _get_alts(self):
        return {'A', 'C', 'G', 'T'} - {self.ref}

    def get_min_edits(self):
        '''
        Returns the lowest edit frequency.
        '''
        return (self.base_frequency(i) for i in self._get_alts())

    def get_max_edits(self):
        '''
        Returns the highest edit frequency.
        '''
        return max(self.base_frequency(i) for i in self._get_alts())

    def get_max_ratio(self):
        '''
        Returns the highest ratio of variant to variant + reference frequency.
        '''
        ref_freq = self.base_frequency(self.ref)
        var_freq = self.get_max_edits()
        return var_freq / (var_freq + ref_freq)

    def get_variants(self):
        '''
        Returns a list of all detected variants.
        '''
        return [f"{self.ref}{i}" for i in self._get_alts() if self.base_frequency(i) > 0]

    def get_strand(self, threshold=0):
        '''
        Determines the mean strandedness of a position.

        Parameters:
            threshold (int) - Confidence minimum for strand identification

        Returns
            '+', '-', or '*'
        '''
        strand_counts = {'+': 0, '-': 0, '*': 0}
        for i in self.strands:
            strand_counts[i] += 1
        total = strand_counts["+"] + strand_counts["-"]
        if total == 0:
            return '*'

        strand = max(strand_counts, key=strand_counts.get)
        if strand_counts[strand] / total >= threshold:
            return strand
        return '*'

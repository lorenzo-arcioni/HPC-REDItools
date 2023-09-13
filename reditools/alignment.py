#!/usr/bin/env python
'''
Wrappers for pysam files.
'''
# https://pysam.readthedocs.io/en/latest/api.html#pysam.HTSFile.parse_region

from itertools import chain
from pysam.libcalignmentfile import AlignmentFile as PysamAlignmentFile
from pysam.libcfaidx import FastaFile as PysamFastaFile

class AlignmentManager:
    '''
    Manages multiple AlignmentFiles with a single fetch
    '''
    def __init__(self, *args, min_quality, min_length, files=None, **kwargs):
        self._bam_args = args
        self._bam_kwargs = kwargs
        self._min_quality = min_quality
        self._min_length = min_length
        self._bams = []

        if files:
            for bam in files:
                self.add_file(bam)

    def add_file(self, fname):
        '''
        Add an alignment file to the manager for analysis.

        Parameters:
            fname (str) - Path to BAM file.
        '''
        new_file = AlignmentFile(
            fname,
            *self._bam_args,
            min_quality=self._min_quality,
            min_length=self._min_length,
            **self._bam_kwargs)
        new_file.check_index()
        self._bams.append(new_file)

    def _clear_none(self, reads):
        return [(i, j) for i, j in reads if j is not None]

    def fetch_by_position(self, *args, **kwargs):
        '''
        Performs the same as reditools.alignment.AlignmentFile, but combines reads
        from all added files.
        '''
        iterators = [bam.fetch_by_position(*args, **kwargs) for bam in self._bams]
        reads = self._clear_none([(i, next(i, None)) for i in iterators])

        while reads:
            position = min(j[0].reference_start for _, j in reads)
            position_reads = [j for _, j in reads if j[0].reference_start == position]
            yield list(chain(*position_reads))
            reads = self._clear_none([
                (i, next(i, None) if j[0].reference_start == position else j) for i, j in reads])

class AlignmentFile(PysamAlignmentFile):
    '''
    Wrapper for pysam.AlignmentFile to provide filtering on fetch.
    '''
    def __new__(cls, *args, **kwargs):
        kwargs.pop("min_quality", None)
        kwargs.pop("min_length", None)
        return PysamAlignmentFile.__new__(cls, *args, **kwargs)

    def __init__(self, *args, min_quality=0, min_length=0, **kwargs): # pylint: disable=unused-argument
        '''
        Creates a wrapper for pysam.AlignmentFile.

        Parameters:
            min_quality (int): Minimum read quality
            min_length (int): Minimum read length
        '''
        PysamAlignmentFile.__init__(self)

        self._checklist = [
            self._check_flags,
            self._check_tags
        ]

        if min_quality > 0:
            self._min_quality = min_quality
            self._checklist.append(self._check_quality)

        if min_length > 0:
            self._min_length = min_length
            self._checklist.append(self._check_length)

    _flag_filters = {
        77: 'NOT_MAPPED',
        141: 'NOT_MAPPED',
        512: 'QC_FAIL',
        256: 'IS_SECONDARY',
        2048: 'IS_SUPPLEMENTARY',
        1024: 'IS_DUPLICATE'}

    _paired_flag_filters = {
        99: 'NOT_PROPER',
        147: 'NOT_PROPER',
        83: 'NOT_PROPER',
        163: 'NOT_PROPER'}

    def _check_flags(self, read):
        '''
        Checks that a read is usable for analysis based on SAM Flags.

        Parameters:
            read (AlignmentRead)

        Returns:
            True if usable, else False
        '''
        return read.flag not in AlignmentFile._flag_filters and (
            not read.is_paired or read.flag in AlignmentFile._paired_flag_filters)

    def _check_quality(self, read):
        '''
        Determines if a read has sufficient quality.

        Parameters:
            read (AlignmentRead)

        Returns:
            True if the mapping quality is >= minimum read quality, else False
        '''
        return read.mapping_quality >= self._min_quality

    def _check_length(self, read):
        '''
        Determines if a read is long enough to use in analysis.

        Parameters:
            read (AlignmentRead)

        Returns:
            True if the read length is >= minimum read length, else False
        '''
        return read.query_length >= self._min_length

    def _check_tags(self, read):
        '''
        Determines whether a read is chimeric.

        Paramters:
            read (AlignmentRead)

        Returns:
            True if the read is NOT chimeric, else False
        '''
        return not read.has_tag('SA')

    def _check_read(self, read):
        '''
        Determines whether a read passes quality checks.

        Parameters:
            read (ALignmentRead)

        Returns:
            True if all checks pass, else False
        '''
        for check in self._checklist:
            if not check(read):
                return False
        return True

    def fetch(self, *args, **kwargs):
        '''fetch reads aligned in a region.'''
        try:
            iterator = super().fetch(*args, **kwargs)
        except ValueError:
            return
        for read in iterator:
            if self._check_read(read):
                yield read

    def fetch_by_position(self, *args, **kwargs):
        '''
        Retrieves groups of reads that all start at the same point on the reference.

        Returns
        -------
        Generator of lists containing reads
        '''
        iterator = self.fetch(*args, **kwargs)

        first_read = next(iterator, None)
        if first_read is None:
            return

        reads = [first_read]
        ref_start = first_read.reference_start

        for read in iterator:
            if read.reference_start == ref_start:
                reads.append(read)
            else:
                yield reads
                reads = [read]
                ref_start = read.reference_start
        yield reads


class FastaFile(PysamFastaFile):
    '''
    Wrapper for pysam.FastaFile to provide sequence cache.
    '''
    def __new__(cls, *args, **kwargs):
        return PysamFastaFile.__new__(cls, *args, **kwargs)

    def __init__(self, *args, **kwargs): # pylint: disable=unused-argument
        '''
        Creates a wrapper for pysam.FastaFile.
        '''
        PysamFastaFile.__init__(self)

        self._contig_name = False
        self._contig_cache = None

    def _update_contig_cache(self, contig):
        keys = (contig, "chr" + contig, contig.replace("chr", ""))
        for ref in keys:
            try:
                self._contig_cache = self.fetch(ref).upper()
                self._contig_name = contig
                return
            except KeyError:
                pass
        raise KeyError(f"Reference name {contig} not found in FASTA file.")


    def get_base(self, contig, *position):
        '''
        Retrieves the base at the given position.

        Parameters:
            contig (string): Chromsome name
            position (int): Zero-indexed position on reference
        Returns:
            Base the position as a string.
        '''
        if contig != self._contig_name:
            self._update_contig_cache(contig)
        try:
            if len(position) == 1:
                return self._contig_cache[position[0]]
            return [self._contig_cache[i] for i in position]
        except IndexError as exc:
            raise IndexError(f"Tried to lookup base at {position} on contig {contig}, " +
                             f"but contig length is {len(self._contig_cache)}. " +
                             "Are you using the correct reference?") from exc

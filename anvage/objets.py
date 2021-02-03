class CdsSeq:
    def __init__(self, seqid, start, end, sequence):
        self.seqid = seqid
        self.start = start
        self.end = end
        self.sequence = sequence
        self.codons = [sequence[i:i+3] for i in range(0, len(sequence), 3)]
    def __eq__(self, other): 
        if not isinstance(other, CdsSeq):
            # don't attempt to compare against unrelated types
            return NotImplemented
        return self.seqid == other.seqid and self.start == other.start and self.end == other.end and self.sequence == other.sequence and self.codons == other.codons
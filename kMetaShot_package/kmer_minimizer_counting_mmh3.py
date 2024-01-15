import re
import numpy as np
import pandas as pd
from gzip import open as gzopen
from mmh3 import hash as mm3hash
from datetime import datetime as dt
from bitarray import bitarray as btx
from numba import njit, int32, uint32
from .reverse_complement import complrev

CONVERSION_TABLE = {
    65: [False, False], #A
    67: [False, True], #C
    71: [True, False], #G
    84: [True, True], #T
    85: [True, True] #U
}


class EmptyError(Exception):
    'When an object is unexpectedly empty '
    pass


@njit
def itou(nr: int32) -> uint32:
    if nr > 0:
        return nr
    elif nr < 0:
        return (2**32) + nr
    else:
        return 0


def return_time(message: str = None):
    t = dt.now()
    s = str(t).split('.')[0]
    if message is None:
        print("[%s] Return time: e' l'ora dei teletubbies" % s)
    else:
        print("[%s] %s" % (s, message))


class SequenceFile:
    def __init__(self, path_of_file, cds_ncrna: bool):
        self.cds_ncrna = cds_ncrna
        self.path = path_of_file
        if cds_ncrna:
            self.path_cds = path_of_file.replace('_genomic.', '_cds_from_genomic.')
            self.path_rna = path_of_file.replace('_genomic.', '_rna_from_genomic.')
            self.buffer_cds_rna = self.open_file()
            self.buffer = []
            try:
                self.buffer = self.buffer_cds_rna[0].readlines()
            except TypeError:
                pass
            try:
                self.buffer.extend(self.buffer_cds_rna[1].readlines())
            except TypeError:
                pass
        else:
            self.buffer = self.open_file()
        self.sequences = list()
        self.nphashtable = pd.DataFrame()

    def open_file(self):
        try:
            if self.cds_ncrna:
                # Creo il buffer delle RNA
                if self.path_cds.endswith('.gz') or self.path_cds.endswith('.gzip'):
                    a = gzopen(self.path_cds,'rb')
                else:
                    a = open(self.path_cds,'rb')
                # Creo il buffer degli RNA
                if self.path_rna.endswith('.gz') or self.path_rna.endswith('.gzip'):
                    b = gzopen(self.path_rna,'rb')
                else:
                    b = open(self.path_rna,'rb')
                return a, b
            else:
                # Creo il Buffer del file genomico
                if self.path.endswith('.gz') or self.path.endswith('.gzip'):
                    return gzopen(self.path,'rb')
                else:
                    return open(self.path,'rb')
        except FileNotFoundError as e:
            return_time(str(format(e)))

    @staticmethod
    def convert(sequence):
        bn = btx('')
        conversion = [bn.extend(CONVERSION_TABLE[digit]) for digit in sequence]
        return bn

    def counter_kmer(self, sequence, kmer_len):
        # sequence, kmer_len = t
        offset = 64 - kmer_len*2 # numero di bit 0 da aggiungere in coda alla bitstring per arrivare a 64 bit
        sequence_seeker = 0
        for i in range(0, (len(sequence)-(kmer_len * 2) + 1), 2):
            kmer = sequence[i: i + (kmer_len*2)]
            kmer.extend('0'*offset)
            yield itou(mm3hash(kmer.tobytes()))
            sequence_seeker += 2

    def counter_kmer_minimizer(self, t):
        ##############################################
        ######### ACHTUNG
        ######### L'OFFSET VA SETTATO SEMPRE ED AGGIUNTO SEMPRE
        ######### IN CODA AL BIT ARRAY ALTRIMENTI LA TRASFORMAZIONE
        ######### IN INTERO A 32 BIT NON E' STABILE
        ##############################################
        sequence, kmer_len, minim_len = t
        kmer_range = range(kmer_len * 2, len(sequence), 2)
        minim_range = range(0, ((kmer_len * 2) - (minim_len * 2) + 2), 2)
        # la sequenza è una stringa di bit {0, 1}
        # offset: numero di bit 0 da aggiungere in coda alla bitstring per arrivare a 32 bit
        offset = 32 - minim_len * 2
        kmer0 = sequence[: kmer_len * 2]
        mers = [kmer0[j: j + (minim_len * 2)] for j in minim_range]
        exe = [m.extend('0' * offset) for m in mers]
        mm = min(mers)
        retn = itou(mm3hash(mm.tobytes()))
        # print(0, retn, mm, [len(m) for m in mers])
        yield retn
        for i in kmer_range:
            minimiz = sequence[i - (minim_len * 2) + 2: i + 2]
            minimiz.extend('0' * offset)
            # print(i - (minim_len * 2) + 2, i+2)
            mers = mers[1:]
            mers.append(minimiz)
            minimum = min(mers)
            # print(i - (kmer_len*2) + 2, minimum)
            # print([len(m) for m in mers])
            yield itou(mm3hash(minimum.tobytes()))


class FastaReference(SequenceFile):
    def __init__(self, path_of_file, cds_ncrna):

        super().__init__(path_of_file, cds_ncrna)

    def parseFasta_no_complrev(self):
        try:
            sequence = list()
            for seq in self.buffer:
                if seq.startswith(b'>') is False:
                    # print(seq)
                    sequence.append(seq.strip())
                else:
                    # print(seq)
                    sequence = b''.join(sequence)
                    sequence = sequence.upper()
                    for s in re.split(b'[RYSWKMBDHVN]', sequence):
                        # print(s)
                        self.sequences.append(self.convert(s))
                    sequence = list()
            sequence = b''.join(sequence)
            sequence = sequence.upper()
            for s in re.split(b'[RYSWKMBDHVN]', sequence):
                # print(s)
                self.sequences.append(self.convert(s))
            #self.buffer.close()
        except EOFError:
            print('An EOFError has just occurred and managed in %s.' % self.path)

    def parseFasta_complrev(self):
        '''It parses Fasta file, converts each sequence and
        its reverse complement in 2 bit strings, stored in self.sequences list'''
        try:
            sequence = list()
            for seq in self.buffer:
                if seq.startswith(b'>') is False:
                    sequence.append(seq.strip())
                else:
                    sequence = b''.join(sequence)
                    sequence = sequence.upper()
                    for s in re.split(b'[RYSWKMBDHVN]', sequence):
                        # print(s)
                        self.sequences.append(self.convert(s))
                        self.sequences.append(self.convert(complrev(s)))
                    sequence = list()
            sequence = b''.join(sequence)
            sequence = sequence.upper()
            for s in re.split(b'[RYSWKMBDHVN]', sequence):
                # print(s)
                self.sequences.append(self.convert(s))
                self.sequences.append(self.convert(complrev(s)))
            #self.buffer.close()
        except EOFError:
            print('An EOFError has just occurred and menaged in %s.' % self.path)
            pass

    def parser(self, complrev=False):
        if complrev is False:
            # NO COMPLREV
            self.parseFasta_no_complrev()
        else:
            # SI COMPLREV
            self.parseFasta_complrev()

    def seq2bit2kmers(self, kmer_len, minimizer=None, complrev=False):#, taxid):

        self.parser(complrev=complrev)
        if minimizer is None:
            # NO MINIMIZER
            iter_exe = (self.counter_kmer(s, kmer_len) for s in self.sequences\
                        if len(s) >= kmer_len*2)
            dt = np.dtype([('kmer', np.uint32()),
                           ('occurr', np.uint8())])
        else:
            # SI MINIMIZER
            iter_exe = (self.counter_kmer_minimizer((s, kmer_len, minimizer))
                        for s in self.sequences if len(s) >= kmer_len*2)
            dt = np.dtype([('kmer', np.uint32()),
                           ('occurr', np.uint8())])
        try:
            hashtable = dict()
            for seq_res in iter_exe:
                for kmer in seq_res:
                    hashtable.setdefault(kmer, 0)
                    hashtable[kmer] += 1
            self.nphashtable = pd.DataFrame(np.array(list(hashtable.items()), dtype=dt))
        except NotImplementedError:
            print(hashtable)

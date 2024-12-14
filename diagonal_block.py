#!/usr/bin/env python3

import collections

NO_RESOURCE = False
try:
    import resource
except:
    NO_RESOURCE = True
import subprocess
import sys
import time
import typing
import gzip
import argparse
import os


TEN_MB: typing.Final = 10_000_000
TEN_KB: typing.Final = 10_000

RUSAGE_ATTRS: typing.Final = [
    "ru_utime",
    "ru_stime",
    "ru_maxrss",
    "ru_minflt",
    "ru_majflt",
    "ru_inblock",
    "ru_oublock",
    "ru_nvcsw",
    "ru_nivcsw",
]

FastaSequence = collections.namedtuple(
    "FastaSequence", ["description", "sequence", "length"], defaults=["", "", 0]
)

if not NO_RESOURCE:
    def debug_start(
        who: int = resource.RUSAGE_SELF, message: str = ""
    ) -> tuple[resource.struct_rusage, int, int]:
        print(f"DEBUG: {message}", file=sys.stderr, flush=True)
        r_beg = resource.getrusage(who)
        beg = time.monotonic_ns()
        return r_beg, beg, who


    def debug_end(
        r_beg: resource.struct_rusage, beg: int, who: int, message: str = ""
    ) -> None:
        ns = time.monotonic_ns() - beg
        r_end = resource.getrusage(who)
        print(f"DEBUG: {message}: {ns} ns", file=sys.stderr, flush=True)
        for rusage_attr in RUSAGE_ATTRS:
            value = getattr(r_end, rusage_attr) - getattr(r_beg, rusage_attr)
            print(f"DEBUG:   {rusage_attr}: {value}", file=sys.stderr, flush=True)
            
class FastaFile:
    """
    From Rico's code
    """
    def __init__(self, pathname: str, debug: bool = False) -> None:
        self.debug = debug
        self.pathname = pathname
        self.sequences: list[FastaSequence] = []
        self._read_fasta()
        self.sequences.sort(key=lambda x: x.length, reverse=True)

    def _read_fasta(self) -> None:
        description = ""
        seqs: list[str] = []

        if self.debug and not NO_RESOURCE:
            debug_r_beg, debug_beg, debug_who = debug_start(
                resource.RUSAGE_SELF, f"loading fasta {self.pathname}"
            )

        with self._get_open_method() as f:
            for line in f:
                line = line.rstrip()

                if line.startswith(">"):
                    if seqs:
                        sequence = "".join(seqs)
                        self.sequences.append(
                            FastaSequence(description, sequence, len(sequence))
                        )
                        seqs.clear()

                    description = line
                else:
                    seqs.append(line)

            if seqs:
                sequence = "".join(seqs)
                self.sequences.append(
                    FastaSequence(description, sequence, len(sequence))
                )

        if self.debug:
            debug_end(
                debug_r_beg,
                debug_beg,
                debug_who,
                f"loaded fasta {self.pathname}",
            )

    def _get_open_method(self) -> typing.TextIO:
        try:
            with open(self.pathname, "rb") as f:
                if f.read(2) == b"\x1f\x8b":
                    return gzip.open(self.pathname, mode="rt")
        except FileNotFoundError:
            sys.exit(f"ERROR: Unable to read file: {self.pathname}")
        except Exception:
            pass

        return open(self.pathname, mode="rt")

    @property
    def total_bases(self) -> int:
        total = 0
        for fasta_sequence in self.sequences:
            total += fasta_sequence.length

        return total

    def __iter__(self) -> typing.Iterator[FastaSequence]:
        for fasta_sequence in self.sequences:
            yield fasta_sequence

    def to_single_seq(self) -> None:
        description = self.sequences[0].description
        sequence = "".join([seq.sequence for seq in self.sequences])

        self.sequences.clear()
        self.sequences.append(FastaSequence(description, sequence, len(sequence)))

    def discard_sequences_after_and_including(
        self, description: str, debug: bool = False
    ) -> None:
        split_index = -1
        for idx, sequence in enumerate(self.sequences):
            if sequence.description == f">{description}":
                split_index = idx
                break

        if split_index == -1:
            sys.exit(f"ERROR: sequence {description} not found")

        if debug:
            print(
                f"DEBUG: discarding sequences after and including {description}",
                file=sys.stdout,
                flush=True,
            )

        if split_index == 0:
            self.sequences.clear()
        else:
            self.sequences = self.sequences[:split_index - 1]
            
    def masking(self, mask_ranges: list[tuple[int, int]], sequence_idx) -> None:
        if not mask_ranges:
            return
        desc, _seq, length = self.sequences[sequence_idx]
        seq = list(_seq)
        for mask_range in mask_ranges:
            start, end = mask_range
            seq[start:end] = [i.lower() for i in seq[start:end]]
        seq = "".join(seq)
        self.sequences[sequence_idx] = FastaSequence(desc, seq, length)
        
    def _write(self, output_file: str, line_size: int = None) -> None:
        with open(output_file, "w") as f:
            for seq in self.sequences:
                f.write(f"{seq.description}\n")
                sequence = seq.sequence
                if line_size is not None:
                    for i in range(0, len(sequence), line_size):
                        f.write(f"{sequence[i:i+line_size]}\n")
                else:
                    f.write(f"{sequence}\n")
    def write(self, output: typing.TextIO , line_size: int = None) -> None:
        for seq in self.sequences:
            print(f"{seq.description}\n", file=output, end="")
            sequence = seq.sequence
            if line_size is not None:
                for i in range(0, len(sequence), line_size):
                    print(f"{sequence[i:i+line_size]}\n", file=output, end="")
            else:
                print(f"{sequence}\n", file=output, end="")
                    
def segment_fasta(sequence: str, description: str, segment_length: int, overlap: int, output_dir: str) -> None:
    start = 0
    end = segment_length

    while start < len(sequence):
        segment = sequence[start:end]
        segment_description = f"{description} {start+1}-{min(end, len(sequence))}"
        segment_filename = os.path.join(
            output_dir, f"segment_{start+1}_{min(end, len(sequence))}.fasta"
        )
        with open(segment_filename, "w") as output_file:
            output_file.write(f"{segment_description}\n")
            output_file.write(f"{segment}\n")

        if start + segment_length >= len(sequence):
            break
        start += segment_length - overlap
        end = start + segment_length
                        
if __name__ == "__main__":
    DEBUG = False
    parser = argparse.ArgumentParser(description="Fasta File Masking")
    parser.add_argument(
        "input",
        type=str,
        help="Path to the input fasta file.",
    )
    parser.add_argument(
        "output_dir",
        type=str,
        help="Directory to output the segmented fasta files.",
    )
    parser.add_argument(
        "--segment_length",
        type=int,
        required=True,
        help="Length of each segment.",
    )
    parser.add_argument(
        "--overlap",
        type=int,
        required=True,
        help="Overlap length between consecutive segments.",
    )
    args = parser.parse_args()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    fasta_file = FastaFile(args.input, debug=DEBUG)
    if len(fasta_file.sequences) != 1:
        sys.exit("ERROR: The input fasta file must contain exactly one chromosome.")

    sequence = fasta_file.sequences[0].sequence
    description = fasta_file.sequences[0].description
    segment_length = args.segment_length
    overlap = args.overlap
    start = 0
    end = segment_length

    segment_fasta(sequence, description, segment_length, overlap, args.output_dir)
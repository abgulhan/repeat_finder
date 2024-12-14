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
from concurrent.futures import ProcessPoolExecutor


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
            self.sequences = self.sequences[: split_index - 1]
            
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
                    
                    
def parse_line(line: str, format: str) -> tuple[str, str, int, int, str, int, int]:
    """
    Returns (direction, query_seq_name, query_start, query_end, target_seq_name, target_start, target_end)
    """
    if format == "paf":
        query_seq_name, query_seq_len, query_start, query_end, direction, \
        target_seq_name, target_seq_len, target_start, target_end = line.split()[:9]    
    elif format == "segment": # target comes first in segment file
        target_seq_name, target_start, target_end, \
        query_seq_name, query_start, query_end, direction, score = line.split()
    else:
        sys.exit(f"Error: invalid format {format}")
    # TODO make custom type for return
    return (direction, query_seq_name, int(query_start), int(query_end), target_seq_name, int(target_start), int(target_end))       

def filter_diagonal(lines: list[str], diagonal_radius: int, format: str, debug=False) -> list[tuple[int, int]]:
    if debug:
        print(f"# Filtering diagonals with radius {diagonal_radius}")
    coordinates = []
    for line in lines:
        #query_seq_name, query_seq_len, query_start, query_end, direction, \
        #target_seq_name, target_seq_len, target_start, target_end = line.split()[:10]
        direction, query_seq_name, query_start, query_end, target_seq_name, target_start, target_end = parse_line(line, format)
        if direction == "-":
            continue
        assert direction == "+", f"Error: invalid direction {direction}"
        assert query_seq_name == target_seq_name, f"Error: query and target sequence names do not match: {query_seq_name}, {target_seq_name}"
        half_dist = int((int(query_end) - int(query_start)) // 2)
        assert int(query_end) > int(query_start)
        assert int(target_end) > int(target_start)
        query_mid = int(query_start) + half_dist
        target_mid = int(target_start) + half_dist
        if abs(query_mid - target_mid) <= diagonal_radius:
            coordinates.append((int(query_start), int(query_end)))
    return coordinates
              
def filter_diagonal_parallel(lines: list[str], diagonal_radius: int, format: str, cores: int = -1, debug=False) -> list[tuple[int, int]]:
    if debug:
        print(f"# Filtering diagonals with radius {diagonal_radius}")

    def process_line(line: str) -> typing.Optional[tuple[int, int]]:
        direction, query_seq_name, query_start, query_end, target_seq_name, target_start, target_end = parse_line(line, format)
        if direction == "-":
            return None
        assert direction == "+", f"Error: invalid direction {direction}"
        assert query_seq_name == target_seq_name, f"Error: query and target sequence names do not match: {query_seq_name}, {target_seq_name}"
        half_dist = int((int(query_end) - int(query_start)) // 2)
        assert int(query_end) > int(query_start)
        assert int(target_end) > int(target_start)
        query_mid = int(query_start) + half_dist
        target_mid = int(target_start) + half_dist
        if abs(query_mid - target_mid) <= diagonal_radius:
            return (int(query_start), int(query_end))
        return None

    chunk_size = max(1, len(lines) // (cores if cores > 0 else 1))
    chunks = [lines[i:i + chunk_size] for i in range(0, len(lines), chunk_size)]

    with ProcessPoolExecutor(max_workers=cores if cores > 0 else None) as executor:
        results = list(executor.map(lambda chunk: [process_line(line) for line in chunk], chunks))
    
    # Flatten the list of results
    results = [item for sublist in results for item in sublist if item is not None]

    # Filter out None values
    coordinates = [result for result in results if result is not None]

    # Filter out None values
    coordinates = [result for result in results if result is not None]

    return coordinates
              
if __name__ == "__main__":
    DEBUG = False
    parser = argparse.ArgumentParser(description="Fasta File Masking")
    parser.add_argument(
        "input",
        type=str,
        help="Path to the input segment file, directory of segment files, .paf file, or 'keg' file. Must only contains HSPs for one self-aligned chromosome.",
    )
    parser.add_argument(
        "output_file", 
        type=str,
        help="Path to the output masked fasta file. If not given, output will be range of regions to mask.",
        nargs='?',
        default=None
    )
    parser.add_argument(
        "-t", "--mask_threshold", 
        type=int, 
        default=10, 
        help="Min number of HSPs in column to begin masking (default: 10)"
    )
    parser.add_argument(
        "--mask", "-m", 
        type=str,
        choices=["soft", "s", "hard", "h"], 
        help="Choose between soft masking or hard masking (default: soft)",
        default="soft"
    )
    
    parser.add_argument(
        "--line_size", 
        type=int,
        help="Select line size for masked fasta file. Set 0 for no line limit (default: 50)",
        default=50
    )
    parser.add_argument(
        "--print_regions", "-r",
        action="store_true",
        help="Print the masked regions to stderr (default: False)",
        default=False
    )
    parser.add_argument(
        "--region_format", "-rf",
        type=str,
        choices=['default', 'bed'],
        help="Format to print masked regions in (default: start:end)",
        default='default'
    )
    parser.add_argument(
        "--filter_diagonal", "-f", 
        type=int,
        help="How many base pairs above and below diagonals to include (default: no filtering)",
        default=None
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Print debug info to stdout (default: False)",
        default=False
    )
    
    parser.add_argument(
        "--multiprocessing", "-mp",
        type=int,
        help="Number of cores to use. -1 to use all cores (default: 1)",
        default=1
    )

    args = parser.parse_args()
    
    ### Read HSP coordinates ###
    input = args.input
    input_type = None
    coordinates: list[tuple[int, int]] = []
    chr_name = None # used only for printing bed formatted regions
    
    if args.debug:
        print(f"# Reading input ...")
    
    if os.path.isdir(input):
        # Input is a directory
        for file_name in os.listdir(input):
            file_path = os.path.join(input, file_name)
            with open(file_path, "r") as file:
                lines = file.readlines()
            coordinates.extend([(int(line.split()[1]), int(line.split()[2])) for line in lines])
    elif input.endswith(".tgz"):
        # Input is a tgz or keg file
        sys.exit("ERROR: keg file support not implemented.")
    elif input.endswith(".paf"):
        # Input is a paf file
        with open(input, "r") as file:
            lines = file.readlines()
        if args.filter_diagonal is not None:
            if args.multiprocessing == 1:
                coordinates = filter_diagonal(lines, args.filter_diagonal, format="paf", debug=args.debug)
            else:
                coordinates = filter_diagonal_parallel(lines, args.filter_diagonal, format="paf", debug=args.debug, cores=args.multiprocessing)
        else:   
            coordinates = [(int(line.split()[2]), int(line.split()[3])) for line in lines]
        chr_name = lines[0].split()[0]
        
    elif os.path.isfile(input):
        # Input is a segment file
        with open(input, "r") as file:
            lines = file.readlines()
        if args.filter_diagonal is not None:
            if args.multiprocessing == 1:
                coordinates = filter_diagonal(lines, args.filter_diagonal, format="segment", debug=args.debug)
            else:
                coordinates = filter_diagonal_parallel(lines, args.filter_diagonal, format="segment", debug=args.debug, cores=args.multiprocessing)
        else:
            coordinates = [(int(line.split()[1]), int(line.split()[2])) for line in lines]
        chr_name = lines[0].split()[0]

    else:
        input_type = "unknown"
        sys.exit(f"ERROR: input {input} not recognized as valid file or directory.")
    
    assert len(coordinates) > 0, "Error: no coordinates found in input file after filtering."
    
    ### Determine masked regions ###
    
    # transform coordinates
    augment_list = []

    for i in range(len(coordinates)):
        augment_list.append((coordinates[i][0], 1)) # start
        augment_list.append((coordinates[i][1], -1)) # end


    # sort by start and end coordinates
    augment_list = sorted(augment_list, key=lambda x: (x[0])) 
    # NOTE: for cases with same coordinates, start should come first (key=lambda x: (x[0], x[-1])) 
    # to prevent cases where end = start coordinate of subsequent masked regions
    # HOWEVER, this increases sorting time by 2.5x
    # and these dont effect the final result


    threshold = args.mask_threshold
    # determine regions with overlap > threshold
    overlaps = []
    current_overlap = 0
    overlap_start = None
    for i in range(len(augment_list)):
        current_overlap += augment_list[i][1]
        if current_overlap >= threshold:
            if overlap_start is None:
                overlap_start = augment_list[i][0]
        else:
            if overlap_start is not None:
                overlap_end = augment_list[i][0]
                overlaps.append((overlap_start, overlap_end))
                overlap_start = None
    assert current_overlap == 0, f"Error in overlap calculation. Final overlap should be 0, but is {current_overlap}"
    
    masked_intervals = overlaps
    
    
    ### Output masked regions ###
    if args.debug:
        print(f"# Total masked base pairs: {sum([j-i for i,j in masked_intervals])}")
    
    if args.debug:
        if len(masked_intervals) > 0:
            sorted_masked = sorted(masked_intervals)
            prev_end = sorted_masked[0][1]
            for idx, (s,e) in enumerate(sorted_masked[1:]):
                #print(f"{idx} - {s}:{e}")
                try:
                    assert prev_end < s, (prev_end, s)
                except:
                    print(f"WARNING: {prev_end}, {s}")
                prev_end = e
    
    if args.output_file is not None:
        target_fasta = FastaFile(args.output_file, debug=False)
        target_fasta.masking(masked_intervals, 0)
        line_size = args.line_size
        if line_size <= 0:
            line_size = None
        target_fasta.write(sys.stdout, line_size=line_size)
    
    if args.print_regions:
        for interval in sorted(masked_intervals):
            if args.region_format == 'bed':
                print(f"{chr_name}\t{interval[0]}\t{interval[1]}\n", file=sys.stderr, end="")
            else:
                print(f"{chr_name} {interval[0]}:{interval[1]}\n", file=sys.stderr, end="")


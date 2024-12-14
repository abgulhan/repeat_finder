# repeat_finder

### Requirements:
- [kegalign](https://github.com/galaxyproject/KegAlign)

## Example usage:

* split input chromosome into overlapping blocks
```
python3 diagonal_block.py HG002_chrY.fasta ./HG002_chrY_split --segment_length 1000000 --overlap 8000
```

* call kegalign and perform self non-gapped alignment for each block
```
mkdir segments
banded_runner.sh HG002_chrY_split ./segments
```

* use the output segment file to determine mask regions and save to bed file

```
python3 hsp_mask.py ./segments/all_segments.plus.segments --mask_threshold 10 --filter_diagonal 4000 --print_regions --region_format bed 2> regions.bed
```

# rumi

Rust UMI based PCR deduplication based on the directional adjacency
as UMI-tools but with a constant time hamming distance implementation.

## Install

```
cargo install --path .
```

## Usage

```bash
rumi input.bam -o output.bam
```

## Performance

I have not sat down and done any serious benchmarking yet. Anecdotally
this is at least 4X faster than umi_tools. There are *A LOT* of low
hanging fruit in terms of optimizations to apply though.

## Known differences from umi_tools

- The criteria for choosing what read to use is different, and in my
  opinion more extensive. For example, if two reads are the same mapq,
  edit distance, etc, it will come down to read length, keeping the
  longer. In cases of an absolute tie, the incumbant read wins. This can
  lead to differences in results, especially among reverse reads.

## TODO

- Clean up the lib and split out into multiple modules
- Add better (any) error handling
- Figure out how to reduce allocations / string clones
- Clarify the workflow / what actually happens in docs
- Allow for selection of number of cores to use

## Prior Art

[UMICollapse](https://www.biorxiv.org/content/10.1101/648683v1)
[UMITools](https://genome.cshlp.org/content/early/2017/01/18/gr.209601.116.abstract)

## Notes

First pass: Collect all reads into a dict that is keyed on position.
Track metrics like umi freq, and extracted umis while building this.
Then iter over that dict and deduplicate at each position.


# Diffs in example.bam (from umi_tools)

- SRR2057595.3354975_CGGGTTGGT: rumi correctly uses the umi starting
  with C since there are two reads with that umi. umi_tools uses the umi
  with only a feq of 1.
- SRR2057595.4915638_TTGGTTAAA: rumi correctly chooses the read with the
  decided upon umi as the best read.
- SRR2057595.5405752_AACGGTTGG: rumi correctly leaves as it's own group.
  umi_tools corrects it dist 3 away to ATTGGTTCG. I expect this to be
  the end source of the 30 extra reads in rumi's output. What causes
  this in umi_tools?

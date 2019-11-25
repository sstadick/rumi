# rumi

Rust UMI based PCR deduplication based on the directional adjacency
as UMI-tools but with a constant time hamming distance implementation.

## Install

```
cargo install --path .
```

## Usage

```bash
$ rumi --help
rumi-dedup 0.1.0
Seth Stadick <sstadick@gmail.com>
Deduplicate reads based on umis

USAGE:
    rumi [FLAGS] [OPTIONS] <INBAM> --output <OUTBAM> --umi_tag <umi_tag>

FLAGS:
        --group_only           Don't deduplicate reads, just group them given them agroup id, and print them. Rules
                                               for filtering out unpaired reads, etc, will still be applied.
    -h, --help                 Prints help information
        --ignore_splice_pos    If two reads have the same start pos, and contain a splice site, they will be
                                                    grouped together, instead of further splitting them based on the
                               splice site
        --is_paired            Input is paired end. Read pairs with unmapped read1 will be ignored.
        --umi_in_read_id       The UMI is located in the read id after the last '_'. Otherwise use the RX tag.
    -V, --version              Prints version information

OPTIONS:
    -o, --output <OUTBAM>                                  Output bam file. Use - if stdout [default: -]
    -c, --allowed_count_factor <allowed_count_factor>
            The factor to multiply the count of a umi by when determining whether or not to group it with other umis
            within allowed_read_dist. include umi_b as adjacent to umi_a if: umi_a.counts >= allowed_count_factor *
            umi_b.counts [default: 2]
    -n, --allowed_network_depth <allowed_network_depth>
            The number of nodes deep to go when creating a group. If allowed_read_dist 1, then allowed_network_depth of
            2 will enable getting all umis with hamming distance of 2 from current umi. [default: 2]
    -d, --allowed_read_dist <allowed_read_dist>
            The distance between umis that will allow them to be counted as adjacent. [default: 1]

    -u, --umi_tag <umi_tag>                                The tag holding the umi information. [default: RX]

ARGS:
    <INBAM>    Input bam file. Use - if stdin [default: -]
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

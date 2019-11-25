use basebits::{hamming_dist_none, BaseBits};
//use rayon::iter::ParBridge;
use rayon::prelude::*;
use rust_htslib::bam::errors::Error;
use rust_htslib::bam::record::{Aux, Cigar, CigarString};
use rust_htslib::bam::{self, Read};
use std::cmp::Ordering;
use std::collections::hash_map::{Entry::Occupied, Entry::Vacant};
use std::collections::{BTreeMap, HashMap, HashSet, VecDeque};
use std::fmt::Display;
use std::process;
use std::sync::mpsc::channel;
use std::sync::{Arc, Mutex};

#[cfg(test)]
mod test;

#[derive(Debug)]
pub struct Config {
    pub allowed_read_dist: u32,
    pub allowed_count_factor: u32,
    pub allowed_network_depth: usize,
    pub umi_tag: String,
    pub input_bam: String,
    pub output_bam: String,
    pub umi_in_read_id: bool,
    pub ignore_splice_pos: bool,
    pub group_only: bool,
    pub is_paired: bool,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Node {
    umi: BaseBits,
    freq: ReadFreq,
    connections: Vec<usize>,
}

#[derive(Hash, PartialEq, Eq, Debug, Clone)]
pub struct Position {
    pos: i32,
    is_spliced: Option<u32>,
    is_rev: bool,
    target: i32,
    tlen: Option<i32>,
}

impl PartialOrd for Position {
    fn partial_cmp(&self, other: &Position) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Position {
    fn cmp(&self, other: &Position) -> Ordering {
        let comp = self.target.cmp(&other.target);
        if comp != Ordering::Equal {
            return comp;
        }

        let comp = self.pos.cmp(&other.pos);
        if comp != Ordering::Equal {
            return comp;
        }

        let comp = self.tlen.cmp(&other.tlen);
        if comp != Ordering::Equal {
            return comp;
        }

        let comp = self.is_spliced.cmp(&other.is_spliced);
        if comp != Ordering::Equal {
            return comp;
        }

        self.is_rev.cmp(&other.is_rev)
    }
}

impl Position {
    /// Takes a read and determins the position to use as a key in the returned group.
    pub fn new(record: &bam::record::Record, ignore_splice_pos: bool, use_tlen: bool) -> Self {
        let mut pos = record.pos();
        let mut is_spliced: Option<u32>;
        let tlen: Option<i32>;
        let cigarview = record.cigar();
        let cigar = &cigarview;

        if record.is_reverse() {
            pos = cigarview.end_pos();
            // if the end of the read was soft clipped, add that amount back to its pos
            if let Cigar::SoftClip(num) = cigar[cigar.len() - 1] {
                pos = pos + num as i32;
            }

            is_spliced = Position::find_splice(&cigar, true);
        } else {
            if let Cigar::SoftClip(num) = cigar[0] {
                pos = pos - num as i32;
            }
            is_spliced = Position::find_splice(&cigar, false);
        }
        if ignore_splice_pos && is_spliced.is_some() {
            is_spliced = Some(0);
        }

        if use_tlen {
            tlen = Some(record.insert_size());
        } else {
            tlen = None;
        }
        Self {
            pos: pos,
            is_rev: record.is_reverse(),
            target: record.tid(),
            is_spliced: is_spliced,
            tlen: tlen,
        }
    }

    /// Takes a cigar string and finds the first splice postion as an offset from the start.
    /// Equivalent of `find_splice` in umi_tools
    fn find_splice(cigar: &CigarString, is_reversed: bool) -> Option<u32> {
        let mut range: Vec<usize> = (0..cigar.len()).collect();
        let mut offset = 0;

        if is_reversed {
            range = (0..cigar.len()).rev().collect();
        }
        if let Cigar::SoftClip(num) = cigar[range[0]] {
            offset = num;
            range.remove(0);
        }
        for i in range {
            match cigar[i] {
                // Found splice
                Cigar::RefSkip(_) | Cigar::SoftClip(_) => return Some(offset),
                // Reference consumeing operations
                Cigar::Match(num) | Cigar::Del(num) | Cigar::Equal(num) | Cigar::Diff(num) => {
                    offset += num
                }
                // Non-reference consuming operations
                Cigar::Ins(_) | Cigar::HardClip(_) | Cigar::Pad(_) => continue,
            }
        }
        None
    }
}

/// Abstraction so ReadFreq can hold a single best read for it's read signature or hold all reads
/// for it's read signature (used for --group_only).
#[derive(Debug, Clone, PartialEq)]
pub enum ReadCollection {
    SingleRead(bam::record::Record),
    ManyReads(Vec<bam::record::Record>),
}

/// A Read or Reads and the number of times that read signature has been seen
/// Read signature meaning Position + UMI
#[derive(Debug, Clone, PartialEq)]
pub struct ReadFreq {
    read: ReadCollection,
    freq: u32,
}

/// A group of reads that have been deduplicated
#[derive(Debug)]
pub struct Group<'a> {
    nodes: Vec<&'a Node>,
    umi: &'a BaseBits,
    master_node: usize,
}

#[derive(Debug)]
pub enum RecordEvent {
    RecordMapped,
    RecordUnmapped,
    RecordUnpaired,
    RecordMateUnmapped,
    RecordChimeric,
}

#[derive(Debug)]
pub struct Stats {
    reads_in: u32,
    reads_out: u32,
    reads_unmapped: u32,
    reads_unpaired: u32,
    mate_unmapped: u32,
    chimeric: u32,
}

impl Stats {
    pub fn new() -> Self {
        Stats {
            reads_in: 0,
            reads_out: 0,
            reads_unmapped: 0,
            reads_unpaired: 0,
            mate_unmapped: 0,
            chimeric: 0,
        }
    }
    pub fn update(&mut self, other: &Self) {
        self.reads_in += other.reads_in;
        self.reads_out += other.reads_out;
        self.reads_unmapped += other.reads_unmapped;
        self.reads_unpaired += other.reads_unpaired;
        self.mate_unmapped += other.mate_unmapped;
        self.chimeric += other.chimeric;
    }
}

impl Display for Stats {
    fn fmt(&self, fmt: &mut std::fmt::Formatter) -> std::result::Result<(), std::fmt::Error> {
        write!(fmt, "Reads In: {}\nReads Out: {}\nReads Unmapped: {}\nReads Unpaired: {}\nMates Unmapped: {}\nReads Chimeric: {}", self.reads_in, self.reads_out, self.reads_unmapped, self.reads_unpaired, self.mate_unmapped, self.chimeric)
    }
}

pub type UmiMap = HashMap<BaseBits, ReadFreq>;
pub type ReadMap = BTreeMap<Position, UmiMap>;

/// The main function to coordinate the deduplication process
pub fn run_dedup(config: &Config) -> Result<(), &'static str> {
    let mut bam = bam::Reader::from_path(&config.input_bam).unwrap();
    let header = bam::Header::from_template(bam.header());
    let mut writer = bam::Writer::from_path(&config.output_bam, &header, bam::Format::BAM).unwrap();
    let mut read_store: HashSet<Vec<u8>> = HashSet::new();
    let (sender, reciever) = channel();
    let global_stats = Arc::new(Mutex::new(Stats::new()));

    let bundler = Bundler {
        records: bam.records(),
        last_chr: None,
        next_bundle: vec![],
    };

    bundler
        .par_bridge()
        .flat_map(|bundle| {
            let (x, stats) = group_reads(bundle, &config);
            let y: Vec<(Position, UmiMap)> = x.into_iter().collect();
            let mut g_stats = global_stats.lock().unwrap();
            g_stats.update(&stats);
            y
        })
        .flat_map(|(_, reads)| dedup(reads, config))
        .for_each_with(sender, |s, x| s.send(x).unwrap());

    let mut reads_out = 0;
    reciever.iter().for_each(|read| {
        reads_out += 1;
        writer.write(&read).unwrap_or_else(|err| {
            eprintln!("Problem writing: {}", err);
            process::exit(1);
        });
        if config.is_paired {
            read_store.insert(read.qname().to_vec());
        }
    });

    if config.is_paired {
        bam.records()
            .map(|read| read.unwrap())
            .filter(|read| read.is_last_in_template() && read_store.contains(read.qname()))
            .for_each(|read| {
                reads_out += 1;
                writer.write(&read).unwrap_or_else(|err| {
                    eprintln!("Problem writing: {}", err);
                    process::exit(1);
                });
            });
    }

    let mut stats = global_stats.lock().unwrap();
    stats.reads_out += reads_out;
    println!("{}", stats);
    Ok(())
}

pub fn run_group(config: &Config) -> Result<(), &'static str> {
    let mut bam = bam::Reader::from_path(&config.input_bam).unwrap();
    let header = bam::Header::from_template(bam.header());
    let mut writer = bam::Writer::from_path(&config.output_bam, &header, bam::Format::BAM).unwrap();
    let mut read_store: HashMap<Vec<u8>, (bam::record::Aux, Vec<u8>)> = HashMap::new();
    let global_stats = Arc::new(Mutex::new(Stats::new()));
    let (sender, reciever) = channel();

    let bundler = Bundler {
        records: bam.records(),
        last_chr: None,
        next_bundle: vec![],
    };

    bundler
        .par_bridge()
        .flat_map(|bundle| {
            let (x, stats) = group_reads(bundle, &config);
            let y: Vec<(Position, UmiMap)> = x.into_iter().collect();
            let mut g_stats = global_stats.lock().unwrap();
            g_stats.update(&stats);
            y
        })
        .flat_map(|(_, reads)| label_groups(reads, config))
        .for_each_with(sender, |s, x| s.send(x).unwrap());

    let mut group_count: i64 = 0;
    let mut reads_out = 0;
    reciever.iter().for_each(|mut group| {
        for read in group.iter_mut() {
            reads_out += 1;
            read.push_aux(b"UG", &bam::record::Aux::Integer(group_count));
            writer.write(&read).unwrap_or_else(|err| {
                eprintln!("Problem writing: {}", err);
                process::exit(1);
            });
            if config.is_paired {
                let umi = read.aux(b"BX").unwrap().string().to_vec();
                read_store.insert(
                    read.qname().to_vec(),
                    (bam::record::Aux::Integer(group_count), umi),
                );
            }
            group_count += 1;
        }
    });

    if config.is_paired {
        bam.records()
            .map(|read| read.unwrap())
            .filter(|read| read.is_last_in_template())
            .for_each(|mut read| {
                if let Some((ug, bx_val)) = read_store.get(read.qname()) {
                    reads_out += 1;
                    read.push_aux(b"UG", ug);
                    read.push_aux(b"BX", &bam::record::Aux::String(&bx_val));
                    writer.write(&read).unwrap_or_else(|err| {
                        eprintln!("Problem writing: {}", err);
                        process::exit(1);
                    });
                }
            });
    }
    let mut stats = global_stats.lock().unwrap();
    stats.reads_out = reads_out;
    println!("{}", stats);
    Ok(())
}

struct Bundler<I>
where
    I: Iterator<Item = Result<rust_htslib::bam::record::Record, Error>>,
{
    records: I,
    last_chr: Option<i32>,
    next_bundle: Vec<rust_htslib::bam::record::Record>,
}

impl<I> Iterator for Bundler<I>
where
    I: Iterator<Item = Result<rust_htslib::bam::record::Record, Error>>,
{
    type Item = Vec<rust_htslib::bam::record::Record>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut bundle = vec![];
        std::mem::swap(&mut self.next_bundle, &mut bundle);
        while let Some(r) = self.records.next() {
            let record = r.unwrap();
            if let Some(tid) = self.last_chr {
                if tid == record.tid() {
                    bundle.push(record);
                } else {
                    self.last_chr = Some(record.tid());
                    self.next_bundle.push(record);
                    break;
                }
            } else {
                self.last_chr = Some(record.tid());
                bundle.push(record);
            }
        }
        if bundle.len() > 0 {
            Some(bundle)
        } else {
            None
        }
    }
}

fn get_tag<'a>(record: &'a bam::record::Record, config: &Config) -> &'a [u8] {
    if config.umi_in_read_id {
        match record.qname().split(|&c| c == b'_').last() {
            Some(tag) => tag,
            None => panic!("No tag in read id"),
        }
    } else {
        match record.aux(config.umi_tag.as_bytes()) {
            Some(tag) => tag.string(),
            None => panic!("No tag on read"),
        }
    }
}

pub fn check_record(record: &bam::record::Record, paired_end: bool) -> RecordEvent {
    if paired_end {
        if record.is_unmapped() {
            return RecordEvent::RecordUnmapped;
        }
        if !record.is_paired() {
            return RecordEvent::RecordUnpaired;
        }
        if record.tid() != record.mtid() {
            return RecordEvent::RecordChimeric;
        }
        if record.is_mate_unmapped() {
            return RecordEvent::RecordMateUnmapped;
        }
    } else {
        if record.is_unmapped() {
            return RecordEvent::RecordUnmapped;
        }
    }
    RecordEvent::RecordMapped
}

/// Group reads together based on their positions.
pub fn group_reads(
    records: Vec<rust_htslib::bam::record::Record>,
    config: &Config,
) -> (ReadMap, Stats) {
    let mut read_map: ReadMap = BTreeMap::new();
    let mut stats = Stats::new();

    for record in records.into_iter() {
        stats.reads_in += 1;

        if config.is_paired && record.is_last_in_template() {
            continue;
        }

        match check_record(&record, config.is_paired) {
            RecordEvent::RecordMapped => {}
            RecordEvent::RecordUnmapped => {
                stats.reads_unmapped += 1;
                continue;
            }
            RecordEvent::RecordUnpaired => {
                stats.reads_unpaired += 1;
                continue;
            }
            RecordEvent::RecordMateUnmapped => {
                stats.mate_unmapped += 1;
            }
            RecordEvent::RecordChimeric => {
                stats.chimeric += 1;
                continue;
            }
        }

        let tag = get_tag(&record, config);
        let position = Position::new(&record, config.ignore_splice_pos, config.is_paired);

        // Add to my reverse lookup
        let bb = BaseBits::new(tag).unwrap();
        let position_map = read_map.entry(position).or_insert(HashMap::new());
        match position_map.entry(bb) {
            Occupied(entry) => {
                let rf = entry.into_mut();
                match &rf.read {
                    ReadCollection::SingleRead(read) => {
                        if !read_a_ge_b(&read, &record) {
                            rf.read = ReadCollection::SingleRead(record);
                        }
                    }
                    ReadCollection::ManyReads(reads) => {
                        // This is stupid but it's only for the group version so....
                        let mut reads = reads.clone();
                        reads.push(record);
                        rf.read = ReadCollection::ManyReads(reads);
                    }
                }
                rf.freq += 1;
            }
            Vacant(entry) => {
                if !config.group_only {
                    entry.insert(ReadFreq {
                        read: ReadCollection::SingleRead(record),
                        freq: 1,
                    });
                } else {
                    entry.insert(ReadFreq {
                        read: ReadCollection::ManyReads(vec![record]),
                        freq: 1,
                    });
                }
            }
        };
    }
    (read_map, stats)
}

/// Create a graph from the UmiMap
/// TODO: Inline?
pub fn build_graph(reads: UmiMap) -> Vec<Node> {
    reads
        .into_iter()
        .map(|(umi, freqs)| Node {
            umi,
            connections: vec![],
            freq: freqs,
        })
        .collect()
}

/// Create the connections between the umis via an all vs all comparison.
/// A Connection will only be formed from a larger node to a smaller node.
/// Larger being defined as node_a >= 2x node_b - 1, the provides the directionality.
/// TODO: Keep a seen list here instead of later? Some connections will be redundant.
pub fn connect_graph(mut graph: Vec<Node>, dist: u32, counts_factor: u32) -> Vec<Node> {
    for i in 0..graph.len() {
        for j in 0..graph.len() {
            if i == j {
                continue;
            }
            if hamming_dist_none(&graph[i].umi, &graph[j].umi) <= dist
                && graph[i].freq.freq >= (counts_factor * graph[j].freq.freq) - 1
            {
                graph[i].connections.push(j);
            }
        }
    }
    graph
}

// TODO: Use proper bk tree for faster lookups
fn determine_umi<'a>(graph: &'a Vec<Node>, allowed_network_depth: usize) -> Vec<Group> {
    // Group the umis by distance
    let mut groups = vec![];
    let mut seen: Vec<usize> = Vec::new();
    // Create a vec of nodes indicies going from highest counts to lowest
    let mut graph_indicies: Vec<usize> = (0..graph.len()).collect();
    &graph_indicies.sort_by(|&a, &b| graph[b].freq.freq.cmp(&graph[a].freq.freq));

    for &x in graph_indicies.iter() {
        if seen.contains(&x) {
            continue;
        }
        seen.push(x);
        let node = &graph[x];
        let mut group_holder: Vec<Vec<&Node>> = Vec::new();

        // Get all the nodes within 1 hamming dist
        let mut group: Vec<&Node> = vec![];
        for &x in node.connections.iter() {
            if !seen.contains(&x) {
                seen.push(x);
                group.push(&graph[x]);
            }
        }

        // Get all the nodes within k hamming dist
        // If two nodes lie equidistant away from a smaller node, it shouldn't matter which node
        // gets the discrepent reads, there would be no real biological way to tell...
        let mut queue: VecDeque<Vec<usize>> = VecDeque::new();
        queue.push_back(group.iter().flat_map(|n| n.connections.clone()).collect());
        for _ in 0..(allowed_network_depth - 1) {
            if let Some(connections) = queue.pop_front() {
                let mut new_group: Vec<&Node> = vec![];
                queue.push_back(
                    connections
                        .iter()
                        .flat_map(|&x| {
                            if !seen.contains(&x) {
                                seen.push(x);
                                new_group.push(&graph[x]);
                            }
                            graph[x].connections.clone()
                        })
                        .collect(),
                );
                group_holder.push(new_group);
            }
        }
        // Must add after, otherwise it will be searched again
        group.push(node);

        // Merge all groups and choose concensus umi
        for g in group_holder {
            group.extend(g.iter());
        }

        let master_node = group.iter().enumerate().fold(0, |max, (i, x)| {
            if x.freq.freq > group[max].freq.freq {
                i
            } else {
                max
            }
        });
        let umi = &group[master_node].umi;
        let group = Group {
            nodes: group,
            umi,
            master_node,
        };

        groups.push(group);
    }
    groups
}

/// Deduplicate a group of reads that all positioned at the same position
fn dedup(reads: UmiMap, config: &Config) -> Vec<bam::record::Record> {
    let graph = build_graph(reads);
    let graph = connect_graph(graph, config.allowed_read_dist, config.allowed_count_factor);
    let groups = determine_umi(&graph, config.allowed_network_depth);
    let mut final_reads = vec![];

    for group in groups.into_iter() {
        let node = group.nodes[group.master_node];
        if let ReadCollection::SingleRead(read) = &node.freq.read {
            let read = read.clone();
            final_reads.push(read);
        } else {
            unreachable!();
        }
    }
    final_reads
}

/// TODO: Don't clone the read :(
fn label_groups(reads: UmiMap, config: &Config) -> Vec<Vec<bam::record::Record>> {
    let graph = build_graph(reads);
    let graph = connect_graph(graph, config.allowed_read_dist, config.allowed_count_factor);
    let groups = determine_umi(&graph, config.allowed_network_depth);
    let mut records = vec![];

    for group in groups.into_iter() {
        let mut group_list = vec![];
        let master_umi = group.nodes[group.master_node].umi;
        for node in group.nodes {
            if let ReadCollection::ManyReads(reads) = &node.freq.read {
                for read in reads.into_iter() {
                    let mut read = read.clone();
                    read.push_aux(b"BX", &bam::record::Aux::String(&master_umi.decode()));
                    group_list.push(read);
                }
            } else {
                unreachable!();
            }
        }
        records.push(group_list);
    }
    records
}

/////////////////////// Helpers
/// Decide wich read is better.
/// For now this uses the simplistic approach of comparing mapq values.
/// Returns true if alpha is better than beta, false otherwise
fn read_a_ge_b(alpha: &bam::record::Record, beta: &bam::record::Record) -> bool {
    // Take the read with the hightest mapq
    match alpha.mapq().cmp(&beta.mapq()) {
        Ordering::Less => false,
        Ordering::Greater => true,
        // Take the read with the lowest number of muli mappings
        Ordering::Equal => match alpha
            .aux(b"NH")
            .unwrap_or(Aux::Integer(0))
            .integer()
            .cmp(&beta.aux(b"NH").unwrap_or(Aux::Integer(0)).integer())
        {
            Ordering::Less => true,
            Ordering::Greater => false,
            // Take the read with the smallest edit distance
            Ordering::Equal => match alpha
                .aux(b"NM")
                .unwrap_or(Aux::Integer(0))
                .integer()
                .cmp(&beta.aux(b"NM").unwrap_or(Aux::Integer(0)).integer())
            {
                Ordering::Less => true,
                Ordering::Greater => false,
                // Take the read with the longer sequence
                Ordering::Equal => match alpha.seq().len().cmp(&beta.seq().len()) {
                    Ordering::Less => false,
                    Ordering::Greater => true,
                    // The incumbant wins if we get this far
                    Ordering::Equal => true,
                },
            },
        },
    }
}

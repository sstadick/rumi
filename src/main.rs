#[macro_use]
extern crate clap;
//extern crate rumi_dedup_lib;
use clap::{App, Arg};
use rumi_lib;
//use basebits::{hamming_dist, BaseBits};
//use rust_htslib::bam;
//use rust_htslib::prelude::*;
//use std::collections::hash_map::{Entry::Occupied, Entry::Vacant};
//use std::collections::HashMap;
use std::process;

fn main() {
    let matches = App::new("rumi")
        .version(crate_version!())
        .author(crate_authors!())
        .about("Deduplicate reads based on umis")
        .arg(
            Arg::with_name("INBAM")
                .help("Input bam file. Use - if stdin")
                .default_value("-")
                .required(true)
                .index(1),
        )
        .arg(
            Arg::with_name("OUTBAM")
                .short("o")
                .long("output")
                .help("Output bam file. Use - if stdout")
                .default_value("-")
                .required(true),
        )
        .arg(
            Arg::with_name("umi_tag")
                .short("u")
                .long("umi_tag")
                .help("The tag holding the umi information.")
                .default_value("RX")
                .required(true),
        )
        .arg(
            Arg::with_name("allowed_read_dist")
                .short("d")
                .long("allowed_read_dist")
                .help("The distance between umis that will allow them to be counted as adjacent.")
                .default_value("1")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("allowed_count_factor")
                .short("c")
                .long("allowed_count_factor")
                .help(
                    "The factor to multiply the count of a umi by when determining \
                     whether or not to group it with other umis within allowed_read_dist. \
                     include umi_b as adjacent to umi_a if: \
                     umi_a.counts >= allowed_count_factor * umi_b.counts",
                )
                .default_value("2")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("allowed_network_depth")
                .short("n")
                .long("allowed_network_depth")
                .help(
                    "The number of nodes deep to go when creating a group. If allowed_read_dist \
                     1, then allowed_network_depth of 2 will enable getting all umis with hamming \
                     distance of 2 from current umi.",
                )
                .default_value("2")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("umi_in_read_id")
                .long("umi_in_read_id")
                .help(
                "The UMI is located in the read id after the last '_'. Otherwise use the RX tag.",
            ),
        )
        .arg(
            Arg::with_name("ignore_splice_pos")
                .long("ignore_splice_pos")
                .help(
                    "If two reads have the same start pos, and contain a splice site, they will be
                     grouped together, instead of further splitting them based on the splice site",
                ),
        )
        .arg(Arg::with_name("group_only").long("group_only").help(
            "Don't deduplicate reads, just group them given them agroup id, and print them. Rules
                for filtering out unpaired reads, etc, will still be applied.",
        ))
        .arg(
            Arg::with_name("is_paired")
                .long("is_paired")
                .help("Input is paired end. Read pairs with unmapped read1 will be ignored."),
        )
        .get_matches();

    // Parse Args
    let config = rumi_lib::Config {
        input_bam: value_t!(matches, "INBAM", String).unwrap(),
        output_bam: value_t!(matches, "OUTBAM", String).unwrap(),
        umi_tag: value_t!(matches, "umi_tag", String).unwrap(),
        allowed_read_dist: value_t!(matches, "allowed_read_dist", u32).unwrap(),
        allowed_count_factor: value_t!(matches, "allowed_count_factor", u32).unwrap(),
        allowed_network_depth: value_t!(matches, "allowed_network_depth", usize).unwrap(),
        umi_in_read_id: matches.is_present("umi_in_read_id"),
        ignore_splice_pos: matches.is_present("ignore_splice_pos"),
        group_only: matches.is_present("group_only"),
        is_paired: matches.is_present("is_paired"),
    };

    if !config.group_only {
        if let Err(e) = rumi_lib::run_dedup(&config) {
            eprintln!("An error occured: {}", e);
            process::exit(1);
        }
    } else if let Err(e) = rumi_lib::run_group(&config) {
        eprintln!("An error occured: {}", e);
        process::exit(1);
    }
}

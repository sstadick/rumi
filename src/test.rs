use super::*;
use rust_htslib::bam;

macro_rules! btreemap {
    ( $b:expr; $($x:expr => $y:expr),* ) => ({
        let mut temp_map = BTreeMap::with_b($b);
        $(
            temp_map.insert($x, $y);
        )*
        temp_map
    });
    ( $($x:expr => $y:expr),* ) => ({
        let mut temp_map = BTreeMap::new();
        $(
            temp_map.insert($x, $y);
        )*
        temp_map
    });
    ( $b:expr; $($x:expr => $y:expr,)* ) => (
        btreemap!{$b; $($x => $y),*}
    );
    ( $($x:expr => $y:expr,)* ) => (
        btreemap!{$($x => $y),*}
    );
}

macro_rules! map(
    { $($key:expr => $value:expr),+ } => {
        {
            let mut m = ::std::collections::HashMap::new();
            $(
                m.insert($key, $value);
            )+
            m
        }
     };
);

fn get_header() -> bam::HeaderView {
    bam::HeaderView::from_bytes(b"@HD	VN:1.0	SO:coordinate
@SQ	SN:chr1	LN:197195432
@SQ	SN:chr10	LN:129993255
@SQ	SN:chr11	LN:121843856
@SQ	SN:chr12	LN:121257530
@SQ	SN:chr13	LN:120284312
@SQ	SN:chr14	LN:125194864
@SQ	SN:chr15	LN:103494974
@SQ	SN:chr16	LN:98319150
@SQ	SN:chr17	LN:95272651
@SQ	SN:chr18	LN:90772031
@SQ	SN:chr19	LN:61342430
@SQ	SN:chr2	LN:181748087
@SQ	SN:chr3	LN:159599783
@SQ	SN:chr4	LN:155630120
@SQ	SN:chr5	LN:152537259
@SQ	SN:chr6	LN:149517037
@SQ	SN:chr7	LN:152524553
@SQ	SN:chr8	LN:131738871
@SQ	SN:chr9	LN:124076172
@SQ	SN:chrM	LN:16299
@SQ	SN:chrX	LN:166650296
@SQ	SN:chrY	LN:15902555
@PG	ID:Bowtie	VN:1.1.2	CL:\"bowtie --wrapper basic-0 --threads 4 -v 2 -m 10 -a /ifs/mirror/genomes/bowtie/mm9 /dev/fd/63 --sam\"
")
}

fn check_readgroups(grouped: ReadMap, expected_group: ReadMap) {
    assert_eq!(grouped.keys().len(), expected_group.keys().len());
    assert_eq!(grouped.values().len(), expected_group.values().len());
    //println!("Test positions: ");
    //for t_position in grouped.keys() {
    //println!("{:#?}", t_position);
    //}
    //println!("Expected positions: ");
    for (e_position, e_umis) in expected_group {
        //println!("{:#?}", e_position);
        let t_umis = grouped.get(&e_position).unwrap();
        assert_eq!(e_umis.keys().len(), t_umis.keys().len());
        assert_eq!(e_umis.values().len(), t_umis.keys().len());
        for (e_umi, e_freq) in e_umis {
            let t_freq = t_umis.get(&e_umi).unwrap();
            assert_eq!(*t_freq, e_freq);
        }
    }
}

fn check_graph(graph: Vec<Node>, expected: Vec<Node>) {
    assert_eq!(graph.len(), expected.len());
    for (i, node) in graph.into_iter().enumerate() {
        assert_eq!(node, expected[i]);
    }
}

fn check_umi(grouping: Vec<&Node>, expected: Vec<&Node>) {
    assert_eq!(grouping.len(), expected.len());
    for (i, node) in grouping.into_iter().enumerate() {
        assert_eq!(node, expected[i]);
    }
}

#[test]
fn test_group_reads_small() {
    let header = get_header();
    let records_raw: Vec<&[u8]> = vec![
        b"SRR2057595.142416_TAGTA	0	chr19	61240266	255	26M	*	0	0	*	*	XA:i:1	MD:Z:12C13	NM:i:1	RX:Z:TAGTA",
        b"SRR2057595.297818_CAGTA	0	chr19	61240266	255	26M	*	0	0	*	*	XA:i:1	MD:Z:12C13	NM:i:1	RX:Z:CAGTA",
        b"SRR2057595.324156_CAGTA	0	chr19	61240266	255	26M	*	0	0	*	*	XA:i:1	MD:Z:12C13	NM:i:1	RX:Z:CAGTA",
        b"SRR2057595.357312_CAGTA	0	chr19	61240266	255	26M	*	0	0	*	*	XA:i:1	MD:Z:12C13	NM:i:1	RX:Z:CAGTA",
        b"SRR2057595.413242_CAGTA	0	chr19	61240266	255	26M	*	0	0	*	*	XA:i:1	MD:Z:12C13	NM:i:1	RX:Z:CAGTA",
        b"SRR2057595.509959_CAGTA	0	chr19	61240266	255	26M	*	0	0	*	*	XA:i:1	MD:Z:12C13	NM:i:1	RX:Z:CAGTA",
        b"SRR2057595.623861_CAGTA	0	chr19	61240266	255	26M	*	0	0	*	*	XA:i:1	MD:Z:12C13	NM:i:1	RX:Z:CAGTA",
    ];
    let records: Vec<bam::record::Record> = records_raw
        .iter()
        .map(|&r| bam::record::Record::from_sam(&header, r).unwrap())
        .collect();
    let expected_group: ReadMap = btreemap![
            Position {pos: 61240265, is_spliced: None, is_rev: false, target: 10, tlen: None} => map![
            BaseBits::new(b"CAGTA").unwrap() => ReadFreq {
                read: ReadCollection::SingleRead(bam::record::Record::from_sam(&header, records_raw[1]).unwrap()),
                freq: 6,
            },
            BaseBits::new(b"TAGTA").unwrap() => ReadFreq {
                read: ReadCollection::SingleRead(bam::record::Record::from_sam(&header, records_raw[0]).unwrap()),
                freq: 1,
            }
        ]
    ];

    let config = Config {
        input_bam: String::from("INPUT"),
        output_bam: String::from("OUTPUT"),
        umi_tag: String::from("RX"),
        allowed_read_dist: 1,
        allowed_count_factor: 2,
        allowed_network_depth: 2,
        umi_in_read_id: false,
        group_only: false,
        ignore_splice_pos: false,
        is_paired: false,
    };

    let (grouped, _) = group_reads(records, &config);

    // Check the read_groups
    check_readgroups(grouped, expected_group);
}

// Test the following:
// - A reverse mapped read is grouped on it's own
// - An unmapped read is not included
// - Different splice sites cause differnt groups to form
// - Different pos causes different groups to form
// - Different tid causes different groups to form
// - First highest mapq read is the representative read for a group
#[test]
fn test_group_reads_complex() {
    let header = get_header();
    let records_raw: Vec<&[u8]> = vec![
        b"SRR2057595.142416_TAGTA	0	chr19	61240266	255	26M	*	0	26	*	*	XA:i:1	MD:Z:12C13	NM:i:1	RX:Z:TAGTA",
        b"SRR2057595.297818_CAGTA	16	chr19	61240266	255	26M	*	0	26	*	*	XA:i:1	MD:Z:12C13	NM:i:1	RX:Z:CAGTA",
        b"SRR2057595.324156_CAGTA	0	chr19	61240266	255	26M	*	0	26	*	*	XA:i:1	MD:Z:12C13	NM:i:1	RX:Z:CAGTA",
        b"SRR2057595.357312_CAGTA	0	chr19	61240266	254	26M	*	0	26	*	*	XA:i:1	MD:Z:12C13	NM:i:1	RX:Z:CAGTA",
        b"SRR2057595.324245_CAGTA	0	chr19	61240266	255	26M	*	0	26	*	*	XA:i:1	MD:Z:12C13	NM:i:1	RX:Z:CAGTA",
        b"SRR2057595.413242_CAGTA	0	chr19	61240266	255	25M	*	0	25	*	*	XA:i:1	MD:Z:12C13	NM:i:1	RX:Z:CAGTA",
        b"SRR2057595.509959_CAGTA	0	chr18	61240266	255	26M	*	0	26	*	*	XA:i:1	MD:Z:12C13	NM:i:1	RX:Z:CAGTA",
        b"SRR2057595.623861_CAGTA	0	chr19	61240265	255	25M	*	0	25	*	*	XA:i:1	MD:Z:12C13	NM:i:1	RX:Z:CAGTA",
    ];
    let records: Vec<bam::record::Record> = records_raw
        .iter()
        .map(|&r| bam::record::Record::from_sam(&header, r).unwrap())
        .collect();
    let expected_group: ReadMap = btreemap![
            Position {pos: 61240265, is_spliced: None, is_rev: false, target: 10, tlen: None} => map![
                BaseBits::new(b"CAGTA").unwrap() => ReadFreq {
                    read: ReadCollection::SingleRead(bam::record::Record::from_sam(&header, records_raw[2]).unwrap()),
                    freq: 4,
                },
                BaseBits::new(b"TAGTA").unwrap() => ReadFreq {
                    read: ReadCollection::SingleRead(bam::record::Record::from_sam(&header, records_raw[0]).unwrap()),
                    freq: 1,
                }
        ],
            Position {pos: 61240291, is_spliced: None, is_rev: true, target: 10, tlen: None} => map![
                BaseBits::new(b"CAGTA").unwrap() => ReadFreq {
                    read: ReadCollection::SingleRead(bam::record::Record::from_sam(&header, records_raw[1]).unwrap()),
                    freq: 1,
                }
        ],
            Position {pos: 61240264, is_spliced: None, is_rev: false, target: 10, tlen: None} => map![
                BaseBits::new(b"CAGTA").unwrap() => ReadFreq {
                    read: ReadCollection::SingleRead(bam::record::Record::from_sam(&header, records_raw[7]).unwrap()),
                    freq: 1,
                }
        ],
            Position {pos: 61240265, is_spliced: None, is_rev: false, target: 9,tlen: None} => map![
                BaseBits::new(b"CAGTA").unwrap() => ReadFreq {
                    read: ReadCollection::SingleRead(bam::record::Record::from_sam(&header, records_raw[6]).unwrap()),
                    freq: 1,
                }
        ]
    ];

    let config = Config {
        input_bam: String::from("INPUT"),
        output_bam: String::from("OUTPUT"),
        umi_tag: String::from("RX"),
        allowed_read_dist: 1,
        allowed_count_factor: 2,
        allowed_network_depth: 2,
        umi_in_read_id: false,
        group_only: false,
        ignore_splice_pos: false,
        is_paired: false,
    };

    let (grouped, _) = group_reads(records, &config);

    // Check the read_groups
    check_readgroups(grouped, expected_group);
}

#[test]
fn test_read_groups_umi_tools() {
    let header = get_header();
    let records_raw: Vec<&[u8]> = vec![
        b"SRR2057595.11597812_ATAAA	16	chr19	4078297	255	38M	*	0	0	*	*	XA:i:1	MD:Z:29A8	NM:i:1	RX:Z:ATAAA	UG:i:52	BX:Z:ATAAA",
        b"SRR2057595.10788_ATAAA	16	chr19	4078298	255	37M	*	0	0	*	*	XA:i:0	MD:Z:37	NM:i:0	RX:Z:ATAAA	UG:i:52	BX:Z:ATAAA",
        b"SRR2057595.42646_ATAAA	16	chr19	4078298	255	37M	*	0	0	*	*	XA:i:0	MD:Z:37	NM:i:0	RX:Z:ATAAA	UG:i:52	BX:Z:ATAAA",
    ];
    let records: Vec<bam::record::Record> = records_raw
        .iter()
        .map(|&r| bam::record::Record::from_sam(&header, r).unwrap())
        .collect();
    let expected_group: ReadMap = btreemap![
            Position {pos: 4078334, is_spliced: None, is_rev: true, target: 10, tlen: None} => map![
                BaseBits::new(b"ATAAA").unwrap() => ReadFreq {
                    read: ReadCollection::SingleRead(bam::record::Record::from_sam(&header, records_raw[1]).unwrap()),
                    freq: 3,
                }
        ]
    ];

    let config = Config {
        input_bam: String::from("INPUT"),
        output_bam: String::from("OUTPUT"),
        umi_tag: String::from("RX"),
        allowed_read_dist: 1,
        allowed_count_factor: 2,
        allowed_network_depth: 2,
        umi_in_read_id: false,
        group_only: false,
        ignore_splice_pos: false,
        is_paired: false,
    };

    let (grouped, _) = group_reads(records, &config);

    // Check the read_groups
    check_readgroups(grouped, expected_group);
}

#[test]
fn test_read_groups_cigars() {
    let header = get_header();
    let records_raw: Vec<&[u8]> = vec![
        b"SRR2057595.11597812_ATAAA	16	chr19	4078297	255	3S35M	*	0	0	*	*	XA:i:1	MD:Z:29A8	NM:i:1	RX:Z:ATAAA	UG:i:52	BX:Z:ATAAA",
        b"SRR2057595.10788_ATAAA	16	chr19	4078294	255	33M4S	*	0	0	*	*	XA:i:0	MD:Z:37	NM:i:0	RX:Z:ATAAA	UG:i:52	BX:Z:ATAAA",
        b"SRR2057595.42646_ATAAA	16	chr19	4078298	255	15M7N15M	*	0	0	*	*	XA:i:0	MD:Z:37	NM:i:0	RX:Z:ATAAA	UG:i:52	BX:Z:ATAAA",
        b"SRR2057595.11597790_ATAAA	0	chr19	4078300	255	3S35M	*	0	0	*	*	XA:i:1	MD:Z:29A8	NM:i:1	RX:Z:ATAAA	UG:i:52	BX:Z:ATAAA",
        b"SRR2057595.10988_ATAAA	0	chr19	4078298	255	33M4S	*	0	0	*	*	XA:i:0	MD:Z:37	NM:i:0	RX:Z:ATAAA	UG:i:52	BX:Z:ATAAA",
        b"SRR2057595.4246_ATAAA	0	chr19	4078298	255	15M7N15M	*	0	0	*	*	XA:i:0	MD:Z:37	NM:i:0	RX:Z:ATAAA	UG:i:52	BX:Z:ATAAA",
    ];
    let records: Vec<bam::record::Record> = records_raw
        .iter()
        .map(|&r| bam::record::Record::from_sam(&header, r).unwrap())
        .collect();
    let expected_group: ReadMap = btreemap![
            Position {pos: 4078330, is_spliced: None, is_rev: true, target: 10, tlen: None} => map![
                BaseBits::new(b"ATAAA").unwrap() => ReadFreq {
                    read: ReadCollection::SingleRead(bam::record::Record::from_sam(&header, records_raw[1]).unwrap()),
                    freq: 1,
                }],
            Position {pos: 4078331, is_spliced: Some(35), is_rev: true, target: 10,tlen: None} => map![
                BaseBits::new(b"ATAAA").unwrap() => ReadFreq {
                    read: ReadCollection::SingleRead(bam::record::Record::from_sam(&header, records_raw[0]).unwrap()),
                    freq: 1,
                }],
            Position {pos: 4078334, is_spliced: Some(15), is_rev: true, target: 10,tlen: None} => map![
                BaseBits::new(b"ATAAA").unwrap() => ReadFreq {
                    read: ReadCollection::SingleRead(bam::record::Record::from_sam(&header, records_raw[2]).unwrap()),
                    freq: 1,
                }],
            Position {pos: 4078296, is_spliced: None, is_rev: false, target: 10,tlen: None} => map![
                BaseBits::new(b"ATAAA").unwrap() => ReadFreq {
                    read: ReadCollection::SingleRead(bam::record::Record::from_sam(&header, records_raw[3]).unwrap()),
                    freq: 1,
                }],
            Position {pos: 4078297, is_spliced: Some(33), is_rev: false, target: 10,tlen: None} => map![
                BaseBits::new(b"ATAAA").unwrap() => ReadFreq {
                    read: ReadCollection::SingleRead(bam::record::Record::from_sam(&header, records_raw[4]).unwrap()),
                    freq: 1,
                }],
            Position {pos: 4078297, is_spliced: Some(15), is_rev: false, target: 10,tlen: None} => map![
                BaseBits::new(b"ATAAA").unwrap() => ReadFreq {
                    read: ReadCollection::SingleRead(bam::record::Record::from_sam(&header, records_raw[5]).unwrap()),
                    freq: 1,
                }]
    ];

    let config = Config {
        input_bam: String::from("INPUT"),
        output_bam: String::from("OUTPUT"),
        umi_tag: String::from("RX"),
        allowed_read_dist: 1,
        allowed_count_factor: 2,
        allowed_network_depth: 2,
        umi_in_read_id: false,
        group_only: false,
        ignore_splice_pos: false,
        is_paired: false,
    };

    let (grouped, _) = group_reads(records, &config);

    // Check the read_groups
    check_readgroups(grouped, expected_group);
}

#[test]
fn test_graph_small() {
    let header = get_header();
    let uncon_graph = vec![
        Node {
            umi: BaseBits::new(b"CAGTA").unwrap(),
            freq: ReadFreq {
            read: ReadCollection::SingleRead(bam::record::Record::from_sam(&header,
                b"SRR2057595.297818_CAGTA	0	chr19	61240266	255	26M	*	0	0	*	*	XA:i:1	MD:Z:12C13	NM:i:1	RX:Z:CAGTA").unwrap()),
                freq: 6,
            },
            connections: vec![],
        },
        Node {
            umi: BaseBits::new(b"TAGTA").unwrap(),
            freq: ReadFreq {
            read: ReadCollection::SingleRead(bam::record::Record::from_sam(&header,
                b"SRR2057595.142416_TAGTA	0	chr19	61240266	255	26M	*	0	0	*	*	XA:i:1	MD:Z:12C13	NM:i:1	RX:Z:TAGTA").unwrap()),
                freq: 1,
            },
            connections: vec![],
        },
    ];

    let expected = vec![
        Node {
            umi: BaseBits::new(b"CAGTA").unwrap(),
            freq: ReadFreq {
            read: ReadCollection::SingleRead(bam::record::Record::from_sam(&header,
                b"SRR2057595.297818_CAGTA	0	chr19	61240266	255	26M	*	0	0	*	*	XA:i:1	MD:Z:12C13	NM:i:1	RX:Z:CAGTA").unwrap()),
                freq: 6,
            },
            connections: vec![1],
        },
        Node {
            umi: BaseBits::new(b"TAGTA").unwrap(),
            freq: ReadFreq {
            read: ReadCollection::SingleRead(bam::record::Record::from_sam(&header,
                b"SRR2057595.142416_TAGTA	0	chr19	61240266	255	26M	*	0	0	*	*	XA:i:1	MD:Z:12C13	NM:i:1	RX:Z:TAGTA").unwrap()),
                freq: 1,
            },
            connections: vec![],
        },
    ];

    let config = Config {
        input_bam: String::from("INPUT"),
        output_bam: String::from("OUTPUT"),
        umi_tag: String::from("RX"),
        allowed_read_dist: 1,
        allowed_count_factor: 2,
        allowed_network_depth: 2,
        umi_in_read_id: false,
        group_only: false,
        ignore_splice_pos: false,
        is_paired: false,
    };

    let graph = connect_graph(
        uncon_graph,
        config.allowed_read_dist,
        config.allowed_count_factor,
    );
    println!("{:#?}", graph);
    check_graph(graph, expected);
}

#[test]
fn test_graph_umi() {
    // Test the umis found in the umi blog post:
    // https://cgatoxford.files.wordpress.com/2015/08/schematic_25-e1443714121688.png
    let header = get_header();
    let uncon_graph = vec![
        Node {
            umi: BaseBits::new(b"ATTG").unwrap(),
            freq: ReadFreq {
                read: ReadCollection::SingleRead(bam::record::Record::from_sam(&header,
                    b"SRR2057595.297818_ATTG	0	chr19	61240266	255	26M	*	0	0	*	*	XA:i:1	MD:Z:12C13	NM:i:1	RX:Z:ATTG").unwrap()),
                freq: 1,
            },
            connections: vec![],
        },
        Node {
            umi: BaseBits::new(b"ATTA").unwrap(),
            freq: ReadFreq {
                read: ReadCollection::SingleRead(bam::record::Record::from_sam(&header,
                    b"SRR2057595.297818_ATTA	0	chr19	61240266	255	26M	*	0	0	*	*	XA:i:1	MD:Z:12C13	NM:i:1	RX:Z:ATTA").unwrap()),
                freq: 456,
            },
            connections: vec![],
        },
        Node {
            umi: BaseBits::new(b"ATTT").unwrap(),
            freq: ReadFreq {
                read: ReadCollection::SingleRead(bam::record::Record::from_sam(&header,
                    b"SRR2057595.297818_ATTT	0	chr19	61240266	255	26M	*	0	0	*	*	XA:i:1	MD:Z:12C13	NM:i:1	RX:Z:ATTT").unwrap()),
                freq: 2,
            },
            connections: vec![],
        },
        Node {
            umi: BaseBits::new(b"AGTA").unwrap(),
            freq: ReadFreq{
                read: ReadCollection::SingleRead(bam::record::Record::from_sam(&header,
                    b"SRR2057595.297818_AGTA	0	chr19	61240266	255	26M	*	0	0	*	*	XA:i:1	MD:Z:12C13	NM:i:1	RX:Z:AGTA").unwrap()),
                freq: 72,
            },
            connections: vec![],
        },
        Node {
            umi: BaseBits::new(b"AGTC").unwrap(),
            freq: ReadFreq {
                read: ReadCollection::SingleRead(bam::record::Record::from_sam(&header,
                    b"SRR2057595.297818_AGTC	0	chr19	61240266	255	26M	*	0	0	*	*	XA:i:1	MD:Z:12C13	NM:i:1	RX:Z:AGTC").unwrap()),
                freq: 1,
            },
            connections: vec![],
        },
        Node {
            umi: BaseBits::new(b"AGGA").unwrap(),
            freq: ReadFreq {
                read: ReadCollection::SingleRead(bam::record::Record::from_sam(&header,
                    b"SRR2057595.142416_AGGA	0	chr19	61240266	255	26M	*	0	0	*	*	XA:i:1	MD:Z:12C13	NM:i:1	RX:Z:AGGA").unwrap()),
                freq: 90,
            },
            connections: vec![],
        },
    ];
    let expected = vec![
        Node {
            umi: BaseBits::new(b"ATTG").unwrap(),
            freq: ReadFreq {
                read: ReadCollection::SingleRead(bam::record::Record::from_sam(&header,
                    b"SRR2057595.297818_ATTG	0	chr19	61240266	255	26M	*	0	0	*	*	XA:i:1	MD:Z:12C13	NM:i:1	RX:Z:ATTG").unwrap()),
                freq: 1,
            },
            connections: vec![],
        },
        Node {
            umi: BaseBits::new(b"ATTA").unwrap(),
            freq: ReadFreq {
                read: ReadCollection::SingleRead(bam::record::Record::from_sam(&header,
                    b"SRR2057595.297818_ATTA	0	chr19	61240266	255	26M	*	0	0	*	*	XA:i:1	MD:Z:12C13	NM:i:1	RX:Z:ATTA").unwrap()),
                freq: 456,
            },
            connections: vec![0, 2, 3],
        },
        Node {
            umi: BaseBits::new(b"ATTT").unwrap(),
            freq: ReadFreq {
                read: ReadCollection::SingleRead(bam::record::Record::from_sam(&header,
                    b"SRR2057595.297818_ATTT	0	chr19	61240266	255	26M	*	0	0	*	*	XA:i:1	MD:Z:12C13	NM:i:1	RX:Z:ATTT").unwrap()),
                freq: 2,
            },
            connections: vec![0],
        },
        Node {
            umi: BaseBits::new(b"AGTA").unwrap(),
            freq: ReadFreq{
                read: ReadCollection::SingleRead(bam::record::Record::from_sam(&header,
                    b"SRR2057595.297818_AGTA	0	chr19	61240266	255	26M	*	0	0	*	*	XA:i:1	MD:Z:12C13	NM:i:1	RX:Z:AGTA").unwrap()),
                freq: 72,
            },
            connections: vec![4],
        },
        Node {
            umi: BaseBits::new(b"AGTC").unwrap(),
            freq: ReadFreq {
                read: ReadCollection::SingleRead(bam::record::Record::from_sam(&header,
                    b"SRR2057595.297818_AGTC	0	chr19	61240266	255	26M	*	0	0	*	*	XA:i:1	MD:Z:12C13	NM:i:1	RX:Z:AGTC").unwrap()),
                freq: 1,
            },
            connections: vec![],
        },
        Node {
            umi: BaseBits::new(b"AGGA").unwrap(),
            freq: ReadFreq {
                read: ReadCollection::SingleRead(bam::record::Record::from_sam(&header,
                    b"SRR2057595.142416_AGGA	0	chr19	61240266	255	26M	*	0	0	*	*	XA:i:1	MD:Z:12C13	NM:i:1	RX:Z:AGGA").unwrap()),
                freq: 90,
            },
            connections: vec![],
        },
    ];

    let config = Config {
        input_bam: String::from("INPUT"),
        output_bam: String::from("OUTPUT"),
        umi_tag: String::from("RX"),
        allowed_read_dist: 1,
        allowed_count_factor: 2,
        allowed_network_depth: 2,
        umi_in_read_id: false,
        group_only: false,
        ignore_splice_pos: false,
        is_paired: false,
    };

    let graph = connect_graph(
        uncon_graph,
        config.allowed_read_dist,
        config.allowed_count_factor,
    );
    println!("{:#?}", graph);
    check_graph(graph, expected);
}

#[test]
fn test_determine_umi() {
    // Test the umis found in the umi blog post:
    // https://cgatoxford.files.wordpress.com/2015/08/schematic_25-e1443714121688.png
    let header = get_header();
    let graph = vec![
        Node {
            umi: BaseBits::new(b"ATTG").unwrap(),
            freq: ReadFreq {
                read: ReadCollection::SingleRead(bam::record::Record::from_sam(&header,
                    b"SRR2057595.297818_ATTG	0	chr19	61240266	255	26M	*	0	0	*	*	XA:i:1	MD:Z:12C13	NM:i:1	RX:Z:ATTG").unwrap()),
                freq: 1,
            },
            connections: vec![],
        },
        Node {
            umi: BaseBits::new(b"ATTA").unwrap(),
            freq: ReadFreq {
                read: ReadCollection::SingleRead(bam::record::Record::from_sam(&header,
                    b"SRR2057595.297818_ATTA	0	chr19	61240266	255	26M	*	0	0	*	*	XA:i:1	MD:Z:12C13	NM:i:1	RX:Z:ATTA").unwrap()),
                freq: 456,
            },
            connections: vec![0, 2, 3],
        },
        Node {
            umi: BaseBits::new(b"ATTT").unwrap(),
            freq: ReadFreq {
                read: ReadCollection::SingleRead(bam::record::Record::from_sam(&header,
                    b"SRR2057595.297818_ATTT	0	chr19	61240266	255	26M	*	0	0	*	*	XA:i:1	MD:Z:12C13	NM:i:1	RX:Z:ATTT").unwrap()),
                freq: 2,
            },
            connections: vec![0],
        },
        Node {
            umi: BaseBits::new(b"AGTA").unwrap(),
            freq: ReadFreq{
                read: ReadCollection::SingleRead(bam::record::Record::from_sam(&header,
                    b"SRR2057595.297818_AGTA	0	chr19	61240266	255	26M	*	0	0	*	*	XA:i:1	MD:Z:12C13	NM:i:1	RX:Z:AGTA").unwrap()),
                freq: 72,
            },
            connections: vec![4, 5],
        },
        Node {
            umi: BaseBits::new(b"AGTC").unwrap(),
            freq: ReadFreq {
                read: ReadCollection::SingleRead(bam::record::Record::from_sam(&header,
                    b"SRR2057595.297818_AGTC	0	chr19	61240266	255	26M	*	0	0	*	*	XA:i:1	MD:Z:12C13	NM:i:1	RX:Z:AGTC").unwrap()),
                freq: 1,
            },
            connections: vec![],
        },
        Node {
            umi: BaseBits::new(b"AGGA").unwrap(),
            freq: ReadFreq {
                read: ReadCollection::SingleRead(bam::record::Record::from_sam(&header,
                    b"SRR2057595.297818_AGTG	0	chr19	61240266	255	26M	*	0	0	*	*	XA:i:1	MD:Z:12C13	NM:i:1	RX:Z:AGTG").unwrap()),
                freq: 5,
            },
            connections: vec![],
        },
        Node {
            umi: BaseBits::new(b"AGGA").unwrap(),
            freq: ReadFreq {
                read: ReadCollection::SingleRead(bam::record::Record::from_sam(&header,
                    b"SRR2057595.142416_AGGA	0	chr19	61240266	255	26M	*	0	0	*	*	XA:i:1	MD:Z:12C13	NM:i:1	RX:Z:AGGA").unwrap()),
                freq: 90,
            },
            connections: vec![5],
        },
    ];

    let node1 = Node {
        umi: BaseBits::new(b"ATTA").unwrap(),
        freq: ReadFreq {
            read: ReadCollection::SingleRead(bam::record::Record::from_sam(
                &header,
                b"SRR2057595.297818_ATTA	0	chr19	61240266	255	26M	*	0	0	*	*	XA:i:1	MD:Z:12C13	NM:i:1	RX:Z:ATTA",
            )
            .unwrap()),
            freq: 456,
        },
        connections: vec![0, 2, 3],
    };
    let node2 = Node {
        umi: BaseBits::new(b"AGGA").unwrap(),
        freq: ReadFreq {
            read: ReadCollection::SingleRead(bam::record::Record::from_sam(
                &header,
                b"SRR2057595.142416_AGGA	0	chr19	61240266	255	26M	*	0	0	*	*	XA:i:1	MD:Z:12C13	NM:i:1	RX:Z:AGGA",
            )
            .unwrap()),
            freq: 90,
        },
        connections: vec![5],
    };
    let expected = vec![&node1, &node2];

    let config = Config {
        input_bam: String::from("INPUT"),
        output_bam: String::from("OUTPUT"),
        umi_tag: String::from("RX"),
        allowed_read_dist: 1,
        allowed_count_factor: 2,
        allowed_network_depth: 2,
        umi_in_read_id: false,
        group_only: false,
        ignore_splice_pos: false,
        is_paired: false,
    };

    let grouping = determine_umi(&graph, config.allowed_network_depth);
    // Test that nodes can't be double added to two different groups. Otherwise the group_only ends
    // up printing them twice
    assert_eq!(grouping[0].nodes.len(), 6);
    assert_eq!(grouping[1].nodes.len(), 1);
    let grouping: Vec<&Node> = grouping.iter().map(|n| n.nodes[n.master_node]).collect();
    check_umi(grouping, expected);
}

use strict;
use warnings;
use Data::Dumper;

my $expected_file = shift;
my $found_file = shift;
my $expected_group_file = shift;
my $found_group_file = shift;

open(my $expected_fh, "<", $expected_file) or die "Couldn't open $expected_file";
open(my $found_fh, "<", $found_file) or die "Couldn't open $found_file";
open(my $expected_groups_fh, "<", $expected_group_file) or die "Couldn't open $expected_group_file";
open(my $found_groups_fh, "<", $found_group_file) or die "Couldn't open $found_group_file";
my $expected_reads = build_read_maps($expected_fh, "id");
my $found_reads = build_read_maps($found_fh, "id");
my $expected_with_groups = build_read_maps($expected_groups_fh, "qname");
my $found_with_groups = build_read_maps($found_groups_fh, "qname");
close($expected_groups_fh);
close($found_groups_fh);
close($expected_fh);
close($found_fh);

compare_reads($expected_reads, $found_reads, $expected_with_groups, $found_with_groups);

sub ham {
	my ($astr, $bstr) = @_;
	my @a = split //, $astr;
	my @b = split //, $bstr;
	my $dist = 0;
	for my $i (0..(scalar(@a) - 1)) {
		$dist ++ if ($a[$i] ne $b[$i]);
	}
	return $dist;
}
sub compare_reads {
	my ($e_reads, $f_reads, $e_w_groups, $f_w_groups) = @_;

	for my $f_id (keys %$f_reads) {
		if (! exists $e_reads->{$f_id}) {
			my $found_read = $found_with_groups->{$f_reads->{$f_id}->{qname}};
			my $expected_read = $expected_with_groups->{$f_reads->{$f_id}->{qname}};
			#print $f_id . "\t" . $f_reads->{$f_id}->{read} . "\n";
			if (ham($expected_read->{nametag}, $expected_read->{BX}) > 2) {
				#print "FOUND:\t " . $found_read->{read} . "\n";
				print STDERR "HAMMMMMM\n";
				#print "EXPECTED:\t " . $expected_read->{read} . "\n";
			} else {
				print "FOUND:\t " . $found_read->{read} . "\n";
				print "EXPECTED ($f_id):\t " . $expected_read->{read} . "\n";
			}
		}
	}
}

sub split_read {
	my ($read) = @_;
	my @vars = split(/\t/, $read);
	my $hash = {
		qname => $vars[0],
		flag => $vars[1],
		rname => $vars[2],
		pos => $vars[3],
		mapq => $vars[4],
		cigar => $vars[5],
		rnext => $vars[6],
		pnext => $vars[7],
		tlen => $vars[8],
		seq => $vars[9],
		qual => $vars[10],
		read => $read,
	};

	my ($name, $nametag) = split(/_/, $hash->{qname});
	$hash->{nametag} = $nametag;
	foreach my $i (11.. scalar(@vars) - 1 ) {
		my ($tag_name, $tag_type, $tag_value) = split(/:/, $vars[$i]);
		$hash->{$tag_name} = $tag_value;
	}
	if ($hash->{flag} == 16) {
		my ($len) = $hash->{cigar} =~ /(\d+)M/;
		my $pos = $hash->{pos} + $len;
		$hash->{id} = join("$;", ($hash->{rname}, $pos, $hash->{nametag}));
	} else {
		$hash->{id} = join("$;", ($hash->{rname}, $hash->{pos}, $hash->{nametag}));
	}
	return $hash;
}

sub build_read_maps{
	my ($fh, $key) = @_;
	my %reads;
	while (my $read = <$fh>) {
		chomp($read);
		my $read_map = split_read($read);	
		$reads{$read_map->{$key}} = $read_map;
	}
	\%reads;
}

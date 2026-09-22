#!/usr/bin/perl

package GPBpath;

use strict;
use warnings;

sub parse_path_name {
	my ($raw, $opts) = @_;
	$opts ||= {};
	my $delim = defined $opts->{delim} ? $opts->{delim} : '.';
	my $coord_base = defined $opts->{coord_base} ? $opts->{coord_base} : 0;

	my %rec = (
		raw        => $raw,
		name       => $raw,
		sample     => undef,
		hap        => undef,
		locus      => undef,
		start      => undef,
		end        => undef,
		coord_base => $coord_base,
		style      => 'raw',
	);

	my $name = $raw;

	if ($name =~ s/#\d+\[(\d+)-(\d+)\]$//) {
		($rec{start}, $rec{end}) = ($1, $2);
	} elsif ($name =~ s/:(\d+)-(\d+)$//) {
		($rec{start}, $rec{end}) = ($1, $2);
	}
	$rec{name} = $name;

	my @f = split /#/, $name, -1;
	if (@f >= 3) {
		$rec{sample} = $f[0];
		$rec{hap}    = $f[1];
		$rec{locus}  = $f[2];
		$rec{style}  = 'pansn';
	} elsif (@f == 2) {
		$rec{sample} = $f[0];
		$rec{locus}  = $f[1];
		$rec{style}  = 'pansn2';
	} elsif ($name =~ /^(.*?)\Q$delim\E((?i:chr)\w*)$/) {
		$rec{sample} = $1;
		$rec{locus}  = $2;
		$rec{style}  = 'delim';
	} else {
		$rec{sample} = $name;
		$rec{style}  = 'raw';
	}

	return \%rec;
}


sub path_identity {
	my ($rec) = @_;
	return join('#', map { defined $_ ? $_ : '' } @$rec{qw(sample hap locus)});
}


sub coord_1based {
	my ($rec) = @_;
	return (undef, undef) unless defined $rec->{start};
	my $start = $rec->{coord_base} == 0 ? $rec->{start} + 1 : $rec->{start};
	return ($start, $rec->{end} + 0);
}


sub resolve_ref_path {
	my ($query, $path_names, $opts) = @_;
	$opts ||= {};

	my $q = parse_path_name($query, $opts);
	my @parsed = map { parse_path_name($_, $opts) } @$path_names;

	my @cand;
	for my $p (@parsed) {
		if ($p->{raw} eq $query) {
			@cand = ($p);
			last;
		}
	}

	if (!@cand && defined $q->{locus}) {
		my $qid = path_identity($q);
		@cand = grep { path_identity($_) eq $qid } @parsed;
	}

	if (!@cand && $q->{style} eq 'raw') {
		@cand = grep { defined $_->{sample} && $_->{sample} eq $q->{sample} } @parsed;
	}

	if (!@cand && $q->{style} eq 'raw') {
		@cand = grep { defined $_->{locus} && $_->{locus} eq $q->{raw} } @parsed;
	}

	if (!@cand) {
		@cand = grep { index($_->{raw}, $query) != -1 } @parsed;
		if (@cand) {
			warn "Warning: --ref-name '$query' matched by substring only: " . join(", ", map { $_->{raw} } @cand) . "\n";
		}
	}

	unless (@cand) {
		my $avail = join("\n    ", map { $_->{raw} } @parsed);
		die "Error: No reference path matches --ref-name '$query'.\n"
		  . "  Available paths (" . scalar(@parsed) . "):\n    $avail\n";
	}
	if (@cand > 1) {
		my $list = join("\n    ", map { $_->{raw} } @cand);
		die "Error: --ref-name '$query' matches " . scalar(@cand) . " paths, but the reference must be a single path:\n"
		  . "    $list\n"
		  . "  Please rerun with the exact full path name.\n";
	}

	my $ref = $cand[0];
	if ($opts->{forbid_coords} && defined $ref->{start}) {
		die "Error: The reference path '$ref->{raw}' carries a coordinate suffix, so it is not a full-length path. Subgraph extraction requires the reference to be a full-length path.\n";
	}
	if ($opts->{require_locus} && !defined $ref->{locus}) {
		die "Error: The reference path '$ref->{raw}' carries no chromosome/locus information, so it cannot be matched against --region or gene annotation coordinates.\n"
		  . "  Subgraph extraction requires a path name like 'P1#0#chr1' or 'P1.chr1'. To visualize this graph directly without extraction, use '--no-extract'.\n";
	}

	return $ref;
}


1;

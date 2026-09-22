#!/usr/bin/perl

package GPBpheno;

use strict;
use warnings;
use List::Util qw(first);


sub _sample_key {
	my ($name) = @_;
	return $1 if $name =~ /^([^#]+)#+/;
	return $1 if $name =~ /^(.*?)(?=\.(?i:chr))/;
	return $name;
}


sub load_pheno {
	my ($pheno_file, $nodes, $MISSING_RE) = @_;
	
	my %node_sample;
	$node_sample{$_} = 1 for split /,/, join ',', map{ $_->[2] } grep { defined $_->[2] && $_->[2] ne '' } @$nodes;
	my @ks = keys %node_sample;
	$node_sample{_sample_key($_)} = 1 for @ks;

	open my $fh_ph, '<', $pheno_file or die "Error: Can't open file '$pheno_file': $!\n";
	my $header = <$fh_ph>;
	$header =~ s/\r?\n$//;
	my @pname = split /\t/, $header;
	$_ =~ s/^"|"$//g for @pname;
	shift @pname;
	
	my %is_cont = map { $_ => 1 } @pname;
	my (%phen, %level_list);
	while (<$fh_ph>) {
		s/\r?\n$//;
		next if /^#/;
		next if /^\s*$/;
		my @arr = split /\t/;
		$_ =~ s/^"|"$//g for @arr;

		my $path = shift @arr;
		next unless $node_sample{$path};
		for my $i (0..$#pname) {
			my $v = $arr[$i];
			$phen{$path}{ $pname[$i] } = $v;
			next unless $is_cont{$pname[$i]};
			next if !defined $v || $v eq '' || $v =~ $MISSING_RE;

			if ($is_cont{$pname[$i]} && $v !~ /^-?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?$/i) {
				$is_cont{$pname[$i]} = 0;
			}
		}
	}

	seek $fh_ph, 0, 0;

	<$fh_ph>; 
	while (<$fh_ph>) {
		s/\r?\n$//;
		next if /^#/;
		next if /^\s*$/;
		my @arr = split /\t/;
		$_ =~ s/^"|"$//g for @arr;

		my $path = shift @arr;
                next unless $node_sample{$path};

		 for my $i (0..$#pname) {
		 	my $v = $arr[$i];
			next if !defined $v || $v eq '' || $v =~ $MISSING_RE;
			unless ($is_cont{$pname[$i]}) {
				$level_list{$pname[$i]}{$v} = 1;
			}
		 }
	}

	close $fh_ph;	

	$level_list{$_} = [ sort keys %{ $level_list{$_} } ] for grep { !$is_cont{$_} } @pname;
	
	return (\%phen, \%is_cont, \@pname, \%level_list);
}


sub add_pheno {
	my ($nodes, $pheno_file) = @_;

	my $MISSING_RE = qr/^(?:-|NA|NaN|null|None|NULL|undefined|unknown|Unknown)$/i;

	my ($pheno_data, $is_cont, $pname, $level_list) = load_pheno($pheno_file, $nodes, $MISSING_RE);
	my $pheno_meta = [ map { [ $_, $is_cont->{$_}, $level_list->{$_} ] } @$pname];
	
	foreach my $node (@$nodes) {
		my @samples = defined $node->[2] ? split /,/, $node->[2] : ();
		my %seen;
		my @sta_uniq;

		for my $i (0..$#$pheno_meta) {
			my ($pheno_name, $cont, $levels) = @{$pheno_meta->[$i]};
			if ($cont) {
				$sta_uniq[$i] = [0, (0)x5, []];
			} else {
				my @empty = (0) x @$levels;
				$sta_uniq[$i] = [0, @empty];
			}
		}

		my @uniq_samples = grep { !$seen{_sample_key($_)}++ } @samples;

		foreach my $s (@uniq_samples) {
			foreach my $i (0..$#$pheno_meta) {
				my ($pheno_name, $cont, $levels) = @{$pheno_meta->[$i]};
				my $val = $pheno_data->{$s}{$pheno_name} // $pheno_data->{_sample_key($s)}{$pheno_name};
				if (!defined $val || $val eq '' || $val =~ $MISSING_RE) {
					$sta_uniq[$i][0]++;
					next;
				}

				if ($cont) {
					push @{$sta_uniq[$i][6]}, $val;
				} else {
					my $idx = (first { $levels->[$_] eq $val } 0 .. $#$levels) // -1;
					next if $idx < 0;
					$sta_uniq[$i][$idx+1]++;		
				}
			}
		}

		for my $i (0..$#$pheno_meta) {
			next unless $pheno_meta->[$i][1];
			my $arr = $sta_uniq[$i][6];
			my $len     = defined $arr ? @$arr : 0;
			$sta_uniq[$i][6] = $len;
			next unless $len > 0;
			my @v = sort { $a <=> $b } @$arr;
			my $n = @v or next;
			$sta_uniq[$i][1] = $v[0] + 0;
			$sta_uniq[$i][2] = _quantile(\@v, 0.25) + 0;
			$sta_uniq[$i][3] = _median(\@v) + 0;
			$sta_uniq[$i][4] = _quantile(\@v, 0.75) + 0;
			$sta_uniq[$i][5] = $v[-1] + 0;   
		}
		push @$node, \@sta_uniq;
	}
	
	return ($nodes, $pheno_meta);
}


sub _median {
    my $arr = shift;
    return _quantile($arr, 0.5);
}


sub _quantile {
    my ($arr, $p) = @_;
    my $n = @$arr;
    my $pos = $p * ($n - 1);
    my $i = int($pos);
    my $w = $pos - $i;
    return $arr->[$i] if $i + 1 >= $n;
    return (1 - $w) * $arr->[$i] + $w * $arr->[$i + 1];
}

1;


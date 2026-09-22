#!/usr/bin/perl

package GPBanno;

use strict;
use warnings;
use List::Util qw(max);

sub get_interval {

        my ($geneid, $genelist, $region, $regionlist, $extend, $anno, $dir_name, $ref_path_name) = @_;

        if ($geneid || $genelist) {
                my @genes;
                if ($geneid) {
                        @genes = ($geneid);
                } else {
                        open my $fh, '<', $genelist or die $!;
                        @genes = map { s/\r?\n$//; $_ } grep { !/^#/ && /\S/ } <$fh>;
                        close $fh;
                        die "Error: No valid gene IDs found in $genelist\n" unless @genes;
                }
                my %seen_gene;
                @genes = grep { !$seen_gene{$_}++ } @genes;
		my $format = detect_anno_format($anno);
                my $interval_arr = extract_gene_annos($anno, \@genes, $dir_name, $format, $extend);
		die "Error: No annotations found for the provided gene ID(s) in $anno.\n" unless @$interval_arr;
                return $interval_arr;
        } elsif ($region || $regionlist) {
                my @regions;
                if ($region) {
                        my ($c, $s, $e) = parse_region_string($region);
                        @regions = (["${c}_${s}-${e}", $c, $s, $e]);
                } else {
                        open my $fh, '<', $regionlist or die "Error: Cannot open $regionlist: $!";
                        my %seen_region;
                        @regions = ();
                        while (<$fh>) {
                                s/\r?\n$//;
                                next if /^#/ || /^\s*$/;
                                my @f = split(/\t/);
                                next unless @f >= 3;
                                my ($chr, $start, $end) = ($f[0], int($f[1]) + 1, int($f[2]));
                                my $key = join("\t", $chr, $start, $end);
                                next if $seen_region{$key};
                                $seen_region{$key} = 1;
                                push @regions, ["${chr}_${start}-${end}", $chr, $start, $end];
                        }
                        close $fh;
                        die "Error: No valid regions found in $regionlist\n" unless @regions;
                }
		if(defined $anno){
			my $format = detect_anno_format($anno);
                	extract_region_annos($anno, \@regions, $dir_name, $format);
		}
                return \@regions;
        } elsif (defined $ref_path_name) {
		my ($c, $s, $e) = parse_region_string($ref_path_name);
		my @regions = (["${c}_${s}-${e}", $c, $s, $e]);
		if(defined $anno){
                        my $format = detect_anno_format($anno);
                        extract_region_annos($anno, \@regions, $dir_name, $format);
                }
		return \@regions;		
	}
}


sub extract_region_annos {
        my ($anno, $regions, $dir_name, $format) = @_;
        my $suffix = $format eq 'gtf' ? 'gtf' : 'gff';

        my %chr_regions;   # chr => [ [start, end, key], ... ]
        my %key_to_coord;  # key => [chr, start, end]（
        for my $r (@$regions) {
                my ($key, $chr, $s, $e) = @$r;
                unless (exists $key_to_coord{$key}) {
                        $key_to_coord{$key} = [$chr, 0+$s, 0+$e];
                        push @{$chr_regions{$chr}}, [0+$s, 0+$e, $key];
                }
        }

        for my $chr (keys %chr_regions) {
                @{$chr_regions{$chr}} = sort { $a->[0] <=> $b->[0] } @{$chr_regions{$chr}};
        }
        my $total_regions = scalar(keys %key_to_coord);
        return unless $total_regions > 0;

        my %fh_cache;
        my @fh_queue;
        my $WINDOW = 100;
        my $get_fh = sub {
                my ($key) = @_;
                (my $key_filename = $key) =~ s/:/_/;
                return get_cached_fh($key_filename, $dir_name, $suffix, $WINDOW, \%fh_cache, \@fh_queue);
        };

	my $sort_cmd = sprintf("grep -v '^#' %s | sort -t \$'\\t' -k1,1 -k4,4n 2>/dev/null", quotemeta($anno));
	open my $in, "-|", $sort_cmd or die "Error: Cannot open sort pipe: $!";

        my $current_chr = '';
        my $region_list = [];
        my $r_idx = 0;
        my @active;

	while (<$in>) {
                next if /^#/ || /^$/;
                my @cols = split /\t/, $_, 9;
                next if @cols < 9;
                my ($chr, $start, $end) = ($cols[0], $cols[3] + 0, $cols[4] + 0);
                if ($chr ne $current_chr) {
                        @active = ();
                        $r_idx = 0;
                        $region_list = $chr_regions{$chr} || [];
                        $current_chr = $chr;
                }
                next unless @$region_list;

                while ($r_idx < @$region_list && $region_list->[$r_idx][0] <= $end) {
                        push @active, $region_list->[$r_idx];
                        $r_idx++;
                }

                my $j = 0;
                for my $i (0..$#active) {
			if ($active[$i][1] >= $start) {
                                $active[$j++] = $active[$i];
                        }
                }
                splice @active, $j;

                for my $reg (@active) {
			next unless $reg->[0] <= $end && $reg->[1] >= $start;
                        my $fh = $get_fh->($reg->[2]);
                        print $fh $_;
                }
        }

        close $in;
        close_all_cached_fhs(\%fh_cache, \@fh_queue);

        return;
}



sub extract_gene_annos {
        my ($anno, $genes, $dir_name, $format, $extend) = @_;

        $extend ||= 0;
        my %target_genes = map { $_ => 1 } @$genes;
        my $suffix = $format eq 'gtf' ? 'gtf' : 'gff';

        my %gene_coords;
        for my $g (@$genes) {
                $gene_coords{$g} = { chr => undef, min => undef, max => undef };
        }

        my $WINDOW = 100;
        my %fh_cache;
        my @fh_queue;
        my $get_fh = sub {
                my ($unit_id) = @_;
                return get_cached_fh($unit_id, $dir_name, $suffix, $WINDOW, \%fh_cache, \@fh_queue);
        };

        if ($format eq 'gtf') {
		#my $sort_cmd = sprintf("grep -v '^#' %s | sort -t \$'\\t' -k1,1 -k4,4n", quotemeta($anno));
		#open my $in, "-|", $sort_cmd or die "Error: Cannot open sort pipe: $!";
		open my $in, '<', $anno or die "Error: Cannot open $anno: $!";

                while (<$in>) {
                        next if /^#/ || /^$/;
                        next unless /gene_id\s+"([^"]+)"/;
                        my $gid = $1;
                        next unless $target_genes{$gid};

                        my $fh = $get_fh->($gid);
                        print $fh $_;

                        my @cols = split /\t/, $_, 9;
                        next if @cols < 8;
                        my ($chr, $s, $e) = @cols[0, 3, 4];
                        my $info = $gene_coords{$gid};
                        $info->{chr} = $chr;
                        $info->{min} = $s if !defined $info->{min} || $s < $info->{min};
                        $info->{max} = $e if !defined $info->{max} || $e > $info->{max};
                }
                close $in;
   	} else {
		#my $sorted_tmp = "$dir_name/.sorted.tmp.$$";
		#my $sort_cmd = sprintf("grep -v '^#' %s | sort -t \$'\\t' -k1,1 -k4,4n > %s", quotemeta($anno), quotemeta($sorted_tmp));
		#system($sort_cmd) == 0 or die "Sort failed: $?";

                my %id_to_parents;
		#open my $f1, '<', $sorted_tmp or die $!;
		open my $f1, '<', $anno or die "Error: Cannot open $anno: $!";
                while (<$f1>) {
                        next if /^#/ || /^$/;
                        my @cols = split /\t/, $_, 9;
                        next if @cols < 9;
                        my $attr = $cols[8];
                        my ($id) = $attr =~ /ID=([^;\s]+)/;
                        my ($parents_str) = $attr =~ /Parent=([^;]+)/;
                        if ($id && $parents_str) {
                                my @parents = split /,/, $parents_str;
                                s/^\s+|\s+$//g for @parents;
                                $id_to_parents{$id} = \@parents if @parents;
                        }
                }
                close $f1;

		#open my $f2, '<', $sorted_tmp or die $!;
		open my $f2, '<', $anno or die "Error: Cannot open $anno: $!";
                while (<$f2>) {
                        next if /^#/ || /^$/;
                        s/\r?\n$//;
                        my @cols = split /\t/, $_, 9;
                        next if @cols < 9;
                        my ($chr, $s, $e, $attr) = @cols[0, 3, 4, 8];
                        my ($id) = $attr =~ /ID=([^;\s]+)/;
                        my ($parents_str) = $attr =~ /Parent=([^;]+)/;
                        my @parents;
                        if ($parents_str) {
                                @parents = split /,/, $parents_str;
                                s/^\s+|\s+$//g for @parents;
                        }

                        my %parent_to_gene;  # key: parent_id, value: target_gene_id
                        if (@parents) {
                                for my $parent (@parents) {
                                        my $gene = trace_to_gene($parent, \%id_to_parents, \%target_genes, 20);
                                        $parent_to_gene{$parent} = $gene if $gene;
                                }
                        } elsif ($id && $target_genes{$id}) {
                                $parent_to_gene{''} = $id;
                        }
                        next unless %parent_to_gene;

                        my %gene_parents;
                        while (my ($par, $gene) = each %parent_to_gene) {
                                push @{$gene_parents{$gene}}, $par;
                        }

                        for my $gene (keys %gene_parents) {
                                my $parents_for_gene = $gene_parents{$gene};
                                my $info = $gene_coords{$gene};
                                $info->{chr} = $chr;
                                $info->{min} = $s if !defined $info->{min} || $s < $info->{min};
                                $info->{max} = $e if !defined $info->{max} || $e > $info->{max};

                                for my $parent (@$parents_for_gene) {
                                        my $new_attr;
                                        if ($parent eq '') {
                                                $new_attr = $attr;
                                        } else {
                                                $new_attr = $attr;
                                                $new_attr =~ s/Parent=[^;]+/Parent=$parent/;
                                        }
                                        my $fh = $get_fh->($gene);
                                        print $fh join("\t", @cols[0..7], $new_attr), "\n";
                                }
                        }
                }
                close $f2;
		#unlink $sorted_tmp;
        }

	close_all_cached_fhs(\%fh_cache, \@fh_queue);

        my @regions;
        for my $gene (@$genes) {
                my $info = $gene_coords{$gene};
                next unless defined $info->{chr};
                my $start = $info->{min} - $extend;
                $start = 1 if $start < 1;
                push @regions, [$gene, $info->{chr}, $start, $info->{max} + $extend];
        }

        return \@regions;
}


sub trace_to_gene {
        my ($start_id, $id_to_parents_ref, $target_genes_ref, $max_depth) = @_;
        my @currents = ($start_id);
        my %visited;
        my $depth = 0;

        while (@currents && $depth < $max_depth) {
                my @next_level;
                for my $current (@currents) {
                        next if $visited{$current}++;
                        return $current if $target_genes_ref->{$current};
                        if (exists $id_to_parents_ref->{$current}) {
                                push @next_level, @{$id_to_parents_ref->{$current}};
                        }
                }
                @currents = @next_level;
                $depth++;
        }
        return undef;
}


sub get_cached_fh {
        my ($unit_id, $dir, $suffix, $window, $cache_ref, $queue_ref) = @_;
        return $cache_ref->{$unit_id} if exists $cache_ref->{$unit_id};
        if (@$queue_ref >= $window) {
                my $old_unit = shift @$queue_ref;
                if (exists $cache_ref->{$old_unit}) {
                        close $cache_ref->{$old_unit};
                        delete $cache_ref->{$old_unit};
                }
        }
        my $file = "$dir/$unit_id.$suffix";
        open my $fh, '>>:utf8', $file or die "Error: Cannot append to $file: $!";
        $cache_ref->{$unit_id} = $fh;
        push @$queue_ref, $unit_id;
        return $fh;
}


sub close_all_cached_fhs {
        my ($cache_ref, $queue_ref) = @_;
        close $_ for values %$cache_ref;
        %$cache_ref = ();
        @$queue_ref = ();
}


sub detect_anno_format {
        my ($file) = @_;
        if ($file =~ /\.gff3?$/i) {
                return 'gff';
        } elsif ($file =~ /\.gtf$/i) {
                return 'gtf';
        } else {
                 die "Error: Please provide the GFF3/GTF file with gene annotation\n";
        }
}


sub parse_region_string {
        my ($str) = @_;
        if ($str =~ /^([^:]+):(\d+)-(\d+)$/) {
                my ($c, $s, $e) = ($1, int($2), int($3));
                die "Error: Invalid region: start > end ($str)\n" if $s > $e;
                return ($c, $s, $e);
        } else {
                die "Error: Invalid region format: $str (expected chr:start-end)\n";
        }
}
	


sub find_gene_anno {
	
	my ($gene_id, $anno) = @_;
	my $format = detect_anno_format($anno);
	
	open my $fh, '<', $anno or die "Error: Can't open file '$anno': $!\n";
	
	my @out;
	my %tx;
	my $gene_hit  = 0;
	my @gene_lines;

	my %is_child = map { $_ => 1 } qw(exon CDS five_prime_UTR three_prime_UTR UTR3 UTR5 UTR);

	if ($format eq 'gtf') {
		while (<$fh>) {
			next if substr($_, 0, 1) eq '#' || /^\s*$/;
			s/\r?\n$//;
			my @f = split /\t/, $_, 9;
			next unless @f >= 9;
			my ($type, $attr) = @f[2,8];

			my $g_id = _attr_gtf($attr, 'gene_id');
			next unless defined $g_id && $g_id eq $gene_id;

			if ($type eq 'gene') {
				$gene_hit++;
				push @gene_lines, $.;
				push @out, $_;
			}

			my $tx_id = _attr_gtf($attr, 'transcript_id');
			if (defined $tx_id) {
				$tx{$tx_id} = 1;
				if ($type eq 'transcript' || $type eq 'mRNA' || $is_child{$type}) {
					push @out, $_;
				}
			}
		}

		die "Error: Gene ID '$gene_id' found $gene_hit times (should be exactly 1) at line " . join(',', @gene_lines) . " in '$anno'\n" if $gene_hit > 1;
	} else {
		while (<$fh>) {
			next if substr($_, 0, 1) eq '#' || /^\s*$/;
			s/\r?\n$//;
			my @f = split /\t/, $_, 9;
			next unless @f >= 9;
			my ($type, $attr) = @f[2,8];

			if ($type eq 'gene') {
				my $id = _attr_gff($attr, 'ID');
				if (defined $id && $id eq $gene_id) {
					$gene_hit++;
					push @gene_lines, $.;
					push @out, $_;
				}
			} elsif ($type eq 'mRNA' || $type eq 'transcript') {
				my $parent = _attr_gff($attr, 'Parent');
				if (defined $parent && $parent eq $gene_id) {
					my $tx_id = _attr_gff($attr, 'ID');
					$tx{$tx_id} = 1 if defined $tx_id;
					push @out, $_;
				}
			} 
		}
		
		die "Error: Gene ID '$gene_id' found $gene_hit times (should be exactly 1) at line " . join(',', @gene_lines) . " in '$anno'\n" if $gene_hit > 1;

		seek $fh, 0, 0;
		while (<$fh>) {
			next if substr($_, 0, 1) eq '#' || /^\s*$/;
			s/\r?\n$//;
			my @f = split /\t/, $_, 9;
			next unless @f >= 9;
			my ($type, $attr) = @f[2,8];
			next unless $is_child{$type};
			my $parent = _attr_gff($attr, 'Parent');
			push @out, $_ if defined $parent && exists $tx{$parent};
		}
	}

	close $fh;
	return \@out;
}


sub _attr_gff {
    my ($attr, $key) = @_;
    return undef if index($attr, "$key=") < 0;
    return $1 if $attr =~ /(?:^|;\s*)$key=([^;]+)/;
    return undef;
}


sub _attr_gtf {
    my ($attr, $key) = @_;
    return undef if index($attr, "$key ") < 0;
    return $1 if $attr =~ /(?:^|;\s*)$key "([^"]+)"/;
    return undef;
}


sub get_gene_pos {

        my ($anno) = @_;
	my $format = detect_anno_format($anno);

	my %rank = (
        	gene            => 1,
        	mRNA            => 2,
        	transcript      => 2,
        	exon            => 3,
        	CDS             => 4,
        	three_prime_UTR => 4,
        	five_prime_UTR  => 4,
        	UTR3            => 4,
        	UTR5            => 4,
		UTR		=> 4
    	);
	
	open my $fh, '<', $anno or die "Error: Can't open file '$anno': $!\n";

	my @rec;
	if ($format eq 'gff') {
		while (<$fh>) {
			s/\r?\n$//;
			next if /^#/;
			next if /^\s*$/;
			my @f = split /\t/;
			next unless $f[2] eq 'gene';
			my $gid;
			for (split /;/, $f[8]) {
				my ($k, $v) = split /=/, $_, 2;
				if ($k eq 'ID') {
					$gid = $v;
					last;
				}
			}
			next unless defined $gid;
			push @rec, [ $f[0], $f[3], $gid ];
		}		
	} else {
		my %gene_info;
		while (<$fh>) {
			s/\r?\n$//;
			next if /^#/;
			next if /^\s*$/;
			my @f = split /\t/;
			next unless @f >= 9;
			my ($gid) = $f[8] =~ /gene_id "([^"]+)"/;
			next unless defined $gid;
			if (!exists $gene_info{$gid} || $f[3] < $gene_info{$gid}[1]) {
				$gene_info{$gid} = [ $f[0], $f[3] ];
			}
		}
		@rec = map { [ $gene_info{$_}[0], $gene_info{$_}[1], $_ ] } keys %gene_info;
	}
	close $fh;

	@rec = sort { $a->[0] cmp $b->[0] || $a->[1] <=> $b->[1] } @rec;

	my %seen;
	my @gene_sorted = grep { !$seen{$_}++ } map { $_->[2] } @rec;

	my @all_gene_pos;
	for my $gene (@gene_sorted) {
		my %gene_pos;
		my @trans_arr;
		my $gene_anno = find_gene_anno($gene, $anno);

		@$gene_anno = sort {
			my ($ta) = (split /\t/, $a)[2];
			my ($tb) = (split /\t/, $b)[2];
			my $ra = $rank{$ta} // 99;
			my $rb = $rank{$tb} // 99;
			$ra <=> $rb || $a cmp $b
		} @$gene_anno;

		my %exons;

		foreach my $line (@$gene_anno) {
			my @f = split /\t/, $line;
			my ($chr, $class, $start, $end, $strand, $attrs) = @f[0,2,3,4,6,8];

			if($format eq 'gff') {
				my %attr = $attrs =~ /([^;=]+)=([^;]*)/g;

				if($class eq "gene"){
					$gene_pos{'gene_id'} = $attr{"ID"};
					$gene_pos{'chr'} = $chr;
					$gene_pos{'start'} = $start;
					$gene_pos{'end'} = $end;
					$gene_pos{'strand'} = $strand;	
					$gene_pos{'gene_name'} = $attr{'gene_name'} // '';
				}elsif($class eq "mRNA" || $class eq "transcript"){
					push(@trans_arr, $attr{"ID"});
					push(@{$gene_pos{'ele'}}, ["transcript", $start, $end, $#trans_arr + 1, $attr{"ID"}]);
				}elsif ($class eq "exon"){
					my @idx = grep { $trans_arr[$_] eq $attr{"Parent"} } 0 .. $#trans_arr;
					next unless @idx;
					my $tx_idx = $idx[0] + 1;
					push(@{$gene_pos{'ele'}}, ["exon", $start, $end, $tx_idx, $attr{"ID"} // '']);
					push @{$exons{$tx_idx}}, [$start, $end, $attr{"ID"} // ''];
				}elsif ($class =~ /^(CDS|three_prime_UTR|five_prime_UTR|UTR3|UTR5|UTR)$/i) {
					my @idx = grep { $trans_arr[$_] eq $attr{"Parent"} } 0 .. $#trans_arr;
					next unless @idx;
					my $tx_idx = $idx[0] + 1;
					$class = 'UTR3' if $class eq 'three_prime_UTR';
					$class = 'UTR5' if $class eq 'five_prime_UTR';
					my $exon_info = _find_overlapping_exon($start, $end, $exons{$tx_idx});
					push(@{$gene_pos{'ele'}}, [$class, $start, $end, $tx_idx, $attr{"ID"} // '', $exon_info]);
				}
			} else {
				my %attr;
				while ($attrs =~ /(\S+)\s+"([^"]*)"/g) {
					$attr{$1} = $2;
				}
				if ($class eq "gene") {
					$gene_pos{'gene_id'} = $attr{"gene_id"} // $gene;
					$gene_pos{'chr'} = $chr;
					$gene_pos{'start'} = $start;
					$gene_pos{'end'} = $end;
					$gene_pos{'strand'}  = $strand;
					$gene_pos{'gene_name'} = $attr{'gene_name'} // '';
				} elsif ($class eq "mRNA" || $class eq "transcript") {
					my $tx_id = $attr{"transcript_id"};
					next unless defined $tx_id;
					push @trans_arr, $tx_id;
					push @{$gene_pos{'ele'}}, ["transcript", $start, $end, $#trans_arr + 1, $tx_id];
				} elsif ($class eq "exon") {
					my $tx_id = $attr{"transcript_id"};
					next unless defined $tx_id;
					my @idx = grep { $trans_arr[$_] eq $tx_id } 0 .. $#trans_arr;
					next unless @idx;
					my $tx_idx = $idx[0] + 1;
					my $exon_info = _find_overlapping_exon($start, $end, $exons{$tx_idx});
					push @{$gene_pos{'ele'}}, ["exon", $start, $end, $tx_idx, $attr{"exon_id"} // '', $exon_info];
					push @{$exons{$tx_idx}}, [$start, $end, $attr{"exon_id"} // ''];
				} elsif ($class =~ /^(CDS|three_prime_UTR|five_prime_UTR|UTR3|UTR5|UTR)$/i) {
					my $tx_id = $attr{"transcript_id"};
					next unless defined $tx_id;
					my @idx = grep { $trans_arr[$_] eq $tx_id } 0 .. $#trans_arr;
					next unless @idx;
					my $tx_idx = $idx[0] + 1;
					$class = 'UTR3' if $class eq 'three_prime_UTR';
					$class = 'UTR5' if $class eq 'five_prime_UTR';
					my $exon_info = _find_overlapping_exon($start, $end, $exons{$tx_idx});
					push @{$gene_pos{'ele'}}, [$class, $start, $end, $tx_idx, $attr{"exon_id"} // '', $exon_info];
				}
			}
		}
		push @all_gene_pos, \%gene_pos;
	}

	return \@all_gene_pos;

}


sub _find_overlapping_exon {
        my ($start, $end, $exon_list) = @_;
        return undef unless defined $exon_list && @$exon_list;
        
        for my $exon (@$exon_list) {
                my ($e_start, $e_end, $e_id) = @$exon;
                if ($start <= $e_end && $end >= $e_start) {
                        return $e_id;
                }
        }
        return undef;
}


1;





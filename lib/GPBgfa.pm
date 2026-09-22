#!/usr/bin/perl

package GPBgfa;

use strict;
use warnings;
use GPBpath;

sub parse_graph_gfa {
	
	my ($gfa_file, $refname) = @_;

	my @nodes;
        my %node_sample;
        my %edge_pair;
        my %ref_info;
	my @node_inv;
	my %node_seq;
	my %dup_map;
	my $max_id = 0;
	my @excluded_paths;
	my @p_names;
	my %warned_raw;

	push(@nodes, ["0+", ""]);
	push(@nodes, ["Inf+", ""]);

	open my $fh_gfa, '<', $gfa_file or die "Error: Can't open file '$gfa_file': $!\n";
        while(<$fh_gfa>){
                s/\r?\n$//;
                s/\r//g;
                my @arr = split(/\t/, $_);
                my $class = $arr[0];
                if($class eq "S"){
			push(@nodes, [$arr[1]."+", $arr[2]]);
			$node_seq{$arr[1]} = $arr[2];
			$max_id = $arr[1] if $arr[1] =~ /^\d+$/ && $arr[1] > $max_id;
                }elsif($class eq "P"){
                        if ($arr[1] =~ /_MINIGRAPH_|^Consensus_/) {
                                push @excluded_paths, $arr[1];
                                next;
                        }
                        push @p_names, $arr[1];
                        my @cur_nodes = split(/,/, $arr[2]);
                        my %seen_in_path;
                        for(my $i = 0; $i < @cur_nodes; $i++){
                                my $cur = $cur_nodes[$i];
                                next unless $seen_in_path{$cur}++;
                                my ($base, $orient) = $cur =~ /^(.+?)([+-])$/ ? ($1, $2) : ($cur, '+');
                                $max_id++;
                                $max_id++ while exists $node_seq{$max_id};
                                my $new_node = $max_id . $orient;
                                $node_seq{$max_id} = $node_seq{$base} // '';
                                push(@nodes, [$max_id."+", $node_seq{$max_id}]);
                                $dup_map{$new_node} = $cur;
                                $cur_nodes[$i] = $new_node;
                        }
                        $arr[2] = join(",", @cur_nodes);
                        unshift(@cur_nodes, "0+");
			push(@cur_nodes, "Inf+");
                        for(my $i = 0; $i < @cur_nodes - 1; $i++){
                                my $from = $cur_nodes[$i];
                                my $to = $cur_nodes[$i + 1];
				push @node_inv, $from =~ s/-$//r if $from =~ /-$/;
				push @node_inv, $to =~ s/-$//r if $to =~ /-$/;
                                if(exists $edge_pair{"${from}_${to}"}){
                                        $edge_pair{"${from}_${to}"} += 1;
                                }else{
                                        $edge_pair{"${from}_${to}"} = 1;
                                }

                        }
                        my $sample_name;
                        if ($arr[1] =~ /^([^#]+)#([^#]+)#/) {
                                $sample_name = "$1#$2";
                        }elsif ($arr[1] =~ /^([^#]+)#+/ || $arr[1] =~ /^(.*?)(?=\.(?i:chr))/) {
                                $sample_name = $1;
                        }else{
                                $sample_name = $arr[1];
                                if (!$warned_raw{$arr[1]}++) {
                                        warn "Warning: cannot extract sample name from path name '$arr[1]'; using the full path name as its sample.\n";
                                }
                        }
                        foreach my $node (@cur_nodes){
                                if(exists $node_sample{$node}){
                                        $node_sample{$node} .= ",$sample_name";
                                }else{
                                        $node_sample{$node} = $sample_name;
                                }
                        }
                        if($arr[1] =~ /^\Q$refname\E:(\d+)-(\d+)$/){
                                 die "Error: multiple paths match the resolved reference name '$refname' in '$gfa_file'\n" if $ref_info{nodes};
                                 $ref_info{chr} = $refname;
                                 $ref_info{start} = $1 + 1;
                                 $ref_info{end} = $2;
                                 $ref_info{nodes} = $arr[2];
                        } elsif ($arr[1] eq $refname) {
                                 die "Error: multiple paths match the resolved reference name '$refname' in '$gfa_file'\n" if $ref_info{nodes};
                                 $ref_info{chr} = $refname;
                                 $ref_info{nodes} = $arr[2];
                                 my $ref_rec = GPBpath::parse_path_name($arr[1]);
                                 if (defined $ref_rec->{start}) {
                                         ($ref_info{start}, $ref_info{end}) = GPBpath::coord_1based($ref_rec);
                                 }
                        }
                }
        }

        close $fh_gfa;

	if (!$ref_info{nodes}) {
		my $msg = "Error: no path matches the resolved reference name '$refname' in '$gfa_file'.\n  Available paths: " . join(", ", @p_names) . "\n";
		$msg .= "  Note: these path(s) were excluded as auxiliary: " . join(", ", @excluded_paths) . "\n" if @excluded_paths;
		die $msg;
	}
	if (!defined $ref_info{start}) {
		my $ref_len = 0;
		foreach my $n (split(/,/, $ref_info{nodes})) {
			(my $b = $n) =~ s/[+-]$//;
			$ref_len += length($node_seq{$b} // '');
		}
		($ref_info{start}, $ref_info{end}) = (1, $ref_len);
		warn "Warning: reference path '$refname' carries no coordinates; displaying 1-$ref_len as local coordinates.\n";
	}

	if (@excluded_paths) {
		warn "Info: excluded " . scalar(@excluded_paths) . " auxiliary path(s) matching '_MINIGRAPH_' or '^Consensus_': " . join(", ", @excluded_paths) . "\n";
	}

	$node_sample{"0+"} = "";
	$node_sample{"Inf+"} = "";

	my %seen;
	@node_inv = grep !$seen{$_}++, @node_inv;
	foreach my $nodei (@node_inv){
		push(@nodes, [$nodei."-", "The reverse complement of ".$nodei."+"]);
	}

	my @edge_out;
        foreach my $key (keys %edge_pair) {
                my ($start, $end) = split /_/, $key;
                push @edge_out, [$start, $end, $edge_pair{$key}];
        }

	my @cycle_nodes = detect_cycle_nodes(\@edge_out);
	if(@cycle_nodes){
		warn "Warning: nodes [".join(", ", @cycle_nodes)."] remain in a cycle that cannot be resolved by path-based duplication. The layout of this region may fail; please check the input graph.\n";
	}

	my @node_add_sample;
        foreach my $node (@nodes){
		push(@node_add_sample, [@$node[0], @$node[1], $node_sample{@$node[0]} // '']);
        }

	return (\@node_add_sample, \@edge_out, \%ref_info, \%dup_map);

}


sub detect_cycle_nodes {

	my ($edge_out) = @_;

	my %in_degree;
	my %graph;
	foreach my $e (@$edge_out){
		my ($from, $to) = @$e[0, 1];
		$in_degree{$from} += 0;
		$in_degree{$to} += 0;
		$graph{$from} ||= [];
		push @{$graph{$from}}, $to;
		$in_degree{$to}++;
	}

	my @queue = grep { $in_degree{$_} == 0 } keys %in_degree;
	my %reachable;
	while(@queue){
		my $cur = shift @queue;
		$reachable{$cur} = 1;
		foreach my $nb (@{$graph{$cur} || []}){
			if(--$in_degree{$nb} == 0){
				push @queue, $nb;
			}
		}
	}

	return grep { !exists $reachable{$_} } sort keys %in_degree;

}


sub parse_variant_gfa {

	my ($gfa_file,  $sample_table_file, $vcf_file, $chr) = @_;

	chomp(my $line = `awk '/^#CHROM/{print; exit}' $vcf_file`);
	 my @h = split /\t/, $line;
	 my $sample_n = @h - 9;
	 my @all_samples = @h[9..$#h];

	my %id2samples;
	if(defined $sample_table_file && -s $sample_table_file){
		open my $st, '<', $sample_table_file or die "Error: Can't open file '$sample_table_file': $!\n";
		while(<$st>){
			s/\r?\n$//;
			my ($id, $samples) = split /\t/, $_, 2;
			$id2samples{$id} = [split /,/, $samples];
		}
		close $st;
	}

	open my $gfa, '<', $gfa_file or die "Error: Can't open file '$gfa_file': $!\n";
        my @nodes;
	my %node_sample;
        my @edges;
        my %edge_pair;
        my %ref_info;
        my %alt_paths;
        my @ref_nodes;
        my @ref_edge_pair;

	push(@nodes, ["0+", ""]);
        push(@nodes, ["Inf+", ""]);

        while(<$gfa>){
                s/\r?\n$//;
                s/\r//g;
                my @arr = split(/\t/, $_);
                my $class = $arr[0];
                if($class eq "S"){
			push(@nodes, [$arr[1]."+", $arr[2]]);
                }elsif($class eq "L"){
			my ($s1, $o1, $s2, $o2) = @arr[1..4];
			($s1,$o1,$s2,$o2) = ($s2,'+',$s1,'+') if $o1 eq '-' && $o2 eq '-';
                        push @edges, ["$s1$o1", "$s2$o2"];
                }elsif($class eq "P"){
                        my $cur_path = $arr[1];
                        my $cur_nodes = $arr[2];
                        my @cur_node_arr = split(/,/, $cur_nodes);
                        if($cur_path =~ /^\Q$chr\E:(\d+)-(\d+)$/){
                                die "Error: multiple paths match the reference name '$chr' in '$gfa_file'\n" if $ref_info{nodes};
                                $ref_info{chr} = $chr;
                                $ref_info{start} = $1 + 1;
                                $ref_info{end} = $2;
                                $ref_info{nodes} = $cur_nodes;
                                @ref_nodes = @cur_node_arr;
                        }elsif($cur_path =~ /_alt_(ID\d+-\d+)_(\d+):(\d+)-(\d+)/){
                                $alt_paths{$1}->{$2} = [@cur_node_arr];
                        }
                }

        }

	if (!$ref_info{nodes}) {
		die "Error: no path matches the reference name '$chr' in '$gfa_file'; the extracted subgraph is expected to contain a '$chr:start-end' reference path\n";
	}

	foreach my $e (@edges){
                my @e = @$e;
                $edge_pair{"$e[0]_$e[1]"} = [];
        }
	unshift @ref_nodes, "0+";
	push @ref_nodes, "Inf+";

        foreach my $x (0..($#ref_nodes-1)){
                $edge_pair{"$ref_nodes[$x]_$ref_nodes[$x+1]"} = [@all_samples];
                push(@ref_edge_pair, "$ref_nodes[$x]_$ref_nodes[$x+1]");
        }

	@node_sample{ @ref_nodes } = ([@all_samples]) x @ref_nodes;

	foreach my $id (sort keys %alt_paths){
        if ($id =~ /(ID\d+-(\d+))$/) {
		my $variant_id = $1;
		my @variant_sample = @{$id2samples{$variant_id}};
                my $sample_count = $2;
                my @ref_node_idx;
                if($alt_paths{$id}{0} && $alt_paths{$id}{1}){
                        my @alt_node_0 = @{$alt_paths{$id}{0}};
                        my @alt_node_1 = @{$alt_paths{$id}{1}};
                        foreach my $node (@alt_node_0){
                                my ($idx) = grep{$ref_nodes[$_] eq $node} 0..$#ref_nodes;
                                if(defined $idx){
                                        push(@ref_node_idx, $idx);
                                }
                        }
                        if(scalar(@ref_node_idx) == scalar(@alt_node_0)){
                                if($ref_node_idx[0] > 0){
                                        unshift(@ref_node_idx, $ref_node_idx[0] - 1);
                                }
                                foreach my $i (@ref_node_idx){
					my $diff_arr = array_setdiff($edge_pair{"$ref_nodes[$i]_$ref_nodes[$i+1]"}, \@variant_sample);
					$edge_pair{"$ref_nodes[$i]_$ref_nodes[$i+1]"} = $diff_arr;
                                }
				foreach my $node (@{$alt_paths{$id}{0}}){
					$node_sample{$node} = array_setdiff($node_sample{$node}, \@variant_sample);
				}
                        }
                        
                        unshift(@alt_node_1, $ref_nodes[$ref_node_idx[0]]);
                        push(@alt_node_1, $ref_nodes[$ref_node_idx[$#ref_node_idx] + 1]);
                        foreach my $j (0..($#alt_node_1-1)){
				my $old_arr = $edge_pair{"$alt_node_1[$j]_$alt_node_1[$j+1]"} || [];
				my $union_arr = array_union($old_arr, \@variant_sample);
				$edge_pair{"$alt_node_1[$j]_$alt_node_1[$j+1]"} = $union_arr;
                        }
			foreach my $node (@{$alt_paths{$id}{1}}){
				next if $node eq '0+' || $node eq 'Inf+';
				$node_sample{$node} = array_union($node_sample{$node}, \@variant_sample);
			}
               }elsif($alt_paths{$id}{0}){
                        my @alt_node_0 = @{$alt_paths{$id}{0}};
                        foreach my $node (@alt_node_0){
                                my ($idx) = grep{$ref_nodes[$_] eq $node} 0..$#ref_nodes;
                                if(defined $idx){
                                        push(@ref_node_idx, $idx);
                                }
                        }
                        if(scalar(@ref_node_idx) == scalar(@alt_node_0)){
                                if($ref_node_idx[0] > 0){
                                        unshift(@ref_node_idx, $ref_node_idx[0] - 1);
                                }
                                foreach my $i (@ref_node_idx){
					my $diff_arr = array_setdiff($edge_pair{"$ref_nodes[$i]_$ref_nodes[$i+1]"}, \@variant_sample);
					$edge_pair{"$ref_nodes[$i]_$ref_nodes[$i+1]"} = $diff_arr;
                                }
				foreach my $node (@{$alt_paths{$id}{0}}){
					$node_sample{$node} = array_setdiff($node_sample{$node}, \@variant_sample);
				}
                        
				my $old_arr = $edge_pair{"$ref_nodes[$ref_node_idx[0]]_$ref_nodes[$ref_node_idx[$#ref_node_idx]+1]"} || [];
				my $new_arr;
				if ($old_arr && @$old_arr){
					$new_arr = array_union($old_arr, array_intersect($node_sample{$ref_nodes[$ref_node_idx[0]]}, \@variant_sample));
				}else{
					$new_arr = array_intersect($node_sample{$ref_nodes[$ref_node_idx[0]]}, \@variant_sample);
				}

				my %in_ref = map { $ref_nodes[$_] => 1 } @ref_node_idx;
				foreach my $k (keys %edge_pair) {
					next unless index($k, "$ref_nodes[$ref_node_idx[0]]_") == 0;
					next if ($k eq "$ref_nodes[$ref_node_idx[0]]_$ref_nodes[$ref_node_idx[$#ref_node_idx]+1]");
					next if (@{$edge_pair{$k}} == 0);

    					my ($s, $e) = split /_/, $k;  
    					next unless exists $in_ref{$s};

    					my ($s_idx) = grep { $ref_nodes[$_] eq $s } 0 .. $#ref_nodes;
    					my ($e_idx) = grep { $ref_nodes[$_] eq $e } 0 .. $#ref_nodes;
    					next if (defined $e_idx && abs($e_idx - $s_idx) == 1);

    					$edge_pair{$k} = array_setdiff($edge_pair{$k}, $new_arr);
				}

				$edge_pair{"$ref_nodes[$ref_node_idx[0]]_$ref_nodes[$ref_node_idx[$#ref_node_idx]+1]"} = $new_arr;

			}
		}elsif($alt_paths{$id}{1}){
                        my @alt_node_1 = @{$alt_paths{$id}{1}};
                        foreach my $e (@ref_edge_pair){
                                my @c = split(/_/, $e);
                                if(exists $edge_pair{"$c[0]_$alt_node_1[0]"} && exists $edge_pair{"$alt_node_1[$#alt_node_1]_$c[1]"}){
                                        unshift(@alt_node_1, $c[0]);
                                        push(@alt_node_1, $c[1]);
					$edge_pair{"$c[0]_$c[1]"} = array_setdiff($edge_pair{"$c[0]_$c[1]"}, \@variant_sample);
					foreach my $j (0..($#alt_node_1-1)){
						 my $old_arr = $edge_pair{"$alt_node_1[$j]_$alt_node_1[$j+1]"} || [];
						 my $union_arr = array_union($old_arr, \@variant_sample);
						 $edge_pair{"$alt_node_1[$j]_$alt_node_1[$j+1]"} = $union_arr;
					}
					foreach my $node (@{$alt_paths{$id}{1}}){
						next if $node eq '0+' || $node eq 'Inf+';
						$node_sample{$node} = array_union($node_sample{$node}, \@variant_sample);
					}
					last;
                                }
                        }
                }
        }

	}

	close $gfa;

	my @edge_out;
	foreach my $key (keys %edge_pair) {
    		my ($start, $end) = split /_/, $key;
		ref($edge_pair{$key}) eq 'ARRAY' or warn "bad key=$key  value=$edge_pair{$key}  ref=" . ref($edge_pair{$key}) . "\n";
    		push @edge_out, [$start, $end, scalar @{$edge_pair{$key}}];
	}

	$node_sample{"0+"} = [];
        $node_sample{"Inf+"} = [];

	my @node_add_sample;
        foreach my $node (@nodes){
		my $node_id = $node->[0];
		my $sample_list = $node_sample{$node_id} || []; 
		my $cur_samples = join ',', @$sample_list; 
                push(@node_add_sample, [$node->[0], $node->[1], $cur_samples]);
        }
	
	return (\@node_add_sample, \@edge_out, \%ref_info);

}


sub array_setdiff {
    my ($arr1_ref, $arr2_ref) = @_;
    my %seen;
    my @difference;

    @seen{@$arr2_ref} = ();

    foreach my $item (@$arr1_ref) {
        push @difference, $item unless exists $seen{$item};
    }

    return \@difference;
}


sub array_union {
    my ($arr1_ref, $arr2_ref) = @_;
    my %seen;
    @seen{ @{ $arr1_ref || [] }, @{ $arr2_ref || [] } } = ();
    return [ keys %seen ];
}


sub array_intersect {
    my ($a_ref, $b_ref) = @_;
    return [] unless $a_ref && $b_ref; 

    my %seen;
    @seen{@$b_ref} = (); 

    my %common;  
    exists $seen{$_} and $common{$_} = 1 for @$a_ref;

    return [ keys %common ];
}

1;


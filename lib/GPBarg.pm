#!/usr/bin/perl

package GPBarg;

use strict;
use warnings;
use GPBpath;


sub check_input {
	my ($vcf, $reffa, $graph, $refname, $anno, $geneid, $genelist, $region, $regionlist, $noextract, $extend, $bed, $pheno, $maxe, $maxd, $threads) = @_;

	my $has_vcf_mode = ($vcf || $reffa);
	my $has_graph_mode = ($graph || $refname);

	if($has_vcf_mode && $has_graph_mode){
		die "Error: Cannot use both Variant mode and Graph mode parameters simultaneously\n";
	}elsif(!$has_vcf_mode && !$has_graph_mode){
		die "Error: Must select one input mode: Variant mode (--variant/--ref-fa) or Graph mode (--graph/--ref-name)\n";
	}else{
		die "Error: --ref-fa is required when using --variant\n" if $vcf && !$reffa;
        	die "Error: --variant is required when using --ref-fa\n" if $reffa && !$vcf;
        	die "Error: --ref-name is required when using --graph\n" if $graph && !$refname;
        	die "Error: --graph is required when using --ref-name\n" if $refname && !$graph;
	}

	die "Error: --no-extract is only available in Graph mode (--graph/--ref-name)\n" if $noextract && $has_vcf_mode;

	if((defined $geneid || defined $genelist) && !defined $anno){
		die "Error: Please provide the GFF/GTF file with gene annotation\n";
	}

	my @range_options = ($geneid, $genelist, $region, $regionlist, $noextract);
	my $range_count = scalar grep { defined $_ } @range_options;
	if ($range_count == 0){
		die "Error: Please specify one analysis target parameter: --geneid, --geneid-list, --region, --region-list, or --no-extract\n";
	}elsif ($range_count > 1) {
		die "Error: Can only specify one analysis target parameter: --geneid, --geneid-list, --region, --region-list, or --no-extract\n";
	}

	if (defined $extend && !$geneid && !$genelist) {
		die "Error: --extend parameter can only be used with --geneid or --geneid-list\n" if $noextract;
    		warn "Warning: --extend parameter can only be used with --geneid or --geneid-list\n";
	}

	die "Error: -e/--extend value must be greater than or equal to 0\n" if defined $extend && $extend < 0;
        die "Error: -d/--max-distance-subpaths value must be greater than or equal to 0\n" if defined $maxd && $maxd < 0;
        die "Error: -m/--max-merging-iterations value must be greater than or equal to 0\n" if defined $maxe && $maxe < 0;
        die "Error: -t/--threads value must be greater than 0\n" if defined $threads && $threads < 1;
        die "Error: Invalid region format: $region. Expected format: 'chr:start-end' (e.g., chr1:1000-2000), and end must be greater than start.\n" if $region && !validate_region_format($region);

	check_file($vcf, "Variants") if defined $vcf;
	check_vcf_index($vcf) if defined $vcf;
	check_file($reffa, "Reference genome") if defined $reffa;
	check_file($graph, "Graph pangenome") if defined $graph;

	check_file($anno, "Gene annotation") if defined $anno;
	check_zip($anno) if defined $anno;
	validate_anno_file($anno) if defined $anno;
	check_file($bed, "Additional annotation") if defined $bed;
	validate_bed_file($bed, 1) if defined $bed;
	check_file($pheno, "Phenotype") if defined $pheno;
	validate_phenotype_file($pheno) if defined $pheno;

	check_file($genelist, "Gene list") if defined $genelist;
	check_file($regionlist, "Region list") if defined $regionlist;
	validate_bed_file($regionlist, 0) if defined $regionlist;

}

sub check_zip {
	my ($file) = @_;

	$file =~ s/'/'\\''/g;	
	my $type = `file -b --mime-type '$file' 2>/dev/null`;
	$? && die "Error: 'file' command failed or not found.\n";
	
	chomp($type);
	if ($type =~ /^application\/(gzip|x-gzip|zip|x-bzip2|xz)$/){
		die "Error: '$file' is compressed ($type). Please provide a plain-text GFF file.\n"
	}

	return 1;
}

sub mk_outdir {
	my ($out) = @_;

	my $dir_name;
	if (defined $out){
        	$dir_name = $out;
	} else {
        	my ($sec, $min, $hour, $mday, $mon, $year) = localtime;
        	$dir_name = sprintf("fgpb_res_%02d-%02d_%02d%02d%02d", $mon+1, $mday, $hour, $min, $sec);
	}
	die "Error: Output directory \"$dir_name\" already exists. To avoid overwriting of existing files, we kindly request that the output directory should not exist.\n" if -e $dir_name;

	mkdir $dir_name or die "Error: Failed to create directory \"$dir_name\": $!\n";
	
	return $dir_name;
}


sub check_file {
    my ($file, $type) = @_;
    return 0 unless defined $file;

    if (!-e $file) {
        die "Error: $type file does not exist: $file\n";
        return 0;
    }
    if (!-r $file) {
        die "Error: $type file is not readable: $file\n";
        return 0;
    }
    if (-z $file) {
        die "Error: $type file is empty: $file\n";
        return 0;
    }

    return 1;
}


sub check_vcf_index {
    my ($vcf) = @_;

    -e $vcf or die "Error: File'$vcf' does not exist\n";

    open my $fh, '<:raw', $vcf or die "Error: Can't open file '$vcf': $!\n";
    my $head;
    read $fh, $head, 16;
    close $fh;
    substr($head, 0, 3) eq "\x1f\x8b\x08" or die "Error: '$vcf' is not bgzip (block-gzipped) format\n";
	
    (-f "$vcf.tbi" || -f "$vcf.csi") or die "Error: file '$vcf' missing index file (.tbi or .csi)\n";

    return 1;
}


sub validate_region_format {
    my ($region) = @_;
    return 0 unless defined $region;

    if ($region =~ /^([^:]+):(\d+)-(\d+)$/) {
        my ($chr, $start, $end) = ($1, $2, $3);
        return 0 if $start >= $end;
        return 0 if $start < 1;
        return 1;
    }

    return 0;
}


sub validate_anno_file {
    my ($anno) = @_;
    open my $fh, '<', $anno or die "Cannot open $anno: $!";

    my ($line_num, $fmt, $has_gene, $has_transcript);

    while (<$fh>) {
    	next if /^#/ || /^\s*$/;
	$line_num++;
	s/\r?\n$//;

	my @cols = split /\t/, $_, -1;
	die "Line $line_num: expected 9 columns, got " . scalar(@cols) . "\n" if @cols != 9;

	unless ($fmt) {
		my $a = $cols[8];
		if ($a =~ /gene_id\s/) { $fmt = 'gtf'; }
		elsif ($a =~ /=/) { $fmt = 'gff'; }
		else {die "Error: Please check the format of gene annotation file '$anno', only GFF3 or GTF format is accepted.\n"; }
	}

	if ($fmt eq 'gtf' && (!$has_gene || !$has_transcript)) {
		$has_gene = 1 if $cols[2] eq 'gene';
		$has_transcript = 1 if $cols[2] eq 'transcript' || $cols[2] eq 'mRNA';
	}

	my ($s, $e) = @cols[3,4];
	die "Error: Line $line_num in file '$anno': start/end must be positive integers\n" unless $s =~ /^\d+$/ && $e =~ /^\d+$/;
	die "Error: Line $line_num in file '$anno': start > end ($s > $e)\n" if $s > $e;
    }

    die "Error: No valid lines found in file '$anno'. Please check the gene annotation file.\n" unless $line_num && $fmt;
    close $fh;

    if ($fmt eq 'gtf') {
    	die "Error: Gene annotation file '$anno' appears to lack gene features, please check the file\n" unless $has_gene;
	die "Error: Gene annotation file '$anno' appears to lack transcript features, please check the file\n" unless $has_transcript;
    }

    return 1;
}


sub validate_bed_file {
    my ($file, $isanno) = @_;
    open my $fh, '<', $file or die "Error: Can't open file '$file': $!\n";
    my $line_no = 0;
    while (my $line = <$fh>) {
        $line_no++;
	next if $line =~ /^#/;
        next if $line =~ /^\s*$/;
        $line =~ s/\s+$//;
        my @f = split /\t/, $line;

	if($isanno){
        	@f == 4 or die "Error: Line $line_no in file '$file': expected 4 columns (sequence, start, end, annotation); got " . scalar(@f) . "\n";
	}else{
		@f == 3 or die "Error: Line $line_no in file '$file': expected 3 columns (sequence, start, end); got " . scalar(@f) . "\n";
	}

        $f[1] =~ /^[0-9]+$/ or die "Error: Line $line_no in file '$file': start '$f[1]' not integer \n";
        $f[2] =~ /^[0-9]+$/ or die "Error: Line $line_no in file '$file': end   '$f[2]' not integer \n";
        $f[1] < $f[2]       or die "Error: Line $line_no in file '$file': start >= end '\n";
    }
    close $fh;
    return 1;
}


sub validate_phenotype_file {
        my ($file) = @_;

        open my $fh, '<', $file or die "Error: Can't open file '$file': $!\n";

        my $header = <$fh>;
        $header =~ s/\r?\n$//;
        my @headers = split(/\t/, $header);

        if (scalar @headers < 2) {
                warn "Warning: Phenotype file should contain at least one phenotype column besides 'Path'\n";
        }

        my $line_num = 1;
        my %path_names;
        while (my $line = <$fh>) {
                $line_num++;
                $line =~ s/\r?\n$//;
		next if $line =~ /^#/;
                next if $line =~ /^\s*$/;
                my @fields = split(/\t/, $line);
                if (scalar @fields < 1) {
                        warn "Warning: Line $line_num in file '$file': Empty line or missing data\n";
                        next;
                }
                my $path_name = $fields[0];
                if (exists $path_names{$path_name}) {
                        warn "Warning: Line $line_num in file '$file': Duplicate path name '$path_name'\n";
                }
                $path_names{$path_name}++;

        }

        close $fh;

	return 1;
}


sub get_ref_path_name {
	my ($og_file, $refname, $thread, $noextract) = @_;
	
	my $cmd = "odgi paths -i '$og_file' -L --threads $thread";
	my $output = `$cmd`;

	if ($? != 0) {
		die "Error: Failed to execute 'odgi paths -i $og_file -L'\n";
	}

	my $q = GPBpath::parse_path_name($refname);
	if (!$noextract && defined $q->{start}) {
		die "Error: --ref-name accepts a path name only, not a coordinate suffix ('$refname'). Specify the analysis range with --geneid/--geneid-list/--region/--region-list instead.\n";
	}

	my @all_paths = split /\n/, $output;
	my %opts = $noextract ? () : (forbid_coords => 1, require_locus => 1);
	my $ref = GPBpath::resolve_ref_path($refname, \@all_paths, \%opts);

	if (@all_paths <= 1) {
		die "Error: Only the reference genome path found. Please provide a pangenome graph containing multiple genome paths.\n"; 
	}

	return $ref->{raw};
}


1;


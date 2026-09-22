# fGPB: Flow-based Graph Pangenome Browser

*A novel flow-based visualization tool for pangenome graphs*

fGPB is a visualization tool that employs a novel flow-based layout to clearly reveal graph topology and population frequency. It integrates the sequence graph with reference genome annotations for direct evaluation of variant impact, and can incorporate phenotypic data to explore genotype-to-phenotype associations. All features are accessible through an interactive web interface.

**Key Features:**

* **1. Flow-Based Layout Algorithm**  
  Visualizes graph structure and population frequency simultaneously.

* **2. Integrated Genome Annotation**  
  Integrates the sequence graph with reference annotations to assess variant effects on genomic features.

* **3. Phenotypic Data Visualization**  
  Enables overlay of phenotypic information for validating variant-to-phenotype links.

* **4. Interactive Web Interface**  
  Provides zooming, panning, and data inspection capabilities for interactive exploration.

## Installation

### Direct Installation

#### Requirements
Basic (Both graph mode and Variant mode):
- Perl
- ODGI (v0.9 or later, https://github.com/pangenome/odgi)

Variant mode:
- VG (bundled in `ext/bin/`)
- Samtools (v1.16 or later, https://github.com/samtools/samtools)
- Bcftools (v1.16 or later, https://github.com/samtools/bcftools)

```
# Download the fGPB from Github
git clone https://github.com/SJTU-CGM/fGPB.git

# Add `fgpb` to `PATH` and add `lib/` to `PERL5LIB`
export PATH="/path/to/fGPB:$PATH"
export PERL5LIB="/path/to/fGPB/lib${PERL5LIB:+:$PERL5LIB}"

fgpb --help
```

### Conda
```
# Download the fGPB from Github
git clone https://github.com/SJTU-CGM/fGPB.git
cd fGPB

# Create conda
conda env create -f fgpb.yml

# Configure PATH and PERL5LIB 
bash scripts/setup_env.sh

# Activate the fGPB
conda activate fGPB

fgpb --help
```

### Docker
```
# Download the fGPB from Github
git clone https://github.com/SJTU-CGM/fGPB.git
cd fGPB

# Build the Docker image from the Dockerfile in the current directory
docker build -f Dockerfile -t fgpb:latest .

# Run fGPB with your data mounted
docker run --rm -v $(pwd):/data fgpb:latest --help

# Or enter the container interactively (run 'fgpb --help' inside)
docker run --rm -it -v $(pwd):/data --entrypoint /bin/bash fgpb:latest
```

### Singularity
```
# Download the pre-built Singularity image
wget https://cgm.sjtu.edu.cn/fGPB/src/fgpb_v1.0.sif

# Run fGPB directly from the Singularity image
singularity exec fgpb_v1.0.sif fgpb --help
```

## Usage 
A listing of all parameters can be obtained with fgpb --help or fgpb -h.
```
Usage:
Two running modes are supported:

Mode 1: Extract gene(s)/region(s) from a pangenome graph, and visualize
    fgpb --graph <graph.og> --ref-name <ref_path> (--geneid <ID> | --geneid-list <list.file> | --region <chr:start-end> | --region-list <list.bed>) [OPTIONS]
    Special case: with '--no-extract', the input graph is visualized as-is (e.g. an already
    extracted subgraph), without subgraph extraction:
    fgpb --graph <graph.og> --ref-name <ref_path> --no-extract [OPTIONS]

Mode 2: Build graph from variant data, extract gene(s)/region(s), and visualize
    fgpb --variant <in.vcf> --ref-fa <ref.fa> (--geneid <ID> | --geneid-list <list.file> | --region <chr:start-end> | --region-list <list.bed>) [OPTIONS]

REQUIRED ARGUMENTS
  Input mode (choose one group):
    -v, --variant       <file>          Variants in VCF format (.vcf). (requires '-r/--ref-fa')
    -r, --ref-fa        <file>          Reference genome in FASTA format (.fa). (requires '-v/--variant')

    -g, --graph         <file>          Variation graph in ODGI format (.og). (requires '-R/--ref-name')
    -R, --ref-name      <string>        Name of the reference path in the graph, e.g. 'P1#0#chr1' or 'P1.chr1'.
                                        Coordinate suffixes are not accepted in extract modes; the analysis range
                                        is given by the target parameters below. (requires '-g/--graph')

  Analysis target (choose one):
    --geneid            <string>        Single gene ID to analyze. (requires '-a/--gene-anno')
    --geneid-list       <file>          File containing list of gene IDs (one per line). (requires '-a/--gene-anno')
    --region            <string>        Single genomic region in 'chr:start-end' format (e.g. chr1:1000-2000).
    --region-list       <file>          BED file with genomic regions (chrom<tab>start<tab>end, 0-based).
    --no-extract                        Visualize the input graph as-is (e.g. an already extracted subgraph),
                                        without subgraph extraction.
                                        (Graph mode only; '-d/-m/-e' do not apply. Gene/extra annotation tracks
                                        require the reference path name to carry locus/coordinate information.)

RECOMMENDED ARGUMENTS (DATA TRACKS)
    -a, --gene-anno     <file>          Gene annotation file in GFF3/GTF format for the reference genome.
    -x, --extra-anno    <file>          Additional reference-genome annotations besides genes in BED format
                                        (e.g. repeats, domains). 
    -p, --pheno         <file>          Tab-delimited phenotype matrix: first column contains path identifiers 
                                        (header: 'Path'), remaining columns represent different phenotype with 
                                        corresponding names as headers.    
                                        Missing values in the file are indicated by '-', 'NA', 'NaN', 'null', 
                                        'None', 'NULL', 'undefined', 'unknown' or 'Unknown'.

OPTIONAL ARGUMENTS
    -o, --out           <string>        Output directory name.

    -e, --extend        <n>             Extend analysis region by N bp upstream and downstream of genes.
                                        Only work with --geneid or --geneid-list.
                                        (Default:10)

    -d,                 <n>             Parameter for 'odgi extract'.
    --max-distance-subpaths             Maximum distance between subpaths allowed for merging them. It reduces
                                        the fragmentation of unspecified paths in the input path ranges.
                                        (Default:10000)

    -m,                 <n>             Parameter for 'odgi extract'.
    --max-merging-iterations            Maximum number of iterations in attempting to merge close subpaths. It
                                        stops early if during an iteration no subpaths were merged.
                                        (Default:3)

    -t, --threads       <n>             Thread number.

    -h, --help                          Print usage page. 

```

## Quick start
### Graph-based pangenome mode
For an already constructed graph pangenome, the following parameters are required:
* Graph pangenome file (`--graph`)
* Reference genome path name (`--ref-name`)
* Analysis target (one of `--geneid`, `--geneid-list`, `--region`, `--region-list`, or `--no-extract` for direct visualization of an already extracted subgraph)
* Reference genome gene annotation (`--gene-anno`, optional)
* Additional reference genome annotations (`--extra-anno`, optional)

**Note on reference path naming**: `--ref-name` must resolve to exactly one path in the graph. PanSN names (`sample#hap#locus`, e.g. `P1#0#chr1`), delimiter-style names (e.g. `P1.chr1`, `chr` case-insensitive) and bare names are all accepted. In extract modes the reference path must be a full-length path carrying locus information (no `:start-end` coordinate suffix); to visualize graphs whose paths carry coordinates or lack locus information, use `--no-extract`. Auxiliary paths such as `_MINIGRAPH_...` and `Consensus_...` are ignored automatically. Sample information is read from P lines only (W lines of GFA 1.1 are not supported yet). Phenotype files are matched at the genome level, so haplotype-specific path names like `P1#0#chr1` and `P1#1#chr1` share the phenotype row of `P1`.

Download demo data:
```bash
# Download chr09_mc.og (9.0G)
wget https://cgm.sjtu.edu.cn/fGPB/demo/rice/chr09_mc.og
# Download rice_gene_anno.gff3 (79M)
wget https://cgm.sjtu.edu.cn/fGPB/demo/rice/rice_gene_anno.gff3
# Download rice_repeat_anno.bed (41M)
wget https://cgm.sjtu.edu.cn/fGPB/demo/rice/rice_repeat_anno.bed
# Download rice_demo_pheno.txt (1.5K)
wget https://cgm.sjtu.edu.cn/fGPB/demo/rice/rice_demo_pheno.txt
```

Visualizing a single gene (Example: *LOC_Os09g28300*):
```
fgpb --graph chr09_mc.og --ref-name IRGSP-1.0 --gene-anno rice_gene_anno.gff3  --extra-anno rice_repeat_anno.bed --geneid LOC_Os09g28300 --out demores_rice_gene --pheno rice_demo_pheno.txt --extend 100
```
Visualizing multiple genes (Example: *LOC_Os09g29820*, *LOC_Os09g26999*, *LOC_Os09g15840*):
```
printf '%s\n' LOC_Os09g29820 LOC_Os09g26999 LOC_Os09g15840 > rice_demo_genelist.txt
fgpb --graph chr09_mc.og --ref-name IRGSP-1.0 --gene-anno rice_gene_anno.gff3  --extra-anno rice_repeat_anno.bed --geneid-list rice_demo_genelist.txt --out demores_rice_genelist --pheno rice_demo_pheno.txt --extend 100
```
Visualizing a specific genomic region (Example: chr09:18669248-18673240):
```
fgpb --graph chr09_mc.og --ref-name IRGSP-1.0 --gene-anno rice_gene_anno.gff3  --extra-anno rice_repeat_anno.bed --region chr09:18669248-18673240 --out demores_rice_region --pheno rice_demo_pheno.txt
```
Visualizing multiple genomic regions (Example: chr09:7231334-7235878, chr09:17324231-17329297, chr09:15385163-15389649):
```
printf 'chr09\t7231333\t7235878\nchr09\t17324230\t17329297\nchr09\t15385162\t15389649\n' > rice_demo_region.bed
fgpb --graph chr09_mc.og --ref-name IRGSP-1.0 --gene-anno rice_gene_anno.gff3 --extra-anno rice_repeat_anno.bed --region-list rice_demo_region.bed --out  demores_rice_regionlist --pheno rice_demo_pheno.txt
```

### Variant data mode
For variant data (VCF file), the following parameters are required:
* Variant data (`--variant`)
* Reference genome sequence (`--ref-fa`)
* Target region (one of `--geneid`, `--geneid-list`, `--region`, `--region-list`)
* Reference genome gene annotation (`--gene-anno`, optional)
* Additional reference genome annotations (`--extra-anno`, optional)

Download demo data:
```bash
# Download human_demo.vcf.gz (634M)
wget https://cgm.sjtu.edu.cn/fGPB/demo/human/human_demo.vcf.gz
wget https://cgm.sjtu.edu.cn/fGPB/demo/human/human_demo.vcf.gz.tbi
# Download hg38.fa (3.1G)
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
# Download human_gene_anno.gff3 (84M)
wget https://cgm.sjtu.edu.cn/fGPB/demo/human/human_gene_anno.gff3
# Download human_repeat_anno.bed (232M)
wget https://cgm.sjtu.edu.cn/fGPB/demo/human/human_repeat_anno.bed
# Download human_demo_pheno.txt (45K)
wget https://cgm.sjtu.edu.cn/fGPB/demo/human/human_demo_pheno.txt
```

Visualizing a single gene (Example: *ENSG00000112695.13*):
```
fgpb --variant human_demo.vcf.gz --ref-fa hg38.fa --gene-anno human_gene_anno.gff3 --extra-anno human_repeat_anno.bed --geneid ENSG00000112695.13 --out demores_human_gene --pheno human_demo_pheno.txt --extend 1000
```

Visualizing multiple genes (Example: *ENSG00000171611.10*, *ENSG00000164430.17*, *ENSG00000205269.6*):
```
printf '%s\n' ENSG00000171611.10 ENSG00000164430.17 ENSG00000205269.6 > human_demo_genelist.txt
fgpb --variant human_demo.vcf.gz --ref-fa hg38.fa --gene-anno human_gene_anno.gff3 --extra-anno human_repeat_anno.bed --geneid-list human_demo_genelist.txt --out demores_human_genelist --pheno human_demo_pheno.txt --extend 1000
```

Visualizing a specific genomic region (Example: chr6:42915066-42940195):
```
fgpb --variant human_demo.vcf.gz --ref-fa hg38.fa --gene-anno human_gene_anno.gff3 --extra-anno human_repeat_anno.bed --region chr6:42915066-42940195 --out demores_human_region --pheno human_demo_pheno.txt
```

Visualizing multiple genomic regions (Example: chr6:2884917-2912669, chr6:18385699-18470573, chr6:73393730-73453504):
```
printf 'chr6\t2884916\t2912669\nchr6\t18385698\t18470573\nchr6\t73393729\t73453504\n' > human_demo_region.bed
fgpb --variant human_demo.vcf.gz --ref-fa hg38.fa --gene-anno human_gene_anno.gff3  --extra-anno human_repeat_anno.bed --region-list human_demo_region.bed --out demores_human_regionlist --pheno human_demo_pheno.txt
```





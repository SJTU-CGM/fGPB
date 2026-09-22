#!/usr/bin/perl

package GPBreport;

use strict;
use warnings;
use JSON;
use MIME::Base64 qw(encode_base64);


my %dna_to_bin = (
    'A' => '000',
    'T' => '001',
    'C' => '010',
    'G' => '011',
    'N' => '100'
);


sub compress_dna {
    my ($dna_seq) = @_;
    $dna_seq = uc($dna_seq);

    $dna_seq =~ s/[^ATCGN]/N/g;
    
    if (length($dna_seq) <= 4) {
        return encode_base64(pack('C', 0) . $dna_seq, '');
    }
    
    my $binary_str = '';
    foreach my $base (split //, $dna_seq) {
        $binary_str .= $dna_to_bin{$base};
    }
    my $padding = (8 - (length($binary_str) % 8)) % 8;
    $binary_str .= '0' x $padding;
    my $bytes = '';
    while (length($binary_str) >= 8) {
        my $byte_bin = substr($binary_str, 0, 8, '');
        $bytes .= pack('C', oct("0b$byte_bin"));
    }
    my $compressed = pack('C', 1) . pack('C', $padding) . $bytes;
    return encode_base64($compressed, '');
}


sub export_figure {
	my ($dir, $pfx, $final_nodes, $edge_out, $ref_info, $all_gene_pos, $bed_data, $pheno_meta, $dup_map) = @_;

	my %all_sample;
	$all_sample{$_} = 1 for map { split /,/, $_->[2] } @$final_nodes;
	my @sample_list = sort keys %all_sample;
	my $has_haplotype = (grep { /#/ } @sample_list) ? 1 : 0;
	my %sample_cnt;

	for my $node (@$final_nodes) {
		%sample_cnt = (); 
		$sample_cnt{$_}++ for split /,/, $node->[2];

		my @freq = map { $sample_cnt{$_} || 0 } @sample_list;
		pop @freq while @freq && $freq[-1] == 0;

		if (@freq == @sample_list && !grep { $_ != 1 } @freq) {
			$node->[2] = 'ALL';
		}else{
			$node->[2] = join ',', map { $_ || '' } @freq;
		}

		if($node->[1] ne ''){
			$node->[1] = compress_dna($node->[1]);
		}
	}
	
	my %out = (
                graphData => {
                        node => $final_nodes,
                        edge => $edge_out,
                        ref => $ref_info
                },
		sampleList => join(',', @sample_list),
		haplotypeNaming => $has_haplotype,
		dupMap => $dup_map,
                geneAnnoData => $all_gene_pos,
                bedAnnoData => $bed_data,
                phenoMeta => $pheno_meta
        );
	my $out_json = encode_json(\%out);
	my $showBed = ($bed_data eq '') ? 'false' : 'true';
	my $showGene = ($all_gene_pos eq '') ? 'false' : 'true';

	open my $fh_fig, '>', "$dir/$pfx.html" or die "Error: Can't write file '$dir/$pfx.html': $!\n";

	my $html = '<!doctype html>
<html>
<head>
<meta charset="utf-8">
<title>fGPB report</title>
<link rel="stylesheet" type="text/css" href="./css/bootstrap_5.1.1_css_bootstrap.min.css">
<link rel="stylesheet" type="text/css" href="./css/colorPicker.css">
<script type="text/javascript" charset="utf8" src="./js/jquery-3.5.1.js"></script> 
<script type="text/javascript" charset="utf8" src="./js/echarts.js"></script> 
<script type="text/javascript" charset="utf8" src="./js/bootstrap_5.1.1_js_bootstrap.bundle.min.js"></script> 
<script type="text/javascript" charset="utf8" src="./js/splitedData.js"></script> 
<script type="text/javascript" charset="utf8" src="./js/nodesLayout.js"></script> 
<script type="text/javascript" charset="utf8" src="./js/plottingSets.js"></script>
<script type="text/javascript" charset="utf8" src="./js/nodesInfo.js"></script>
<script type="text/javascript" charset="utf8" src="./js/colorPicker.js"></script>
</head>

<body>
<div class="container" style="margin-top:20px;width:90%">
  <div class="row">
    <table width="350px">
      <tr>
        <td><b>Chromosome:</b></td>
        <td id = "region_chr"></td>
      </tr>
      <tr>
        <td><b>Start:</b></td>
        <td id = "region_start"></td>
      </tr>
      <tr>
        <td><b>End:</b></td>
        <td id = "region_end"></td>
      </tr>
    </table>
  </div>
  <br>
  <hr>
  <br>
  <div class="row align-items-center">
    <div class="col-1  d-flex align-items-center gap-2">
      <label class="label-inline">Colors:</label>
    </div>
    <div class="col-4 d-flex align-items-center gap-2">
      <label class="label-inline">Reference node</label>
      <div class="custom-select flex-grow-1" id="refWrapper">
        <div class="select-trigger" onclick="toggleDropdown(\'ref\')">
          <div class="trigger-colors" id="refTriggerColors"></div>
          <span class="trigger-arrow">▼</span> </div>
        <div class="select-dropdown" id="refDropdown"></div>
      </div>
    </div>
    <div class="col-4  d-flex align-items-center gap-2">
      <label class="label-inline">Non-reference node</label>
      <div class="custom-select flex-grow-1" id="altWrapper">
        <div class="select-trigger" onclick="toggleDropdown(\'alt\')">
          <div class="trigger-colors" id="altTriggerColors"></div>
          <span class="trigger-arrow">▼</span> </div>
        <div class="select-dropdown" id="altDropdown"></div>
      </div>
    </div>
    <div class="col-1  d-flex align-items-center gap-2">
      <label class="label-inline">Edge</label>
      <div class="color-trigger" onclick="document.getElementById(\'arrowInput\').click()">
        <div class="color-box" id="arrowBox">
          <input type="color" id="arrowInput" value="#bababa" onchange="updatePreview()">
        </div>
      </div>
    </div>
    <div class="col-1  d-flex align-items-center gap-2">
      <label class="label-inline">Reference region</label>
      <div class="color-trigger" onclick="document.getElementById(\'refAreaInput\').click()">
        <div class="color-box" id="refAreaBox">
          <input type="color" id="refAreaInput" value="#ebf0e4" onchange="updatePreview()">
        </div>
      </div>
    </div>
  </div>
  <div id="phenoPanel" class="mt-4">
    <div class="d-flex align-items-center gap-3 flex-wrap">
      <div id="phenoControls" class="d-flex align-items-center gap-3">
        <label class="label-inline">Phenotype:</label>
        <select id="phenSelected" class="form-select form-select-sm w-auto">
        </select>
        <div class="d-flex align-items-center gap-2">
          <input class="form-check-input" type="checkbox" id="monoCheck" style="margin-top: 0;">
          <label class="form-check-label text-nowrap mb-0" for="monoCheck">Use node colors</label>
        </div>
      </div>
      <button id="refreshBtn" class="btn btn-primary btn-sm"> <i class="bi bi-arrow-clockwise" id="refreshIcon"></i> Refresh </button>
    </div>
  </div>
  <br>
  <div class="row">
    <div id="main" style="width:100%;height:800px"></div>
  </div>
  <p id="nodePanelText"> <small style="color: gray">Please click on the sequence blocks in the figure to view the corresponding sequence, path list.</small> </p>
  <div class="row" id="nodePanel" style="width:90%; margin:0 auto;">
    <div id="nodePanelSeq"></div>
    <div class="col-12 mb-3"></div>
    <div id="nodePanelSample"></div>
  </div>
  <br>
</div>
<script>
	
$(document).ready(function () {
  const chart = echarts.init(document.getElementById("main"));
    const jsonText = `'.$out_json.'`;
    const data = JSON.parse(jsonText);
    const graphData = data.graphData;
    const refInfo = graphData.ref;
    const geneData = data.geneAnnoData;
    const bedData = data.bedAnnoData;
    const phenoMetaData = data.phenoMeta;
    const dupMap = data.dupMap || {};
    const sampleUnit = data.haplotypeNaming ? "haplotype" : "genome";
    const nodeData = getNodeData(graphData.node, data.sampleList);

    var strucColors = {
      "gene": "#36648B",
      "transcript": "#A1A1A1",
      "exon": "#A1A1A1",
      "CDS": "#1368BD",
      "UTR3": "#68A3DE",
      "UTR5": "#68A3DE",
      "UTR": "#68A3DE"
    };

    let refNodeColors, altNodePalette, refAreaColor, arrowColor;
    const axisPointerColor = "#e38e28";
    const bedColor = "gray";

    function updateColorsFromPicker() {
      const config = getColorConfig();
      refNodeColors = config.refNodeColors;
      altNodePalette = config.altNodePalette;
      arrowColor = config.arrowColor;
      refAreaColor = config.refAreaColor;
    }
    updateColorsFromPicker();

    const arrowWidth = 0.5;
    const yMax = graphData.edge.filter(r=>r[0]==="0+").reduce((s,r)=>s+r[2],0);

    document.getElementById("region_chr").innerHTML = refInfo.chr
    document.getElementById("region_start").innerHTML = refInfo.start
    document.getElementById("region_end").innerHTML = refInfo.end

    let phenoNameList = [];
    if (phenoMetaData !== "") {
      phenoNameList = phenoMetaData.map(item => item[0]);
      const sel = document.getElementById("phenSelected");
      sel.innerHTML = "";
      phenoNameList.forEach(v => {
        sel.add(new Option(v, v));
      });
    } else {
      document.getElementById("phenoControls").style.setProperty("display", "none", "important");
    }

    const graphStart = Number(graphData.ref.start);
    const nodeSeq = getNodeSeq(nodeData);
    const nodeSample = new Map(nodeData.map(row => [row[0], row[2]]));
    const nodeSampleN = statNodeSampleN(graphData.edge);

    const refNodes = graphData.ref.nodes.split(",");
    const allEdges = graphData.edge;
    const allNodes = new Set(allEdges.flatMap(edge => [edge[0], edge[1]]));
    const sortedNodes = sortNodes(allNodes, allEdges, refNodes);
    const groupedEdge = getGroupedEdge(graphData, sortedNodes);
    mergeReverse(sortedNodes, groupedEdge, nodeSeq, nodeSample, nodeData, nodeSampleN);

    const nodeXPos = getNodeXPos(nodeSeq, sortedNodes);
    const xMin = 0.5;
    const xMax = nodeXPos[sortedNodes[sortedNodes.length - 1]].xend;
    const refMergedRange = getRefRange(nodeXPos, refNodes);
    const refMarkedArea = convertRangesToMarkedAreas(refMergedRange);

    const {
      nodeYRange,
      blockArrow,
      lineData
    } = layoutAll(sortedNodes, groupedEdge, nodeXPos, yMax);

    const newXLabel = splitXAxis(refMergedRange, graphStart);

    const {
      strucData,
      strucYtext
    } = getStrucData(geneData, strucColors);
    const strucSplitedData = splitStruc(strucData, refMergedRange, graphStart);
    const seriesGeneStrand = getStrandSeries(strucSplitedData);

    const {
      bedSplitedData,
      bedYtext
    } = splitBed(bedData, graphStart, refMergedRange, bedColor);

    let phenGroupName = "";
    let seriesPhenTrackList = [];
    let phenoYmax = yMax;
    let phenTrackMaxN = 0;
    let blockNode = null;

    function refreshChart() {
      updateColorsFromPicker(); 
      const {
        seriesBlockList,
        blockNode: bn,
        blockColors
      } = getBlockSeries(nodeXPos, nodeYRange, blockArrow, sortedNodes, nodeSampleN, arrowWidth, refNodes, refNodeColors, altNodePalette, dupMap);
      blockNode = bn;
      const seriesLineList = getLineSeries(lineData, arrowWidth, arrowColor);

      const {
        graphNodePos,
        seriesGraphNodeList
      } = getGraphNodeSeries(nodeXPos, blockNode, blockColors, 0.9, dupMap);
      const seriesGraphEdgeList = getGraphEdgeSeries(graphNodePos, groupedEdge);

      const seriesGeneStruc = getStrucSeries(strucSplitedData, refAreaColor, refMarkedArea);
      const seriesBed = getBedSeries(bedSplitedData, refAreaColor, refMarkedArea);

      const withColor = document.getElementById("monoCheck").checked;

      if (phenoMetaData !== "") {
        const curPhen = document.getElementById("phenSelected").value;
        const curPhenIdx = phenoNameList.indexOf(curPhen);
        const curPhenInfo = phenoMetaData.find(item => item[0] === curPhen);
        const phenData = Object.fromEntries(
          nodeData.map(a => {
            const arr = a[3][curPhenIdx];
            return [a[0], arr.slice(1).concat(arr[0])];
          })
        );
        if (curPhenInfo[1]) {
          phenGroupName = [phenoMetaData[curPhenIdx][0]];
          phenoYmax = Math.max(...Object.values(phenData).flatMap(a => a.slice(2, 5)));
          seriesPhenTrackList = getConPhenTrackSeries(phenData, phenGroupName, nodeXPos, 4, blockNode, blockColors, refAreaColor, refMarkedArea, withColor);
        } else {
          phenGroupName = [...phenoMetaData[curPhenIdx][2]];
          const hasLastGt0 = Object.values(phenData)
            .some(arr => Array.isArray(arr) && arr.length > 0 && arr.at(-1) > 0);
          if (hasLastGt0) {
            phenGroupName.push("MISSING");
          }
          phenoYmax = Math.max(...Object.values(phenData).flat());
          seriesPhenTrackList = getDisPhenTrackSeries(phenData, phenGroupName, nodeXPos, 4, blockNode, blockColors, refAreaColor, refMarkedArea, withColor);
        }
        phenTrackMaxN = phenGroupName.length;

      }
      const {
        option,
        totalHeight
      } = getChartOption(axisPointerColor, xMin, xMax, newXLabel, yMax, strucYtext, bedYtext, phenTrackMaxN, phenGroupName, phenoYmax, seriesBlockList, seriesGraphNodeList, seriesGraphEdgeList, seriesLineList, seriesGeneStruc, seriesGeneStrand, seriesBed, seriesPhenTrackList, '.$showGene.', '.$showBed.');

      const dom = document.getElementById("main");
      dom.style.height = totalHeight + "px";

      chart.resize();
      chart.clear();
      chart.setOption(option);
    }

    refreshChart();
    document.getElementById("refreshBtn").addEventListener("click", refreshChart);

    chart.off("click").on("click", function (params) {
      if (params.componentSubType == "custom" && params.value.length == 2) {
        const curNode = blockNode[params.seriesIndex];
        if (curNode) {
          document.getElementById("nodePanelText").hidden = true;
          document.getElementById("nodePanel").hidden = false;
          fillNodePanel(curNode, nodeSeq, nodeSample, refInfo, dupMap, sampleUnit);
        }
      }
    });
})

	
</script>
</body>
</html>';

	print $fh_fig $html;
	close $fh_fig;

}


1;

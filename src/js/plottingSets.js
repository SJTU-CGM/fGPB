function getGraphSeries(nodeXPos, groupedEdge, blockNode, blockColors) {

  const nodeHeight = 10;
  const nodeWidhtPer = 0.2;
  const nodeData = [];
  const linkData = [];

  Object.keys(nodeXPos).forEach(key => {
    const curNode = nodeXPos[key];
    const curWidth = curNode.xend - curNode.xstart;
    const curX = curNode.xstart + curWidth / 2;
    const curY = 0;
    nodeData.push({
      name: key,
      value: [curX, curY],
      symbolSize: [curWidth * nodeWidhtPer, nodeHeight],
      itemStyle: {
        color: blockColors[blockNode.indexOf(key)]
      }
    })
    if (groupedEdge[key]) {
      groupedEdge[key].forEach(item => {
        linkData.push({
          source: item[0],
          target: item[1],
          lineStyle: {
            width: 3,
            curveness: 0.1
          }
        });
      });
    }
  });

  return {
    type: 'graph',
    layout: 'none',
    coordinateSystem: 'cartesian2d',
    symbol: 'rect',
    draggable: true,
    symbolSize: 1,
    roam: true,
    label: {
      show: false
    },
    edgeSymbolSize: [4, 10],
    data: nodeData,
    links: linkData,
    lineStyle: {
      opacity: 0.9,
      width: 1,
      curveness: 0
    },
    xAxisIndex: 0,
    yAxisIndex: 0,
    clip: true,
    large: true,
    emphasis: {
      disabled: true
    },
    tooltip: {
      show: false
    }
  };

}


function getGraphNodeSeries(nodeXPos, blockNode, blockColors, blockWidth, dupMap) {

  const graphDataList = [];
  const graphNodePos = {};
  const blockHeight = 1;

  Object.keys(nodeXPos).forEach(key => {
    const curWidth = nodeXPos[key].xend - nodeXPos[key].xstart;
    const curP = 0.9 + 0.1 * (1 - 1 / (1 + Math.log2(curWidth + 1)));
    const curW = curWidth * (1 - curP)
    if (curWidth > 0) {
      graphDataList.push([
        [nodeXPos[key].xstart + curW / 2, 0 - blockHeight / 2],
        [nodeXPos[key].xstart + curW / 2, 0 + blockHeight / 2],
        [nodeXPos[key].xend - curW / 2, 0 + blockHeight / 2],
        [nodeXPos[key].xend - curW / 2, 0 - blockHeight / 2],
        blockColors[blockNode.indexOf(key)],
        key
      ]);
      graphNodePos[key] = {
        xstart: nodeXPos[key].xstart + curW / 2,
        xend: nodeXPos[key].xend - curW / 2
      }
    }
  });

  const gKeys = Object.keys(graphNodePos).sort((a, b) => graphNodePos[a].xstart - graphNodePos[b].xstart);
  graphDataList.forEach(entry => {
    const i = gKeys.indexOf(entry[5]);
    const gapB = i > 0 ? graphNodePos[entry[5]].xstart - graphNodePos[gKeys[i - 1]].xend : Infinity;
    const gapA = i < gKeys.length - 1 ? graphNodePos[gKeys[i + 1]].xstart - graphNodePos[entry[5]].xend : Infinity;
    entry.push(gapB, gapA);
  });

  const seriesGraphNodeList = graphDataList.map(data => {
    return {
      type: 'custom',
      renderItem: function (params, api) {
        if (params.context.rendered) {
          return;
        }
        params.context.rendered = true;
        let points = [];
        for (let i = 0; i < 4; i++) {
          points.push(api.coord(data[i]));
        }
        const x0 = points[0][0], x1 = points[2][0];
        const y0 = Math.min(points[0][1], points[2][1]);
        const y1 = Math.max(points[0][1], points[2][1]);
        const w = x1 - x0, h = y1 - y0;
        const cy = y0 + h / 2;
        const pxPerUnit = api.coord([1, 0])[0] - api.coord([0, 0])[0];
        const gapPx = (data[5].endsWith('-') ? data[6] : data[7]) * pxPerUnit;
        const ext = Math.min(h * 0.8, 6, gapPx);
        const ts = Math.min(h * 0.8, w * 0.3, 8);
        let shapePoints;
        if (w < 4) {
          shapePoints = points;
        } else if (data[5].endsWith('-')) {
          shapePoints = ext >= 1.5 ?
            [[x1, y0], [x1, y1], [x0, y1], [x0 - ext, cy], [x0, y0]] :
            [[x1, y0], [x1, y1], [x0 + ts, y1], [x0, cy], [x0 + ts, y0]];
        } else {
          shapePoints = ext >= 1.5 ?
            [[x0, y0], [x0, y1], [x1, y1], [x1 + ext, cy], [x1, y0]] :
            [[x0, y0], [x0, y1], [x1 - ts, y1], [x1, cy], [x1 - ts, y0]];
        }
        const children = [{
          type: 'polygon',
          transition: ['shape'],
          shape: {
            points: shapePoints
          },
          style: {
            fill: data[4]
          }
        }];
        if (data[5].endsWith('-') || (dupMap && (dupMap[data[5]] || data[5].split(',').some(p => dupMap[p])))) {
          const dots = [];
          for (let yy = y0 + 2; yy < y1; yy += 4) {
            for (let xx = x0 + 2; xx < x1; xx += 4) {
              dots.push({
                type: 'circle',
                shape: {
                  cx: xx,
                  cy: yy,
                  r: 0.8
                },
                style: {
                  fill: 'rgba(0,0,0,0.35)'
                }
              });
            }
          }
          children.push({
            type: 'group',
            children: dots,
            clipPath: {
              type: 'polygon',
              shape: {
                points: shapePoints
              }
            }
          });
        }
        return {
          type: 'group',
          children: children
        };
      },
      xAxisIndex: 0,
      yAxisIndex: 0,
      encode: {
        x: 0,
        y: 0
      },
      clip: true,
      data: data[0, 1, 2, 3],
      large: true,
      emphasis: {
        disabled: true
      },
      tooltip: {
        show: false
      }

    }
  });

  return {
    graphNodePos,
    seriesGraphNodeList
  };

}


function getGraphEdgeSeries(graphNodePos, groupedEdge) {

  const seriesGraphEdgeList = Object.keys(groupedEdge).map(key => {
    const graphNodes = [];
    const graphLinks = [];

    if (graphNodePos[key]) {
      graphNodes.push({
        name: key,
        value: [graphNodePos[key].xend, 0.5]
      });
      groupedEdge[key].forEach(item => {
        if (item[1] != "Inf+") {
          graphNodes.push({
            name: item[1],
            value: [graphNodePos[item[1]].xstart, 0.5]
          })
          graphLinks.push({
            source: item[0],
            target: item[1],
            lineStyle: {
              width: 3,
              curveness: 0.4
            }
          });
        }
      });
    }

    return {
      type: 'graph',
      layout: 'none',
      coordinateSystem: 'cartesian2d',
      symbol: 'none',
      symbolSize: 1,
      roam: true,
      label: {
        show: false
      },
      edgeSymbolSize: [4, 10],
      data: graphNodes,
      links: graphLinks,
      lineStyle: {
        opacity: 0.9,
        width: 1,
        curveness: 0
      },
      xAxisIndex: 0,
      yAxisIndex: 0,
      clip: true,
      large: true,
      emphasis: {
        disabled: true
      },
      tooltip: {
        show: false
      }
    }
  });

  return seriesGraphEdgeList;

}


function getConPhenTrackSeries(phenData, phenGroupName, nodeXPos, axisIdx, blockNode, blockColors, refAreaColor, refMarkedArea, withColor = false) {

  const phenTrackSta = [];
  Object.keys(nodeXPos).forEach(node => {
    if (node != "0+" && node != "Inf+") {
      for (let x = nodeXPos[node].xstart; x <= nodeXPos[node].xend; x += 0.5) {
        phenTrackSta.unshift({
          x: x,
          l: phenData[node][1],
          value: phenData[node][2],
          u: phenData[node][3],
          n: phenData[node][5],
          missingn: phenData[node][6],
        });
      }
    }
  });

  let medianSeries;
  if (!withColor) {
    const monoData = [];
    const lineColors = [];

    Object.keys(nodeXPos).forEach(node => {
      if (node != "0+" && node != "Inf+") {
        if (monoData.length > 0) {
          monoData.push([null, null]);
          lineColors.push('transparent');
        }
        const nodeColor = '#36648B';
        for (let x = nodeXPos[node].xstart; x <= nodeXPos[node].xend; x += 0.5) {
          monoData.push([x, phenData[node][5] === 0 ? null : phenData[node][2]]);
          lineColors.push(nodeColor);
        }
      }
    });

    medianSeries = [{
      type: 'line',
      name: 'Median:',
      xAxisIndex: axisIdx,
      yAxisIndex: axisIdx,
      data: monoData,
      showSymbol: false,
      large: true,
      clip: true,
      animation: false,
      connectNulls: false,
      lineStyle: {
        color: '#36648B',
        width: 2
      },
      tooltip: {
        trigger: 'axis'
      }
    }];
  } else {
    medianSeries = [];
    Object.keys(nodeXPos).forEach((node, idx) => {
      if (node != "0+" && node != "Inf+") {
        const nodeData = [];
        for (let x = nodeXPos[node].xstart; x <= nodeXPos[node].xend; x += 0.5) {
          nodeData.push([x, phenData[node][5] === 0 ? null : phenData[node][2]]);
        }
        medianSeries.push({
          type: 'line',
          name: 'Median:',
          xAxisIndex: axisIdx,
          yAxisIndex: axisIdx,
          data: nodeData,
          showSymbol: false,
          large: true,
          clip: true,
          animation: false,
          lineStyle: {
            color: blockColors[blockNode.indexOf(node)],
            width: 2
          },
          tooltip: {
            trigger: 'axis',
            formatter: function (params) {
              return 'Node ' + node + '<br/>Median: ' + params[0].value[1];
            }
          },
          legendHoverLink: false
        });
      }
    });
  }

  const seriesPhenTrackList = [{
      name: 'Sample number',
      type: 'line',
      xAxisIndex: axisIdx,
      yAxisIndex: axisIdx,
      data: phenTrackSta.map(item => [item.x, item.n]),
      lineStyle: {
        opacity: 0
      },
      itemStyle: {
        opacity: 0
      },
      symbol: 'none'
    },
    {
      name: 'Missing',
      type: 'line',
      xAxisIndex: axisIdx,
      yAxisIndex: axisIdx,
      data: phenTrackSta.map(item => [item.x, item.missingn]),
      lineStyle: {
        opacity: 0
      },
      itemStyle: {
        opacity: 0
      },
      symbol: 'none'
    },
    {
      type: 'line',
      name: 'Upper quartile:',
      xAxisIndex: axisIdx,
      yAxisIndex: axisIdx,
      data: phenTrackSta.map(item => [item.x, item.n === 0 ? null : item.u]),
      showSymbol: false,
      large: true,
      clip: true,
      animation: false,
      lineStyle: {
        opacity: 0
      },
      areaStyle: {
        color: '#ccc',
        opacity: 1
      },
      symbol: 'none',
      emphasis: {
        focus: 'none',
        scale: false
      },
      silent: true
    },
    ...medianSeries,
    {
      type: 'line',
      name: 'Lower quartile:',
      xAxisIndex: axisIdx,
      yAxisIndex: axisIdx,
      data: phenTrackSta.map(item => [item.x, item.n === 0 ? null : item.l]),
      showSymbol: false,
      large: true,
      clip: true,
      animation: false,
      lineStyle: {
        opacity: 0
      },
      areaStyle: {
        color: 'white',
        opacity: 1
      },
      symbol: 'none',
      silent: true
    }
  ];

  return seriesPhenTrackList;
}


function getDisPhenTrackSeries(phenData, phenGroupName, nodeXPos, axisIdx, blockNode, blockColors, refAreaColor, refMarkedArea, withColor = false) {

  const phenTrackData = [];
  phenGroupName.forEach((pheno, phenIdx) => {
    phenTrackData[phenIdx] = [];
    Object.keys(nodeXPos).forEach((node, idx) => {
      if (node != "0+" && node != "Inf+") {
        phenTrackData[phenIdx].push({
          name: node,
          value: [nodeXPos[node].xstart, nodeXPos[node].xend, phenData[node][phenIdx], ''],
          itemStyle: {
            color: withColor ? blockColors[blockNode.indexOf(node)] : 'gray'
          }
        });
      }
    });
  });

  const seriesPhenTrackList = phenTrackData.map((data, idx) => {
    return {
      type: 'custom',
      renderItem: function (params, api) {
        var yValue = api.value(2);
        var start = api.coord([api.value(0), yValue]);
        var size = api.size([api.value(1) - api.value(0), yValue]);
        var rectShape = echarts.graphic.clipRectByRect({
          x: start[0],
          y: start[1],
          width: size[0],
          height: size[1]
        }, {
          x: params.coordSys.x,
          y: params.coordSys.y,
          width: params.coordSys.width,
          height: params.coordSys.height
        });
        return (rectShape && {
          type: 'rect',
          transition: ['shape'],
          shape: rectShape,
          style: api.style()
        });
      },
      xAxisIndex: axisIdx + idx,
      yAxisIndex: axisIdx + idx,
      large: true,
      encode: {
        x: 3,
        y: 3
      },
      data: data,
      markArea: {
        itemStyle: {
          color: 'white',
        },
        data: refMarkedArea,
        emphasis: {
          disabled: true
        }
      },
      tooltip: {
        trigger: "item",
        formatter: function (param) {
          return `
           <div style="font-weight:bold">Node ${param.name}</div>
           <div>Value: <b>${param.value[2]}</b></div>
         `;
        }
      }
    }
  });

  return seriesPhenTrackList;
}


function getBlockSeries(nodeXPos, nodeYRange, blockArrow, sortedNodes, nodeSampleN, arrowWidth, refNodes, refNodeColors, altNodePalette, dupMap) {

  const blockDataList = [];
  const blockColors = [];
  const blockNode = [];
  let altNodeNum = 0;
  let refNodeNum = 0;
  const maxSampleN = Math.max(...Array.from(nodeSampleN.values()));

  function colorInterpolate(color1, color2, ratio) {
    function parseColor(color) {
      color = color.trim();
      if (color.startsWith('#')) {
        let hex = color.slice(1);
        if (hex.length === 3) {
          hex = hex.split('').map(c => c + c).join('');
        }
        return [
          parseInt(hex.slice(0, 2), 16),
          parseInt(hex.slice(2, 4), 16),
          parseInt(hex.slice(4, 6), 16)
        ];
      } else if (color.startsWith('rgb')) {
        const nums = color.match(/\d+/g);
        return [parseInt(nums[0]), parseInt(nums[1]), parseInt(nums[2])];
      }
      return [128, 128, 128];
    }
    const c1 = parseColor(color1);
    const c2 = parseColor(color2);
    const r = Math.round(c1[0] + (c2[0] - c1[0]) * ratio);
    const g = Math.round(c1[1] + (c2[1] - c1[1]) * ratio);
    const b = Math.round(c1[2] + (c2[2] - c1[2]) * ratio);
    return 'rgb(' + r + ',' + g + ',' + b + ')';
  }

  sortedNodes.forEach((node, idx) => {
    if (node != "0+" && node != "Inf+") {
      const ratio = nodeSampleN.get(node) / maxSampleN;
      let curColor;

      if (refNodes.includes(node)) {
        if (refNodeColors && refNodeColors.length === 2) {
          curColor = colorInterpolate(refNodeColors[0], refNodeColors[1], ratio);
        } else if (refNodeColors && refNodeColors.length > 0) {
          curColor = refNodeColors[refNodeNum % refNodeColors.length];
        } else {
          const gray = 255 - 60 - ratio * 100;
          curColor = 'rgb(' + gray + ',' + gray + ',' + gray + ')';
        }
        refNodeNum++;
      } else {
        if (altNodePalette && altNodePalette.length === 2) {
          curColor = colorInterpolate(altNodePalette[0], altNodePalette[1], ratio);
        } else if (altNodePalette && altNodePalette.length > 0) {
          curColor = altNodePalette[altNodeNum % altNodePalette.length];
        } else {
          curColor = '#cccccc';
        }
        altNodeNum++;
      }

      if (blockArrow[node]) {
        arrPos = blockArrow[node];
        nodeYRange[node].forEach(block => {
          const points = [
            [nodeXPos[node].xstart, -block.ystart],
            [nodeXPos[node].xstart, -block.yend],
            [nodeXPos[node].xend, -block.yend]
          ];
          arrPos.forEach(item => {
            if (item[0] > block.ystart && item[1] < block.yend) {
              points.push([nodeXPos[node].xend, -item[0]]);
              points.push([nodeXPos[node].xend + arrowWidth, -(item[0] + item[1]) / 2]);
              points.push([nodeXPos[node].xend, -item[1]]);
            }
          });
          points.push([nodeXPos[node].xend, -block.ystart]);
          blockDataList.unshift(points);
          blockNode.unshift(node);
          blockColors.unshift(curColor);
        });
      } else {
        nodeYRange[node].forEach(block => {
          blockDataList.unshift([
            [nodeXPos[node].xstart, -block.ystart],
            [nodeXPos[node].xstart, -block.yend],
            [nodeXPos[node].xend, -block.yend],
            [nodeXPos[node].xend, -block.ystart]
          ]);
          blockNode.unshift(node);
          blockColors.unshift(curColor);
        });
      }
    }
  });

  seriesBlockList = blockDataList.map((data, idx) => {
    return {
      type: 'custom',
      renderItem: function (params, api) {
        if (params.context.rendered) {
          return;
        }
        params.context.rendered = true;
        let points = [];
        for (let i = 0; i < data.length; i++) {
          points.push(api.coord(data[i]));
        }
        let color = blockColors[idx];
        return {
          type: 'polygon',
          transition: ['shape'],
          shape: {
            points: points
          },
          style: {
            fill: color,
            stroke: echarts.color.lift(color, 0.1)
          }
        };
      },
      xAxisIndex: 1,
      yAxisIndex: 1,
      encode: {
        x: 1,
        y: 1
      },
      clip: true,
      data: data,
      large: true,
      tooltip: {
        trigger: "item",
        formatter: function (params) {
          curVal = params.value
          const bn = blockNode[idx];
          const dupNote = (dupMap && dupMap[bn]) ? ' <span style="color:#999">(duplicated from node ' + dupMap[bn] + ')</span>' : '';
          return '<b>Node ' + bn + '</b>' + dupNote;
        }
      }
    }
  })

  return {
    seriesBlockList,
    blockNode,
    blockColors
  };

}


function getLineSeries(lineData, arrowWidth, arrowColor, stripeWidth, stripeGap) {

  stripeWidth = 2;
  stripeGap = 1;
  const lineDataList = [];
  lineData.forEach(item => {
    lineDataList.push([
      [item[0], -item[2]],
      [item[0] + arrowWidth, -(item[3] + item[2]) / 2],
      [item[0], -item[3]],
      [item[1], -item[3]],
      [item[1] + arrowWidth, -(item[3] + item[2]) / 2],
      [item[1], -item[2]]
    ]);
  });

  const seriesLineList = lineDataList.map(data => ({
    type: 'custom',
    xAxisIndex: 1,
    yAxisIndex: 1,
    clip: true,
    data: data,
    large: true,
    emphasis: {
      disabled: true
    },
    tooltip: {
      show: false
    },
    renderItem: function (params, api) {
      if (params.context.rendered) return;
      params.context.rendered = true;
      const points = data.map(pt => api.coord(pt));
      return {
        type: 'polygon',
        shape: {
          points: points
        },
        style: {
          decal: {
            dashArrayX: [1, stripeWidth + stripeGap],
            dashArrayY: [stripeWidth, stripeGap],
            rotation: Math.PI / 4,
            color: arrowColor,
            backgroundColor: 'white'
          },
          fill: 'transparent'
        }
      }
    }
  }));

  return seriesLineList;

}


function renderItem(params, api) {

  var categoryIndex = api.value(0);
  var start = api.coord([api.value(1), categoryIndex]);
  var end = api.coord([api.value(2), categoryIndex]);
  var type = api.value(3);

  if (type == 'gene') {
    var height = api.size([0, 1])[1] * 0.2;
  } else if (type == 'transcript') {
    var height = api.size([0, 1])[1] * 0.1;
  } else if (type == 'exon') {
    var height = api.size([0, 1])[1] * 0.3;
  } else if (type == 'CDS') {
    var height = api.size([0, 1])[1] * 0.6;
  } else if (type == 'UTR3') {
    var height = api.size([0, 1])[1] * 0.3;
  } else if (type == 'UTR5') {
    var height = api.size([0, 1])[1] * 0.3;
  } else if (type == 'UTR') {
    var height = api.size([0, 1])[1] * 0.3;
  } else if (type == 'bed') {
    var height = api.size([0, 1])[1] * 0.5;
  } else {
    var height = api.size([0, 1])[1] * 0;
  }

  var rectShape = echarts.graphic.clipRectByRect({
    x: start[0],
    y: start[1] - height / 2,
    width: end[0] - start[0],
    height: height
  }, {
    x: params.coordSys.x,
    y: params.coordSys.y,
    width: params.coordSys.width,
    height: params.coordSys.height
  });

  return (
    rectShape && {
      type: 'rect',
      transition: ['shape'],
      shape: rectShape,
      style: api.style()
    }
  );

}


function convertRangesToMarkedAreas(ranges) {
  return ranges.map(item => [{
      xAxis: item[0]
    },
    {
      xAxis: item[1]
    }
  ]);
}


function getStrucSeries(strucSplitedData, refAreaColor, refMarkedArea) {
  return {
    type: 'custom',
    renderItem: renderItem,
    itemStyle: {
      opacity: 1
    },
    large: true,
    xAxisIndex: 2,
    yAxisIndex: 2,
    encode: {
      x: 4,
      y: 4
    },
    data: strucSplitedData,
    markArea: {
      itemStyle: {
        color: refAreaColor
      },
      data: refMarkedArea,
      emphasis: {
        disabled: true
      }
    },
    tooltip: {
      trigger: 'item',
      backgroundColor: 'rgba(255,255,255,.95)',
      borderColor: '#ddd',
      borderWidth: 1,
      padding: [10, 14],
      extraCssText: 'border-radius:4px;box-shadow:0 1px 4px rgba(0,0,0,.15);',
      formatter: (params) => {
        const curVal = params.value;
        const part = curVal[3];
        const fixed = part === part.toLowerCase()
          ? part.charAt(0).toUpperCase() + part.slice(1)
          : part;
        const hasExonInfo = curVal.length > 10 && curVal[10];
        const typeColorMap = {
          'Gene': '#e74c3c',
          'Transcript': '#e74c3c',
          'Exon': '#A1A1A1',
          'CDS': '#1368BD',
          'UTR5': '#68A3DE',
          'UTR3': '#68A3DE',
	  'UTR': '#68A3DE'
        };
        const tagColor = typeColorMap[fixed] || '#666';

        return `
          <div style="line-height:1.7;min-width:220px;">
            <div style="display:flex;align-items:center;gap:8px;margin-bottom:6px;">
              <span style="background:${tagColor};color:#fff;padding:1px 8px;border-radius:3px;font-size:12px;font-weight:bold;">${fixed}</span>
              <span style="font-weight:bold;color:#333;">${curVal[7]}</span>
            </div>
            <div style="background:#f8f9fa;padding:8px 10px;border-radius:4px;">
              <div style="display:flex;gap:12px;margin-bottom:2px;">
                <span style="color:#666;min-width:45px;">ID</span><b>${curVal[4]}</b>
              </div>
              <div style="display:flex;gap:12px;margin-bottom:2px;">
                <span style="color:#666;min-width:45px;">Chr</span><b>${curVal[8]}</b>
              </div>
              <div style="display:flex;gap:12px;margin-bottom:2px;">
                <span style="color:#666;min-width:45px;">Start</span><b>${curVal[5]}</b>
              </div>
              <div style="display:flex;gap:12px;margin-bottom:2px;">
                <span style="color:#666;min-width:45px;">End</span><b>${curVal[6]}</b>
              </div>
              <div style="display:flex;gap:12px;margin-bottom:2px;">
                <span style="color:#666;min-width:45px;">Strand</span><b>${decodeURIComponent(curVal[9])}</b>
              </div>
              ${hasExonInfo ? `<div style="display:flex;gap:12px;margin-bottom:2px;">
                <span style="color:#666;min-width:45px;">Exon</span><b>${curVal[10]}</b>
              </div>` : ''}
            </div>
          </div>`;
      }
    }
  };
}

const PATH_RIGHT = 'M2 16L12 16M12 16L8 12M12 16L8 20';
const PATH_LEFT = 'M30 16L20 16M20 16L24 12M20 16L24 20';

function renderArrow(param, api) {
  const point = api.coord([
    api.value(2),
    api.value(1)
  ]);
  const strand = api.value(3);
  const arrowSize = 10;
  return {
    type: 'path',
    shape: {
      pathData: strand === '+' ? PATH_RIGHT : PATH_LEFT,
      x: strand === '+' ? 0 : -arrowSize,
      y: -arrowSize / 2,
      width: arrowSize,
      height: arrowSize
    },
    rotation: 0,
    position: point,
    style: api.style({
      stroke: 'black',
      lineWidth: 1,
      fill: 'none'
    })
  };
}


function getStrandSeries(strucSplitedData) {
  const arrowPos = Array.from(
    strucSplitedData.reduce((m, {
      value: v
    }) => {
      if (v[3] === 'gene') {
        const gene = v[7];
        const old = m.get(gene);
        const y = Number(v[0]);
        let x;
        if (v[9] === '+') {
          x = Number(v[2]);
          m.set(gene, {
            y,
            x: old === undefined ? x : Math.max(old.x, x),
            strand: v[9]
          });
        } else if (v[9] === '-') {
          x = Number(v[1]);
          m.set(gene, {
            y,
            x: old === undefined ? x : Math.min(old.x, x),
            strand: v[9]
          });
        }
      }
      return m;
    }, new Map())
  ).map(([gene, {
    y,
    x,
    strand
  }]) => [gene, y, x, strand]);
  return {
    type: 'custom',
    renderItem: renderArrow,
    data: arrowPos,
    xAxisIndex: 2,
    yAxisIndex: 2,
    encode: {
      x: 2,
      y: 1
    },
    tooltip: {
      show: false
    }
  };
}


function getBedSeries(bedSplitedData, refAreaColor, refMarkedArea) {
  return {
    type: 'custom',
    renderItem: renderItem,
    itemStyle: {
      opacity: 1
    },
    large: true,
    xAxisIndex: 3,
    yAxisIndex: 3,
    encode: {
      x: 4,
      y: 4
    },
    data: bedSplitedData,
    markArea: {
      itemStyle: {
        color: refAreaColor
      },
      data: refMarkedArea,
      emphasis: {
        disabled: true
      }
    },
    tooltip: {
      trigger: 'item',
      formatter: function (params) {
        curVal = params.value;
        return '<b>' + curVal[4] + '</b></br>Start: <b>' + curVal[5] + '</b></br>End: <b>' + curVal[6] + '</b>';

      }
    }
  };
}


function generateGridConfig(phenTrackCount, showGene, showBed) {

  const sliderHeight = 80;
  const barHeight = 80;
  const blockHeight = 200;
  const geneHeight = showGene === false ? 0 : 200; 
  const bedHeight = showBed === false ? 0 : 30;
  const phenHeight = 50;
  const gap = 10;

  const totalHeight = sliderHeight + barHeight + blockHeight + geneHeight + bedHeight
    + (phenTrackCount * phenHeight)
    + ((4 + phenTrackCount) * gap);

  const grids = [];
  let currentTop = sliderHeight;

  grids.push({
    left: '15%',
    top: currentTop + 'px',
    height: barHeight + 'px',
    width: '80%'
  });
  currentTop += barHeight + gap;

  grids.push({
    left: '15%',
    top: currentTop + 'px',
    height: blockHeight + 'px',
    width: '80%'
  });
  currentTop += blockHeight + gap;

  grids.push({
    left: '15%',
    top: currentTop + 'px',
    height: geneHeight + 'px',
    width: '80%'
  });
  currentTop += geneHeight + gap;

  grids.push({
    left: '15%',
    top: currentTop + 'px',
    height: bedHeight + 'px',
    width: '80%'
  });
  currentTop += bedHeight + gap;

  for (let i = 0; i < phenTrackCount; i++) {
    grids.push({
      left: '15%',
      top: currentTop + 'px',
      height: phenHeight + 'px',
      width: '80%'
    });
    currentTop += phenHeight + gap;
  }

  return {
    totalHeight: totalHeight,
    grid: grids
  };

}


function generateXAxisConfig(phenTrackCount, axisIdx, xMin, xMax, newXLabel) {

  const xAxis = [{
      gridIndex: 0,
      type: 'value',
      scale: true,
      position: 'top',
      min: xMin,
      max: xMax,
      minInterval: 1,
      show: false,
      splitLine: {
        show: false
      }

    }, {
      gridIndex: 1,
      scale: true,
      position: 'top',
      type: "value",
      min: xMin,
      max: xMax,
      minInterval: 1,
      splitLine: {
        show: false
      },
      axisPointer: {
        label: {
          show: true
        }
      },
      axisTick: {
        show: true
      },
      axisLabel: {
        show: true,
        fontWeight: 'bold',
        color: 'black',
        fontSize: 10,
        formatter: function (value) {
          return newXLabel[value - 1];
        }
      },
    }, {
      gridIndex: 2,
      scale: true,
      min: xMin,
      max: xMax,
      minInterval: 1,
      show: false,
      splitLine: {
        show: false
      }
    },
    {
      gridIndex: 3,
      scale: true,
      min: xMin,
      max: xMax,
      minInterval: 1,
      show: false,
      splitLine: {
        show: false
      }
    }
  ];

  const phenTrackAxis = idx => ({
    gridIndex: idx,
    type: 'value',
    scale: true,
    min: xMin,
    max: xMax,
    minInterval: 1,
    show: false,
    splitLine: {
      show: false
    }
  });

  for (let i = 0; i < phenTrackCount; i++) xAxis.push(phenTrackAxis(axisIdx + i));

  return xAxis;

}


function generateYAxisConfig(phenTrackCount, axisIdx, yMax, strucYtext, bedYtext, phenGroupName, phenoYmax, showGene, showBed) {

  const yAxis = [{
      gridIndex: 0,
      type: 'value',
      min: -3,
      max: 3,
      show: false,
      splitLine: {
        show: false
      }
    }, {
      gridIndex: 1,
      name: 'Path Number',
      nameLocation: 'middle',
      nameGap: 30,
      min: -yMax,
      max: 0,
      splitLine: {
        show: false
      },
      nameTextStyle: {
        fontWeight: 'bold',
        color: 'black'
      },
      axisLabel: {
        show: true,
        fontWeight: 'bold',
        color: 'black',
        fontSize: 10,
        formatter: function (value) {
          return value * (-1);
        }
      }
    },
    {
      gridIndex: 2,
      data: strucYtext,
      show: showGene,
      inverse: true,
      type: 'category',
      name: "Gene annotation",
      nameRotate: 0,
      nameTextStyle: {
        fontWeight: 'bold',
        color: 'black'
      },
      nameLocation: 'middle',
      nameGap: 30,
      splitLine: {
        show: false
      },
      axisTick: {
        alignWithLabel: true,
        interval: 0
      },
      axisLabel: {
        show: true,
        fontWeight: 'bold',
        color: 'black',
        fontSize: 10
      }
    },
    {
      gridIndex: 3,
      data: bedYtext,
      show: showBed,
      inverse: true,
      type: 'category',
      name: "Additional annotation",
      nameRotate: 0,
      nameTextStyle: {
        fontWeight: 'bold',
        color: 'black'
      },
      nameLocation: 'middle',
      nameGap: 30,
      splitLine: {
        show: false
      },
      axisTick: {
        alignWithLabel: true,
        interval: 0
      },
      axisLabel: {
        show: true,
        fontWeight: 'bold',
        color: 'black',
        fontSize: 10
      }
    }
  ];

  const phenTrackAxis = idx => ({
    gridIndex: idx,
    type: 'value',
    name: phenGroupName[idx - axisIdx],
    min: 0,
    max: phenoYmax,
    show: true,
    nameRotate: 0,
    nameTextStyle: {
      fontWeight: 'bold',
      color: 'black'
    },
    nameLocation: 'middle',
    nameGap: 30,
    splitLine: {
      show: false
    }
  });

  for (let i = 0; i < phenGroupName.length; i++) yAxis.push(phenTrackAxis(axisIdx + i));

  return yAxis;

}


function getChartOption(axisPointerColor, xMin, xMax, newXLabel, yMax, strucYtext, bedYtext, phenTrackMaxN, phenGroupName, phenoYmax, seriesBlockList, seriesGraphNodeList, seriesGraphEdgeList, seriesLineList, seriesGeneStruc, seriesGeneStrand, seriesBed, seriesPhenTrackList, showGene, showBed) {

  const gridHeight = generateGridConfig(phenTrackMaxN, showGene, showBed);
  const totalHeight = gridHeight.totalHeight;
  const xAxisConfig = generateXAxisConfig(phenTrackMaxN, 4, xMin, xMax, newXLabel);
  const yAxisConfig = generateYAxisConfig(phenTrackMaxN, 4, yMax, strucYtext, bedYtext, phenGroupName, phenoYmax, showGene, showBed);
  const option = {
    grid: gridHeight.grid,
    tooltip: {
      trigger: 'axis',
      axisPointer: {
        axis: "x",
        type: "line",
        lineStyle: {
          color: "red"
        }
      }
    },
    axisPointer: {
      link: {
        xAxisIndex: 'all'
      },
      label: {
        show: true,
        backgroundColor: axisPointerColor,
        formatter: function (value) {
          return newXLabel[Math.round(value.value) - 1];
        }
      }
    },
    xAxis: xAxisConfig,
    yAxis: yAxisConfig,
    dataZoom: [{
        type: 'slider',
        xAxisIndex: [...Array(phenTrackMaxN + 4).keys()],
        top: 5,
        showDataShadow: false
      },
      {
        type: 'inside',
        xAxisIndex: [...Array(phenTrackMaxN + 4).keys()],
        filterMode: 'none'
      },
      {
        type: 'slider',
        yAxisIndex: [2],
        width: 20,
        start: 0,
        end: 100,
        handleSize: '80%',
        handleIcon: 'path://M30.9,53.2C16.8,53.2,5.3,41.7,5.3,27.6S16.8,2,30.9,2C45,2,56.4,13.5,56.4,27.6S45,53.2,30.9,53.2z M30.9,3.5C17.6,3.5,6.8,14.4,6.8,27.6c0,13.3,10.8,24.1,24.101,24.1C44.2,51.7,55,40.9,55,27.6C54.9,14.4,44.1,3.5,30.9,3.5z M36.9,35.8c0,0.601-0.4,1-0.9,1h-1.3c-0.5,0-0.9-0.399-0.9-1V19.5c0-0.6,0.4-1,0.9-1H36c0.5,0,0.9,0.4,0.9,1V35.8z M27.8,35.8 c0,0.601-0.4,1-0.9,1h-1.3c-0.5,0-0.9-0.399-0.9-1V19.5c0-0.6,0.4-1,0.9-1H27c0.5,0,0.9,0.4,0.9,1L27.8,35.8L27.8,35.8z',
        showDetail: false
      },
      {
        type: 'slider',
        yAxisIndex: [3],
        width: 20,
        start: 0,
        end: 100,
        handleSize: '80%',
        handleIcon: 'path://M30.9,53.2C16.8,53.2,5.3,41.7,5.3,27.6S16.8,2,30.9,2C45,2,56.4,13.5,56.4,27.6S45,53.2,30.9,53.2z M30.9,3.5C17.6,3.5,6.8,14.4,6.8,27.6c0,13.3,10.8,24.1,24.101,24.1C44.2,51.7,55,40.9,55,27.6C54.9,14.4,44.1,3.5,30.9,3.5z M36.9,35.8c0,0.601-0.4,1-0.9,1h-1.3c-0.5,0-0.9-0.399-0.9-1V19.5c0-0.6,0.4-1,0.9-1H36c0.5,0,0.9,0.4,0.9,1V35.8z M27.8,35.8 c0,0.601-0.4,1-0.9,1h-1.3c-0.5,0-0.9-0.399-0.9-1V19.5c0-0.6,0.4-1,0.9-1H27c0.5,0,0.9,0.4,0.9,1L27.8,35.8L27.8,35.8z',
        showDetail: false
      }
    ],
    toolbox: {
      show: true,
      feature: {
        saveAsImage: {
          show: true,
          title: "Download",
          name: "gpb_screenshot"
        }
      }
    },
    series: [
      ...seriesBlockList,
      ...seriesGraphNodeList,
      ...seriesGraphEdgeList,
      ...seriesLineList,
      seriesGeneStruc,
      seriesGeneStrand,
      seriesBed,
      ...seriesPhenTrackList
    ]
  }

  return {
    option,
    totalHeight
  };

}


function fillNodePanel(curNode, nodeSeq, nodeSample, refInfo, dupMap) {
	const fileName = refInfo.chr + '_' + refInfo.start + '-' + refInfo.end
	
  const seq = nodeSeq.get(curNode) || '';
  const dupNote = (dupMap && dupMap[curNode]) ? ` (duplicated from node ${dupMap[curNode]})` : '';
  const revNote = curNode.endsWith('-') ? ` (reverse complement of ${curNode.split(',').map(s => s.slice(0, -1) + '+').reverse().join(',')})` : '';
  document.getElementById("nodePanelSeq").innerHTML =
    `<br><div style="display:flex;justify-content:space-between;align-items:center;margin-bottom:4px;">
      <div style="white-space:pre-wrap;word-break:break-all;font-weight:bold;">Sequence of node ${curNode}${dupNote}${revNote}:</div>
      <button class="btn btn-sm btn-success" onclick="downloadFasta('${curNode}', \`${seq.replace(/\\/g, '\\\\').replace(/`/g, '\\`').replace(/\n/g, '\\n')}\`, '${fileName}')" style="font-size:12px;padding:2px 10px;">
        <span style="margin-right:4px;">⬇</span>Download FASTA
      </button>
    </div>`
    + '<textarea class="form-control" readonly style="width:100%;font-family:Courier New,monospace;font-size:13px;resize:none;margin-left:0%" rows="5">'
    + seq
    + "</textarea>";

  const curSampleAll = nodeSample.get(curNode).split(",");
  const countMap = {};
  curSampleAll.forEach(g => countMap[g] = (countMap[g] || 0) + 1);
  const uniqueGenomes = Object.keys(countMap).sort((a, b) => a.localeCompare(b));
  const textList = uniqueGenomes.map(g => countMap[g] > 1 ? `${g}(${countMap[g]})` : g);
  document.getElementById("nodePanelSample").innerHTML =
    `<div style="display:flex;justify-content:space-between;align-items:center;margin-bottom:4px;">
      <b>Detected ${curSampleAll.length} path(s) from ${uniqueGenomes.length} genome(s):</b>
      <button class="btn btn-sm btn-success" onclick="downloadNodeCsv('${curNode}', '${uniqueGenomes.join(',')}', '${uniqueGenomes.map(g => countMap[g]).join(',')}', '${fileName}')" style="font-size:12px;padding:2px 10px;">
        <span style="margin-right:4px;">⬇</span>Download CSV
      </button>
    </div>`
    + '<textarea class="form-control" readonly style="width:100%;font-family:Courier New,monospace;font-size:13px;resize:none;margin-left:0%;height:180px">'
    + textList.join(", ")
    + "</textarea>";
}


function downloadFasta(nodeId, sequence, fileName) {
  const fastaContent = `>${nodeId}\n${sequence.match(/.{1,60}/g).join('\n')}`;
  const blob = new Blob([fastaContent], {
    type: 'text/plain'
  });
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = `${fileName}_node${nodeId}.fasta`;
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
  URL.revokeObjectURL(url);
}


function downloadNodeCsv(nodeId, genomesStr, countsStr, fileName) {
  const genomes = genomesStr.split(',');
  const counts = countsStr.split(',').map(Number);
  let csvContent = 'Genome,Count\n';
  genomes.forEach((g, i) => {
    csvContent += `${g},${counts[i]}\n`;
  });
  const blob = new Blob([csvContent], {
    type: 'text/csv;charset=utf-8;'
  });
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = `${fileName}_node${nodeId}_genomes.csv`;
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
  URL.revokeObjectURL(url);
}

 const saveSvg=require("saveSvg")
var color = d3.scaleOrdinal(d3.schemeCategory10);
var radius = 5;

//var attractForce = d3.forceManyBody().strength(-10).distanceMax(400).distanceMin(60);
var repelForce = d3.forceManyBody().strength(20).distanceMax(500).distanceMin(15);

var simulation = d3.forceSimulation()
    .force("link", d3.forceLink().id(function(d) { return d.id; }).strength(0.1))
//   .force("charge", d3.forceManyBody().strength(-100).distanceMax(300).distanceMin(30))
//    .force("attractForce",attractForce)
    .force("repelForce",repelForce)
    .force("center", d3.forceCenter(width / 2, height / 2))
    .force("collide", d3.forceCollide().radius(radius*3+1.5).iterations(2));






r2d3.onRender(function(graph, svg, width, height, options) {
    var g = svg.append("g");
    
    var zoomFunc=d3.zoom()    
        .scaleExtent([1 / 2, 4])
        .on("zoom", zoomed);
      
    function zoomed() {
    g.attr("transform", d3.event.transform);
  }

    
    g.call(zoomFunc);

  var link = g.append("g")
      .attr("class", "links")
      .selectAll("line")
    .data(graph.links)
    .enter().append("line")
    .attr('marker-end','url(#end)')
      .attr("stroke-width", 2);
      

  var nodes = g.selectAll(".nodes")
            .data(graph.nodes).enter()
            .append("g")
        .attr("class", "nodes")
        .on("mouseover", moveoverNode)
        .on("mouseout", moveoutNode)
      .call(d3.drag()
          .on("start", dragstarted)
          .on("drag", dragged)
          .on("end", dragended));
          
   //get legend type and text       
var keys = d3.map(graph.nodes,function(d){return d.parent}).keys();
// Add one dot in the legend for each name.
var legend = g.selectAll("legend")
  .data(keys).enter()
  .append("circle")
    .attr("cx", 110)
    .attr("cy", function(d,i){ return 50 + i*25}) // 100 is where the first dot appears. 25 is the distance between dots
    .attr("r", 7)
    .style("fill", function(d){ return color(d)});

  
// Add one dot in the legend for each name.
var legendText = g.selectAll("mylabels")
  .data(keys).enter()
  .append("text")
    .attr("x", 20)
    .attr("y", function(d,i){ return 50 + i*25}) // 100 is where the first dot appears. 25 is the distance between dots
    .style("fill", function(d){ return color(d)})
    .text(function(d){ return d})
    .attr("text-anchor", "left")
    .style("alignment-baseline", "middle");
    

    
  
  // group each type of nodes
    var LigandNode=nodes.filter(function(d){return d.parent=="ligand";}).attr("class","ligand");
    var ReceptorNode=nodes.filter(function(d){return d.parent=="receptor";}).attr("class","receptor");
    var LigandCellNode=nodes.filter(function(d){return d.parent=="ligand cell";}).attr("class","ligand_cell");
    var ReceptorCellNode=nodes.filter(function(d){return d.parent=="receptor cell";}).attr("class","receptor_cell");


    LigandNode.append("circle")
    .attr("r", radius)
    .attr("data-legend",function(d) { return d.parent})
    .attr("fill", function(d) { return color(d.parent); });
    
    ReceptorNode.append("circle")
    .attr("r", radius)
    .attr("data-legend",function(d) { return d.parent})
    .attr("fill", function(d) { return color(d.parent); });
    
    LigandCellNode.append("circle")
    .attr("r", radius*3)
    .attr("data-legend",function(d) { return d.parent})
    .attr("fill", function(d) { return color(d.parent); });
    
    ReceptorCellNode.append("circle")
    .attr("r", radius*3)
    .attr("data-legend",function(d) { return d.parent})
    .attr("fill", function(d) { return color(d.parent); });
    

    LigandNode.append("text")
    .attr("dx", ".60em")
    .attr("dy", ".35em")
    .attr("font-size",10)
    .text(function(d) { return d.name; });
    
    ReceptorNode.append("text")
    .attr("dx", ".60em")
    .attr("dy", ".35em")
    .attr("font-size",10)
    .text(function(d) { return d.name; });
    
    LigandCellNode.append("text")
    .attr("dx", "-0.80em")
    .attr("dy", ".10em")
    .attr("font-size",15)
    .text(function(d) { return d.name; });
    
    ReceptorCellNode.append("text")
    .attr("dx", "-0.80em")
    .attr("dy", ".10em")
    .attr("font-size",15)
    .text(function(d) { return d.name; });
    
    
  simulation
    .nodes(graph.nodes)
    .on("tick", ticked);

    simulation.force("link")
      .links(graph.links);
 
  function ticked() {
    //var nodesValue=d3.map(graph.nodes,function(d){return path.centroid(d);});
     // var q = d3.geom.quadtree(nodesValue),
     // i = 0,
    //  n = nodesValue.length;
    //  while (++i < n) q.visit(collide(nodesValue[i]));
      
      
    link
        .attr("x1", function(d) { return d.source.x; })
        .attr("y1", function(d) { return d.source.y; })
        .attr("x2", function(d) { return d.target.x; })
        .attr("y2", function(d) { return d.target.y; });
        
            nodes
        //constrains the nodes to be within a box
    //.attr("cx", function(d) { return d.x = Math.max(radius, Math.min(width - radius, d.x)); })
  //.attr("cy", function(d) { return d.y = Math.max(radius, Math.min(height - radius, d.y)); })
      .attr("transform", function(d) { return "translate(" + d.x + "," + d.y + ")"; });

    LigandCellNode
    .attr("cx",function(d){return d.x=Math.max(radius,Math.min(150-radius,d.x));})
    .attr("cy",function(d){return d.y=Math.max(50+radius,Math.min(height-radius,d.y));})
    
      LigandNode
    .attr("cx",function(d){return d.x=Math.max(150+radius,Math.min(350-radius,d.x));})
    .attr("cy",function(d){return d.y=Math.max(50+radius,Math.min(height-radius,d.y));})
        ReceptorNode
    .attr("cx",function(d){return d.x=Math.max(450+radius,Math.min(650-radius,d.x));})
    .attr("cy",function(d){return d.y=Math.max(50+radius,Math.min(height-radius,d.y));})
        ReceptorCellNode
    .attr("cx",function(d){return d.x=Math.max(500+radius,Math.min(850-radius,d.x));})
    .attr("cy",function(d){return d.y=Math.max(50+radius,Math.min(height-radius,d.y));})

    
  }
  
  function dragstarted(d) {
  if (!d3.event.active) simulation.alphaTarget(0.3).restart();
  d.fx = d.x;
  d.fy = d.y;
}

function dragged(d) {
  d.fx = d3.event.x;
  d.fy = d3.event.y;
}

function dragended(d) {
  if (!d3.event.active) simulation.alphaTarget(0);
  d.fx = null;
  d.fy = null;
}
     
function getNeighbors(node) {
  return graph.links.reduce((neighbors, link) => {
    if (link.target.id ===node.id) {
      neighbors.push(link.source.id);
    } else if (link.source.id === node.id) {
      neighbors.push(link.target.id);
    }
    return neighbors;
  }, [node.id]);
}

function isNeighborLink(node, link) {
  return link.target.id === node.id || link.source.id === node.id;
}
function getOpacity(node,link){
  return isNeighborLink(node,link)? 0.9:0.2;
}

function moveoverNode(moveoveredNode){
       //   d3.select(this).select("circle").transition()
        //  .duration(500)
      //    .attr("r", 10);
      //    d3.select(this).select("text").transition()
      //    .duration(500)
      //    .attr("font-size",14 );
  const neighbors=getNeighbors(moveoveredNode);

  nodes.classed("inactiveNode",function(d){
    return neighbors.indexOf(d.id) ===-1;
  });
  link.attr("opacity",link=>getOpacity(moveoveredNode,link));
}

function moveoutNode(moveouttedNode){
     //     d3.select(this).select("circle").transition()
      //    .duration(500)
      //    .attr("r", radius);
    //      d3.select(this).select("text").transition()
    //      .duration(500)
    //      .attr("font-size",10 );
          nodes.classed("inactiveNode",false);
          link.attr("opacity",0.9);
}
          
  
Shiny.addCustomMessageHandler("refreshUpExp", refresh);
function refresh(message){
  link.remove();
  nodes.remove();
  LigandNode.remove();
  ReceptorNode.remove();
  LigandCellNode.remove();
  ReceptorCellNode.remove();
legend.remove();
legendText.remove();
simulation.alpha(1).restart();
}    

Shiny.addCustomMessageHandler("downloadUpStreamNetwork1", downloadPlot);
// Set-up the export button
function downloadPlot(dir){
      saveSvg.saveSvgAsPng(svg.node(),"upstreamNetworkPlot.png");
}
  
});


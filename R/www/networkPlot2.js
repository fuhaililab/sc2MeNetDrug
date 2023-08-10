const saveSvg=require("saveSvg")
var color = d3.scaleOrdinal(d3.schemeCategory20);
var radius =5;

//var attractForce = d3.forceManyBody().strength(-10).distanceMax(400).distanceMin(60);
var repelForce = d3.forceManyBody().strength(40).distanceMax(1000).distanceMin(500);

var simulation = d3.forceSimulation()
    .force("link", d3.forceLink().id(function(d) { return d.id; }).strength(0.1))
//   .force("charge", d3.forceManyBody().strength(-100).distanceMax(300).distanceMin(30))
//    .force("attractForce",attractForce)
      .force("charge", d3.forceManyBody())
    .force("center", d3.forceCenter(width / 2, height / 2))
    .force("collide", d3.forceCollide().radius(function(d){return d.size*radius*3}).iterations(2));






r2d3.onRender(function(graph, svg, width, height, options) {
        var g = svg.append("g");
    
    var zoomFunc=d3.zoom()    
        .scaleExtent([1 / 2, 4])
        .on("zoom", zoomed);
      
    function zoomed() {
    g.attr("transform", d3.event.transform);
  }

    
    g.call(zoomFunc);

  g.append('defs').append('marker')
                .attr('id','end')
                .attr('viewBox','-0 -5 10 10')
                .attr('refX',17)
                .attr('refY',0)
                .attr('orient','auto')
                .attr('markerWidth',5)
                .attr('markerHeight',5)
                .attr('xoverflow','visible')
                .append('svg:path')
                    .attr('d', 'M 0,-5 L 10 ,0 L 0,5')
                    .attr('fill', '#ccc')
                    .attr('stroke','#ccc');
                    
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
var legend =g.selectAll("legend")
  .data(keys).enter()
  .append("circle")
    .attr("cx", 160)
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
    var LigandNode=nodes.filter(function(d){return d.parent=="Ligand";}).attr("class","Ligand");
    var ReceptorNode=nodes.filter(function(d){return d.parent=="Receptor";}).attr("class","Receptor");
    var TFNode=nodes.filter(function(d){return d.parent=="TranscriptionFactor";}).attr("class","TF");
    var TargetsNode=nodes.filter(function(d){return d.parent=="Target";}).attr("class","Targets");
    var linkNode=nodes.filter(function(d){return d.parent=="LinkNode";}).attr("class","LinkNode");

    LigandNode.append("circle")
    .attr("r", function(d){return d.size*radius*0.8})
    .attr("data-legend",function(d) { return d.parent})
    .attr("fill", function(d) { return color(d.parent); });
    
    ReceptorNode.append("circle")
    .attr("r", function(d){return d.size*radius*0.8})
    .attr("data-legend",function(d) { return d.parent})
    .attr("fill", function(d) { return color(d.parent); });
    
    TFNode.append("circle")
    .attr("r", function(d){return d.size*radius*0.8})
    .attr("data-legend",function(d) { return d.parent})
    .attr("fill", function(d) { return color(d.parent); });
    
    TargetsNode.append("circle")
    .attr("r", function(d){return d.size*radius*0.8})
    .attr("data-legend",function(d) { return d.parent})
    .attr("fill", function(d) { return color(d.parent); });
    
    linkNode.append("circle")
    .attr("r", function(d){return d.size*radius*0.8})
    .attr("data-legend",function(d) { return d.parent})
    .attr("fill", function(d) { return color(d.parent); });
    
    

    nodes.append("text")
    .attr("dx", ".35em")
    .attr("dy", ".35em")
    .attr("font-size",10)
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
          d3.select(this).select("circle").transition()
          .duration(500)
          .attr("r", function(d){return (d.size*radius)});
          d3.select(this).select("text").transition()
          .duration(500)
          .attr("font-size",14 );
  const neighbors=getNeighbors(moveoveredNode);

  nodes.classed("inactiveNode",function(d){
    return neighbors.indexOf(d.id) ===-1;
  });
  link.attr("opacity",link=>getOpacity(moveoveredNode,link));
}

function moveoutNode(moveouttedNode){
          d3.select(this).select("circle").transition()
          .duration(500)
          .attr("r", function(d){return(d.size*radius*0.8)});
          d3.select(this).select("text").transition()
          .duration(500)
          .attr("font-size",10 );
          nodes.classed("inactiveNode",false);
          link.attr("opacity",0.9);
}
          
  
Shiny.addCustomMessageHandler("refreshNetwork2", refresh);
function refresh(message){
  link.remove();
  nodes.remove();
  LigandNode.remove();
legend.remove();
legendText.remove();
simulation.alpha(1).restart();
}

Shiny.addCustomMessageHandler("downloadNetwork2", downloadPlot);
// Set-up the export button
function downloadPlot(dir){
      saveSvg.saveSvgAsPng(svg.node(),"networkPlot.png");
}

});



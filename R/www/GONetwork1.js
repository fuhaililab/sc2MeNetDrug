//const saveSvg=require("saveSvg")
var color = d3.scaleOrdinal(d3.schemeCategory10);
var radius = 7;

//var attractForce = d3.forceManyBody().strength(-10).distanceMax(400).distanceMin(60);
var repelForce = d3.forceManyBody().strength(-500).distanceMax(70).distanceMin(15);

var simulation = d3.forceSimulation()
    .force("link", d3.forceLink().id(function(d) { return d.node; }).strength(0.1))
//   .force("charge", d3.forceManyBody().strength(-100).distanceMax(300).distanceMin(30))
//    .force("attractForce",attractForce)
    .force("repelForce",repelForce)
    .force("center", d3.forceCenter(width / 2, height / 2))
    .force("collide", d3.forceCollide().radius(radius*3+2).iterations(2));






r2d3.onRender(function(graph, svg, width, height, options) {

  var link = svg.append("g")
      .attr("class", "links")
      .selectAll("line")
    .data(graph.links)
    .enter().append("line")
    .attr('marker-end','url(#end)')
      .attr("stroke-width", 2);
      

  var nodes = svg.selectAll(".nodes")
            .data(graph.nodes).enter()
            .append("g")
        .attr("class", "nodes")
        .on("mouseover", moveoverNode)
        .on("mouseout", moveoutNode)
      .call(d3.drag()
          .on("start", dragstarted)
          .on("drag", dragged)
          .on("end", dragended));
      

    
  
  // group each type of nodes
    var GONode=nodes.filter(function(d){return d.type=="GO";}).attr("class","GO");
    var geneNode=nodes.filter(function(d){return d.type=="gene";}).attr("class","gene");

    GONode.append("circle")
    .attr("r", radius*2)
    .attr("data-legend",function(d) { return d.type})
    .attr("fill", function(d) { return color(d.type); });
    
    geneNode.append("circle")
    .attr("r", radius)
    .attr("data-legend",function(d) { return d.type})
    .attr("fill", function(d) { return color(d.type); });
    

    GONode.append("text")
    .attr("dx", ".60em")
    .attr("dy", ".35em")
    .attr("font-size",15)
    .text(function(d) { return d.name; });
    
    geneNode.append("text")
    .attr("dx", "-0.80em")
    .attr("dy", ".10em")
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
    .attr("cx", function(d) { return d.x = Math.max(radius, Math.min(width - radius, d.x)); })
   .attr("cy", function(d) { return d.y = Math.max(radius, Math.min(height - radius, d.y)); })
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
    if (link.target.node ===node.node) {
      neighbors.push(link.source.node);
    } else if (link.source.node === node.node) {
      neighbors.push(link.target.node);
    }
    return neighbors;
  }, [node.node]);
}

function isNeighborLink(node, link) {
  return link.target.node=== node.node || link.source.node === node.node;
}
function getOpacity(node,link){
  return isNeighborLink(node,link)? 0.9:0.2;
}

function moveoverNode(moveoveredNode){
  const neighbors=getNeighbors(moveoveredNode);

  nodes.classed("inactiveNode",function(d){
    return neighbors.indexOf(d.node) ===-1;
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
          
  
Shiny.addCustomMessageHandler("refreshGONetwork1", refresh);
function refresh(message){
  link.remove();
  nodes.remove();
  GONode.remove();
  geneNode.remove();
simulation.alpha(1).restart();
}    
});
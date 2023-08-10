var color = d3.scaleOrdinal()
    .domain(d3.range(4))
    .range(["#d6b0fe","#d6d8fe"]);
var radius = 5;

//var attractForce = d3.forceManyBody().strength(-10).distanceMax(400).distanceMin(60);
var repelForce = d3.forceManyBody().strength(-500).distanceMax(200).distanceMin(100);

var simulation = d3.forceSimulation()
    .force("link", d3.forceLink().id(function(d) { return d.node; }).strength(0.1))
//   .force("charge", d3.forceManyBody().strength(-100).distanceMax(300).distanceMin(30))
//    .force("attractForce",attractForce)
    .force("repelForce",repelForce)
    .force("center", d3.forceCenter(width / 2, height / 2))
    .force("collide", d3.forceCollide().radius(radius+2).iterations(2));


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
      .attr("stroke-width", 2);
      

  var nodes = g.selectAll(".nodes")
            .data(graph.nodes).enter()
            .append("g")
        .attr("class", "nodes")
      .call(d3.drag()
          .on("start", dragstarted)
          .on("drag", dragged)
          .on("end", dragended));
          
   //get legend type and text       

  
  
  // group each type of nodes
    var normalNode=nodes.filter(function(d){return d.type=="normal";}).attr("class","normal");
    var exemNode=nodes.filter(function(d){return d.type=="exemplar";}).attr("class","exemplar");

    normalNode.append("circle")
    .attr("r", radius)
    .attr("data-legend",function(d) { return d.type})
    .attr("fill", function(d) { return color(d.type); });
    
    exemNode.append("circle")
    .attr("r", radius*2)
    .attr("data-legend",function(d) { return d.type})
    .attr("fill", function(d) { return color(d.type); });
    


    nodes.append("text")
    .attr("dx", ".40em")
    .attr("dy", ".40em")
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
      .attr("transform", function(d) { return "translate(" + d.x + "," + d.y + ")"; })
      .attr("cx",function(d){return d.x=Math.max(radius,Math.min(width-radius,d.x));})
    .attr("cy",function(d){return d.y=Math.max(radius,Math.min(height-radius,d.y));});
      
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
     

          
  
Shiny.addCustomMessageHandler("refreshCluster2", refresh);
function refresh(message){
  link.remove();
  nodes.remove();
  normalNode.remove();
  exemNode.remove();

simulation.alpha(1).restart()
}          
  
  
});


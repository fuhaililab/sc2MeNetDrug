
// Initialization
svg.attr("font-family", "sans-serif")
  .attr("font-size", "8")
  .attr("text-anchor", "middle")
  .attr("class","graph");
    
    
var format = d3.format(",d");
var color = d3.scaleOrdinal(d3.schemeCategory20c);

var group = svg.append("g");

let centerX = width * 0.5;
let centerY = height * 0.5;
let strength = 0.03;


// use the force
let simulation = d3.forceSimulation()
	// .force('link', d3.forceLink().id(d => d.id))
	.force('charge', d3.forceManyBody())
	.force('x', d3.forceX(centerX ).strength(strength))
	.force('y', d3.forceY(centerY ).strength(strength))
  .force("collide",d3.forceCollide().radius(d => Math.min(80,Math.max(40,d.length*2))).iterations(2));

// Rendering
r2d3.onRender(function(data, svg, width, height, options) {
    //svg.on("click", function(e){if(e.target==this){hideMarkerGene}})
    svg.on("click", hideMarkerGene)
    var node = group.selectAll(".node")
    .data(data)
    .enter().append("g")
      .attr("class", "node")
      .attr("transform", function(d) { return "translate(" + d.x + "," + d.y + ")"; })
      .on("click",showMarkerGene);

  node.append("circle")
      .attr("id", function(d) { return d.cellType; })
      .attr("r", function(d) { return Math.min(80,Math.max(40,d.length*2)); })
      .attr("fill", function(d) { return color(d.cellType); })
      .call(d3.drag()
      .on("start", dragstarted)
      .on("drag", dragged)
      .on("end", dragended))
      .transition().duration(2000).ease(d3.easeElasticOut)
				.tween('circleIn', (d) => {
					let i = d3.interpolateNumber(0, d.radius);
					return (t) => {
						d.r = i(t);
						simulation.force('collide', forceCollide);
					}
				})
				;

    var text=node.append("text")
    .attr("font-size",10)
    .attr("cy",-2)
    
    text.selectAll("tspan")
    .data(d=>d.cellType.split("="))
    .enter()
    .append("tspan")
    .text(function(d) { return d ; })
    .attr("dy","1em")
    .attr("x","0");
    

      

  simulation.nodes(data).on('tick', ticked);
  
  		function ticked() {
			node
				.attr("transform", function(d) { return "translate(" + d.x + "," + d.y + ")"; })
        .attr("cy",function(d){return d.y=Math.max(30,Math.min(height-30,d.y));});
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

function showMarkerGene(cell){

  
    var g=svg.selectAll(".geneContainer").data(cell.cellType).enter().append("g");
    g.attr("id","geneContainer")
    g.classed(".container")
    var container=g.append("rect")
            .attr('width', width)
            .attr('height', height/3)
            .attr("y",height-height/3)
            .attr("x",0)
            .attr("id","containerRect")
            .attr("fill","white");
           // .attr("style","stroke:gray;stroke-width:2px 0px 2px 0px;transition: opacity 250ms;-webkit-transition: opacity 250ms;-moz-transition: opacity 250ms;box-shadow: 2px 0 0 0;");
        g.append("line")
          .attr("x1",0)
          .attr("y1",height-height/3)
          .attr("x2",width)
          .attr("y2",height-height/3)
          .attr("style","stroke:gray;stroke-width:2px;stroke-opacity:0.5;box-shadow:2px")

      g.append("text")
       .attr("font-size",24)
       .attr("x",125)
       .attr("y",height-height/3+30)
       .text("Biomarker gene list:");
       
    var genes=g.append("text")
      .attr("font-size",16)
    .attr("y",height-height/3+30+20);
      
    genes.selectAll("tspan")
    .data(d=>cell.markerGene.split(","))
    .enter()
    .append("tspan")
    .text(function(d) { return d ; })
    .attr("x",35)
    .attr("dy","1em");
    
                    
}

function hideMarkerGene(d){
        var event = d3.event,
          target =  event.srcElement, 
          data = d3.select(target).datum(); 
      if(data == undefined) {
          var c=svg.selectAll("#geneContainer");
          c.remove();
          //event.stopPropagation();
      } 
      
}

Shiny.addCustomMessageHandler("refreshBubble", refresh);
function refresh(message){
  node.remove();
simulation.alpha(1).restart();
}

});
/**
 * @module
 * @description  Print coloured circles into the cell output
 * @file  d3-simple-example.js
 */

 // Inspired by: https://www.d3-graph-gallery.com/intro_d3js.html

// create SVG element in the output area
// the ``element`` is a contextual binding to the output of the current cell
let svg = d3.select(element.get(0)).append('svg');

// create group
let g = svg.append('g');

g.append("circle")
  .attr("cx", 2).attr("cy", 2).attr("r", 40).style("fill", "blue");
g.append("circle")
  .attr("cx", 140).attr("cy", 70).attr("r", 40).style("fill", "red");
g.append("circle")
  .attr("cx", 300).attr("cy", 100).attr("r", 40).style("fill", "green");
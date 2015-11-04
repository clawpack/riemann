# Interactive apps to show solutions of Riemann problems.

# Left and right states in the phase plane can be dragged
# and dropped, as well as the time bar on the x-t plane.
# MPLD3: D3-javascript plugin for matplotlib http://mpld3.github.io/
# Questions contact Mauricio Del Razo S. @maojrs: maojrs@gmail.com
"""
In order to run the following code in python, try running the tests at 
the end of the code. If you are using the jupyter notebooks, remember 
switching mpld3.show() to mpld3.display()
"""


import numpy as np
import jinja2
import json
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.path as mpath
import matplotlib.patches as mpatches

import mpld3
from mpld3 import plugins, utils

# Shallow water interactive phase plane plugin in mpld3
class PPlaneNLPlugin(plugins.PluginBase):
    JAVASCRIPT = r"""
    // init custom PPLaneNL plugin
    mpld3.register_plugin("drag", PPlaneNLPlugin);
    PPlaneNLPlugin.prototype = Object.create(mpld3.Plugin.prototype);
    PPlaneNLPlugin.prototype.constructor = PPlaneNLPlugin;
    PPlaneNLPlugin.prototype.requiredProps = ["id", "idmpoint", 
                                              "idg", "iditer", "iditer_charac", "idtime", "idhmax", "idoffm",
                                              "idhugol", "idhugor", "idintcl", "idintcr",
                                              "idhugol2", "idhugor2", "idintcl2", "idintcr2", 
                                              "idqlm", "idqmm", "idqrm", 
                                              "idshock1", "idshock2", "idrar1", "idrar2",
                                              "idrar1a", "idrar2a", "idrar1b", "idrar2b",
                                              "idrar1c", "idrar2c", "idrar1d", "idrar2d",
                                              "idtimedot", "idtimeline",
                                              "idq1", "idq2"];
    PPlaneNLPlugin.prototype.defaultProps = {}
    function PPlaneNLPlugin(fig, props){
        mpld3.Plugin.call(this, fig, props);
        mpld3.insert_css("#" + fig.figid + " path.dragging",
                         {"fill-opacity": "1.0 !important",
                          "stroke-opacity": "1.0 !important"});
    };

    // Call draw function, this function is being looped all the time
    PPlaneNLPlugin.prototype.draw = function(){
        // Get elements into script variables
        var obj = mpld3.get_element(this.props.id);
        var midpoint = mpld3.get_element(this.props.idmpoint);
        var g = this.props.idg;
        var iter = this.props.iditer;
        var iter_charac = this.props.iditer_charac;
        var time = this.props.idtime;
        var hmax = this.props.idhmax;
        var offm = this.props.idoffm;
        var hugol = mpld3.get_element(this.props.idhugol);
        var hugor = mpld3.get_element(this.props.idhugor);
        var intcl = mpld3.get_element(this.props.idintcl);
        var intcr = mpld3.get_element(this.props.idintcr);
        var hugol2 = mpld3.get_element(this.props.idhugol2);
        var hugor2 = mpld3.get_element(this.props.idhugor2);
        var intcl2 = mpld3.get_element(this.props.idintcl2);
        var intcr2 = mpld3.get_element(this.props.idintcr2);
        var qlm = mpld3.get_element(this.props.idqlm);
        var qmm = mpld3.get_element(this.props.idqmm);
        var qrm = mpld3.get_element(this.props.idqrm);
        var shock1 = mpld3.get_element(this.props.idshock1);
        var shock2 = mpld3.get_element(this.props.idshock2);
        var rar1 = mpld3.get_element(this.props.idrar1);
        var rar2 = mpld3.get_element(this.props.idrar2);
        var rar1a = mpld3.get_element(this.props.idrar1a);
        var rar2a = mpld3.get_element(this.props.idrar2a);
        var rar1b = mpld3.get_element(this.props.idrar1b);
        var rar2b = mpld3.get_element(this.props.idrar2b);
        var rar1c = mpld3.get_element(this.props.idrar1c);
        var rar2c = mpld3.get_element(this.props.idrar2c);
        var rar1d = mpld3.get_element(this.props.idrar1d);
        var rar2d = mpld3.get_element(this.props.idrar2d);
        var timedot = mpld3.get_element(this.props.idtimedot);
        var timeline = mpld3.get_element(this.props.idtimeline);
        var q1 = mpld3.get_element(this.props.idq1);
        var q2 = mpld3.get_element(this.props.idq2);
        
        // Set initial conditions for javascript calculations
        //var qleft = obj.offsets[0];
        //var qright = obj.offsets[1];
        //var qmid = midpoint.offsets[0];
        var offx = obj.ax.x(offm);
        var offy = offx;
        // Define initial values of hl, hr, hul and hur
        var hl = obj.offsets[0][0];
        var hul = obj.offsets[0][1];
        var hr = obj.offsets[1][0];
        var hur = obj.offsets[1][1];
        hmfinal = 0;
        humfinal = 0;
        
        // Main d3 drag function
        var drag = d3.behavior.drag()
            .origin(function(d) { return {x:obj.ax.x(d[0]),
                                          y:obj.ax.y(d[1])}; }) 
            .on("dragstart", dragstarted)
            .on("drag", dragged)
            .on("dragend", dragended);
        
        // Dragtime function for draggable time-dot (analog to main)
        var dragtime = d3.behavior.drag()
            .origin(function(d) { return {x:timedot.ax.x(d[0]),
                                          y:timedot.ax.y(d[1])}; })                              
            .on("dragstart", dragstarted)
            .on("drag", dragged)
            .on("dragend", dragended);
    
        // Set elements of ql and qr draggable points and call main drag function
        obj.elements()
           .data(obj.offsets)
           .style("cursor", "default")
           .call(drag);
            
        // Set elements for timedot draggable point and call dragtime function  
        timedot.elements()
           .data(timedot.offsets)
           .style("cursor", "default")
           .call(dragtime);

        // Begin phi and phi prime function
        function phi(h,hl,hr,hul,hur) {
            var ul = hul/hl;
            var ur = hur/hr;
            if (h < hl) {
                var termIC1 = 2.0*Math.sqrt(g*hl);
                var termIC2 = 2.0*Math.sqrt(g*h);
                var phil = ul + termIC1 - termIC2;
                var philp = -Math.sqrt(g/h);
            }
            if (h >= hl) {
                var termHL1 = (h - hl);
                var termHL2 = Math.sqrt(0.5*g*(1.0/h + 1.0/hl));
                var termHL3 = termHL1*termHL2;
                var phil = ul - termHL3 ;
                var nume = termHL2*4*h*h;
                var philp = -termHL2 + termHL1*g/nume;
            }
            if (h < hr) {
                var termIC1 = 2.0*Math.sqrt(g*hr);
                var termIC2 = 2.0*Math.sqrt(g*h);
                var phir = ur - termIC1 + termIC2;
                var phirp = Math.sqrt(g/h);
            }
            if (h >= hr) {
                var termHL1 = (h - hr);
                var termHL2 = Math.sqrt(0.5*g*(1.0/h + 1.0/hr));
                var termHL3 = termHL1*termHL2;
                var phir = ur + termHL3 ;
                var nume = termHL2*4*h*h;
                var phirp = termHL2 - termHL1*g/nume;
            }
            var phi = phil-phir;
            var phiprime = philp - phirp
            var um = 0.5*(phil + phir);
            return [phi,phiprime,um];
        }
        
        // Begin Newton iteration
        function newton(hstar,hl,hr,hul,hur) {
            var hn = hstar;
            var error = 1;
            while (error > 0.005) {
                var hold = hn;
                var phin = phi(hn,hl,hr,hul,hur);
                var hn = hold - phin[0]/phin[1];
                var um = phin[2];
                error = Math.abs(hold - hn);
            }
            return [hn,um];
        }
        
        // Calculate hugoniot loci
        function hugoloci(hm,hside,huside,sign) {
                var termHL1 = hm*(hm - hside);
                if (hm == 0) {
                    var termHL2 = 1.0; }
                else {
                   var termHL2 = Math.sqrt(0.5*g*(1.0/hm + 1.0/hside));}
                   var termHL3 = termHL1*termHL2;
                   var hloci = hm*huside/hside + sign*termHL3;
            return hloci;} 
            
        // Calculate integral curve
        function integralcurve(hm,hside,huside,sign) {
                var termIC1 = 2.0*hm*Math.sqrt(g*hside);
                var termIC2 = 2.0*hm*Math.sqrt(g*hm);
                var intcurve = hm*huside/hside - sign*(termIC1 - termIC2);
            return intcurve;}
            
        // Calculate solution plot of h and hu as function of x
        function solplot(hl,hm,hr,hul,hum,hur) {
            var lam = [];
            lam.push([-1000000, 0,0])
            var ul = hul/hl;
            var um = hum/hm;
            var ur = hur/hr;
            var shock1 = false;
            var shock2 = false;
            var qsol1 = d3.range(iter);
            var qsol2 = d3.range(iter);
            if (hm >= hl) {
                var s1 = time*(hum - hul)/(hm - hl);
                shock1 = true;
            } else {
                var s1b = [0,0];
                s1b[0] = time*(ul - Math.sqrt(g*hl));
                s1b[1] = time*(um - Math.sqrt(g*hm));
                s1b.sort(d3.ascending);
            }
            if (hm >= hr) {
                var s2 = time*(hur - hum)/(hr - hm);
                shock2 = true;
            } else {
                var s2b = [0,0];
                s2b[0] = time*(um + Math.sqrt(g*hm));
                s2b[1] = time*(ur + Math.sqrt(g*hr));
                s2b.sort(d3.ascending);
            }
            
            for (var ii=0; ii <2*iter; ii++){
                var xx = q1.data[ii][0];
                // Calculate plot solution for 1 charactersitic (shock or rarefaction)
                if (shock1) {
                    if (xx <= s1){ qsol1[ii] = hl; qsol2[ii] = hul;}
                    else { qsol1[ii] = hm; qsol2[ii] = hum;}
                }
                else {
                    if (xx <= s1b[0]) { qsol1[ii] = hl; qsol2[ii] = hul;}
                    else if (xx <= s1b[1]) {
                        var m1 = (hm - hl)/(s1b[1] - s1b[0]);
                        var m2 = (hum - hul)/(s1b[1] - s1b[0]);
                        qsol1[ii] = m1*(xx - s1b[0]) + hl;
                        qsol2[ii] = m2*(xx - s1b[0]) + hul;}
                    else { qsol1[ii] = hm; qsol2[ii] = hum;}
                }
                // Calculate plot solution for 2 charactersitic (shock or rarefaction)
                if (shock2) {
                    if (xx > s2) {qsol1[ii] = hr; qsol2[ii] = hur;}
                }
                else{
                    if (xx > s2b[0] && xx < s2b[1]) {
                        var m1 = (hr - hm)/(s2b[1] - s2b[0]);
                        var m2 = (hur - hum)/(s2b[1] - s2b[0]);
                        qsol1[ii] = m1*(xx - s2b[0]) + hm;
                        qsol2[ii] = m2*(xx - s2b[0]) + hum;}
                    else if (xx > s2b[1]) { qsol1[ii] = hr; qsol2[ii] = hur;}
                }
            }
            var solution = [qsol1,qsol2];
            return solution;
        }

        // Function to update middle state   
        function update_midstate(){
            var hstar = 0.05*(hl + hr);
            var solution = newton(hstar,hl,hr,hul,hur);
            hmfinal = solution[0];
            humfinal = solution[1]*hmfinal;
            var xx = obj.ax.x(hmfinal);
            var yy = obj.ax.y(humfinal);
            // Update middle state point and marker position
            midpoint.elements().transition().duration(5)
            .attr("transform", "translate(" + [xx, yy] + ")");
            // Move marker
            qmm.elements().transition().duration(1)
              .attr("transform", "translate(" + [xx + 0.7*offx, yy + 0.7*offy] + ")");
        }
        
        // Functon to update xt-plane
        function update_xtplane() {
         // Calculate shcok speeds from R-H conditions
         var lam1 = (humfinal - hul)/(hmfinal - hl);   
         var lam2 = (hur - humfinal)/(hr - hmfinal);   
         var lam1m = humfinal/hmfinal - Math.sqrt(g*hmfinal);
         var lam2m = humfinal/hmfinal + Math.sqrt(g*hmfinal);
         var lam1l = hul/hl - Math.sqrt(g*hl);
         var lam2r = hur/hr + Math.sqrt(g*hr);
         var color1 = "red";
         var color2 = "red";
         var thick1 = 4;
         var thick2 = 4;
         var fan = 0;
         for (var ii=0; ii<iter_charac; ii++) {
             if (hmfinal >= hl) {
                 shock1.data[ii][1] = shock1.data[ii][0]/lam1;
                 rar1.data[ii][1] = rar1.data[ii][0]/lam1;
                 rar1a.data[ii][1] = rar1a.data[ii][0]/lam1;
                 rar1b.data[ii][1] = rar1b.data[ii][0]/lam1;
                 rar1c.data[ii][1] = rar1c.data[ii][0]/lam1;
                 rar1d.data[ii][1] = rar1d.data[ii][0]/lam1;
                 console.log(rar2a.data);
                 color1 = "red";
                 thick1 = 4;
             } else {
                 shock1.data[ii][1] = shock1.data[ii][0]/lam1l;
                 rar1.data[ii][1] = rar1.data[ii][0]/lam1m;
                 fan = Math.abs(lam1m - lam1l)/5;
                 rar1a.data[ii][1] = rar1a.data[ii][0]/(lam1l + fan);
                 rar1b.data[ii][1] = rar1b.data[ii][0]/(lam1l + 2*fan);
                 rar1c.data[ii][1] = rar1c.data[ii][0]/(lam1l + 3*fan);
                 rar1d.data[ii][1] = rar1d.data[ii][0]/(lam1l + 4*fan);
                 color1 = "blue";
                 thick1 = 1;
             }
             if (hmfinal >= hr) {
                 shock2.data[ii][2] = shock2.data[ii][0]/lam2;
                 rar2.data[ii][2] = rar2.data[ii][0]/lam2;
                 rar2a.data[ii][2] = rar2a.data[ii][0]/lam2;
                 rar2b.data[ii][2] = rar2b.data[ii][0]/lam2;
                 rar2c.data[ii][2] = rar2c.data[ii][0]/lam2;
                 rar2d.data[ii][2] = rar2d.data[ii][0]/lam2;
                 color2 = "red";
                 thick2 = 4;
            } else {
                 shock2.data[ii][2] = shock2.data[ii][0]/lam2r;
                 rar2.data[ii][2] = rar2.data[ii][0]/lam2m;
                 fan = Math.abs(lam2r - lam2m)/5;
                 rar2a.data[ii][2] = rar2a.data[ii][0]/(lam2m + fan);
                 rar2b.data[ii][2] = rar2b.data[ii][0]/(lam2m + 2*fan);
                 rar2c.data[ii][2] = rar2c.data[ii][0]/(lam2m + 3*fan);
                 rar2d.data[ii][2] = rar2d.data[ii][0]/(lam2m + 4*fan);
                 color2 = "blue";
                 thick2 = 1;
            }                     
         }
         // Do transitions
         shock1.elements().transition().duration(5)
             .attr("d", shock1.datafunc(shock1.data))
             .style("stroke", color1)
             .style("stroke-width", thick1);
         shock2.elements().transition().duration(5)
             .attr("d", shock2.datafunc(shock2.data))
             .style("stroke", color2)
             .style("stroke-width", thick2);
         rar1.elements().transition().duration(5)
             .attr("d", rar1.datafunc(rar1.data))
             .style("stroke", color1)
             .style("stroke-width", thick1);
         rar2.elements().transition().duration(5)
             .attr("d", rar2.datafunc(rar2.data))
             .style("stroke", color2)
             .style("stroke-width", thick2);
         rar1a.elements().transition().duration(5)
             .attr("d", rar1a.datafunc(rar1a.data))
             .style("stroke", color1)
             .style("stroke-width", thick1);
         rar2a.elements().transition().duration(5)
             .attr("d", rar2a.datafunc(rar2a.data))
             .style("stroke", color2)
             .style("stroke-width", thick2);
         rar1b.elements().transition().duration(5)
             .attr("d", rar1b.datafunc(rar1b.data))
             .style("stroke", color1)
             .style("stroke-width", thick1);
         rar2b.elements().transition().duration(5)
             .attr("d", rar2b.datafunc(rar2b.data))
             .style("stroke", color2)
             .style("stroke-width", thick2);
         rar1c.elements().transition().duration(5)
             .attr("d", rar1c.datafunc(rar1c.data))
             .style("stroke", color1)
             .style("stroke-width", thick1);
         rar2c.elements().transition().duration(5)
             .attr("d", rar2c.datafunc(rar2c.data))
             .style("stroke", color2)
             .style("stroke-width", thick2);
         rar1d.elements().transition().duration(5)
             .attr("d", rar1d.datafunc(rar1d.data))
             .style("stroke", color1)
             .style("stroke-width", thick1);
         rar2d.elements().transition().duration(5)
             .attr("d", rar2d.datafunc(rar2d.data))
             .style("stroke", color2)
             .style("stroke-width", thick2);
        }
        
        // Function to update solution plots
        function update_solplots() {
            var qsol = solplot(hl,hmfinal,hr,hul,humfinal,hur);
                for (var ii=0; ii<2*iter; ii++){
                    q1.data[ii][1] = qsol[0][ii]; 
                    q2.data[ii][2] = qsol[1][ii]; 
                }
                // Do transitions
                 q1.elements().transition().duration(5)
                     .attr("d", q1.datafunc(q1.data));
                 q2.elements().transition().duration(5)
                     .attr("d", q2.datafunc(q2.data));
        }
        
        // Initialize solution with given initial states before interacting
        update_midstate();
        update_xtplane();
        update_solplots();
        
        // Begin drag function
        function dragstarted(d) {
          d3.event.sourceEvent.stopPropagation();
          d3.select(this).classed("dragging", true);
        }
   
        // The drag function called while dragging is happening (meat of code here)
        function dragged(d,i) {
          if (i == 0  || i ==1) {
              // Convert mouse coordinates in drag event (d3.event) to python coordinates d
              d[0] = obj.ax.x.invert(d3.event.x);
              d[1] = obj.ax.y.invert(d3.event.y);             
              // Move ql and qr stored in obj (they have been selected in drag)
              d3.select(this)
                .attr("transform", "translate(" + [d3.event.x,d3.event.y] + ")");
              // If obj corresponds to ql, move all the other left elements 
              }
          if (i==0){
              // Move marker
              qlm.elements().transition().duration(1)
                  .attr("transform", "translate(" + [d3.event.x + offx, d3.event.y + offy] + ")");     
              // Re-calculate inital left variables when dragging 
              hl = obj.offsets[0][0];
              hul = obj.offsets[0][1];
              
              // Draw Hugoniot loci through left state    
              for (var ii=0; ii<iter; ii++) {
                //Left from hl
                hugol.data[ii][0] = ii*d[0]/(1.0*iter);
                var hm = hugol.data[ii][0];
                hugol.data[ii][1] = hugoloci(hm,hl,hul,-1) ;
                //Right from hl
                hugol2.data[ii][0] = d[0] + ii*(hmax-d[0])/(1.0*iter);
                var hm2 = hugol2.data[ii][0];
                hugol2.data[ii][1] = hugoloci(hm2,hl,hul,-1) ;
                }
              // Do transitions      
              hugol.elements().transition().duration(5)
                  .attr("d", hugol.datafunc(hugol.data))
                  .style("stroke-dasharray", ("7,7"));
              hugol2.elements().transition().duration(5)
                  .attr("d", hugol2.datafunc(hugol2.data));
                               
                  
              // Draw integral curve through left state (note index bug for intcl.data[ii][j])
              // and intcl2 it should be j =0,1 not j=1,2. Arrays somehow grow of sixe in mpld3
              for (var ii=0; ii<iter; ii++) {
                //Left from hl
                intcl.data[ii][1] = ii*d[0]/(1.0*iter);
                var hm = intcl.data[ii][1];
                intcl.data[ii][2] = integralcurve(hm,hl,hul,-1);
                //Right from hl
                intcl2.data[ii][1] = d[0] + ii*(hmax-d[0])/(1.0*iter);
                var hm2 = intcl2.data[ii][1];
                intcl2.data[ii][2] = integralcurve(hm2,hl,hul,-1);
              }
              // Do transitions      
              intcl.elements().transition().duration(5)
                  .attr("d", intcl.datafunc(intcl.data));
              intcl2.elements().transition().duration(5)
                  .attr("d", intcl2.datafunc(intcl2.data))
                  .style("stroke-dasharray", ("7,7"));
                  
          }      
          // if element corresponds to qr
          else if (i==1) {
              // Move marker
              qrm.elements().transition().duration(1)
                  .attr("transform", "translate(" + [d3.event.x + offx, d3.event.y + offy] + ")");
              // Re-calculate inital right variables when dragging
              hr = obj.offsets[1][0];
              hur = obj.offsets[1][1];
              
              // Draw Hugoniot loci through right state    
              for (var ii=0; ii<iter; ii++) {
                //Left from hr
                hugor.data[ii][0] = ii*d[0]/(1.0*iter);
                var hm = hugor.data[ii][0];
                hugor.data[ii][1] = hugoloci(hm,hr,hur,1) ;
                //Right from hr
                hugor2.data[ii][0] = d[0] + ii*(hmax-d[0])/(1.0*iter);
                var hm2 = hugor2.data[ii][0];
                hugor2.data[ii][1] = hugoloci(hm2,hr,hur,1) ;
                }
              // Do transitions      
              hugor.elements().transition().duration(5)
                  .attr("d", hugor.datafunc(hugor.data))
                  .style("stroke-dasharray", ("7,7"));
              hugor2.elements().transition().duration(5)
                  .attr("d", hugor2.datafunc(hugor2.data));
                    
              // Draw integral curve through right state (note index bug for intcr.data[ii][j])
              // and intcr2 it should be j =0,1 not j=1,2. Arrays somehow grow of sixe in mpld3   
              for (var ii=0; ii<iter; ii++) {
                //Left from hr
                intcr.data[ii][1] = ii*d[0]/(1.0*iter);
                var hm = intcr.data[ii][1];
                intcr.data[ii][2] = integralcurve(hm,hr,hur,1);
                //Right from hr
                intcr2.data[ii][1] = d[0] + ii*(hmax-d[0])/(1.0*iter);
                var hm2 = intcr2.data[ii][1];
                intcr2.data[ii][2] = integralcurve(hm2,hr,hur,1);
              }
              // Do transitions      
              intcr.elements().transition().duration(5)
                  .attr("d", intcl.datafunc(intcr.data));
              intcr2.elements().transition().duration(5)
                  .attr("d", intcr2.datafunc(intcr2.data))
                  .style("stroke-dasharray", ("7,7"));
            }
            // If time marker is moved
            else if (i==2) {
                // Convert mouse coordinates in drag event (d3.event) to python coordinates d
                d[0] = timedot.ax.x.invert(d3.event.x);
                d[1] = timedot.ax.y.invert(d3.event.y);      
                d3.select(this)
                    .attr("transform", "translate(" + [d3.event.x,d3.event.y] + ")");
                // Calculate timedot position and assign it
                var ty = timedot.ax.y.invert(d3.event.y);
                if (ty >=0) {
                    timeline.data[0][1] = ty;
                    timeline.data[1][1] = ty;
                    timeline.data[2][1] = ty;
                    time = ty;
                }
                // Do transitions
                timeline.elements().transition().duration(5)
                    .attr("d", timeline.datafunc(timeline.data));
            }
            
            // Calculate middle state    
            update_midstate();
                  
            // Calculate solution plots of h and hu
            update_solplots();
            
            // Update characteristic in x-t plane
            update_xtplane();
             
        }
        // End dragging
        function dragended(d) {
          d3.select(this).classed("dragging", false);
        }
    }
    mpld3.register_plugin("drag", PPlaneNLPlugin); 
    """

    def __init__(self, points, midpoint, 
                 g, iters, iter_charac, time, hmax, offm,
                 hugol, hugor, intcl, intcr,
                 hugol2, hugor2, intcl2, intcr2, 
                 qlmarker, qmmarker ,qrmarker, 
                 shock1, shock2, rar1, rar2,
                 rar1a,rar2a,rar1b,rar2b,
                 rar1c,rar2c,rar1d,rar2d,
                 timedot, timeline,
                 q1, q2):
        if isinstance(points, mpl.lines.Line2D):
            suffix = "pts"
        else:
            suffix = None

        self.dict_ = {"type": "drag",
                      "id": utils.get_id(points, suffix),
                      "idmpoint": utils.get_id(midpoint,suffix),
                      "idg" : g,
                      "iditer" : iters,
                      "iditer_charac" : iter_charac,
                      "idtime" : time,
                      "idhmax" : hmax,
                      "idoffm": offm,
                      "idhugol" : utils.get_id(hugol),
                      "idhugor" : utils.get_id(hugor),
                      "idintcl" : utils.get_id(intcl),
                      "idintcr" : utils.get_id(intcr),
                      "idhugol2" : utils.get_id(hugol2),
                      "idhugor2" : utils.get_id(hugor2),
                      "idintcl2" : utils.get_id(intcl2),
                      "idintcr2" : utils.get_id(intcr2),
                      "idqlm": utils.get_id(qlmarker,suffix),
                      "idqmm": utils.get_id(qmmarker,suffix),
                      "idqrm": utils.get_id(qrmarker,suffix),
                      "idshock1" : utils.get_id(shock1),
                      "idshock2" : utils.get_id(shock2),
                      "idrar1" : utils.get_id(rar1),
                      "idrar2" : utils.get_id(rar2),
                      "idrar1a" : utils.get_id(rar1a),
                      "idrar2a" : utils.get_id(rar2a),
                      "idrar1b" : utils.get_id(rar1b),
                      "idrar2b" : utils.get_id(rar2b),
                      "idrar1c" : utils.get_id(rar1c),
                      "idrar2c" : utils.get_id(rar2c),
                      "idrar1d" : utils.get_id(rar1d),
                      "idrar2d" : utils.get_id(rar2d),
                      "idtimedot" : utils.get_id(timedot,suffix),
                      "idtimeline" : utils.get_id(timeline),
                      "idq1" : utils.get_id(q1),
                      "idq2" : utils.get_id(q2)
                      }

# Python part for shallow water interactive part
# Define hugoloci and integral curve functions
def hugoloci(g,hm,hside,huside,sign):
    term1 = hm*(hm - hside)
    hm[hm==0] = hm[hm==0] + 0.00000001
    term2 = np.sqrt(0.5*g*(1.0/hm + 1.0/hside))
    hloci = hm*huside/hside + sign*term1*term2
    return hloci

def intcurve(g,hm,hside,huside,sign):
    term1 = 2.0*hm*np.sqrt(g*hside);
    term2 = 2.0*hm*np.sqrt(g*hm);
    intcurve = hm*huside/hside - sign*(term1 - term2);
    return intcurve

# Define interactive plot routine in python (calls mpld3 plugin)
# Loaded with default values, can be called without any argument
def shallow_water(ql=np.array([3.0, 5.0]), qr=np.array([3.0, -5.0]), g=1.0, time=2.0, tmax=5.0, 
                  hmax=None, humin=None, humax=None):
    # Create figure
    # Create a figure
    #fig, axfull = plt.subplots(2,2, figsize=(13, 12))
    fig, axfull = plt.subplots(2,2, figsize=(10, 8))
    fig.subplots_adjust(left=0.05, right=0.9, bottom=0.1, top=0.9,
                        hspace=0.3, wspace=0.15)

    # Have to do this to avoid issue with mpld3
    ax = axfull[0] # First row of plots
    axsol = axfull[1] # Second row of plost 

    # Calculate dq with ql and qr and other parameters
    dq = np.array([qr[0]-ql[0], qr[1]-ql[1]])
    iters = 100
    iter_charac = 2
    # eps required for bug in mpld3 (ql[0] cannot be the same than qr[0])
    eps = 0.00000001
    
    # Calculate plotting boundaries if not specified
    if hmax is None:
        hmax = max(ql[0],qr[0]) + 7.0
    if humax is None:
        humax = 3*max(abs(ql[1]),abs(qr[1]))
    if humin is None:    
        humin = -3*max(abs(ql[1]),abs(qr[1]))

    # PLOT PHASE PLANE
    xxL = np.linspace(0,ql[0],iters)
    xxL2 = np.linspace(ql[0],hmax,iters)
    xxR = np.linspace(0,qr[0]+eps,iters)
    xxR2 = np.linspace(qr[0]+eps,hmax,iters)

    #Plot midpoint
    qm = -1 +0.0*ql
    midpoint = ax[0].plot(qm[0],qm[1],'ok', alpha=0.9, markersize=8, markeredgewidth=1)

    # Plot hugoloci initial state
    yyL = hugoloci(g,xxL,ql[0],ql[1],-1)
    yyL2 = hugoloci(g,xxL2,ql[0],ql[1],-1)
    yyR = hugoloci(g,xxR,qr[0],qr[1],1)
    yyR2 = hugoloci(g,xxR2,qr[0],qr[1],1)
    hugol = ax[0].plot(xxL,yyL, '--r', linewidth=1.5)
    hugol2 = ax[0].plot(xxL2,yyL2, '-r', linewidth=2, label = 'Hugoniot Loci')
    hugor = ax[0].plot(xxR,yyR, '--r', linewidth=1.5, label = 'Hugoniot Loci (unphysical)')
    hugor2 = ax[0].plot(xxR2,yyR2, '-r', linewidth=2)

    # Plot integral curve initial state
    yyL = intcurve(g,xxL,ql[0],ql[1],-1)
    yyL2 = intcurve(g,xxL2,ql[0],ql[1],-1)
    yyR = intcurve(g,xxR,qr[0],qr[1],1)
    yyR2 = intcurve(g,xxR2,qr[0],qr[1],1)
    intcl = ax[0].plot(xxL,yyL, '-b', linewidth=2, label = 'Integral Curves')
    intcl2 = ax[0].plot(xxL2,yyL2, '--b', linewidth=1.5, label = 'Integral Curves (unphysical)')
    intcr = ax[0].plot(xxR,yyR, '-b', linewidth=2)
    intcr2 = ax[0].plot(xxR2,yyR2, '--b', linewidth=1.5)

    # Plot ql and qr
    points = ax[0].plot([ql[0],qr[0]], [ql[1], qr[1]], 'ok', alpha=0.7, markersize=10, markeredgewidth=1)
    #data = ["q_l", "q_r"] 

    # Plot markers
    offsetx = 0.3*hmax/10
    offsety = -3*offsetx
    qlmarker = ax[0].plot(ql[0]+offsetx, ql[1]+offsety, 'ok', marker=(r"$ q_l $"), markersize=15)
    qmmarker = ax[0].plot(qm[0]+offsetx, qm[1]+offsety, 'ok', marker=(r"$ q_m $"), markersize=15)
    qrmarker = ax[0].plot(qr[0]+offsetx, qr[1]+offsety, 'ok', marker=(r"$ q_r $"), markersize=15)

    # Set axis 1 properties
    ax[0].set_title("Phase plane", fontsize=18)
    ax[0].axis([0,hmax,humin,humax])
    ax[0].set_xlabel('h', fontsize=17)
    ax[0].set_ylabel('hu', fontsize=17)
    #ax[0].set_aspect('equal')
    ax[0].grid(alpha=0.1,color='k', linestyle='--')
    legend = ax[0].legend(loc='upper left', shadow=True, fontsize = 8)


    # PLOT x-t PLANE
    x_xtp = np.linspace(-10,10,iter_charac)
    x_xtp2 = np.linspace(-11,11,iter_charac)
    # Shock speeds lam1 and lam2
    lam2 = ql[1]/ql[0] - np.sqrt(g*ql[0])
    lam1 = qr[1]/qr[0] + np.sqrt(g*qr[0])
    char1 = x_xtp/lam1
    char2 = x_xtp/lam2
    char3 = x_xtp2/lam1
    char4 = x_xtp2/lam2
    shock1 = ax[1].plot(x_xtp, char1, '-k', linewidth=4, label="1 or 2 shock")
    shock2 = ax[1].plot(x_xtp, char2, '-k', linewidth=4)
    rar1 = ax[1].plot(x_xtp2, char3, '-k', linewidth=1, label="1 or 2 rarefaction")
    rar2 = ax[1].plot(x_xtp2, char4, '-k', linewidth=1)
    
    # For other 4 rarefaction lines on each side
    x_xtp3 = np.linspace(-11.1,11.1,iter_charac)
    x_xtp4 = np.linspace(-11.2,11.2,iter_charac)
    x_xtp5 = np.linspace(-11.3,11.3,iter_charac)
    x_xtp6 = np.linspace(-11.4,11.4,iter_charac)
    char1a = x_xtp3/lam1; char2a = x_xtp3/lam2 
    char1b = x_xtp4/lam1; char2b = x_xtp4/lam2 
    char1c = x_xtp5/lam1; char2c = x_xtp5/lam2 
    char1d = x_xtp6/lam1; char2d = x_xtp6/lam2 
    rar1a = ax[1].plot(x_xtp3, char1a, '-k', linewidth=1)
    rar2a = ax[1].plot(x_xtp3, char2a, '-k', linewidth=1)
    rar1b = ax[1].plot(x_xtp4, char1b, '-k', linewidth=1)
    rar2b = ax[1].plot(x_xtp4, char2b, '-k', linewidth=1)
    rar1c = ax[1].plot(x_xtp5, char1c, '-k', linewidth=1)
    rar2c = ax[1].plot(x_xtp5, char2c, '-k', linewidth=1)
    rar1d = ax[1].plot(x_xtp6, char1d, '-k', linewidth=1)
    rar2d = ax[1].plot(x_xtp6, char2d, '-k', linewidth=1)
    
    timedot = ax[1].plot([100000,1000000,9], [-10,-10,time], 'ok' , alpha=0.7, markersize=10)
    timeline = ax[1].plot([-12,0,12], [time, time, time], '--k', linewidth = 3, label="time")

    # Set axis 2 properties
    ax[1].set_title("x-t plane", fontsize=18)
    ax[1].set_xlabel('x', fontsize=17)
    ax[1].set_ylabel('t', fontsize=17)
    ax[1].axis([-10,10,-0,tmax])
    ax[1].grid(alpha=0.1,color='k', linestyle='--')
    legend = ax[1].legend(loc='upper center', shadow=True, fontsize = 8)

    # PLOT SOLUTIONS
    xsol = np.linspace(-10,10,2*iters)
    hsol = 0*xsol + ql[0]
    husol = 0*xsol + ql[1]
    q1 = axsol[0].plot(xsol,hsol, '-k', linewidth = 4, alpha = 1.0)
    q2 = axsol[1].plot(xsol,husol, '-k', linewidth = 4, alpha = 1.0)

    def solplot(xsol,ql,qr,qm,g):
        hl = ql[0];    hm = qm[0];    hr = qr[0]
        ul = ql[1]/hl; um = qm[1]/hm; ur = qr[1]/hr     
        lam = np.empty(4, dtype=float)

    # Set axis 3 properties
    axsol[0].set_title("Depth h at time = t", fontsize=18)
    axsol[0].set_xlabel('x', fontsize=17)
    axsol[0].set_ylabel('h', fontsize=17)
    axsol[0].axis([-10,10,0,hmax])
    axsol[0].grid(alpha=0.1,color='k', linestyle='--')

    # Set axis 4 properties
    axsol[1].set_title("Momentum hu at time = t", fontsize=18)
    axsol[1].set_xlabel('x', fontsize=17)
    axsol[1].set_ylabel('hu', fontsize=17)
    axsol[1].axis([-10,10,humin,humax])
    axsol[1].grid(alpha=0.1,color='k', linestyle='--')


    # Remove defult mpld3 plugins
    plugins.clear(fig)

    # Call mpld3 custom PPLane plugin to interact with plot
    plugins.connect(fig, PPlaneNLPlugin(points[0],midpoint[0], 
                                        g,iters,iter_charac,time, hmax, offsetx,
                                        hugol[0],hugor[0],intcl[0],intcr[0],
                                        hugol2[0],hugor2[0],intcl2[0],intcr2[0],
                                        qlmarker[0],qmmarker[0],qrmarker[0],
                                        shock1[0],shock2[0],rar1[0],rar2[0],
                                        rar1a[0],rar2a[0],rar1b[0],rar2b[0],
                                        rar1c[0],rar2c[0],rar1d[0],rar2d[0],
                                        timedot[0],timeline[0],
                                        q1[0],q2[0]))
    return fig


##########################################

# Plugin for interactive phase plane plot for 2D linear case 
class PPlanePlugin(plugins.PluginBase):
    JAVASCRIPT = r"""
    // init custom PPLane plugin
    mpld3.register_plugin("drag", PPlanePlugin);
    PPlanePlugin.prototype = Object.create(mpld3.Plugin.prototype);
    PPlanePlugin.prototype.constructor = PPlanePlugin;
    PPlanePlugin.prototype.requiredProps = ["id", "idmpoint", 
                                            "idlinesla", "idlineslb", "idlinesra", "idlinesrb",  
                                            "idqlm", "idqmm", "idqrm", "idqone", "idqtwo", 
                                            "idrl0", "idrl1", "idrr0", "idrr1"];
    PPlanePlugin.prototype.defaultProps = {}
    function PPlanePlugin(fig, props){
        mpld3.Plugin.call(this, fig, props);
        mpld3.insert_css("#" + fig.figid + " path.dragging",
                         {"fill-opacity": "1.0 !important",
                          "stroke-opacity": "1.0 !important"});
    };

    // Call draw function, this function is being looped all the time
    PPlanePlugin.prototype.draw = function(){
        // Get elements into script variables
        var obj = mpld3.get_element(this.props.id);
        var midpoint = mpld3.get_element(this.props.idmpoint);
        var linesla = mpld3.get_element(this.props.idlinesla);
        var lineslb = mpld3.get_element(this.props.idlineslb);
        var linesra = mpld3.get_element(this.props.idlinesra);
        var linesrb = mpld3.get_element(this.props.idlinesrb);
        var qlm = mpld3.get_element(this.props.idqlm);
        var qmm = mpld3.get_element(this.props.idqmm);
        var qrm = mpld3.get_element(this.props.idqrm);
        var qone = mpld3.get_element(this.props.idqone);
        var qtwo = mpld3.get_element(this.props.idqtwo);
        var rl0 = this.props.idrl0;
        var rl1 = this.props.idrl1;
        var rr0 = this.props.idrr0;
        var rr1 = this.props.idrr1;
        
        // Set initial conditions for javascript calculations
        var qleft = obj.offsets[0];
        var qright = obj.offsets[1];
        var qmid = midpoint.offsets[0];
        var off = 13;
        // Calculate slopes for eigenlines
        var ml = rl1/rl0;
        var mr = rr1/rr0;
        
        // Main d3 drag function
        var drag = d3.behavior.drag()
            .origin(function(d) { return {x:obj.ax.x(d[0]),
                                          y:obj.ax.y(d[1])}; })                              
            .on("dragstart", dragstarted)
            .on("drag", dragged)
            .on("dragend", dragended);
    
        // Set elements of ql and qr points and call main drag function
        obj.elements()
           .data(obj.offsets)
           .style("cursor", "default")
           .call(drag);      

        // Begin drag function
        function dragstarted(d) {
          d3.event.sourceEvent.stopPropagation();
          d3.select(this).classed("dragging", true);
        }
    
        // The drag function called while dragging is happening (meat of code here)
        function dragged(d,i) {
          // Convert mouse coordinates in drag event (d3.event) to python coordinates d
          d[0] = obj.ax.x.invert(d3.event.x);
          d[1] = obj.ax.y.invert(d3.event.y);             
          // Move ql and qr stored in obj (they have been selected in drag)
          d3.select(this)
            .attr("transform", "translate(" + [d3.event.x, d3.event.y] + ")");
          // If obj corresponds to ql, move all the other left elements 
          if (i==0){
              // Move text marker
              qlm.elements().transition().duration(1)
                  .attr("transform", "translate(" + [d3.event.x + off, d3.event.y + off] + ")");
                  
              // Move eigenlines
              for (var ii=0; ii<2; ii++) {
                linesla.data[ii][1] = ml*(linesla.data[ii][0] - d[0]) + d[1];
                lineslb.data[ii][1] = mr*(lineslb.data[ii][0] - d[0]) + d[1]; }
                
              // In script calculations of middle state
              qleft = [d[0], d[1]];
              var det = rl0*rr1 - rr0*rl1;
              var alphal = (rr1*(qright[0] - qleft[0]) - rr0*(qright[1] - qleft[1]))/det
              qmid[0] = qleft[0] + alphal*rl0;
              qmid[1] = qleft[1] + alphal*rl1;
              var xx = obj.ax.x(qmid[0]);
              var yy = obj.ax.y(qmid[1]); 
             
                }
          // if element corresponds to qr
          else {
              // Move text marker
              qrm.elements().transition().duration(1)
                  .attr("transform", "translate(" + [d3.event.x + off, d3.event.y + off] + ")");
                  
              // Move eigenlines
              for (var ii=0; ii<2; ii++) {
                linesra.data[ii][1] = ml*(linesra.data[ii][0] - d[0]) + d[1];
                linesrb.data[ii][1] = mr*(linesrb.data[ii][0] - d[0]) + d[1]; }
                
              // In script calculations of middle state    
              qright = [d[0], d[1]];
              var det = rl0*rr1 - rr0*rl1;
              var alphal = (rr1*(qright[0] - qleft[0]) - rr0*(qright[1] - qleft[1]))/det
              qmid[0] = qleft[0] + alphal*rl0;
              qmid[1] = qleft[1] + alphal*rl1;
              var xx = obj.ax.x(qmid[0]);
              var yy = obj.ax.y(qmid[1]);
                }
              // Update middle state point and marker position
              midpoint.elements().transition().duration(5)
                .attr("transform", "translate(" + [xx, yy] + ")");
              qmm.elements().transition().duration(5)
                .attr("transform", "translate(" + [xx + 0.7*off, yy + 0.7*off] + ")"); 
                
              // Update eigenlines
              linesla.elements().transition().duration(5)
                  .attr("d", linesla.datafunc(linesla.data));
              lineslb.elements().transition().duration(5)
                  .attr("d", lineslb.datafunc(lineslb.data));
              linesra.elements().transition().duration(5)
                  .attr("d", linesra.datafunc(linesra.data));
              linesrb.elements().transition().duration(5)
                  .attr("d", linesrb.datafunc(linesrb.data));
                
              // Update subplots of q1 and q2
              qone.data[0][1] = qleft[0];
              qone.data[1][1] = qleft[0];
              qone.data[2][1] = qmid[0];
              qone.data[3][1] = qmid[0];
              qone.data[4][1] = qright[0];
              qone.data[5][1] = qright[0];
              qtwo.data[0][2] = qleft[1];
              qtwo.data[1][2] = qleft[1];
              qtwo.data[2][2] = qmid[1];
              qtwo.data[3][2] = qmid[1];
              qtwo.data[4][2] = qright[1];
              qtwo.data[5][2] = qright[1];
              qone.elements().transition().duration(5)
                .attr("d", qone.datafunc(qone.data));
              qtwo.elements().transition().duration(5)
                .attr("d", qtwo.datafunc(qtwo.data));  
        }
        // End dragging
        function dragended(d) {
          d3.select(this).classed("dragging", false);
        }
    }
    mpld3.register_plugin("drag", PPlanePlugin); 
    """

    def __init__(self, points, midpoint, 
                 linesla, lineslb, linesra, linesrb, 
                 qlmarker, qmmarker, qrmarker, 
                 qone, qtwo, rl0, rl1, rr0, rr1):
        if isinstance(points, mpl.lines.Line2D):
            suffix = "pts"
        else:
            suffix = None

        self.dict_ = {"type": "drag",
                      "id": utils.get_id(points, suffix),
                      "idmpoint": utils.get_id(midpoint,suffix),
                      "idlinesla": utils.get_id(linesla),
                      "idlineslb": utils.get_id(lineslb),
                      "idlinesra": utils.get_id(linesra),
                      "idlinesrb": utils.get_id(linesrb),
                      "idqlm": utils.get_id(qlmarker,suffix),
                      "idqmm": utils.get_id(qmmarker,suffix),
                      "idqrm": utils.get_id(qrmarker,suffix),
                      "idqone": utils.get_id(qone),
                      "idqtwo": utils.get_id(qtwo),
                      "idrl0": rl0,
                      "idrl1": rl1,
                      "idrr0": rr0,
                      "idrr1": rr1}


def linear_phase_plane(ql=np.array([-2.0, 2.0]),qr=np.array([0.0, 3.0]), 
                       r1=None, r2=None, lam1=None, lam2=None, 
                       q1min=-5, q1max=5, q2min =-5, q2max=5, domain=5, time=1.0, 
                       title1=None, title2=None):
    # Create figure
    # Create a figure
    fig, ax = plt.subplots(1,3, figsize=(13.5, 5.0))
    fig.subplots_adjust(left=0.05, right=0.9, bottom=0.2, top=0.9,
                        hspace=0.3, wspace=0.15)

    #Create subfigure ax1 for phaseplane
    #ax[0] = fig.add_subplot(2,2,1)

    # Calculate dq from ql and qr
    dq = np.array([qr[0]-ql[0], qr[1]-ql[1]])

    # If eigenvectors not specified, default to Acoustics
    rho0 = 1.0
    zz = 2.0
    if r1 is None:
        r1 = np.array([-zz, 1])
    if r2 is None:
        r2 = np.array([zz,  1])
    if lam1 is None:
        lam1 = -zz/rho0
    if lam2 is None:
        lam2 = zz/rho0
    
    # Choose left and right eigenvectors
    if (lam1 < lam2):
        laml = lam1; lamr = lam2
        rl = r1; rr = r2 
    else:
        laml = lam2; lamr = lam1
        rl = r2; rr = r1      
        
    # Assign titles to solution plots
    if title1 is None:
        title1 = "q1"
    if title2 is None:
        title2 = "q2"
    
    #Plot eigenlines across ql and qr
    eps = 0.00000000001 # To avoid bugs in D3
    ml = rl[1]/rl[0]
    mr = rr[1]/rr[0] + eps
    linesla = ax[0].plot([q1min,q1max],[ml*(q1min - ql[0]) + ql[1],ml*(q1max - ql[0]) + ql[1]], '-k')
    lineslb = ax[0].plot([q1min-eps,q1max+eps],[mr*(q1min - ql[0]) + ql[1],mr*(q1max - ql[0]) + ql[1]], '-k')
    linesra = ax[0].plot([q1min-2*eps,q1max+2*eps],[mr*(q1min - qr[0]) + qr[1],mr*(q1max - qr[0]) + qr[1]], '-k')
    linesrb = ax[0].plot([q1min-3*eps,q1max+3*eps],[ml*(q1min - qr[0]) + qr[1],ml*(q1max - qr[0]) + qr[1]], '-k')

    # Plot ql and qr
    points = ax[0].plot([ql[0],qr[0]], [ql[1], qr[1]], 'ok', alpha=0.7, markersize=10, markeredgewidth=1)
    data = ["q_l", "q_r"]
    offset = 0.4*0.5*(q1max-q1min)/5.0
    qlmarker = ax[0].plot(ql[0] + offset, ql[1] - offset, 'ok', marker=(r"$ q_l $"), markersize=15)
    qrmarker = ax[0].plot(qr[0] + offset, qr[1] - offset, 'ok', marker=(r"$ q_r $"), markersize=15)

    #Plot midpoint
    det = rl[0]*rr[1] - rr[0]*rl[1]
    alL = (rr[1]*dq[0] - rr[0]*dq[1])/det
    qm = ql + alL*rl 
    midpoint = ax[0].plot(qm[0],qm[1],'ok', alpha=0.9, markersize=8, markeredgewidth=1)
    qmmarker = ax[0].plot(qm[0]+offset,qm[1]-0.7*offset, 'k',marker=(r"$ q_m $"),markersize=20)

    # Set axis 1 properties
    ax[0].set_title("Phase Plane", fontsize=18)
    ax[0].axis([q1min,q1max,q2min,q2max])
    ax[0].grid(alpha=0.1,color='k', linestyle='--')
    ax[0].set_xlabel(title1)
    ax[0].set_ylabel(title2)

    # Remove defult mpld3 plugins
    plugins.clear(fig)

    # Create solutionl line with six points where solution should be
    ctleft = laml*time
    ctright = lamr*time
    domain = max(domain,abs(ctleft)+1,abs(ctright)+1) # Readjust xdomain if necessarry
    xsol = np.array([-domain,ctleft,ctleft,ctright,ctright,domain])
    qsol1 = 1*xsol
    qsol1[0:2] = ql[0]
    qsol1[2:4] = qm[0]
    qsol1[4:6] = qr[0]

    # Set axis 2 properties
    ax[1].set_title(title1, fontsize=18)
    ax[1].axis([-domain,domain,q1min,q1max])
    ax[1].grid(alpha=0.1,color='k', linestyle='--')
    ax[1].set_xlabel('x')
    #ax[1].set_ylabel(title1)
    
    # Plot solution
    qone = ax[1].plot(xsol, qsol1, '-k', linewidth = 4, alpha = 1.0)

    #Create subfigure ax2 for solution plane
    #ax[2] = fig.add_subplot(2,2,2)
    # Create solutionl line with six points
    xsol2 = np.array([-domain,ctleft,ctleft,ctright,ctright,domain])
    qsol2 = 1*xsol2
    qsol2[0:2] = ql[1]
    qsol2[2:4] = qm[1]
    qsol2[4:6] = qr[1]

    # Set axis 2 properties
    ax[2].set_title(title2, fontsize=18)
    ax[2].axis([-domain,domain,q2min,q2max])
    ax[2].grid(alpha=0.1,color='k', linestyle='--')
    ax[2].set_xlabel('x')
    #ax[2].set_ylabel(title2)

    # Plot solution
    qtwo = ax[2].plot(xsol2, qsol2, '-k', linewidth = 4, alpha = 1.0)

    # Call mpld3 custom PPLane plugin to interact with plot
    plugins.connect(fig, PPlanePlugin(points[0],midpoint[0],
                                      linesla[0],lineslb[0],linesra[0],linesrb[0],
                                      qlmarker[0],qmmarker[0],qrmarker[0],qone[0],qtwo[0],
                                      rl[0],rl[1],rr[0],rr[1]))
    
    return fig


### Tests for interactive apps
def test_interactive_shallowPP(save_to_html=False):
    # Define left and right state (h,hu)
    ql = np.array([3.0, 5.0])
    qr = np.array([3.0, -5.0])
    # Define optional parameters (otherwise chooses default values)
    plotopts = {'g':1.0, 'time':2.0, 'tmax':5, 'hmax':10, 'humin':-15, 'humax':15}
    # Call interactive function (can be called without any argument)
    pt = shallow_water(ql,qr,**plotopts)
    if save_to_html:
        mpld3.save_html(pt, "test_shallow.html")
    mpld3.show()

def test_interactive_linearPP(save_to_html=False):
    ## Define left and right state 
    ql = np.array([-2.0, 2.0]) 
    qr = np.array([0.0, -3.0])
    # Define two eigenvectors and eigenvalues (acoustics)
    zz = 2.0
    rho0 = 1.0
    r1 = np.array([zz,1.0])
    r2 = np.array([-zz,1.0])
    lam1 = zz/rho0
    lam2 = -zz/rho0
    plotopts={'q1min':-5, 'q1max':5, 'q2min':-5, 'q2max':5, 'domain':5, 'time':1, 
            'title1':"Pressure", 'title2':"Velocity"}
    pt = linear_phase_plane(ql,qr,r1,r2,lam1,lam2,**plotopts)
    if save_to_html:
        mpld3.save_html(pt, "test_linearPP.html")
    mpld3.show()

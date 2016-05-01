(function() {
    'use strict';

    angular
        .module('oncoscape')
        .directive('osSunburst', sunburst);
//

    /** @ngInject */
    function sunburst() {

        var directive = {
            restrict: 'E',
            templateUrl: 'app/components/sunburst/sunburst.html',
            controller: SunburstController,
            controllerAs: 'vm',
            bindToController: true
        };

        return directive;

        /** @ngInject */
        function SunburstController(osApi, $state, $timeout, $scope, $stateParams, d3, $window) {

            if (angular.isUndefined($stateParams.datasource)){
                $state.go("datasource");
                return;
            }

            // Elements
            var d3Chart = d3.select("#sunburst-chart").append("svg").attr("id", "chart");
            var ptRow = 1;
            var svg;
            var g, path, text;
            var arc;
            var partition = d3.layout.partition()
                .sort(null)
                .value(function(d) { return d.size; });

           // Properties
            var cohortPatient = osApi.getCohortPatient();
            var rawPatientData, tbl;
            var width, height, radius;
            var x, y;
            var color = d3.scale.category20c();
            var root = {"name": "live", "size": 1, "children": [
                        {"name": "CD45", "size": 1, "children": [
                            {"name": "B cell", "size": 1},
                            {"name": "NK", "size": 1},
                            {"name": "PMN-MDSC","size": 1},
                            {"name": "Macrophage", "size": 1},
                            {"name": "Monocyte", "size": 1},
                            {"name": "M-MDSC", "size": 1},
                            {"name": "CD14-CD33-", "size": 1},
                            {"name": "CD3","size": 1, "children": [
                    			    {"name": "Tgd", "size": 1, "children": [
                                {"name": "Tgd_INFG", "size": 1},
                        				{"name": "Tgd_IL22", "size": 1},
                        				{"name": "Tgd_IL17", "size": 1}				]},
                      			{"name": "NKT", "size": 1},
                      			{"name": "CD8", "size": 1, "children": [
                      				{"name": "CD8_INFG", "size": 1}			]},
                      			{"name": "CD4", "size": 1, "children": [
                              {"name": "Th22", "size": 1},
                              {"name": "Th17", "size": 1},
                              {"name": "Th1", "size": 1},
                              {"name": "Treg", "size": 1}  ]}  ]} ]},
                        {"name": "epCAM","size": 1}   ]}

            // View Model
            var vm = this;
            vm.datasource = $stateParams.datasource;
            vm.optCohortPatients = cohortPatient.get();
            vm.optCohortPatient = vm.optCohortPatients[0];

            // Event Handlers
            vm.resize = draw;


            osApi.setBusy(true)("Loading Dataset");
            osApi.setDataset(vm.datasource).then(function(response) {

                // Percent Population Data
              osApi.getFlowData().then(function(response) {

                    tbl = response.payload;
                    osApi.setBusy(false);

                    setScale();
                    draw();

                    // set SVG size and location
                    function setScale(){
                      width = $window.innerWidth - 100;
                      height = $window.innerHeight - 190;
                      if (angular.element(".tray").attr("locked")=="true") width -= 300;

                       radius = Math.min(width, height) / 2;
                       x = d3.scale.linear().range([0, 2 * Math.PI]);
                       y = d3.scale.sqrt().range([0, radius]);

                      d3Chart.attr("width", width)
                             .attr("height", height)

                      svg = d3Chart.append("g")
                              .attr("transform", "translate(" + width / 2 + "," + (height / 2 + 10) + ")");

                      arc = d3.svg.arc()
                              .startAngle(function(d)  { return Math.max(0, Math.min(2 * Math.PI, x(d.x))); })
                              .endAngle(  function(d)  { return Math.max(0, Math.min(2 * Math.PI, x(d.x + d.dx))); })
                              .innerRadius(function(d) { return Math.max(0, y(d.y)); })
                              .outerRadius(function(d) { return Math.max(0, y(d.y + d.dy)); });

                    }

                }); //flow
            }); //dataset

            d3.select(self.frameElement).style("height", height + "px");

            // update zoom level of sunburst to chosen node
            function click(d) {

                text.transition().attr("opacity", 0);

                path.transition()
                  .duration(1000)
                  .attrTween("d", arcTweenZoom(d))
                  .each("end", function(e, i) {
                      // check if the animated element's data e lies within the visible angle span given in d
                      if (e.x >= d.x && e.x < (d.x + d.dx)) {
                      // get a selection of the associated text element
                      var arcText = d3.select(this.parentNode).select("text");
                      // fade in the text element and recalculate positions
                      arcText.transition().duration(750)
                        .attr("opacity", 1)
                        .attr("transform", function() { return "rotate(" + computeTextRotation(e) + ")" })
                        .attr("x", function(d) { return y(d.y); });
                      }
                  });
            } //click

            // add/update data to path and text
            function draw(){

                g = svg.selectAll("g")
                  .data(partition.nodes(root))
                  .enter().append("g");

                path = g.append("path")
                  .attr("d", arc)
                  .style("fill", function(d) { return color((d.children ? d : d.parent).name); })
                  .on("click", click);

                text = g.append("text")
                  .attr("transform", function(d) { return "rotate(" + computeTextRotation(d) + ")"; })
                  .attr("x", function(d) { return y(d.y); })
                  .attr("dx", "6") // margin
                  .attr("dy", ".35em") // vertical-align
                  .text(function(d) { return d.name; });

              d3.selectAll("input").on("change", function change() {
                var value = this.value === "size"
                  ? function() { return 1; }
                  : function(d){ return getRelativeSize(d) };

                path
                  .data(partition.value(value).nodes)
                  .transition()
                  .duration(1000)
                  .attrTween("d", arcTweenData);
              }); // change input path

            }

            //----------------------------------
            function getRelativeSize(d){
                if(d.depth == 0){ return +tbl[ptRow][d.name];}
                else return setRelSize(d.parent) * +tbl[ptRow][d.name] /100;
            }


            // Setup for switching data: stash the old values for transition.
            function stash(d) {
              d.x0 = d.x;
              d.dx0 = d.dx;
            }

            // When switching data: interpolate the arcs in data space.
            function arcTweenData(a, i) {
              var oi = d3.interpolate({x: a.x0, dx: a.dx0}, a);
              function tween(t) {
                var b = oi(t);
                a.x0 = b.x;
                a.dx0 = b.dx;
                return arc(b);
              }
              if (i == 0) {
               // If we are on the first arc, adjust the x domain to match the root node
               // at the current zoom level. (We only need to do this once.)
                var xd = d3.interpolate(x.domain(), [node.x, node.x + node.dx]);
                return function(t) {
                  x.domain(xd(t));
                  return tween(t);
                };
              } else {
                return tween;
              }
            }

            // When zooming: interpolate the scales.
            function arcTweenZoom(d) {
              var xd = d3.interpolate(x.domain(), [d.x, d.x + d.dx]),
                  yd = d3.interpolate(y.domain(), [d.y, 1]),
                  yr = d3.interpolate(y.range(), [d.y ? 20 : 0, radius]);
              return function(d, i) {
                return i
                    ? function(t) { return arc(d); }
                    : function(t) { x.domain(xd(t)); y.domain(yd(t)).range(yr(t)); return arc(d); };
              };
            }

            function computeTextRotation(d) {
              return (x(d.x + d.dx / 2) - Math.PI / 2) / Math.PI * 180;
            }


         } //controller
    } //sunburst
})();

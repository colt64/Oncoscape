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
            sunburstUrl: 'app/components/sunburst/sunburst.html',
            controller: SunburstController,
            controllerAs: 'vm',
            bindToController: true
        };

        return directive;

        /** @ngInject */
        function SunburstController(osApi, $state, $timeout, $scope, $stateParams) {

            if (angular.isUndefined($stateParams.datasource)){
                $state.go("datasource");
                return;
            }

           // Properties
            var cohortPatient = osApi.getCohortPatient();

            // View Model
            var vm = this;
            vm.datasource = $stateParams.datasource;
            vm.filter;
            vm.diagnosisMin = vm.diagnosisMinValue = 1;
            vm.diagnosisMax = vm.diagnosisMaxValue = 99;
            vm.optCohortPatients = cohortPatient.get();
            vm.optCohortPatient = vm.optCohortPatients[0];


         }

           // Load Datasets
            osApi.setBusy(true);
            osApi.setDataset(vm.datasource).then(function() {
                osApi.getPatientHistoryTable(vm.datasource).then(function(response) {
                    rawData = response.payload;


                 });
            });



          var width = 960,
              height = 700,
              radius = Math.min(width, height) / 2;

          var x = d3.scale.linear()
              .range([0, 2 * Math.PI]);

          var y = d3.scale.sqrt()
              .range([0, radius]);

          var color = d3.scale.category20c();

          var svg = d3.select("body").append("svg")
              .attr("width", width)
              .attr("height", height)
            .append("g")
              .attr("transform", "translate(" + width / 2 + "," + (height / 2 + 10) + ")");

          var partition = d3.layout.partition()
              .sort(null)
              .value(function(d) { return d.size; });

          var arc = d3.svg.arc()
              .startAngle(function(d) { return Math.max(0, Math.min(2 * Math.PI, x(d.x))); })
              .endAngle(function(d) { return Math.max(0, Math.min(2 * Math.PI, x(d.x + d.dx))); })
              .innerRadius(function(d) { return Math.max(0, y(d.y)); })
              .outerRadius(function(d) { return Math.max(0, y(d.y + d.dy)); });


          // Keep track of the node that is currently being displayed as the root.
          var node;
          var ptRow = 1

          var root = tree =  d3.json("ImmuneBiomarker_Tree_Lung.json")
          var tbl =  d3.tsv("Lung_FlowPercent_2-29-16.txt")

          var g = svg.selectAll("g")
            .data(partition.nodes(root))
          .enter().append("g");

          var path = g.append("path")
          .attr("d", arc)
          .style("fill", function(d) { return color((d.children ? d : d.parent).name); })
          .on("click", click);

          var text = g.append("text")
          .attr("transform", function(d) { return "rotate(" + computeTextRotation(d) + ")"; })
          .attr("x", function(d) { return y(d.y); })
          .attr("dx", "6") // margin
          .attr("dy", ".35em") // vertical-align
          .text(function(d) { return d.name; });

          d3.selectAll("input").on("change", function change() {
          var value = this.value === "size"
            ? function() { return 1; }
            : function(d) {
              return d.size; };

          path
            .data(
            partition.value(value).nodes)
            .transition()
            .duration(1000)
            .attrTween("d", arcTweenData);
          });


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

          d3.select(self.frameElement).style("height", height + "px");

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

    }
})();

(function() {
    'use strict';

    angular
        .module('oncoscape')
        .directive('osBarplot', barplot);
//

    /** @ngInject */
    function barplot() {

        var directive = {
            restrict: 'E',
            templateUrl: 'app/components/barplot/barplot.html',
            controller: BarplotController,
            controllerAs: 'vm',
            bindToController: true
        };

        return directive;

        /** @ngInject */
        function BarplotController(osApi, $state, $timeout, $scope, $stateParams, d3) {

            if (angular.isUndefined($stateParams.datasource)){
                $state.go("datasource");
                return;
            }

            // Elements
            var d3Chart = d3.select("#barplot-chart").append("svg").attr("id", "chart");
            var tbl;

           // Properties
            var cohortPatient = osApi.getCohortPatient();
            var margin, width, height, x, y, xAxis, yAxis, svg;
            var data = []
            var orderedColNames = ["sample", "a3ss_NumGenes","a5ss_NumGenes","se_NumGenes","ri_NumGenes"]
            var genericColNames = ["sample", "a3ss","a5ss","se","ri"]

            // View Model
            var vm = this;
            vm.datasource = $stateParams.datasource;
            vm.optCohortPatients = cohortPatient.get();
            vm.optCohortPatient = vm.optCohortPatients[0];

            // Event Listeners
          	d3.selectAll("input").on("change", change, this.value);

            osApi.setBusy(true)("Loading Dataset");
            osApi.setDataset(vm.datasource).then(function() {

                // Patient Data
                osApi.getCountsData().then(function(response) {

                    tbl = response.payload;
                    osApi.setBusy(false);

                    setScale()
                    draw()
                });  //ptHistoryTable
            }); //dataset

          function setScale(){

              margin = {top: 20, right: 20, bottom: 100, left: 40},
                  width = 660 - margin.left - margin.right,
                  height = 300 - margin.top - margin.bottom;

              x = d3.scale.ordinal()
                  .rangeRoundBands([0, width], .1);

              y = d3.scale.linear()
                  .range([height, 0]);

              xAxis = d3.svg.axis()
                  .scale(x)
                  .orient("bottom");

              yAxis = d3.svg.axis()
                  .scale(y)
                  .orient("left")
                  .ticks(10);

              svg = d3Chart
                  .attr("width", 2*(width + margin.left + margin.right))
                  .attr("height", 2*(height + margin.top + margin.bottom))

            } //setScale

          function draw(){

                for (var i=0; i < tbl.length; i++) {
                  data[i] = cloneAndPluck(tbl[i], orderedColNames, genericColNames)
              //		data[i].frequency = parseInt(data[i].frequency)
              //		if(isNaN(data[i].frequency)) data[i].frequency = 0;
                }

                x.domain(data.map(function(d) { return d.sample; }));
                y.domain([0, d3.max(data, function(d){ return d3.max([parseInt(d.a3ss), parseInt(d.a5ss), parseInt(d.ri), parseInt(d.se)]) })]);

              //	for(var i=0;i<genericColNames.length-1;i++){
                for(var i=0;i<4;i++){

                  var graph = svg.append("g").attr("id", genericColNames[i+1])
                  graph.attr("transform", "translate(" + (margin.left+i%2*(margin.left+width+margin.right)) + "," + (margin.top+Math.floor(i/2)*(margin.top+height+margin.bottom))+ ")");

                  graph.append("g")
                    .attr("class", "x axis")
                    .attr("transform", "translate(0," + height + ")")
                    .call(xAxis)
                    .selectAll("text")
                  .style("text-anchor", "end")
                  .attr("dx", "-.8em")
                  .attr("dy", ".15em")
                  .attr("transform", function(d) {
                    return "rotate(-65)"
                  });;

                  graph.append("g")
                    .attr("class", "y axis")
                    .call(yAxis)
                  .append("text")
                    .attr("transform", "rotate(-90)")
                    .attr("y", 6)
                    .attr("dy", ".71em")
                    .style("text-anchor", "end")
                    .text(genericColNames[i+1]);

                  graph.selectAll(".bar")
                    .data(data)
                  .enter().append("rect")
                    .attr("class", "bar")
                    .attr("x", function(d) { return x(d.sample); })
                    .attr("width", x.rangeBand())
                    .attr("y", function(d) {
                      return y(parseInt(d[genericColNames[i+1]])); })
                    .attr("height", function(d) {
                      return height - y(parseInt(d[genericColNames[i+1]])); })
                    .attr("fill", function(d){
                    if(d.sample.match(/tumor/) !== null) return "#ff9900";
                    return "steelblue";
                    });
                }


            }  //draw

          function cloneAndPluck(sourceObject, keys, newKeys) {
                var newObject = {};
                for(var i =0;i<keys.length;i++){
                  newObject[newKeys[i]] = sourceObject[keys[i]];
                }
                return newObject;
            }; //cloneAndPluck

          function change() {

                var Valtype = this.value;

                  // Copy-on-write since tweens are evaluated after a delay.
                  var x0 = x.domain(data.sort(this.checked
                      ? function(a, b) {
                        return parseInt(b[Valtype]) - parseInt(a[Valtype]); }
                      : function(a, b) { return d3.ascending(a.sample, b.sample); })
                      .map(function(d) { return d.sample; }))
                      .copy();

                for(var i=0;i<4;i++){

                  var graph = svg.select("#"+genericColNames[i+1])
                      graph.selectAll(".bar")
                           .sort(function(a, b) { return x0(a.sample) - x0(b.sample); });

                  var transition = graph.transition().duration(5),
                      delay = function(d, i) { return i * 50; };

                  transition.selectAll(".bar")
              //        .delay(delay)
                      .attr("x", function(d) { return x0(d.sample); });

                  transition.select(".x.axis")
                   .attr("transform", "translate(0," + height + ")")
                      .call(xAxis)
                      .selectAll("text")
                   .style("text-anchor", "end")
                   .attr("dx", "-.8em")
                   .attr("dy", ".15em")
                   .attr("transform", function(d) {
                      return "rotate(-65)"
                  })
              //      .selectAll("g")
              //        .delay(delay);
                }
          } //change

          function type(d) {
            d.frequency = +d.frequency;
            return d;
          }

         } //BarplotController
    } //barplot
})();

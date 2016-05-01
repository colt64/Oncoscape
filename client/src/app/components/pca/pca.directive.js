(function() {
    'use strict';

    angular
        .module('oncoscape')
        .directive('osPca', explore);

    /** @ngInject */
    function explore() {

        var directive = {
            restrict: 'E',
            templateUrl: 'app/components/pca/pca.html',
            controller: PcaController,
            controllerAs: 'vm',
            bindToController: true
        };

        return directive;

        /** @ngInject */
        function PcaController(osApi, $state, $stateParams, $timeout, $scope, d3, moment, $window, _) {

            if (angular.isUndefined($stateParams.datasource)) {
                $state.go("datasource");
                return;
            }

            // Elements
            var d3Chart = d3.select("#pca-chart").append("svg").attr("id", "chart");

            // Properties
            var cohortPatient = osApi.getCohortPatient();
            var width, height, xScale, yScale, xMax, yMax, xAxis, yAxis;
            var rawData, rawPatientData;

            // View Model
            var vm = this;
            vm.datasource = $stateParams.datasource;
            vm.geneSets = [];
            vm.geneSet = null;
            vm.optNodeColors = [{name: 'Default'}, {name: 'Gender'}, {name: 'Age At Diagnosis'}];
            vm.optNodeColor = vm.optNodeColors[0];
            vm.optCohortPatients = cohortPatient.get();
            vm.optCohortPatient = vm.optCohortPatients[0];

            // Initalizae
            osApi.setBusy(true)("Loading Dataset");
            osApi.setDataset(vm.datasource).then(function(response) {
                var mtx = response.payload.rownames.filter(function(v) {
                    return v.indexOf("mtx.mrna") >= 0
                });

                // Patient Data
                osApi.getPatientHistoryTable(vm.datasource).then(function(response) {

                    rawPatientData = response.payload.tbl;
                    mtx = mtx[mtx.length - 1].replace(".RData", "");
                    osApi.getPCA(vm.datasource, mtx).then(function() {

                        osApi.getGeneSetNames().then(function(response) {

                            // Load Gene Sets
                            vm.geneSets = response.payload;
                            vm.geneSet = vm.geneSets[0];
                            update();

                        }); //geneset
                    }); //pca
                });  //ptHistoryTable
            }); //dataset

            // API Call To Calculate PCA
            var update = function() {
                osApi.setBusyMessage("Calculating PCA");
                osApi.getCalculatedPCA(vm.geneSet).then(function(response) {
                    var payload = response.payload;
                    var scores = payload.scores;
                    var ids = payload.ids;
                    rawData = scores.map(function(d, i) {
                        d.id = ids[i];
                        return d;
                    }, payload.ids);
                    draw()
                    osApi.setBusy(false);
                });
            };

            function setScale() {
                width = $window.innerWidth - 100;
                height = $window.innerHeight - 190;
                if (angular.element(".tray").attr("locked")=="true") width -= 300;

                d3Chart
                    .attr("width", '100%')
                    .attr("height", height);
                xScale = d3.scale.linear()
                    .domain([-xMax, xMax])
                    .range([0, width]).nice();

                yScale = d3.scale.linear()
                    .domain([-yMax, yMax])
                    .range([height, 0]).nice();
            }

            // Render
            function draw() {

                var dataset = rawData;

                var max, min;
                max = Math.abs(d3.max(dataset, function(d) {return +d[0];                }));
                min = Math.abs(d3.min(dataset, function(d) {return +d[0];                }));
                xMax = ((max > min) ? max : min) * 1.2;
                max = Math.abs(d3.max(dataset, function(d) {return +d[1]; }));
                min = Math.abs(d3.min(dataset, function(d) {return +d[1]; }));
                yMax = ((max > min) ? max : min) * 1.2;

                setScale();

                xAxis = d3.svg.axis().scale(xScale).orient("top").ticks(5);

                yAxis = d3.svg.axis().scale(yScale).orient("left").ticks(5);

                var circles = d3Chart.selectAll("circle").data(rawData, function(d) { return d; });

                circles.enter()
                    .append("circle")
                    .attr({"class": "pca-node","cx":  width * .5, "cy": height * .5,"r": 3})
                    .style("fill-opacity", "0")
                    .attr("cx", function(d) { return xScale(d[0]); })
                    .attr("cy", function(d) { return yScale(d[1]); })
                    .style("fill-opacity", 1);

                circles.exit()
                    .attr("cx", width * .5)
                    .attr("cy", height * .5)
                    .style("fill-opacity", "0")
                    .remove();

            // Listen For Resize
            angular.element($window).bind('resize',
                _.debounce(vm.resize, 300)
            );

        } //draw
      } //controller
    } //explore
})();

(function() {
    'use strict';

    angular
        .module('oncoscape')
        .directive('osExample', example);
//

    /** @ngInject */
    function example() {

        var directive = {
            restrict: 'E',
            templateUrl: 'app/components/example/example.html',
            controller: ExampleController,
            controllerAs: 'vm',
            bindToController: true
        };

        return directive;

        /** @ngInject */
        function ExampleController(osApi, $state, $timeout, $scope, $stateParams, d3) {

            if (angular.isUndefined($stateParams.datasource)){
                $state.go("datasource");
                return;
            }

            // Elements
            var d3Chart = d3.select("#example-chart").append("svg").attr("id", "chart");

           // Properties
            var cohortPatient = osApi.getCohortPatient();
            var rawPatientData;

            // View Model
            var vm = this;
            vm.datasource = $stateParams.datasource;
            vm.diagnosisMin = vm.diagnosisMinValue = 1;
            vm.diagnosisMax = vm.diagnosisMaxValue = 99;
            vm.optCohortPatients = cohortPatient.get();
            vm.optCohortPatient = vm.optCohortPatients[0];


            osApi.setBusy(true)("Loading Dataset");
            osApi.setDataset(vm.datasource).then(function(response) {

                // Patient Data
                osApi.getPatientHistoryTable(vm.datasource).then(function(response) {

                    rawPatientData = response.payload.tbl;
                    osApi.setBusy(false);

                    var circles = d3Chart.selectAll("circle").data(rawPatientData, function(d) { return d; });
                    circles.enter()
                        .append("circle")
                        .attr("cx", function(d) { return 10; })
                        .attr("cy", function(d) { return 10; })
                        .attr("r",  function(d) { return 5;})

                });  //ptHistoryTable
            }); //dataset


         }
    }
})();

(function() {
    'use strict';

    angular
        .module('oncoscape')
        .directive('osTemplate', template);
//

    /** @ngInject */
    function template() {

        var directive = {
            restrict: 'E',
            templateUrl: 'app/components/template/template.html',
            controller: TemplateController,
            controllerAs: 'vm',
            bindToController: true
        };

        return directive;

        /** @ngInject */
        function TemplateController(osApi, $state, $timeout, $scope, $stateParams) {

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
    }
})();

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

            // View Model
            var vm = this;
         }
    }
})();

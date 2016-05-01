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
            var ptRow = 2;
            var svg, node;
            var Fullg, path, text;
            var arc;
            var partition = d3.layout.partition()
                .sort(null)
                .value(function(d) { return d.size; });

           // Properties
            var cohortPatient = osApi.getCohortPatient();
            var tbl;
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
            osApi.setDataset(vm.datasource).then(function() {

                // Percent Population Data
              osApi.getFlowData().then(function(response) {

                    tbl = response.payload;
                    node = root;
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
              	node = d;

                path.transition()
                  .duration(1000)
                  .attrTween("d", arcTweenZoom(d))
                  .each("end", function(e) {
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

              //----------------------------------
              Fullg = svg.datum(root).selectAll("path")
                  .data(partition.nodes)
                .enter().append("g")

              path = Fullg.append("path")
                  .attr("d", arc)
                  .attr("id", function(d) { return (d.name)})
                  .style("fill", function(d) { return color((d.children ? d : d.parent).name); })
                  .on("click", click)
                  .each(stash);

             $('g').tipsy({
                  gravity: 'w',
                  title: function() {
                    var d = this.__data__;
                    return d.name + ': '+ tbl[ptRow][d.name] ;
                  }
                });

              text = Fullg.append("text")
                  .attr("transform", function(d) { return "rotate(" + computeTextRotation(d) + ")"; })
                  .attr("x", function(d) { return y(d.y); })
                  .attr("dx", "6") // margin
                  .attr("dy", ".35em") // vertical-align
                  .text(function(d) { return d.name; });

              d3.selectAll("input").on("change", function change() {
                var value = this.value === "size"
                  ? function() { return 1; }
                  : function(d){ return getRelativeSize(d) };

                path.data(partition.value(value).nodes)
                  .transition()
                  .duration(1000)
                  .attrTween("d", arcTweenData);

                var arcText = Fullg.selectAll("text");
                      // fade in the text element and recalculate positions
                      arcText.transition().duration(750)
                        .attr("opacity", 1)
                        .attr("transform", function(d) { return "rotate(" + computeTextRotation(d) + ")" })
                        .attr("x", function(d) { return y(d.y); });

              }); // change input path

            }

            //----------------------------------
            function getRelativeSize(d){
                if(d.depth == 0){ return +tbl[ptRow][d.name];}
                else return getRelativeSize(d.parent) * +tbl[ptRow][d.name] /100;
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


(function($) {

    function maybeCall(thing, ctx) {
        return (typeof thing == 'function') ? (thing.call(ctx)) : thing;
    }

    function Tipsy(element, options) {
        this.$element = $(element);
        this.options = options;
        this.enabled = true;
        this.fixTitle();
    }

    Tipsy.prototype = {
        show: function() {
            var title = this.getTitle();
            if (title && this.enabled) {
                var $tip = this.tip();

                $tip.find('.tipsy-inner')[this.options.html ? 'html' : 'text'](title);
                $tip[0].className = 'tipsy'; // reset classname in case of dynamic gravity
                $tip.remove().css({top: 0, left: 0, visibility: 'hidden', display: 'block'}).prependTo(document.body);

                var pos = $.extend({}, this.$element.offset(), {
                    width: this.$element[0].offsetWidth || 0,
                    height: this.$element[0].offsetHeight || 0
                });

                if (typeof this.$element[0].nearestViewportElement == 'object') {
                    // SVG
					var el = this.$element[0];
                    var rect = el.getBoundingClientRect();
					pos.width = rect.width;
					pos.height = rect.height;
                }


                var actualWidth = $tip[0].offsetWidth,
                    actualHeight = $tip[0].offsetHeight,
                    gravity = maybeCall(this.options.gravity, this.$element[0]);

                var tp;
                switch (gravity.charAt(0)) {
                    case 'n':
                        tp = {top: pos.top + pos.height + this.options.offset, left: pos.left + pos.width / 2 - actualWidth / 2};
                        break;
                    case 's':
                        tp = {top: pos.top - actualHeight - this.options.offset, left: pos.left + pos.width / 2 - actualWidth / 2};
                        break;
                    case 'e':
                        tp = {top: pos.top + pos.height / 2 - actualHeight / 2, left: pos.left - actualWidth - this.options.offset};
                        break;
                    case 'w':
                        tp = {top: pos.top + pos.height / 2 - actualHeight / 2, left: pos.left + pos.width + this.options.offset};
                        break;
                }

                if (gravity.length == 2) {
                    if (gravity.charAt(1) == 'w') {
                        tp.left = pos.left + pos.width / 2 - 15;
                    } else {
                        tp.left = pos.left + pos.width / 2 - actualWidth + 15;
                    }
                }

                $tip.css(tp).addClass('tipsy-' + gravity);
                $tip.find('.tipsy-arrow')[0].className = 'tipsy-arrow tipsy-arrow-' + gravity.charAt(0);
                if (this.options.className) {
                    $tip.addClass(maybeCall(this.options.className, this.$element[0]));
                }

                if (this.options.fade) {
                    $tip.stop().css({opacity: 0, display: 'block', visibility: 'visible'}).animate({opacity: this.options.opacity});
                } else {
                    $tip.css({visibility: 'visible', opacity: this.options.opacity});
                }

                var t = this;
                var set_hovered  = function(set_hover){
                    return function(){
                        t.$tip.stop();
                        t.tipHovered = set_hover;
                        if (!set_hover){
                            if (t.options.delayOut === 0) {
                                t.hide();
                            } else {
                                setTimeout(function() {
                                    if (t.hoverState == 'out') t.hide(); }, t.options.delayOut);
                            }
                        }
                    };
                };
               $tip.hover(set_hovered(true), set_hovered(false));
            }
        },

        hide: function() {
            if (this.options.fade) {
                this.tip().stop().fadeOut(function() { $(this).remove(); });
            } else {
                this.tip().remove();
            }
        },

        fixTitle: function() {
            var $e = this.$element;

            if ($e.attr('title') || typeof($e.attr('original-title')) != 'string') {
                $e.attr('original-title', $e.attr('title') || '').removeAttr('title');
            }
            if (typeof $e.context.nearestViewportElement == 'object'){
                if ($e.children('title').length){
                    $e.append('<original-title>' + ($e.children('title').text() || '') + '</original-title>')
                        .children('title').remove();
                }
            }
        },

        getTitle: function() {

            var title, $e = this.$element, o = this.options;
            this.fixTitle();

            if (typeof o.title == 'string') {
                var title_name = o.title == 'title' ? 'original-title' : o.title;
                if ($e.children(title_name).length){
                    title = $e.children(title_name).html();
                } else{
                    title = $e.attr(title_name);
                }

            } else if (typeof o.title == 'function') {
                title = o.title.call($e[0]);
            }
            title = ('' + title).replace(/(^\s*|\s*$)/, "");
            return title || o.fallback;
        },

        tip: function() {
            if (!this.$tip) {
                this.$tip = $('<div class="tipsy"></div>').html('<div class="tipsy-arrow"></div><div class="tipsy-inner"></div>');
            }
            return this.$tip;
        },

        validate: function() {
            if (!this.$element[0].parentNode) {
                this.hide();
                this.$element = null;
                this.options = null;
            }
        },

        enable: function() { this.enabled = true; },
        disable: function() { this.enabled = false; },
        toggleEnabled: function() { this.enabled = !this.enabled; }
    };

    $.fn.tipsy = function(options) {

        if (options === true) {
            return this.data('tipsy');
        } else if (typeof options == 'string') {
            var tipsy = this.data('tipsy');
            if (tipsy) tipsy[options]();
            return this;
        }

        options = $.extend({}, $.fn.tipsy.defaults, options);

        if (options.hoverlock && options.delayOut === 0) {
	    options.delayOut = 100;
	}

        function get(ele) {
            var tipsy = $.data(ele, 'tipsy');
            if (!tipsy) {
                tipsy = new Tipsy(ele, $.fn.tipsy.elementOptions(ele, options));
                $.data(ele, 'tipsy', tipsy);
            }
            return tipsy;
        }

        function enter() {
            var tipsy = get(this);
            tipsy.hoverState = 'in';
            if (options.delayIn === 0) {
                tipsy.show();
            } else {
                tipsy.fixTitle();
                setTimeout(function() { if (tipsy.hoverState == 'in') tipsy.show(); }, options.delayIn);
            }
        }

        function leave() {
            var tipsy = get(this);
            tipsy.hoverState = 'out';
            if (options.delayOut === 0) {
                tipsy.hide();
            } else {
                var to = function() {
                    if (!tipsy.tipHovered || !options.hoverlock){
                        if (tipsy.hoverState == 'out') tipsy.hide();
                    }
                };
                setTimeout(to, options.delayOut);
            }
        }

        if (options.trigger != 'manual') {
            var binder = options.live ? 'live' : 'bind',
                eventIn = options.trigger == 'hover' ? 'mouseenter' : 'focus',
                eventOut = options.trigger == 'hover' ? 'mouseleave' : 'blur';
            this[binder](eventIn, enter)[binder](eventOut, leave);
        }

        return this;

    };

    $.fn.tipsy.defaults = {
        className: null,
        delayIn: 0,
        delayOut: 0,
        fade: false,
        fallback: '',
        gravity: 'n',
        html: false,
        live: false,
        offset: 0,
        opacity: 0.8,
        title: 'title',
        trigger: 'hover',
        hoverlock: false
    };

    // Overwrite this method to provide options on a per-element basis.
    // For example, you could store the gravity in a 'tipsy-gravity' attribute:
    // return $.extend({}, options, {gravity: $(ele).attr('tipsy-gravity') || 'n' });
    // (remember - do not modify 'options' in place!)
    $.fn.tipsy.elementOptions = function(ele, options) {
        return $.metadata ? $.extend({}, options, $(ele).metadata()) : options;
    };

    $.fn.tipsy.autoNS = function() {
        return $(this).offset().top > ($(document).scrollTop() + $(window).height() / 2) ? 's' : 'n';
    };

    $.fn.tipsy.autoWE = function() {
        return $(this).offset().left > ($(document).scrollLeft() + $(window).width() / 2) ? 'e' : 'w';
    };

    /**
     * yields a closure of the supplied parameters, producing a function that takes
     * no arguments and is suitable for use as an autogravity function like so:
     *
     * @param margin (int) - distance from the viewable region edge that an
     *        element should be before setting its tooltip's gravity to be away
     *        from that edge.
     * @param prefer (string, e.g. 'n', 'sw', 'w') - the direction to prefer
     *        if there are no viewable region edges effecting the tooltip's
     *        gravity. It will try to vary from this minimally, for example,
     *        if 'sw' is preferred and an element is near the right viewable
     *        region edge, but not the top edge, it will set the gravity for
     *        that element's tooltip to be 'se', preserving the southern
     *        component.
     */
     $.fn.tipsy.autoBounds = function(margin, prefer) {
		return function() {
			var dir = {ns: prefer[0], ew: (prefer.length > 1 ? prefer[1] : false)},
			    boundTop = $(document).scrollTop() + margin,
			    boundLeft = $(document).scrollLeft() + margin,
			    $this = $(this);

			if ($this.offset().top < boundTop) dir.ns = 'n';
			if ($this.offset().left < boundLeft) dir.ew = 'w';
			if ($(window).width() + $(document).scrollLeft() - $this.offset().left < margin) dir.ew = 'e';
			if ($(window).height() + $(document).scrollTop() - $this.offset().top < margin) dir.ns = 's';

			return dir.ns + (dir.ew ? dir.ew : '');
		};
    };
})(jQuery);

// template

(function() {

	var ContigData = function() {

	var parseData = function (data) {
		chart = { assemblies: {} };

		for (var i = 0; i < data.length - 1; i++) {
			if (!chart.assemblies[data[i].assembly])
				chart.assemblies[data[i].assembly] = []

			var sublane = 0;

			while (isOverlapping(data[i], chart.assemblies[data[i].assembly][sublane]))
				sublane++;

			if (!chart.assemblies[data[i].assembly][sublane])
				chart.assemblies[data[i].assembly][sublane] = [];

			chart.assemblies[data[i].assembly][sublane].push(data[i]);
		}		

		return collapseLanes(chart);
	};

	var isOverlapping = function(item, lane) {
		if (lane)
			for (var i = 0; i < lane.length; i++)
				if (item.start <= lane[i].end && lane[i].start <= item.end)
					return true;

		return false;
	};


	var collapseLanes = function (chart) {
		var lanes = [], items = [], laneId = 0;

		for (var assemblyName in chart.assemblies) {
			var lane = chart.assemblies[assemblyName];

			for (var i = 0; i < lane.length; i++) {
				var subLane = lane[i];

				lanes.push({
					id: laneId, 
					label: i === 0 ? assemblyName : ''
				});

				for (var j = 0; j < subLane.length; j++) {
					var item = subLane[j];
					item.lane = laneId
					items.push(item);
				}

				laneId++;
			}
		}

		return {lanes: lanes, items: items};
	}

		var randomNumber = function(min, max) {
			return Math.floor(Math.random(0, 1) * (max - min)) + min;
		};

		// return parseData(generateRandomWorkItems());
		return parseData(contig_data);
	};

/**
* Allow library to be used within both the browser and node.js
*/
var root = typeof exports !== "undefined" && exports !== null ? exports : window;
root.contigData = ContigData;

}).call(this);

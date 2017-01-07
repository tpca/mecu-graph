/*
 * mecu-viz
 * https://github.com/sacdallago/mecu-viz
 *
 * Copyright (c) 2016 Christian Dallago
 * Licensed under the Apache-2.0 license.
 */

const cytoscape = require("cytoscape");
const util = require(__dirname + "/util");

/**
 @class MecuGraph
 */

class MecuGraph {

    constructor(opts) {

        if (typeof(opts) === 'object' && opts.constructor !== Array) {
            this.options = opts || {};
        } else {
            this.options = {};
        }

        this.base = this.options.element || "#mecugraph";
        this.data = [];

        this.graph = cytoscape({container: document.getElementById('cy')});

        return;
    }

    add(data) {
        if (typeof(data) === 'undefined' || data === null || typeof(data) !== 'object') {
            return;
        }

        if (data.constructor !== Array) {
            data = [data];
        }

        this.data = this.data.concat(data);

        for (let i = 0; i < data.length; i++) {
            // Add the valueline path.
            this.lineSvg.append("path")
                .attr("class", "line MECU"+this.data[i].experiment)
                .attr("fill", "none")
                .attr("transform", "translate(" + this.margin.left + ",0)")
                .attr("stroke", this.options.strokeColor || util.sequentialColor(this.data.length, i) || "steelblue")
                .attr("stroke-width", this.options.strokeWidth || '1em')
                .attr("d", this.valueline(this.data[i].reads))
                //.data(data, function(d) { return d._id; })
                .style("stroke-linecap", "round");
        }
    }

    remove(data) {
        if (typeof(data) === 'undefined' || data === null) {
            return;
        }

        if (data.constructor !== Array && typeof(data) === 'object') {
            data = [data];
        } else {
            return;
        }

        // TODO: Data can be stored in array in binary fashion, so removal is speed up, but insertion is slower. --> Are there more inserts or removes?
        this.data = this.data.filter(function (listItem) {
            return data.find(function (match) {
                return match.experiment === listItem.experiment;
            })
        });

        return this.add(this.data);
    }
}

module.exports = Mecu;
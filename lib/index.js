/*
 * mecu-graph
 * https://github.com/sacdallago/mecu-graph
 *
 * Copyright (c) 2016 Christian Dallago
 * Licensed under the Apache-2.0 license.
 */

const cytoscape = require('cytoscape');
const cycola = require('cytoscape-cola');
const cola = require('webcola');
cycola(cytoscape, cola); // register extension
const everpolate = require('everpolate');

// Go from the API response format to a Protein+Experiment+TempReads format
// Optionally get the UNSORTED!! range of temperatures
let extrapolate = function(data, getTemps = false) {

    let result = [];
    let temperatures;

    // Put if outside loops to reduce amounts of checks. Downside: Duplicated code
    if(getTemps){
        temperatures = new Set();
        data.forEach(function(protein) {
            protein.reads.forEach(function(read) {
                let sortedReads = read.reads.sort(function(a,b) {

                    temperatures.add(a.t);
                    temperatures.add(b.t);

                    return a.t < b.t;
                });
                result.push({
                    p: protein.uniprotId,
                    e: read.experiment,
                    r: sortedReads
                });
            })
        });
    } else {
        data.forEach(function(protein) {
            protein.reads.forEach(function(read) {
                let sortedReads = read.reads.sort(function(a,b) {
                    return a.t < b.t;
                });
                result.push({
                    p: protein.uniprotId,
                    e: read.experiment,
                    r: sortedReads
                });
            })
        });
    }

    return [result, temperatures];
};

// Get two vectors containing temperatures and corresponding ratios
let extractTemperaturesAndRatio = function(extrapolatedProtein, temperatures = undefined) {
    if(temperatures !== undefined){

        let [originTemp, originRatio] = [extrapolatedProtein.r.map(function(d){
            return d.t;
        }), extrapolatedProtein.r.map(function(d){
            return d.r;
        })];

        // TODO - make the interpolation more efficient by storing the function globally and calling it only on undefined
        let ratios = [];
        temperatures.forEach(function(temp){

            let ratio = extrapolatedProtein.r.filter(function(read) {
                return read.t == temp;
            }).r;

            if(ratio !== undefined){
                ratios.push(r);
            } else {
                ratios.push(everpolate.polynomial(temp, originTemp, originRatio)[0]);
            }
        });

        return [temperatures, ratios];
    } else {
        return [extrapolatedProtein.r.map(function(d){
            return d.t;
        }), extrapolatedProtein.r.map(function(d){
            return d.r;
        })];
    }
};

// Simple Euclidian distance between two vectors
// Assumes equal length of vectors and matching array positions
let euclidian = function(A,B) {
    let sum = 0;

    for(let i=0; i<A.length; i++){
        sum += Math.pow(A[i] - B[i], 2);
    }
    return Math.sqrt(sum);
};


// Format a collection of proteins in the vectorized format (uniprot+experiment+tempReads) into a cytoscape digestable format
let getGraphDataFromExtractedProteins = function(proteins, distance = euclidian){
    // Extrapolate vectors
    let vectors = [];
    let nodes = [];
    let edges = [];

    // Create a mapping between proteins and vectors
    for(let i=0; i<proteins.length; i++){
        let [temp, ratios] = extractTemperaturesAndRatio(proteins[i]);
        vectors.push(ratios);

        nodes.push({
            data: {
                id: proteins[i].p + "-E" + proteins[i].e
            }
        });
    }

    // Find distance between item and everyone else
    for(let i=0; i<proteins.length; i++){
        for(let j=i+1; j<proteins.length; j++){
            edges.push({
                data: {
                    id: proteins[i].p + "-E" + proteins[i].e + "<->" + proteins[j].p + "-E" + proteins[j].e,
                    source: proteins[i].p + "-E" + proteins[i].e,
                    target: proteins[j].p + "-E" + proteins[j].e,
                    weight: distance(vectors[i], vectors[j])
                }
            });
        }
    }

    return [nodes, edges]
};

// Extend string prototypes to allow color hashing
String.prototype.getHashCode = function() {
    var hash = 0;
    if (this.length == 0) return hash;
    for (var i = 0; i < this.length; i++) {
        hash = this.charCodeAt(i) + ((hash << 5) - hash);
        hash = hash & hash; // Convert to 32bit integer
    }
    return hash;
};
Number.prototype.intToHSL = function() {
    var shortened = this % 360;
    return "hsl(" + shortened + ",100%,40%)";
};

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

        this.base = this.options.element || "#mecuGraph";
        this.nodes = [];
        this.edges = [];
        this.data = [];

        this.options.edgeLength = this.options.edgeLength || 100;

        let self = this;

        this.graph = cytoscape({
            container: document.getElementById(this.base.substr(1,this.base.length)),
            elements: [],
            style: cytoscape.stylesheet()
                .selector('node')
                .style({
                    'background-color': function (e) {
                        return e.data('id').getHashCode().intToHSL();
                    },
                    'label': 'data(id)',
                    'text-opacity': '0.4'
                })
                .selector('edge')
                .style({
                    'opacity': '0.2'
                })
            ,
            layout: {
                name: 'cola',
                fit: false,
                infinite: true,
                edgeLength: function(e){ return (e.data('weight') || 0.001)*self.options.edgeLength; }
            }
        });
    }

    add(data){
        if (typeof(data) === 'undefined' || data === null || typeof(data) !== 'object') {
            return;
        }

        if (data.constructor !== Array) {
            data = [data];
        }

        let [extrapolatedProteins,] = extrapolate(data);

        let self = this;

        let realNewOnes = extrapolatedProteins.filter(function(newElement) {
            return self.data.filter(function(oldElement) {
                return oldElement.p == newElement.p && oldElement.e == newElement.e;
            })
        });

        this.data = this.data.concat(realNewOnes);

        // Use all the data because you need all the edges!
        let [nodes, edges] = getGraphDataFromExtractedProteins(this.data, this.options.distanceMetric);

        this.nodes = nodes;
        this.edges = edges;

        let result = [].concat(this.nodes).concat(this.edges);

        this.graph.add(result);
        this.setEdgeLength(this.options.edgeLength);
        return result;
    }

    remove(data){
        if (typeof(data) === 'undefined' || data === null || typeof(data) !== 'object') {
            return;
        }

        if (data.constructor !== Array) {
            data = [data];
        }

        let [extrapolatedProteins,] = extrapolate(data);
        let self = this;

        extrapolatedProteins.forEach(function(element) {
            let item = self.graph.elements("node[id = '"+ element.p + "-E" + element.e + "']");
            self.graph.remove(item);
        });

        this.data = this.data.filter(function(element){
            return !extrapolatedProteins.find(function(exProt) {
                return element.p == exProt.p && element.e == exProt.e;
            });
        });

        let [nodes, edges] = getGraphDataFromExtractedProteins(this.data, this.options.distanceMetric);

        this.nodes = nodes;
        this.edges = edges;

        let result = [].concat(nodes).concat(edges);

        return this.data;
    }

    setEdgeLength(length = 100){
        this.options.edgeLength = length;
        let self = this;
        this.graph.layout({
            name: 'cola',
            fit: false,
            infinite: true,
            edgeLength: function(e){ return e.data('weight')*self.options.edgeLength; }
        });
    }

    applyMetric(distance){
        this.options.distanceMetric = distance;
    }
}

module.exports = MecuGraph;
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

/**
 * Go from the API response format to a Protein+Experiment+TempReads format
 *  -- Optionally get the UNSORTED!! range of temperatures
 *
 * @param data - Data as from API call
 * @param getTemps - optionally normalize temperatures across dataset and return
 * @returns compoundElement - the first returned element represents digestible data for the internal methods,
 * the second returned element will be non-empty array of temperatures only if parameter getTemps has been passed as true.
 */
let extrapolate = function(data, getTemps = false) {

    let result = [];
    let temperatures;

    // Put if outside loops to reduce amounts of checks. Downside: Duplicated code
    if(getTemps === true){
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

/**
 * Get two vectors containing temperatures and corresponding ratios
 *
 * @param extrapolatedProtein - digestible melting reads data
 * @param temperatures - Optional array of temperatures for which one wants the ratios. Otherwise inferred from the data.
 * @returns compoundElement - compoundElement = [temperatures, ratio]
 */
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

/**
 * Simple implementation of Euclidian distance between two vectors.
 * Assumes equal length of vectors and matching array positions.
 *
 * @param A - A vector of length N
 * @param B - Another vector of length N
 * @returns {number} - Euclidian distance between two vectors
 */
let euclidian = function(A,B) {
    let sum = 0;

    for(let i=0; i<A.length; i++){
        sum += Math.pow(A[i] - B[i], 2);
    }
    return Math.sqrt(sum);
};


/**
 * Format a collection of proteins in the vectorized format (uniprot+experiment+tempReads) into a cytoscape digestible format
 *
 * @param proteins
 * @param distance
 * @returns {Number}
 */
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

/**
 * Extend string prototypes to allow color hashing
 */
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

    /**
     * Constructor for a Mecu Graph element. This will initialize the element in the DOM, but not add any data to it.
     *
     *
     * @param opts - an object containing customization parameters. See README for more information
     */
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

        this.options.edgeLength = this.options.edgeLength || 200;

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
                    'text-opacity': '0.4',
                    'line-color': '#e2dfe1',
                })
                .selector('edge')
                .style({
                    'label': function(e){
                        return Math.round(e.data('weight')*100)/100;
                    },
                    'text-opacity': '0',
                    'opacity': '0.2',
                    'line-opacity': '0',
                    'line-color': '#e2dfe1',
                })
                .selector('edge:selected')
                .css({
                    'opacity': '1',

                    'text-opacity': '1'
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

    /**
     * Method to add nodes (and edges) to the initialized graph.
     *
     *
     * @param proteins - and object or array of objects representing a protein melting read. The object must comply with the format {uniprotId:..., reads: [ { experiment:..., reads: [ { t:.., r:... }]}]}
     */
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

    /**
     * Method to remove nodes (and edges) from the initialized graph.
     *
     *
     * @param proteins - and object or array of objects representing a protein melting read. The object must comply with the format {uniprotId:..., reads: [ { experiment:..., reads: [ { t:.., r:... }]}]}
     */
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

        return result;
    }

    /**
     * Redefine the distance between the edges (weights) by an arbitrary function
     *
     * @param distanceMetric - An arbitrary function accepting two vector parameters
     * @returns {Array.<*>} - An array of object being represented in the graph
     */
    changeDistanceMetric(distanceMetric){
        this.options.distanceMetric = distanceMetric;

        let [nodes, edges] = getGraphDataFromExtractedProteins(this.data, this.options.distanceMetric);

        this.nodes = nodes;
        this.edges = edges;

        let result = [].concat(this.nodes).concat(this.edges);

        this.graph.remove(this.graph.elements("node"));
        this.graph.add(result);
        this.setEdgeLength(this.options.edgeLength);
        return result;
    }

    /**
     * An arbitrary value being multiplied to the weight. Useful to minimize effect of big distances or small similarities
     *
     * @param length - a Number
     */
    setEdgeLength(length){
        this.options.edgeLength = length || this.options.edgeLength;
        let self = this;
        this.graph.layout({
            name: 'cola',
            fit: false,
            infinite: true,
            edgeLength: function(e){ return e.data('weight')*self.options.edgeLength; }
        });
    }
}

module.exports = MecuGraph;
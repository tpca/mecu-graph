/*
 * mecu-viz
 * https://github.com/sacdallago/mecu-viz
 *
 * Copyright (c) 2016 Christian Dallago
 * Licensed under the Apache-2.0 license.
 */

const cytoscape = require('cytoscape');
const cycola = require('cytoscape-cola');
const cola = require('webcola');
cycola(cytoscape, cola); // register extension
const everpolate = require('everpolate');
const util = require(__dirname + '/util');

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

let euclidian = function(A,B) {
    let sum = 0;

    for(let i=0; i<A.length; i++){
        sum += Math.pow(A[i] - B[i], 2);
    }
    return Math.sqrt(sum);
};

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
    }

    add(data) {
        if (typeof(data) === 'undefined' || data === null || typeof(data) !== 'object') {
            return;
        }

        if (data.constructor !== Array) {
            data = [data];
        }

        let [extrapolatedProteins,] = extrapolate(data);

        let [nodes, edges] = getGraphDataFromExtractedProteins(extrapolatedProteins);

        this.nodes = this.nodes.concat(nodes);
        this.edges = this.edges.concat(edges);

        return [].concat(nodes).concat(edges);
    }

    render(){
        let data = [].concat(this.edges).concat(this.nodes);
        this.graph = cytoscape({
            container: document.getElementById(this.base.substr(1,this.base.length)),
            elements: data,
            style: cytoscape.stylesheet()
                .selector('node')
                .style({
                    'background-color': 'blue'
                }),
        });

        this.graph.layout({
            name: 'cola',
            fit: false,
            infinite: true,
            edgeLength: function(e){ return 100 / e.data('weight'); }
        });

        // cytoscape({
        //     container: document.getElementById(this.base.substr(1,this.base.length)),
        //     elements: [{"data":{"id":"P56192-E2"}},{"data":{"id":"P56192-E3"}},{"data":{"id":"P56537-E2"}},{"data":{"id":"P56537-E3"}},{"data":{"id":"P56192-E2<->P56192-E3","source":"P56192-E2","target":"P56192-E3","weight":0.2767654304708042}},{"data":{"id":"P56192-E2<->P56537-E2","source":"P56192-E2","target":"P56537-E2","weight":1.0218655049960128}},{"data":{"id":"P56192-E2<->P56537-E3","source":"P56192-E2","target":"P56537-E3","weight":1.2633673450314178}},{"data":{"id":"P56192-E3<->P56537-E2","source":"P56192-E3","target":"P56537-E2","weight":0.7996353502444211}},{"data":{"id":"P56192-E3<->P56537-E3","source":"P56192-E3","target":"P56537-E3","weight":1.0413743923456298}},{"data":{"id":"P56537-E2<->P56537-E3","source":"P56537-E2","target":"P56537-E3","weight":0.24508149183520306}}]
        // });
    }

    // remove(data) {
    //     if (typeof(data) === 'undefined' || data === null) {
    //         return;
    //     }
    //
    //     if (data.constructor !== Array && typeof(data) === 'object') {
    //         data = [data];
    //     } else {
    //         return;
    //     }
    //
    //     // TODO: Data can be stored in array in binary fashion, so removal is speed up, but insertion is slower. --> Are there more inserts or removes?
    //     this.data = this.data.filter(function (listItem) {
    //         return data.find(function (match) {
    //             return match.experiment === listItem.experiment;
    //         })
    //     });
    //
    //     return this.add(this.data);
    // }
}

module.exports = MecuGraph;
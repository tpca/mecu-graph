module.exports = {
    entry: __dirname + '/lib/index.js',
    output: {
         path: __dirname + '/build',
         filename: 'mecu-graph.js',
         libraryTarget: 'var',
         library: 'MecuGraph'
    }
 };
module.exports = {
    entry: __dirname + '/lib/index.js',
    output: {
         path: __dirname + '/build',
         filename: 'mecu.js',
         libraryTarget: 'var',
         library: 'Mecu'
    }
 };
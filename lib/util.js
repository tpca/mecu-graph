module.exports = {
    sequentialColor(numOfSteps, step) {
        // This function generates vibrant, "evenly spaced" colours (i.e. no clustering). This is ideal for creating easily distinguishable vibrant markers in Google Maps and other apps.
        // Adam Cole, 2011-Sept-14
        // HSV to RBG adapted from: http://mjijackson.com/2008/02/rgb-to-hsl-and-rgb-to-hsv-color-model-conversion-algorithms-in-javascript
        var r, g, b;
        var h = step / numOfSteps;
        var i = ~~(h * 6);
        var f = h * 6 - i;
        var q = 1 - f;
        switch(i % 6){
            case 0: r = 1; g = f; b = 0; break;
            case 1: r = q; g = 1; b = 0; break;
            case 2: r = 0; g = 1; b = f; break;
            case 3: r = 0; g = q; b = 1; break;
            case 4: r = f; g = 0; b = 1; break;
            case 5: r = 1; g = 0; b = q; break;
        }
        var c = "#" + ("00" + (~ ~(r * 255)).toString(16)).slice(-2) + ("00" + (~ ~(g * 255)).toString(16)).slice(-2) + ("00" + (~ ~(b * 255)).toString(16)).slice(-2);
        return (c);
    },
    stringToColour(str) {
        var hash = 0;
        for (var i = 0; i < str.length; i++) {
            hash = str.charCodeAt(i) + ((hash << 5) - hash);
        }
        var colour = '#';
        for (var i = 0; i < 3; i++) {
            var value = (hash >> (i * 8)) & 0xFF;
            colour += ('00' + value.toString(16)).substr(-2);
        }
        return colour;
    }
}



/*
 function colorByHashCode(value) {
 return "<span style='color:" + value.getHashCode().intToHSL() + "'>" + value + "</span>";
 }
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
 return "hsl(" + shortened + ",100%,30%)";
 };


 In HSL its Hue, Saturation, Lightness. So the hue between 0-359 will get all colors, saturation is how rich you want the color, 100% works for me. And Lightness determines the deepness, 50% is normal, 25% is dark colors, 75% is pastel. I have 30% because it fit with my color scheme best.





 */
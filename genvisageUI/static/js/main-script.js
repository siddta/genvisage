/**
 * Created by an5ra on 6/5/2016.
 */
var currentOptions = {
    xAxis: 1,
    yAxis: 2,
    zAxis: 3, 
    algorithm: 'Naive Algorithm'
};
var resultsWindow;
currentRepresentativePlot = {}

var currentData = {};

function updatePaneHeight() {
    var width = (window.innerWidth > 0) ? window.innerWidth : screen.width;
    if (width <= 990) {
        return;
    }
    $('.colored-pane').each(function (index) {
        $(this).css('height', $('.main-row').height());
    });
}

$(document).ready(function () {
    $('select').material_select();
    setTimeout(updatePaneHeight, 300);
    setDataset();
    updateRepresentativePlot();

});

/**
 *
 * @param dataset
 * @param xAxis
 * @param yAxis
 * @param zAxis
 */
 function receiveGeneIdentifiers(){
    var txtName = document.getElementById("txtName");
    var txtOutput = document.getElementById("txtOutput");

    var genes = txtName.value;
    var genesArray = new Array();
    genesArray = genes.split(",");

    txtOutput.value = "Positive Set Received: " + genesArray;

    return genesArray;
  } 


function formatDataset(renderTo, dataset, dataset2, xAxis, yAxis, zAxis) {
    var data = []; //If inputted Gene Identifiers is empty, this will contain all data points. Otherwise data contains postive genes
    var data2 = []; //If there are inputted gene identifers, tis data set contains all other points

    var geneIdentifiers = receiveGeneIdentifiers(); 

    for (var i = 0; i < dataset.data.length; i++) {

            if(geneIdentifiers.length == 1){
                data.push({x: dataset.data[i][xAxis], y: dataset.data[i][yAxis]}); 
            }
            else{
                if(checkArray(geneIdentifiers, dataset.data[i][0])){
                    data.push({x: dataset.data[i][xAxis], y: dataset.data[i][yAxis]});  
                }    
                else{
                    data2.push({x: dataset.data[i][xAxis], y: dataset.data[i][yAxis]});
                }
            }
    }

//    window.alert("Data 2: " + data2);

    var options = {
        'renderTo': renderTo,
        'data': data,
        'data2': data2,
        'x-label': dataset.column_names[xAxis] ,
        'y-label': dataset.column_names[yAxis] ,
        'title': dataset.dataset_name
    };
    return options;
}

function checkArray(arr, val){
    var check = false; 
    for (var i = 0; i < arr.length; i++){
        if(arr[i].valueOf() == val.valueOf()){
            check = true; 
            break; 
        }
    }
    return check; 
}

function updateOptions(index, value) {
    enableButton('update-dataset-settings');
    if (index == 0) currentOptions.xAxis = parseInt(value);
    if (index == 1) currentOptions.yAxis = parseInt(value);
    if (index == 2) currentOptions.zAxis = parseInt(value);
    if (index == 3) currentOptions.algorithm = value;

}

function updateRepresentativePlot() {
    var options = formatDataset("#main-chart", currentData, currentOptions.xAxis, currentOptions.yAxis, currentOptions.zAxis);
    //window.alert("Data2 length: " + data2.length);
    //if (data2.length in options > 1)
    //    var options 2 = formatDataset("#main-chart2", currentData, currentOptions.xAxis, currentOptions.yAxis, currentOptions.zAxis);

    currentRepresentativePlot = drawScatter(options);
    //currentRepresentativePlot2 = drawScatter(options2);

}

function updateResults(results, xMax, yMax) {
    for (var i = 0; i < results.result.length; i++) {
        var resultID = "#result-svg-" + i;
        var options = formatDataset(resultID, results.result[i], results.x, results.y, currentOptions.zAxis);
        options.xMax = xMax;
        options.yMax = yMax;
        drawScatter(options);
    }
}

function enableButton(buttonId) {
    $('#' + buttonId).removeClass('disabled');
}

function disableButton(buttonId) {
    $('#' + buttonId).addClass('disabled');
}

// function getPolygons() {
//     var result = [];
//     var polygons = $('.userPolygon');
//     for (var i = 0; i < polygons.length; i++) {
//         var polygon = polygons[i];
//         var newPolygon = {};
//         // points
//         newPolygon.points = [];
//         for (var j = 0; j < polygon.points.length; j++) {
//             var xCoordinate = polygon.points[j].x;
//             var scaledXCoordinate = currentRepresentativePlot.xScale.invert(xCoordinate - currentRepresentativePlot.leftMargin);
//             var yCoordinate = polygon.points[j].y;
//             var scaledYCoordinate = currentRepresentativePlot.yScale.invert(yCoordinate - currentRepresentativePlot.topMargin);
//             var currentCoordinates = [scaledXCoordinate, scaledYCoordinate];
//             // checking for duplicate last two points
//             //if (!newPolygon.points[newPolygon.points.length - 1] == currentCoordinates)
//             newPolygon.points.push(currentCoordinates);
//         }
//         // type
//         var polygonClasses = $(polygon).attr('class');
//         if (polygonClasses.indexOf('red') > -1) {
//             newPolygon.type = 'red';
//         }
//         else newPolygon.type = "green";
//         result.push(newPolygon);
//     }
//     return result;
// }

function getResults() {
    var data = {polygons: getPolygons(), options: currentOptions};
    $.ajax({
        type: 'POST',
        contentType: 'application/json',
        url: '/getResults',
        dataType: 'html',
        data: JSON.stringify(data),
        success: function (results) {
            if (!resultsWindow) {
                resultsWindow = window.open();
            }
            else {
                resultsWindow.close();
                resultsWindow = window.open();
            }

            var w = resultsWindow;

            w.document.write(results);
            w.focus();
            setTimeout(function () {
                var results = w.getResults();
                w.updateResults(results, currentRepresentativePlot.xMax, currentRepresentativePlot.yMax);
            }, 1000);

            //var updateResultPage = function () {
            //    if (!w.getResults) {
            //        setTimeout(function () {
            //            console.log('still not loaded');
            //            updateResultPage();
            //        }, 2000);
            //    } else {
            //        var results = w.getResults();
            //        w.updateResults(results);
            //    }
            //};
            //updateResultPage();
        }, error: function (result) {
            console.log(result);
        }
    });
}

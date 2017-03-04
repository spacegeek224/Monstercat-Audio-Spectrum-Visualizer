//Heavily modified from https://github.com/caseif/vis.js/blob/gh-pages/js/analysis/spectrum_algorithms.js


var barWidth = (SpectrumBarCount + Bar1080pSeperation) / SpectrumBarCount - Bar1080pSeperation;
var spectrumDimensionScalar = 4.5
var spectrumMaxExponent = 5
var spectrumMinExponent = 3
var spectrumExponentScale = 2

var SpectrumStart = 4
var SpectrumEnd = 1200
var SpectrumLogScale = 2.55
var marginDecay = 1.6;

var resRatio = (window.innerWidth/window.innerHeight)
var spectrumWidth = 1568 * resRatio;
spectrumSpacing = 7 * resRatio;
var spectrumSize = SpectrumBarCount;
spectrumWidth = (Bar1080pWidth + Bar1080pSeperation) * SpectrumBarCount - Bar1080pSeperation;

var spectrumHeight = 255

function SpectrumEase(Value) {
  return Math.pow(Value, SpectrumLogScale)
}
/*
function GetVisualBins(Array) {
  var SamplePoints = []
  var NewArray = []
  var LastSpot = 0
  for (var i = 0; i < SpectrumBarCount; i++) {
    var Bin = Math.round(SpectrumEase(i / SpectrumBarCount) * (SpectrumEnd - SpectrumStart) + SpectrumStart)
    if (Bin <= LastSpot) {
      Bin = LastSpot + 1
    }
    LastSpot = Bin
    SamplePoints[i] = Bin
  }

  for (var i = 0; i < SpectrumBarCount; i++) {
    var CurSpot = SamplePoints[i]
    var NextSpot = SamplePoints[i + 1]
    if (NextSpot == null) {
      NextSpot = SpectrumEnd
    }

    var CurMax = Array[CurSpot]
    var Dif = NextSpot - CurSpot
    for (var j = 1; j < Dif; j++) {
      CurMax = Math.max(Array[CurSpot + j],CurMax)
    }
    NewArray[i] = CurMax
  }

  UpdateParticleAttributes(NewArray)
  return NewArray
}*/

function GetVisualBins(Array) {
  var SamplePoints = []
  var NewArray = []
  var LastSpot = 0
  for (var i = 0; i < SpectrumBarCount; i++) {
    var Bin = Math.round(SpectrumEase(i / SpectrumBarCount) * (SpectrumEnd - SpectrumStart) + SpectrumStart)
    if (Bin <= LastSpot) {
      Bin = LastSpot + 1
    }
    LastSpot = Bin
    SamplePoints[i] = Bin
  }

  var MaxSamplePoints = []
  for (var i = 0; i < SpectrumBarCount; i++) {
    var CurSpot = SamplePoints[i]
    var NextSpot = SamplePoints[i + 1]
    if (NextSpot == null) {
      NextSpot = SpectrumEnd
    }

    var CurMax = Array[CurSpot]
    var MaxSpot = CurSpot
    var Dif = NextSpot - CurSpot
    for (var j = 1; j < Dif; j++) {
      var NewSpot = CurSpot + j
      if (Array[NewSpot] > CurMax) {
        CurMax = Array[NewSpot]
        MaxSpot = NewSpot
      }
    }
    MaxSamplePoints[i] = MaxSpot
  }

  for (var i = 0; i < SpectrumBarCount; i++) {
    var CurSpot = SamplePoints[i]
    var NextMaxSpot = MaxSamplePoints[i]
    var LastMaxSpot = MaxSamplePoints[i - 1]
    if (LastMaxSpot == null) {
      LastMaxSpot = SpectrumStart
    }
    var LastMax = Array[LastMaxSpot]
    var NextMax = Array[NextMaxSpot]

    NewArray[i] = (LastMax + NextMax)/2
  }

  UpdateParticleAttributes(NewArray)
  return NewArray
}

function TransformToVisualBins(Array) {
  //Array = AverageTransform(Array)
  Array = exponentialTransform(Array)

  return Array;
}

function AverageTransform(Array) {
    var Length = Array.length


    var Values = []
    for (var i = 0; i < Length; i++) {
        var Value = 0
        if (i == 0) {
            Value = Array[i];
        } else {
            var PrevValue = Array[i - 1]
            var NextValue = Array[i + 1]
            var CurValue = Array[i]

            Value = ((CurValue * 4) + ((NextValue + PrevValue)/2 * 2))/6
        }
        Value = Math.min(Value + 1, spectrumHeight)

        Values[i] = Value;
    }

    return Values
/*
    var SamplePoints = []
    for (var i = 0; i < Length; i = i + 2) {
      SamplePoints[SamplePoints.length] = i
    }

    function Interpolate(S,E,A) {
      return S + (E-S)*A
    }

    for (var i = 0; i < SamplePoints.length; i++) {
      var CurSamplePoint = SamplePoints[i]
      var NextSamplePoint = SamplePoints[i + 1]
      if (NextSamplePoint) {
        var Dif = NextSamplePoint - CurSamplePoint
        for (var j = 1; j < Dif; j++) {
          Array[CurSamplePoint + j] = Interpolate(Array[CurSamplePoint],Array[NextSamplePoint],j/Dif)
        }
      }
    }

    return Array*/
}

function exponentialTransform(array) {
    var newArr = [];
    for (var i = 0; i < array.length; i++) {
        var exp = spectrumMaxExponent + (spectrumMinExponent - spectrumMaxExponent) * (i/array.length)
        newArr[i] = Math.max(Math.pow(array[i] / spectrumHeight, exp) * spectrumHeight, 1);
    }
    return newArr;
}

function tailTransform(array) {
	var values = [];
	for (var i = 0; i < spectrumSize; i++) {
		var value = array[i];
		if (i < headMargin) {
			value *= headMarginSlope * Math.pow(i + 1, marginDecay) + minMarginWeight;
		} else if (spectrumSize - i <= tailMargin) {
			value *= tailMarginSlope * Math.pow(spectrumSize - i, marginDecay) + minMarginWeight;
		}
		values[i] = value;
	}
	return values;
}

function exponentialTransform(array) {
	var newArr = [];
	for (var i = 0; i < array.length; i++) {
		var exp = (spectrumMaxExponent - spectrumMinExponent) * (1 - Math.pow(i / spectrumSize, spectrumExponentScale)) + spectrumMinExponent;
		newArr[i] = Math.max(Math.pow(array[i] / spectrumHeight, exp) * spectrumHeight, 1);
	}
	return newArr;
}

// top secret bleeding-edge shit in here
function experimentalTransform(array) {
	var resistance = 3; // magic constant
	var newArr = [];
	for (var i = 0; i < array.length; i++) {
		var sum = 0;
		var divisor = 0;
		for (var j = 0; j < array.length; j++) {
			var dist = Math.abs(i - j);
			var weight = 1 / Math.pow(2, dist);
			if (weight == 1) weight = resistance;
			sum += array[j] * weight;
			divisor += weight;
		}
		newArr[i] = sum / divisor;
	}
	return newArr;
}

function peak1(A, c, d) {
	var m = Math.floor((c + d) / 2);
	if (A[m - 1] <= A[m] && A[m] >= A[m + 1]) {
		return m;
	} else if (A[m - 1] > A[m]) {
		return peak1(A, c, m - 1);
	} else if (A[m] < A[m + 1]) {
		return peak1(A, m + 1, d);
	}
}

function doPeak(array) {
	var newArr = [];
	for (var i in array) {
		newArr[i] = peak1(array, opt.i, opt.j);
	}
	return newArr;
}

function ms(array) {
	var newArr = [];
	var m = math.mean(array);
	var s = math.std(array);
	for (var i in array) {
		if (array[i] - m > 2 * s) {
			newArr[i] = array[i];
		} else {
			newArr[i] = array[i] / 2;
		}
	}
	return newArr;
}

function powTransform(array) {
	var newArr = array.map(v => {
		return Math.pow(v = v / 255, 1 - v) * 255
	});

	return newArr;
}

function normalize(value, max, min, dmax, dmin) {
	return (dmax - dmin) / (max - min) * (value - max) + dmax
}

var base = Math.pow(2, 1 / 3);

function spike(array) {
	var newArr = []
	newArr = array.map(v => {
		var _v = normalize(v, 255, 0, 1, 0);
		return getValFromX(_v, 30, 0);
	});
	return newArr;
}

function compute(x) {
	return base ** x;
}

function getValFromX(x, max, min) {
	return compute(x * (max - min) + min);
}

//Heavily modified from https://github.com/caseif/vis.js/blob/gh-pages/js/analysis/spectrum_algorithms.js


// NOTE: Not all config values may necessarily be changed by the user at
// runtime. Some are mutated internally after intialization, meaning changing
// them at runtime may not work as expected.

/* *********************** */
/* * Audio node settings * */
/* *********************** */
var volumeStep = 0.05; // the step for each volume notch as a fraction of 1

/* *************************** */
/* * Basic spectrum settings * */
/* *************************** */
// BASIC ATTRIBUTES
var spectrumSize = 63; // number of bars in the spectrum
var spectrumDimensionScalar = 4.5; // the ratio of the spectrum width to its height
var spectrumSpacing = 7; // the separation of each spectrum bar in pixels at width=1920
var maxFftSize = 16384; // the preferred fftSize to use for the audio node (actual fftSize may be lower)
var audioDelay = 0.4; // audio will lag behind the rendered spectrum by this amount of time (in seconds)
// BASIC TRANSFORMATION
var spectrumStart = 4; // the first bin rendered in the spectrum
var spectrumEnd = 1200; // the last bin rendered in the spectrum
var spectrumScale = 2.5; // the logarithmic scale to adjust spectrum values to
// EXPONENTIAL TRANSFORMATION
var spectrumMaxExponent = 6; // the max exponent to raise spectrum values to
var spectrumMinExponent = 3; // the min exponent to raise spectrum values to
var spectrumExponentScale = 2; // the scale for spectrum exponents
// DROP SHADOW
var spectrumShadowBlur = 6; // the blur radius of the spectrum's drop shadow
var spectrumShadowOffsetX = 0; // the x-offset of the spectrum's drop shadow
var spectrumShadowOffsetY = 0; // the y-offset of the spectrum's drop shadow

/* ********************** */
/* * Smoothing settings * */
/* ********************** */
var smoothingPoints = 3; // points to use for algorithmic smoothing. Must be an odd number.
var smoothingPasses = 1; // number of smoothing passes to execute
var temporalSmoothing = 0.2; // passed directly to the JS analyzer node

/* ************************************ */
/* * Spectrum margin dropoff settings * */
/* ************************************ */
var headMargin = 7; // the size of the head margin dropoff zone
var tailMargin = 0; // the size of the tail margin dropoff zone
var minMarginWeight = 0.7; // the minimum weight applied to bars in the dropoff zone

/* *************************** */
/* * Basic particle settings * */
/* *************************** */
// COUNT
var baseParticleCount = 2000; // the number of particles at 1080p
var fleckCount = 50; // total fleck count
var bokehCount = 250; // total bokeh count
// OPACITY
var particleOpacity = 0.7; // opacity of primary particles
var bokehOpacity = 0.5; // opacity of bokeh (raising this above 0.5 results in weird behavior)
// SIZE
var minParticleSize = 4; // the minimum size scalar for particle systems
var maxParticleSize = 7; // the maximum size scalar for particle systems
var particleSizeExponent = 2; // the exponent to apply during dynamic particle scaling (similar to spectrum exponents)
// POSITIONING
var yVelRange = 3; // the range for particle y-velocities
var xPosBias = 4.5; // bias for particle x-positions (higher values = more center-biased)
var zPosRange = 450; // the range of z-particles
var zModifier = -250; // the amount to add to z-positions
var zPosBias = 2.3; // bias for particle z-positions (higher values = more far-biased)
var leftChance = 0.88; // the chance for a particle to spawn along the left edge of the screen
var rightChance = 0.03; // the chance for a particle to spawn along the right edge of the screen
var topBottomChance = 0.09; // the chance for a particle to spawn along the top/bottom edges of the screen
// VELOCITY
var velBias = 1.8; // bias for particle velocities (higher values = more center-biased)
var minParticleVelocity = 2; // the minimum scalar for particle velocity
var maxParticleVelocity = 5; // the maximum scalar for particle velocity
var absMinParticleVelocity = 0.001; // the absolute lowest speed for particles
var fleckVelocityScalar = 1.75; // velocity of flecks relative to normal particles
var fleckYVelScalar = 0.75; // y-velocity range of flecks relative to x-velocity
var bokehMinVelocity = maxParticleVelocity * 0.15; // the minimum velocity of bokeh
var bokehMaxVelocity = maxParticleVelocity * 0.3; // the maximum velocity of bokeh

/* ****************************** */
/* * Particle analysis settings * */
/* ****************************** */
var ampLower = 7; // the lower bound for amplitude analysis (inclusive)
var ampUpper = 30; // the upper bound for amplitude analysis (exclusive)
var particleExponent = 4.5; // the power to raise velMult to after initial computation

/* ***************** */
/* * Misc settings * */
/* ***************** */
var cycleSpeed = 4; // the (arbitrary) scalar for cycling rainbow spectrums
var blockWidthRatio = 0.63; // the width of the Monstercat logo relative to its containing block
var blockHeightRatio = 0.73; // the height of the Monstercat logo relative to its containing block
var mouseSleepTime = 1000; // inactivity period in milliseconds before the bottom text is hidden

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

// mostly for debugging purposes
function smooth(array) {
	return savitskyGolaySmooth(array);
}

/**
 * Applies a Savitsky-Golay smoothing algorithm to the given array.
 *
 * See {@link http://www.wire.tu-bs.de/OLDWEB/mameyer/cmr/savgol.pdf} for more
 * info.
 *
 * @param array The array to apply the algorithm to
 *
 * @return The smoothed array
 */
function savitskyGolaySmooth(array) {
	var lastArray = array;
	for (var pass = 0; pass < smoothingPasses; pass++) {
		var sidePoints = Math.floor(smoothingPoints / 2); // our window is centered so this is both nL and nR
		var cn = 1 / (2 * sidePoints + 1); // constant
		var newArr = [];
		for (var i = 0; i < sidePoints; i++) {
			newArr[i] = lastArray[i];
			newArr[lastArray.length - i - 1] = lastArray[lastArray.length - i - 1];
		}
		for (var i = sidePoints; i < lastArray.length - sidePoints; i++) {
			var sum = 0;
			for (var n = -sidePoints; n <= sidePoints; n++) {
				sum += cn * lastArray[i + n] + n;
			}
			newArr[i] = sum;
		}
		lastArray = newArr;
	}
	return newArr;
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
function normalizeAmplitude(array) {
	var values = [];
	for (var i = 0; i < spectrumSize; i++) {
		if (begun) {
			values[i] = array[i] / 255 * spectrumHeight;
		} else {
			value = 1;
		}
	}
	return values;
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

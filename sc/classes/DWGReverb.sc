EarlyRef : MultiOutUGen {

    *ar { arg input = 0, source = [0,0,0], receiver = [0,0,0], roomsize = [1,1,1], hw = 0.2, b = 0.97, n = 0, p = 0, allp_lens = [347,113,37], allp_c = 0.7;
        (source.size == 3).not.if { ^Error("Length of source in EarlyRef must be 3").throw; };
        (receiver.size == 3).not.if { ^Error("Length of receiver in EarlyRef must be 3").throw; };
        (roomsize.size == 3).not.if { ^Error("Length of roomsize in EarlyRef must be 3").throw; };
        (allp_lens.size == 3).not.if { ^Error("Length of allp_lens in EarlyRef must be 3").throw; };
        ^this.multiNewList(['audio', input] ++ source ++ receiver ++ roomsize ++ [hw, b, n, p] ++ allp_lens ++ [allp_c]);
    }

    checkInputs {
        if (inputs.at(0).rate != 'audio', {
            ^("Input at index 0 must be audio rate");
        });
        ^this.checkValidInputs;
    }

    init { |...theInputs|
        inputs = theInputs;
        ^this.initOutputs(2, 'audio')
    }

}

EarlyRefGen : UGen {

    *kr { arg bufL = 0, bufR = 0, source = [0,0,0], receiver = [0,0,0], roomsize = [1,1,1], hw = 0.2, b = 0.97, n = 0, hangle=0 ;
        (source.size == 3).not.if { ^Error("Length of source in EarlyRefGen must be 3").throw; };
        (receiver.size == 3).not.if { ^Error("Length of receiver in EarlyRefGen must be 3").throw; };
        (roomsize.size == 3).not.if { ^Error("Length of roomsize in EarlyRefGen must be 3").throw; };
        ^this.multiNewList(['control', bufL, bufR] ++ source ++ receiver ++ roomsize ++ [hw, b, n, hangle]);
    }

}

EarlyRefAtkGen : UGen {

    *kr { arg bufW = 0, bufX = 0, bufY = 0, bufZ = 0, source = [0,0,0], receiver = [0,0,0], roomsize = [1,1,1], hw = 0.2, b = 0.97, n = 0;
        (source.size == 3).not.if { ^Error("Length of source in EarlyRefAtkGen must be 3").throw; };
        (receiver.size == 3).not.if { ^Error("Length of receiver in EarlyRefAtkGen must be 3").throw; };
        (roomsize.size == 3).not.if { ^Error("Length of roomsize in EarlyRefAtkGen must be 3").throw; };
        ^this.multiNewList(['control', bufW, bufX, bufY, bufZ] ++ source ++ receiver ++ roomsize ++ [hw, b, n]);
    }

}

DWGReverbC1C3 : MultiOutUGen {

    *ar {
        arg
            in = 0,
            len = 2000,
            c1 = 1,
            c3 = 10,
            mix = 1,
            coefs = [1,0.9464,0.87352,0.83,0.8123,0.7398,0.69346,0.6349],
			perm = [1,2,3,4,5,6,7,0],
            doprime = 0;
        (coefs.size == 8).not.if { ^Error("Length of coefs in DWGReverbC1C3 must be 8").throw; };
		(perm.size == 8).not.if { ^Error("Length of perm in DWGReverbC1C3 must be 8").throw; };
        ^this.multiNewList(['audio', in, len, c1, c3, mix, doprime] ++ coefs ++ perm);
    }

    checkInputs {
        if (inputs.at(0).rate != 'audio', {
            ^("Input at index 0 must be audio rate");
        });
        ^this.checkValidInputs;
    }

    init { |...theInputs|
        inputs = theInputs;
        ^this.initOutputs(2, 'audio')
    }

}

DWGReverbC1C3_16 : MultiOutUGen {

    *ar {
        arg
            in = 0,
            len = 2000,
            c1 = 1,
            c3 = 10,
            mix = 1,
            coefs = [1, 0.97498666666667, 0.94997333333333, 0.917248, 0.88323733333333, 0.85901333333333, 0.838704, 0.82528, 0.81702, 0.7978, 0.76396666666667, 0.73362133333333, 0.711996, 0.689556, 0.662228, 0.6349],
			perm = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,0],
            doprime = 0;
        (coefs.size == 16).not.if { ^Error("Length of coefs in DWGReverbC1C3_16 must be 16").throw; };
		(perm.size == 16).not.if { ^Error("Length of perm in DWGReverbC1C3_16 must be 16").throw; };
        ^this.multiNewList(['audio', in, len, c1, c3, mix, doprime] ++ coefs ++ perm);
    }

    checkInputs {
        if (inputs.at(0).rate != 'audio', {
            ^("Input at index 0 must be audio rate");
        });
        ^this.checkValidInputs;
    }

    init { |...theInputs|
        inputs = theInputs;
        ^this.initOutputs(2, 'audio')
    }

}

DWGReverb3Band : MultiOutUGen {

    *ar {
        arg
            in = 0,
            len = 2000,
            xover = 200,
            rtlow = 3,
			rtmid = 2,
            fdamp = 6000,
            coefs = [1,0.9464,0.87352,0.83,0.8123,0.7398,0.69346,0.6349],
			perm = [1,2,3,4,5,6,7,0],
            doprime = 0;
        (coefs.size == 8).not.if { ^Error("Length of coefs in DWGReverb3Band must be 8").throw; };
		(perm.size == 8).not.if { ^Error("Length of perm in DWGReverb3Band must be 8").throw; };
        ^this.multiNewList(['audio', in, len, xover, rtlow, rtmid, fdamp, doprime] ++ coefs ++ perm);
    }

    checkInputs {
        if (inputs.at(0).rate != 'audio', {
            ^("Input at index 0 must be audio rate");
        });
        ^this.checkValidInputs;
    }

    init { |...theInputs|
        inputs = theInputs;
        ^this.initOutputs(2, 'audio')
    }

}

DWGReverb3Band_16 : MultiOutUGen {

    *ar {
        arg
            in = 0,
            len = 2000,
            xover = 200,
            rtlow = 3,
			rtmid = 2,
            fdamp = 6000,
            coefs = [1, 0.97498666666667, 0.94997333333333, 0.917248, 0.88323733333333, 0.85901333333333, 0.838704, 0.82528, 0.81702, 0.7978, 0.76396666666667, 0.73362133333333, 0.711996, 0.689556, 0.662228, 0.6349],
			perm = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,0],
            doprime = 0;
        (coefs.size == 16).not.if { ^Error("Length of coefs in DWGReverb3Band_16 must be 16").throw; };
		(perm.size == 16).not.if { ^Error("Length of perm in DWGReverb3Band_16 must be 16").throw; };
        ^this.multiNewList(['audio', in, len, xover, rtlow, rtmid, fdamp, doprime] ++ coefs ++ perm);
    }

    checkInputs {
        if (inputs.at(0).rate != 'audio', {
            ^("Input at index 0 must be audio rate");
        });
        ^this.checkValidInputs;
    }

    init { |...theInputs|
        inputs = theInputs;
        ^this.initOutputs(2, 'audio')
    }

}

PartConvT : UGen {
	*ar { arg in, fftsize, irbufnum, trig;
		^this.multiNew('audio', in, fftsize, irbufnum, trig);
	}
}
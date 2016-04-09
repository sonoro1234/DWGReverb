EarlyRef27Gen : UGen {

    *kr { arg bufL = 0, bufR = 0, source = [0,0,0], receiver = [0,0,0], roomsize = [1,1,1], hw = 0.2, b = 0.97, n = 0;
        (source.size == 3).not.if { ^Error("Length of source in EarlyRef27Gen must be 3").throw; };
        (receiver.size == 3).not.if { ^Error("Length of receiver in EarlyRef27Gen must be 3").throw; };
        (roomsize.size == 3).not.if { ^Error("Length of roomsize in EarlyRef27Gen must be 3").throw; };
        ^this.multiNewList(['control', bufL, bufR] ++ source ++ receiver ++ roomsize ++ [hw, b, n]);
    }

}

EarlyRefAtkGen : UGen {

    *kr { arg bufW = 0, bufX = 0, bufY = 0, bufZ = 0, source = [0,0,0], receiver = [0,0,0], roomsize = [1,1,1], hw = 0.2, b = 0.97, n = 0;
        (source.size == 3).not.if { ^Error("Length of source in EarlyRef27Gen must be 3").throw; };
        (receiver.size == 3).not.if { ^Error("Length of receiver in EarlyRef27Gen must be 3").throw; };
        (roomsize.size == 3).not.if { ^Error("Length of roomsize in EarlyRef27Gen must be 3").throw; };
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
            doprime = 0;
        (coefs.size == 8).not.if { ^Error("Length of coefs in DWGReverbC1C3 must be 8").throw; };
        ^this.multiNewList(['audio', in, len, c1, c3, mix, doprime] ++ coefs);
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
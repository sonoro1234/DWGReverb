--Generates Room impulse response (RIR) with the source image method
--and outputs trigger when buffer changes
--this is the ambisonics version (slightly less efficient)
--bW,bX,bY,bZ: W,X,Y,Z buffer numbers
--L : room sizes
--Ps: source coordinates
--Pr: receiver coordinates
--HW: ear separation
--B: reflection coef (-1,1)
--N: value 0 generates 27 mirror images, 1 to 5 generates 8*(N*2+1)^3 mirror images

EarlyRefAtkGen = UGen:new{name='EarlyRefAtkGen'}
function EarlyRefAtkGen.kr(bW,bX,bY,bZ,Ps,Pr,L,HW,B,N)
	bW = bW or 0;bX = bX or 0;bY = bY or 0;bZ = bZ or 0;
	Ps = Ps or {0,0,0};Pr = Pr or {0,0,0};L = L or {1,1,1};HW = HW or 0.2;B= B or 0.97;N = N or 0
	return EarlyRefAtkGen:MultiNew(concatTables({1,bW,bX,bY,bZ},Ps,Pr,L,{HW,B,N}))
end

--late reverb stereo FDN8
--inp: input
--len: max. fdn length
--c1: inverse of decay time
--c3: inverse of hight frequency decay time
--mix
DWGReverb = MultiOutUGen:new{name="DWGReverb"}
function DWGReverb.ar(inp,len,c1,c3,mix,coefs,doprime)
	inp = inp or 0;c1 = c1 or 1;c3 = c3 or 10;len = len or 2000;mix = mix or 1
	coefs = coefs or {1,0.9464,0.87352,0.83,0.8123,0.7398,0.69346,0.6349}
	doprime = doprime or 0
	assert(#coefs==8)
	return DWGReverb:MultiNew{2,2,inp,len,c1,c3,mix,doprime,unpack(coefs)}
end

--The synthdef

--CIPIC kernel
require"sc.atk"
dec = FoaDecoderKernel.newCIPIC()
--buffer must be power of two for Convolution2
--choose size according to generated RIR
Nbuf = 2048*4
SynthDef("testearly",{out=Master.busin,len=2000,L=Ref{6.2,8,2.7},Ps = Ref{2,3.1,1.2},Pr = Ref{3,3,1.2},B=0.74,HW=0.4,N = 3,revl=1,c1=8,c3=10,bypass=0},function()

	local input = Impulse.ar(2)
	input = HPF.ar(LPF.ar(input,3000),200) *3
	input = LeakDC.ar(input)

	--for moving source
	Ps[1] = SinOsc.kr(0.1):range(2,4)

	local bw = LocalBuf(Nbuf)
	local bx = LocalBuf(Nbuf)
	local by = LocalBuf(Nbuf)
	local bz = LocalBuf(Nbuf)
	local trig = EarlyRefAtkGen.kr(bw,bx,by,bz,Psmod,Pr,L,HW,-B,N)
	local sigw = Convolution2L.ar(input,bw,trig,Nbuf)
	local sigx = Convolution2L.ar(input,bx,trig,Nbuf)
	local sigy = Convolution2L.ar(input,by,trig,Nbuf)
	local sigz = Convolution2L.ar(input,bz,trig,Nbuf)
	local ambis = {sigw,sigx,sigy,sigz}
	local early = AtkKernelConv.ar(ambis,dec.kernel)
	local sig = DWGReverb.ar(Mix(early)*0.5,len,c1,c3)*revl + early  
	sig = Select.ar(bypass,{sig,input})
	Out.ar(out,sig);
end):guiplay()

Master.inserts = {{"to_mono",{bypass=1}}}

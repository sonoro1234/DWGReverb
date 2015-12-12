EarlyRef27Gen = UGen:new{name='EarlyRef27Gen'}
function EarlyRef27Gen.kr(bufL,bufR,Ps,Pr,L,HW,B,N)
	bufL = bufL or 0;bufR = bufR or 0
	Ps = Ps or {0,0,0};Pr = Pr or {0,0,0};L = L or {1,1,1};HW = HW or 0.2;B= B or 0.97;N = N or 0
	return EarlyRef27Gen:MultiNew(concatTables({1,bufL,bufR},Ps,Pr,L,{HW,B,N}))
end
DWGReverb = MultiOutUGen:new{name="DWGReverb"}
function DWGReverb.ar(inp,len,c1,c3,mix,coefs,doprime)
	inp = inp or 0;c1 = c1 or 1;c3 = c3 or 1;len = len or 32000;mix = mix or 1
	coefs = coefs or {1,0.9464,0.87352,0.83,0.8123,0.7398,0.69346,0.6349}
	doprime = doprime or 0
	assert(#coefs==8)
	return DWGReverb:MultiNew{2,2,inp,len,c1,c3,mix,doprime,unpack(coefs)}
end

Nbuf = 2048*4

SynthDef("testearly",{out=Master.busin,len=2000,L=Ref{6.2,8,2.7},Ps = Ref{2,3.1,1.2},Pr = Ref{3,3,1.2},B=0.74,HW=0.4,N = 3,revl=0,c1=8,c3=10,bypass=0},function()

	local input = Impulse.ar(2)
	input = HPF.ar(LPF.ar(input,3000),200) *3
	input = LeakDC.ar(input)

	--for moving source
	Ps[1] = SinOsc.kr(0.1):range(2,4)

	local bL = LocalBuf(Nbuf)
	local bR = LocalBuf(Nbuf)
	local trig = EarlyRef27Gen.kr(bL,bR,Ps,Pr,L,HW,-B,N)
	local sigL = Convolution2.ar(input,bL,trig,Nbuf)
	local sigR = Convolution2.ar(input,bR,trig,Nbuf)
	local early = {sigL,sigR}
	local sig = DWGReverb.ar(Mix(early)*0.5,len,c1,c3)*revl + early  
	sig = Select.ar(bypass,{sig,input})
	Out.ar(out,sig);
end):guiplay()

Master.inserts = {{"to_mono",{bypass=1}}}
--DiskOutBuffer("earlyref27.wav")
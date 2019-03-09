EarlyRef27Gen = UGen:new{name='EarlyRef27Gen'}
function EarlyRef27Gen.kr(bufL,bufR,Ps,Pr,L,HW,B,N,Hangle)
	bufL = bufL or 0;bufR = bufR or 0
	Ps = Ps or {0,0,0};Pr = Pr or {0,0,0};L = L or {1,1,1};HW = HW or 0.2;B= B or 0.97;N = N or 0; Hangle = Hangle or 0
	return EarlyRef27Gen:MultiNew(concatTables({1,bufL,bufR},Ps,Pr,L,{HW,B,N,Hangle}))
end


function SetBuf.create(...)
	local   buf, values, offset   = assign({ 'buf', 'values', 'offset' },{ nil, nil, 0 },...)
	return SetBuf:MultiNew{0,buf,offset,#values,unpack(values)}
end
PartConvT = UGen:new{name='PartConvT'}
function PartConvT.ar(inp,fftsize,irbufnum,trig)
	assert(inp)
	trig = trig or 1
	return PartConvT:MultiNew{2,inp,fftsize,irbufnum,trig}
end

Nbuf = 1024*16*4
--function initCb()
print"zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz"
sdef = SynthDef("testearly",{out=Master.busin,len=3500,L=Ref{32,45.4,7.4},Ps = Ref{8,8,1.2},Pr = Ref{15,10,1.2},B=0.7,HW=0.2,dist=5,angle=0,N = 6,revl=0,c1=2.3,c3=5,bypass=0,size=1,Hangle=0},function()
	--local input=Pluck.ar(PinkNoise.ar(),Impulse.ar(1),0.1,0.005)
	L = L*size
	Ps = Ps*size
	Pr = Pr*size
	local input = WhiteNoise.ar()*Dust.ar(10) --Impulse.ar(1);
    --input = HPF.ar(LPF.ar(input, 8000), 200) * 3;
	input = LPF.ar(input, 8000) * 3;
	--local input = Impulse.ar(4)
	input = LeakDC.ar(input)

	--for moving source
	--Ps[1] = SinOsc.kr(0.1):range(2,4)
	--local angle = LFSaw.kr(0.1):range(-math.pi,math.pi) --SinOsc.kr(0.1):range(-math.pi,math.pi)
	--local dist = 1
	
angle = angle*math.pi
	local Psmod = {Pr[1] + angle:sin()*dist,Pr[2] + angle:cos()*dist,Ps[3]}

	local bL = LocalBuf(Nbuf)
	local bR = LocalBuf(Nbuf)

	local trig = EarlyRef27Gen.kr(bL,bR,Psmod,Pr,L,HW,-B,N,Hangle*math.pi)

	local sigL = PartConvT.ar(input,1024*4,bL,trig)
	local sigR = PartConvT.ar(input,1024*4,bR,trig)
	local early = {sigL,sigR}*dist

	local sig = DWGReverb.ar(Mix(early)*0.5,len,c1,c3)*revl + early  --+ input
	sig = Select.ar(bypass,{sig,input})

	Out.ar(out,sig);
end)
--sdef:dumpInputs() 
sdef:guiplay()

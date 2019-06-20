

Nbuf = 1024*16*4
--function initCb()
print"zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz"
sdef = SynthDef("testearly",{out=Master.busin,len=3500,L=Ref{6.2,8,2.7},Ps = Ref{2,3.1,1.2},Pr = Ref{3,3,1.2},B=0.7,HW=0.4,dist=1.5,angle=0,N = 6,revl=0.5,c1=8,c3=10,bypass=0,size=1,Hangle=0},function()
	--local input=Pluck.ar(PinkNoise.ar(),Impulse.ar(1),0.1,0.005)
	L = L*size
	Ps = Ps*size
	Pr = Pr*size

	local input = Impulse.ar(2)
	input = HPF.ar(LPF.ar(input,3000),200) *3
	input = LeakDC.ar(input)

	--for moving source
	Ps[1] = SinOsc.kr(0.1):range(2,4)
	local angle = LFSaw.kr(0.1):range(-math.pi,math.pi) --SinOsc.kr(0.1):range(-math.pi,math.pi)
	
	local Psmod = Ps--{Pr[1] + angle:sin()*dist,Pr[2] + angle:cos()*dist,Ps[3]}

	local bL = LocalBuf(Nbuf)
	local bR = LocalBuf(Nbuf)

	local trig = EarlyRefGen.kr(bL,bR,Psmod,Pr,L,HW,-B,N,Hangle*math.pi)

	local sigL = PartConvT.ar(input,1024*4,bL,trig)
	local sigR = PartConvT.ar(input,1024*4,bR,trig)
	local early = {sigL,sigR}*dist

	local sig = DWGReverb.ar(Mix(early)*0.5,len,c1,c3)*revl + early  --+ input
	sig = Select.ar(bypass,{sig,input})

	Out.ar(out,sig);
end)

sdef:guiplay()

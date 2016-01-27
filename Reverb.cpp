/*
 by Victor Bombi
*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "SC_PlugIn.h"
InterfaceTable *ft;

///////////////////////////////////////////////////////
#include "dwglib/DWG.cpp"
/*
struct FilterC1C3 : public LTITv<1,1>
{
	float freq;
	float c1;
	float c3;
	FilterC1C3(){c1 = 0;c3 = 0; freq = 0;};
	void setcoeffs(float freq,float c1,float c3){
		if(this->freq != freq || this->c1 != c1 || this->c3 != c3){
			float g = 1.0 - c1/freq;
			float b = 4.0*c3+freq;
			float a1 = (-b+sqrt(b*b-16.0*c3*c3))/(4.0*c3);
			KernelB[0] = g*(1+a1);
			KernelA[0] = a1;
			this->freq = freq;
			this->c1 = c1;
			this->c3 = c3;
		}
	}
};
struct Biquad : public LTITv<3,2>
{
	enum biquadtype {
		pass = 0,
		low,
		high,
		notch
	};
	void setcoeffs(float f0, float fs, float Q, int type){
		float a = 1/(2*tan(M_PI*f0/fs));
		float a2 = a*a;
		float aoQ = a/Q;
		float d = (4*a2+2*aoQ+1);

		KernelA[0] = -(8*a2-2) / d;
		KernelA[1] = (4*a2 - 2*aoQ + 1) / d;

		switch(type) {
		case pass:
			KernelB[0] = 2*aoQ/d;
			KernelB[1] = 0;
			KernelB[2] = -2*aoQ/d;
			break;
		case low:
			KernelB[0] = 1/d;
			KernelB[1] = 2/d;
			KernelB[2] = 1/d;
			break;
		case high:
			KernelB[0] = 4*a2/d;
			KernelB[1] = -8*a2/d;
			KernelB[2] = 4*a2/d;
			break;
		case notch:
			KernelB[0] = (1+4*a2)/d;
			KernelB[1] = (2-8*a2)/d;
			KernelB[2] = (1+4*a2)/d;
			break;
		}
	}
};
*/
struct DWGReverb : public Unit
{
	FilterC1C3 decay[8];
	//LagrangeT<1024*128> delay[8];
	CircularBuffer2POWSizedT<1024*32> delay[8];
	void setcoeffs(float c1, float c3, float mix, float Fs);
	void setlengths(float len);
	void setlengths2(float len);
	void go(float *in,float *out[2],int N) ;
	float mix;
	float len;
	float A[8][8];
	float o[8];
	float b[8];
	float c[8];
	float cR[8];
	float factors[8];
	float lengths[8] = {37,87,181,271,359,593,688,721};//{37,87,181,271,359,592,687,721};
	float doprime;
	DWGReverb(DWGReverb *unit);
};

extern "C"
{
	void DWGReverb_next(DWGReverb *unit, int inNumSamples);
	void DWGReverb_Ctor(DWGReverb* unit);
	void DWGReverb_Dtor(DWGReverb* unit);
}
void DWGReverb_Ctor(DWGReverb* unit)
{
	new(unit) DWGReverb(unit);
}
void DWGReverb_Dtor(DWGReverb* unit)
{
	unit->~DWGReverb();
}
void DWGReverb_next(DWGReverb  *unit, int inNumSamples)
{
	float *out[2];
	out[0] = OUT(0);
	out[1] = OUT(1);
	float *in = IN(0);
	float len = ZIN0(1);
	float c1 = ZIN0(2);
	float c3 = ZIN0(3);
	float mix = ZIN0(4);
	unit->setlengths(len);
	unit->setcoeffs(c1,c3,mix,SAMPLERATE);
	//for(int i = 0; i < inNumSamples; i++)
	//	out[i] = unit->reverb(in[i]);
	unit->go(in,out,inNumSamples);
}

DWGReverb::DWGReverb(DWGReverb *unit){
	float a = -1.0/4.0;
	float aa[8] = {a,1+a,a,a,a,a,a,a};
	float tcR[8] = {-1,1,1,-1,-1,1,1,-1};
	for(int k=0;k<8;k++){
			o[k] = 0;
			b[k] = 1;
			c[k] = k<8?((k%2==0)?1.0/8.0:-1.0/8.0):0.0;
			cR[k] = tcR[k]*1.0/8.0;
	}
	
	for(int j=0;j<8;j++)
		for(int k=0;k<8;k++)
			A[j][k] = aa[(8+(k-j))%8];
			
	this->len = 0;
	float len = ZIN0(1);
	float c1 = ZIN0(2);
	float c3 = ZIN0(3);
	float mix = ZIN0(4);
	doprime = IN0(5);
	for(int i=0;i<8;i++)
		this->factors[i] = IN0(6+i);
	
	unit->setlengths(len);
	unit->setcoeffs(c1,c3,mix,SAMPLERATE);
	SETCALC(DWGReverb_next);
};
bool isprime(int n)
{
  unsigned int i;
  const unsigned int lim = (int)sqrtf((float)n);

  if (n == 2) return(true);
  if ((n & 1) == 0) return(false);
  for(i = 3; i <= lim; i += 2)
    if ((n % i) == 0) return(false);
  return(true);
}

int nearestprime(int n, float rerror)
{
  int bound,k;

  if (isprime(n)) return(n);
  /* assume n is large enough and n*rerror enough smaller than n */
  bound = (int)(n * rerror);
  for(k = 1; k <= bound; k++) {
    if (isprime(n+k)) return(n+k);
    if (isprime(n-k)) return(n-k);
  }
  return(-1);
}
void DWGReverb::setlengths(float len){
	//float factors[8] = {1,0.9464,0.87352,0.83,0.8123,0.7398,0.69346,0.6349};
	if(this->len == len)
		return;
	this->len = len;
	float lent;
	for(int i = 0;i < 8 ;i++){
		if(doprime==1)
			lent = nearestprime(len * factors[i],0.5);
		else
			lent = len * factors[i];
		if(lent == -1){
			Print("DWGReverb error nearestprime\n");
			lent = len * factors[i];
		}else{
			//Print("nearest %g",lent);
		}
		//Print("\n");
		lengths[i] = lent;
	}
	//Print("\n");
}
void DWGReverb::setlengths2(float len){
	int primes[8] = {19,17,13,11,7,5,3,2};
	//float factors[8] = {1,0.9464,0.87352,0.83,0.8123,0.7398,0.69346,0.6349};
	if(this->len == len)
		return;
	this->len = len;
	float lent;
	int lent2;
	for(int i = 0;i < 8 ;i++){
		lent = len * factors[i];
		int mi = floor(0.5 + (log(lent)/log(primes[i])));
		lent2 = pow(primes[i],mi);
		Print("nearest %g %d %d\n",lent,lent2,mi);
		lengths[i] = lent2;
	}
	Print("\n");
}

void DWGReverb::setcoeffs(float c1, float c3, float mix, float Fs){

	this->mix = mix;
	for(int k=0;k<8;k++) {
		decay[k].setcoeffs(Fs/lengths[k],c1,c3);
	}

}
//Householder Feedback Matrix
void DWGReverb :: go(float *inA,float *outA[2],int inNumSamples)
{
  float i[8];
/*
  for(int j=0;j<8;j++) {
    i[j] = in; //b[j] * in;
    for(int k=0;k<8;k++) {
      i[j] += A[j][k] * o[k];
    }
  }
*/
	for(int k=0; k< inNumSamples ;k++){
		float sumo = 0;
		for(int j=0;j<8;j++)
			sumo += o[j];
		sumo *= 0.25;
		sumo -= inA[k];
		for(int j=0;j<8;j++)
			i[j] = o[j] - sumo;
			
		float out[2];out[0]=0;out[1]=0;
		//float out = 0;
		for(int j=0;j<8;j++) {
			delay[j].push(i[j]);
			o[j] = decay[j].filter(delay[j].delay(lengths[j]-1));
            //kill_denormals(o[j]);
			//out += c[j] * o[j];//*.5;
			//out[j%2] += c[j] * o[j];
			out[0] += c[j] * o[j];
			out[1] += cR[j] * o[j];
		}
		//out[1] = out[0];
		
		outA[0][k] =  mix*out[0] + (1.0-mix)*inA[k];
		outA[1][k] =  mix*out[1] + (1.0-mix)*inA[k];
	}
}
////////////////////////////////////////////EarlyRef

float dist(float im[3],float r[3]){
    float x = im[0] - r[0];
    float y = im[1] - r[1];
    float z = im[2] - r[2];
    return sqrt(x*x+y*y+z*z);
}

const int MaxNits = 2;
//const int Nits = MaxNits;
//const int Nrefs = pow((2*Nits + 1),3)*8;
const int MaxNrefs = pow((2*MaxNits + 1),3);
//const int Nrefs = MaxNrefs;
struct EarlyRef:public Unit
{
	EarlyRef(Unit* unit);
    float CalcOne(int n,float exp,float ux,float uy,float uz,float lx,float ly,float lz);
    void refsCalculation();
    void filters_init();
    void filters_tick(float in);
    void filters_tickLR(float *inL,float *inR,float &outL,float &outR);
    void findImage(float ufx,float ufy,float ufz,float lfx,float lfy,float lfz,float * res);
    void getargs(Unit * unit);
    CircularBuffer2POWSizedT<4096*16> MtapDel;
    CircularBuffer2POWSizedT<4096*16> DelR;
    CircularBuffer2POWSizedT<4096*16> DelL;
    LTITv<1,1> filters[4];
    LTITv<1,1> filtersL[4];
    LTITv<1,1> filtersR[4];
    float filters_out[5];
	float L[3];
    float Ps[3];
    float Pr[3];
    float Ps_[3];//center room is 0,0
    float Pr_[3];
    float B,HW,d0,p;
    int N,Nrefs;
    float samprate;
    float delL[MaxNrefs],delR[MaxNrefs],ampL[MaxNrefs],ampR[MaxNrefs];
    int RefN[MaxNrefs];
};
SCWrapClass(EarlyRef);
EarlyRef::EarlyRef(Unit *unit)
{
    Ps[0] = ZIN0(1);
    Ps[1] = ZIN0(2);
    Ps[2] = ZIN0(3);
    Pr[0] = ZIN0(4);
    Pr[1] = ZIN0(5);
    Pr[2] = ZIN0(6);
    L[0] = ZIN0(7);
    L[1] = ZIN0(8);
    L[2] = ZIN0(9);
    HW = ZIN0(10);
    B = ZIN0(11);
    N = sc_clip(ZIN0(12),0.0,(float)MaxNits);
    Nrefs = pow((2*N + 1),3);
    p = sc_clip(ZIN0(13),0.0,0.9);
    samprate = SAMPLERATE;
    //printf("%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n",Ps[0],Ps[1],Ps[2],Pr[0],Pr[1],Pr[2],L[0],L[1],L[2],HW,B);
    filters_init();
    refsCalculation();
    SETCALC(EarlyRef_next);
}
void EarlyRef::filters_init(){
    printf("filter_init %g\n",p);
    for(int i= 0;i<4;i++){
        filtersL[i].KernelA = filtersR[i].KernelA = filters[i].KernelA = -p;
        filtersL[i].KernelB = filtersR[i].KernelB = filters[i].KernelB = 1.0 - p;
    }
}
void EarlyRef::filters_tick(float in)
{
    filters_out[0] = in;
    for(int i= 0;i<4;i++){
        filters_out[i+1] = filters[i].filter(filters_out[i]);
    }
}
//filtering for first version
void EarlyRef::filters_tickLR(float *inL,float *inR,float &outL,float &outR)
{
    outL = outR = 0.0;
    for(int i= 0;i<4;i++){
        outL = filtersL[i].filter(outL + inL[4-i]);
        outR = filtersR[i].filter(outR + inR[4-i]);
    }
    outL += inL[0];
    outR += inR[0];
}
float EarlyRef::CalcOne(int n,float exp,float ux,float uy,float uz,float lx,float ly,float lz)
{
    float image[3];
    findImage(ux,uy,uz,lx,ly,lz,image);
    float rec[3];
    rec[0] = Pr_[0] - HW;
    rec[1] = Pr_[1];
    rec[2] = Pr_[2];
    float distL = dist(image,rec);
    rec[0] = Pr_[0] + HW;
    float distR = dist(image,rec);
    float preA = pow(B,exp);
    ampL[n] = preA/(distL + 0.001);//d0*preA/(distL);
    ampR[n] = preA/(distR + 0.001);//d0*preA/(distR);
    delL[n] = samprate*distL/340.0;
    delR[n] = samprate*distR/340.0;
    RefN[n] = sc_min(exp,4);//only have until reflection 4 filtered
    //printf("%d, %g, %g, %g, %g, %g, %g, %g\n",n,exp,ux,uy,uz,lx,ly,lz);
    return distL;
}
void EarlyRef::refsCalculation(){
    /*
    for(int i=-10;i<=10;i++){
        int a = i%2;
        int b = i & -1;
        printf("i :%d a is %d, b is %d\n",i,a,b);
    }
    */
    //translate coordinates to origin in center room
    Ps_[0] = Ps[0] - L[0]*0.5;
    Ps_[1] = Ps[1] - L[1]*0.5;
    Ps_[2] = Ps[2] - L[2]*0.5;
    Pr_[0] = Pr[0] - L[0]*0.5;
    Pr_[1] = Pr[1] - L[1]*0.5;
    Pr_[2] = Pr[2] - L[2]*0.5;
    unsigned long num = 0;
    int u,v,w;
    int absl,absm,absn;
    for(int l = -N;l <=N;l++){
            //u = 1.0 - 2.0*(l%2);
            u = 1.0 - 2.0*(l & 1);
            absl = abs(l);
            for(int m = -N;m <=N;m++){
                //v = 1.0 - 2.0*(m%2);
                v = 1.0 - 2.0*(m & 1);
                absm = abs(m);
                for(int n = -N;n <=N;n++){
                    //w = 1.0 - 2.0*(n%2);
                    w = 1.0 - 2.0*(n & 1);
                    absn = abs(n);
                    float exp = absl + absm + absn;
                    CalcOne(num,exp,u,v,w,l,m,n);
                    num ++;
                    }
            }
    }
    printf("num is %d MaNrefs is %d Nrefs is %d\n",num,MaxNrefs,Nrefs);
}
void EarlyRef::findImage(float ufx,float ufy,float ufz,float lfx,float lfy,float lfz,float * res){
    res[0] = ufx*Ps_[0] + lfx*L[0];
    res[1] = ufy*Ps_[1] + lfy*L[1];
    res[2] = ufz*Ps_[2] + lfz*L[2];
}
void EarlyRef::getargs(Unit *unit){
    float Pst[3],Prt[3],Lt[3],HWt,Bt,pt;
    int Nt;
    Pst[0] = ZIN0(1);
    Pst[1] = ZIN0(2);
    Pst[2] = ZIN0(3);
    Prt[0] = ZIN0(4);
    Prt[1] = ZIN0(5);
    Prt[2] = ZIN0(6);
    Lt[0] = ZIN0(7);
    Lt[1] = ZIN0(8);
    Lt[2] = ZIN0(9);
    HWt = ZIN0(10);
    Bt = ZIN0(11);
    Nt = sc_clip(ZIN0(12),0.0,(float)MaxNits);

    bool changed = (Pst[0] != Ps[0] || Pst[1] != Ps[1] || Pst[2] != Ps[2] || Prt[0] != Pr[0] || Prt[1] != Pr[1] || Prt[2] != Pr[2] || Lt[0] != L[0] || Lt[1] != L[1] || Lt[2] != L[2] || HWt != HW || Bt != B || Nt != N);
    if (changed) {
        Ps[0] = Pst[0]; Ps[1] = Pst[1]; Ps[2] = Pst[2]; Pr[0] = Prt[0]; Pr[1] = Prt[1];Pr[2] = Prt[2]; L[0] = Lt[0]; L[1] = Lt[1] ; L[2] = Lt[2] ; HW = HWt ; B = Bt;N = Nt;
        //printf("%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n",Ps[0],Ps[1],Ps[2],Pr[0],Pr[1],Pr[2],L[0],L[1],L[2],HW,B);
        Nrefs = pow((2*N + 1),3);
        refsCalculation();
    }
    
    pt = sc_clip(ZIN0(13),0.0,0.9);
    if(p != pt){
        p = pt;
        filters_init();
    }
}
void EarlyRef_next(EarlyRef* unit,int inNumSamples)
{
    unit->getargs(unit);
    float * outL = OUT(0);
    float * outR = OUT(1);
    float * in = IN(0);
    float * delL = unit->delL;
    float * delR = unit->delR;
    float * ampL = unit->ampL;
    float * ampR = unit->ampR;
    int Nrefs = unit->Nrefs;
    /*/version 1 quicker
    for(int i=0;i<inNumSamples; i++){
        unit->MtapDel.push(in[i]);
        float oL = 0.0;
        float oR = 0.0;
        for(int j=0;j<Nrefs;j++){
            oL += unit->MtapDel.delay(unit->delL[j])*unit->ampL[j];
            oR += unit->MtapDel.delay(unit->delR[j])*unit->ampR[j];
        }
        outL[i] = oL;
        outR[i] = oR;
    }
    //*/
    //*/version 1 filtered
    float oL[5];
    float oR[5];
    for(int i=0;i<inNumSamples; i++){
        unit->MtapDel.push(in[i]);

        //for(int j=0;j<5;j++)
            //oR[j] = oL[j] = 0.0;
        memset(oR,0,sizeof(float)*5);
        memset(oL,0,sizeof(float)*5);
    
        for(int j=0;j<Nrefs;j++){
            oL[unit->RefN[j]] += unit->MtapDel.delay(unit->delL[j])*unit->ampL[j];
            oR[unit->RefN[j]] += unit->MtapDel.delay(unit->delR[j])*unit->ampR[j];
        }
        
        unit->filters_tickLR(oL,oR,outL[i],outR[i]);
        //outL[i] = oL;
        //outR[i] = oR;
    }
    //*/
    /*// version 2:Mas lento no se bien por que!!
    CircularBuffer2POWSizedT<4096*16> *DelR = &(unit->DelR);
    CircularBuffer2POWSizedT<4096*16> *DelL = &(unit->DelL);
    float *BufferL = DelL->Buffer;
    float *BufferR = DelR->Buffer;
    for(int i=0;i<inNumSamples; i++){
        for(int j=0;j<Nrefs;j++){
            DelL->add(in[i]*ampL[j], -(int)delL[j]);//neg is future
            DelR->add(in[i]*ampR[j], -(int)delR[j]);
            //BufferL[(DelL->pointer - (int)delL[j]) & DelL->mask] += in[i]*ampL[j];
            //BufferR[(DelR->pointer - (int)delR[j]) & DelR->mask] += in[i]*ampR[j];
        }
        outL[i] = unit->DelL.get_tick();
        outR[i] = unit->DelR.get_tick();        
        
    }
    //*/
    /*// version 2 filtered Mas lento no se bien por que!!
    CircularBuffer2POWSizedT<4096*16> *DelR = &(unit->DelR);
    CircularBuffer2POWSizedT<4096*16> *DelL = &(unit->DelL);
    float *BufferL = DelL->Buffer;
    float *BufferR = DelR->Buffer;
    for(int i=0;i<inNumSamples; i++){
        unit->filters_tick(in[i]);
        for(int j=0;j<Nrefs;j++){
            float filt_in = unit->filters_out[unit->RefN[j]];
            DelL->add(filt_in*ampL[j], -(int)delL[j]);//neg is future
            DelR->add(filt_in*ampR[j], -(int)delR[j]);
            //BufferL[(DelL->pointer - (int)delL[j]) & DelL->mask] += in[i]*ampL[j];
            //BufferR[(DelR->pointer - (int)delR[j]) & DelR->mask] += in[i]*ampR[j];
        }
        outL[i] = unit->DelL.get_tick();
        outR[i] = unit->DelR.get_tick();        
        
    }
    //*/
    
}
/////////////////////////////////////////////////
struct EarlyRef27:public Unit
{
	EarlyRef27(Unit* unit);
    float CalcOne(int n,float exp,float ux,float uy,float uz,float lx,float ly,float lz);
    void refsCalculation();
    void findImage(float ufx,float ufy,float ufz,float lfx,float lfy,float lfz,float * res);
    void getargs(Unit * unit);
    LagrangeMtapT<4096*16,27*2> LagDel;
    //CircularBuffer2POWSizedT<4096*16> LagDel;
	float L[3];
    float Ps[3];
    float Pr[3];
    float B,HW,d0;
    float samprate;
    float delL[27],delR[27],ampL[27],ampR[27];
};
SCWrapClass(EarlyRef27);
EarlyRef27::EarlyRef27(Unit *unit)
{
    Ps[0] = ZIN0(1);
    Ps[1] = ZIN0(2);
    Ps[2] = ZIN0(3);
    Pr[0] = ZIN0(4);
    Pr[1] = ZIN0(5);
    Pr[2] = ZIN0(6);
    L[0] = ZIN0(7);
    L[1] = ZIN0(8);
    L[2] = ZIN0(9);
    HW = ZIN0(10);
    B = ZIN0(11);
    samprate = SAMPLERATE;
    //printf("%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n",Ps[0],Ps[1],Ps[2],Pr[0],Pr[1],Pr[2],L[0],L[1],L[2],HW,B);
    refsCalculation();
    SETCALC(EarlyRef27_next);
}
float EarlyRef27::CalcOne(int n,float exp,float ux,float uy,float uz,float lx,float ly,float lz)
{
    float image[3];
    findImage(ux,uy,uz,lx,ly,lz,image);
    float rec[3];
    rec[0] = Pr[0] - HW;
    rec[1] = Pr[1];
    rec[2] = Pr[2];
    float distL = dist(image,rec);
    rec[0] = Pr[0] + HW;
    float distR = dist(image,rec);
    float preA = pow(B,exp);
    ampL[n] = preA/(distL + 0.001);//d0*preA/(distL);
    ampR[n] = preA/(distR + 0.001);//d0*preA/(distR);
    delL[n] = samprate*distL/340.0;
    delR[n] = samprate*distR/340.0;
    //printf("%g, %g, %g, %g, %g\n",d0,delL[n],delR[n],ampL[n],ampR[n]);
    return distL;
}
void EarlyRef27::refsCalculation(){
    //printf("calc\n");
    d0 = 1.0;
   d0 = CalcOne(0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0);	
CalcOne(1, 1.0, 1.0, 1.0,-1.0, 0.0, 0.0, 0.0);	
CalcOne(2, 1.0, 1.0,-1.0, 1.0, 0.0, 0.0, 0.0);	
CalcOne(3, 2.0, 1.0,-1.0,-1.0, 0.0, 0.0, 0.0);	
CalcOne(4, 1.0,-1.0, 1.0, 1.0, 0.0, 0.0, 0.0);	
CalcOne(5, 2.0,-1.0, 1.0,-1.0, 0.0, 0.0, 0.0);	
CalcOne(6, 2.0,-1.0,-1.0, 1.0, 0.0, 0.0, 0.0);	
CalcOne(7, 3.0,-1.0,-1.0,-1.0, 0.0, 0.0, 0.0);	
CalcOne(8, 1.0, 1.0, 1.0,-1.0, 0.0, 0.0, 2.0);	
CalcOne(9, 2.0, 1.0,-1.0,-1.0, 0.0, 0.0, 2.0);	
CalcOne(10, 2.0,-1.0, 1.0,-1.0, 0.0, 0.0, 2.0);	
CalcOne(11, 3.0,-1.0,-1.0,-1.0, 0.0, 0.0, 2.0);	
CalcOne(12, 1.0, 1.0,-1.0, 1.0, 0.0, 2.0, 0.0);	
CalcOne(13, 2.0, 1.0,-1.0,-1.0, 0.0, 2.0, 0.0);	
CalcOne(14, 2.0,-1.0,-1.0, 1.0, 0.0, 2.0, 0.0);	
CalcOne(15, 3.0,-1.0,-1.0,-1.0, 0.0, 2.0, 0.0);	
CalcOne(16, 2.0, 1.0,-1.0,-1.0, 0.0, 2.0, 2.0);	
CalcOne(17, 3.0,-1.0,-1.0,-1.0, 0.0, 2.0, 2.0);	
CalcOne(18, 1.0,-1.0, 1.0, 1.0, 2.0, 0.0, 0.0);	
CalcOne(19, 2.0,-1.0, 1.0,-1.0, 2.0, 0.0, 0.0);	
CalcOne(20, 2.0,-1.0,-1.0, 1.0, 2.0, 0.0, 0.0);	
CalcOne(21, 3.0,-1.0,-1.0,-1.0, 2.0, 0.0, 0.0);	
CalcOne(22, 2.0,-1.0, 1.0,-1.0, 2.0, 0.0, 2.0);	
CalcOne(23, 3.0,-1.0,-1.0,-1.0, 2.0, 0.0, 2.0);	
CalcOne(24, 2.0,-1.0,-1.0, 1.0, 2.0, 2.0, 0.0);	
CalcOne(25, 3.0,-1.0,-1.0,-1.0, 2.0, 2.0, 0.0);	
CalcOne(26, 3.0,-1.0,-1.0,-1.0, 2.0, 2.0, 2.0);	
}
void EarlyRef27::findImage(float ufx,float ufy,float ufz,float lfx,float lfy,float lfz,float * res){
    res[0] = ufx*Ps[0] + lfx*L[0];
    res[1] = ufy*Ps[1] + lfy*L[1];
    res[2] = ufz*Ps[2] + lfz*L[2];
}
void EarlyRef27::getargs(Unit *unit){
    float Pst[3],Prt[3],Lt[3],HWt,Bt;
    Pst[0] = ZIN0(1);
    Pst[1] = ZIN0(2);
    Pst[2] = ZIN0(3);
    Prt[0] = ZIN0(4);
    Prt[1] = ZIN0(5);
    Prt[2] = ZIN0(6);
    Lt[0] = ZIN0(7);
    Lt[1] = ZIN0(8);
    Lt[2] = ZIN0(9);
    HWt = ZIN0(10);
    Bt = ZIN0(11);
    bool changed = (Pst[0] != Ps[0] || Pst[1] != Ps[1] || Pst[2] != Ps[2] || Prt[0] != Pr[0] || Prt[1] != Pr[1] || Prt[2] != Pr[2] || Lt[0] != L[0] || Lt[1] != L[1] || Lt[2] != L[2] || HWt != HW || Bt != B);
    if (changed) {
        Ps[0] = Pst[0]; Ps[1] = Pst[1]; Ps[2] = Pst[2]; Pr[0] = Prt[0]; Pr[1] = Prt[1];Pr[2] = Prt[2]; L[0] = Lt[0]; L[1] = Lt[1] ; L[2] = Lt[2] ; HW = HWt ; B = Bt;
        //printf("%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n",Ps[0],Ps[1],Ps[2],Pr[0],Pr[1],Pr[2],L[0],L[1],L[2],HW,B);
        refsCalculation();
    }
}
void EarlyRef27_next(EarlyRef27* unit,int inNumSamples)
{
    unit->getargs(unit);
    float * outL = OUT(0);
    float * outR = OUT(1);
    float * in = IN(0);
    for(int i=0;i<inNumSamples; i++){
        unit->LagDel.push(in[i]);
        float oL = 0.0;
        float oR = 0.0;
        for(int j=0;j<27;j++){
            oL += unit->LagDel.delay(unit->delL[j],j*2)*unit->ampL[j];
            oR += unit->LagDel.delay(unit->delR[j],j*2+1)*unit->ampR[j];
            //oL += unit->LagDel.delay(unit->delL[j])*unit->ampL[j];
            //oR += unit->LagDel.delay(unit->delR[j])*unit->ampR[j];
        }
        outL[i] = oL;
        outR[i] = oR;
    }
}
////////////////////////////////////////////////
//include local buffer test in one place
static SndBuf * GetBuffer(Unit * unit, uint32 bufnum)
{
	SndBuf *buf;
	World *world = unit->mWorld;

	if (bufnum >= world->mNumSndBufs) {
		int localBufNum = bufnum - world->mNumSndBufs;
		Graph *parent = unit->mParent;
		if (localBufNum <= parent->localMaxBufNum) {
			buf = parent->mLocalSndBufs + localBufNum;
		} else {
			if (unit->mWorld->mVerbosity > -1)
				Print("invalid buffer number (%d).\n", bufnum);
			return NULL;
		}
	} else {
		buf = world->mSndBufs + bufnum;
	}

	if (buf->data == NULL) {
		if (unit->mWorld->mVerbosity > -1)
			Print("uninitialized buffer (%i).\n", bufnum);
		return NULL;
	}

	return buf;

}
struct EarlyRef27Gen:public Unit
{
	EarlyRef27Gen(Unit* unit);
    float CalcOne(int n,float exp,float ux,float uy,float uz,float lx,float ly,float lz);
    void refsCalculation();
    void findImage(float ufx,float ufy,float ufz,float lfx,float lfy,float lfz,float * res);
    bool getargs(Unit * unit, bool force=false);
	float L[3];
    float Ps[3];
    float Pr[3];
    float B,HW,d0;
    int N;
    float samprate;
    float bufnumL,bufnumR;
    bool let_trig;
    SndBuf *sndbufL, *sndbufR;
    unsigned int framesize,framesize_1, m_pos;
};
SCWrapClass(EarlyRef27Gen);
EarlyRef27Gen::EarlyRef27Gen(Unit *unit)
{
    int ins = 0;
    bufnumL = IN0(ins++);
    bufnumR = IN0(ins++);
    Ps[0] = IN0(ins++) - 1.0;//hack to force getargs
    samprate = FULLRATE;
    sndbufL = GetBuffer(unit,bufnumL);
    sndbufR = GetBuffer(unit,bufnumR);
    if (!sndbufL || !sndbufR){
        SETCALC(*ClearUnitOutputs);
        unit->mDone = true;
    }
    Clear(sndbufL->samples, sndbufL->data);
    Clear(sndbufR->samples, sndbufR->data);
    framesize = sc_min(sndbufL->frames,sndbufR->frames);
    framesize_1 = framesize - 1;
    let_trig = true;
    m_pos = 0;
    //printf("%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n",Ps[0],Ps[1],Ps[2],Pr[0],Pr[1],Pr[2],L[0],L[1],L[2],HW,B);
    //refsCalculation();
    ZOUT0(0) = 0.0;
    SETCALC(EarlyRef27Gen_next);
}
inline void fracdel2(float fsample,float val,float *buf, unsigned int size){
    if(fsample >=size)
        return;
	int pos1 = (int)floor(fsample);
	float frac = fsample - (float)pos1;
	buf[pos1] += val*(1.0-frac);
	buf[pos1 + 1] += val*frac;
}
float EarlyRef27Gen::CalcOne(int n,float exp,float ux,float uy,float uz,float lx,float ly,float lz)
{
    float image[3];
    float dell,delr,ampl,ampr;
    findImage(ux,uy,uz,lx,ly,lz,image);
    float rec[3];
    rec[0] = Pr[0] - HW;
    rec[1] = Pr[1];
    rec[2] = Pr[2];
    float distL = dist(image,rec);
    rec[0] = Pr[0] + HW;
    float distR = dist(image,rec);
    float preA = pow(B,exp);
    ampl = preA/(distL + 0.001);//d0*preA/(distL);
    ampr = preA/(distR + 0.001);//d0*preA/(distR);
    dell = samprate*distL/340.0;
    delr = samprate*distR/340.0;
    //printf("%g, %g, %g, %g, %g\n",d0,delL[n],delR[n],ampL[n],ampR[n]);
    /*
    if (dell < sndbufL->samples)
        sndbufL->data[(int)dell] += ampl;
    if (delr < sndbufR->samples)
        sndbufR->data[(int)delr] += ampr;
    */
    fracdel2(dell,ampl,sndbufL->data,sndbufL->frames);
    fracdel2(delr,ampr,sndbufR->data,sndbufR->frames);
    return distL;
}
/*
void EarlyRef27Gen::predist(u,v,w,l,m,n,L,Ps){
	local x = (2*u-1)*Ps[1] - 2*l*L[1] 
	local y = (2*v-1)*Ps[2] - 2*m*L[2] 
	local z = (2*w-1)*Ps[3] - 2*n*L[3]
	return {x,y,z},-TA{2*u-1,2*v-1,2*w-1},-TA{-2*l,-2*m,-2*n}
}
*/
void EarlyRef27Gen::refsCalculation(){
    
    Clear(sndbufL->samples,sndbufL->data);
    Clear(sndbufR->samples,sndbufR->data);
    d0 = 1.0;
    if(N == 0){
        //printf("EarlyRef27Gen::refsCalculation\n");
        CalcOne(0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0);	
        CalcOne(1, 1.0, 1.0, 1.0,-1.0, 0.0, 0.0, 0.0);	
        CalcOne(2, 1.0, 1.0,-1.0, 1.0, 0.0, 0.0, 0.0);	
        CalcOne(3, 2.0, 1.0,-1.0,-1.0, 0.0, 0.0, 0.0);	
        CalcOne(4, 1.0,-1.0, 1.0, 1.0, 0.0, 0.0, 0.0);	
        CalcOne(5, 2.0,-1.0, 1.0,-1.0, 0.0, 0.0, 0.0);	
        CalcOne(6, 2.0,-1.0,-1.0, 1.0, 0.0, 0.0, 0.0);	
        CalcOne(7, 3.0,-1.0,-1.0,-1.0, 0.0, 0.0, 0.0);	
        CalcOne(8, 1.0, 1.0, 1.0,-1.0, 0.0, 0.0, 2.0);	
        CalcOne(9, 2.0, 1.0,-1.0,-1.0, 0.0, 0.0, 2.0);	
        CalcOne(10, 2.0,-1.0, 1.0,-1.0, 0.0, 0.0, 2.0);	
        CalcOne(11, 3.0,-1.0,-1.0,-1.0, 0.0, 0.0, 2.0);	
        CalcOne(12, 1.0, 1.0,-1.0, 1.0, 0.0, 2.0, 0.0);	
        CalcOne(13, 2.0, 1.0,-1.0,-1.0, 0.0, 2.0, 0.0);	
        CalcOne(14, 2.0,-1.0,-1.0, 1.0, 0.0, 2.0, 0.0);	
        CalcOne(15, 3.0,-1.0,-1.0,-1.0, 0.0, 2.0, 0.0);	
        CalcOne(16, 2.0, 1.0,-1.0,-1.0, 0.0, 2.0, 2.0);	
        CalcOne(17, 3.0,-1.0,-1.0,-1.0, 0.0, 2.0, 2.0);	
        CalcOne(18, 1.0,-1.0, 1.0, 1.0, 2.0, 0.0, 0.0);	
        CalcOne(19, 2.0,-1.0, 1.0,-1.0, 2.0, 0.0, 0.0);	
        CalcOne(20, 2.0,-1.0,-1.0, 1.0, 2.0, 0.0, 0.0);	
        CalcOne(21, 3.0,-1.0,-1.0,-1.0, 2.0, 0.0, 0.0);	
        CalcOne(22, 2.0,-1.0, 1.0,-1.0, 2.0, 0.0, 2.0);	
        CalcOne(23, 3.0,-1.0,-1.0,-1.0, 2.0, 0.0, 2.0);	
        CalcOne(24, 2.0,-1.0,-1.0, 1.0, 2.0, 2.0, 0.0);	
        CalcOne(25, 3.0,-1.0,-1.0,-1.0, 2.0, 2.0, 0.0);	
        CalcOne(26, 3.0,-1.0,-1.0,-1.0, 2.0, 2.0, 2.0);	
    }else{
        //printf("EarlyRef27Gen::refsCalculation N:%d\n",N);
        for(int l = -N;l <=N;l++)
            for(int m = -N;m <=N;m++)
                for(int n = -N;n <=N;n++)
                    for(int u = 0;u <=1;u++)
                        for(int v = 0;v <=1;v++)
                            for(int w = 0;w <=1;w++){
                                float exp = abs(l-u)+abs(l)+abs(m-v)+abs(m)+abs(n-w)+abs(n);
                                float ux = 2*u-1;
                                float uy = 2*v-1;
                                float uz = 2*w-1;
                                float lx = 2*l;
                                float ly = 2*m;
                                float lz = 2*n;
                                CalcOne(0,exp,-ux,-uy,-uz,lx,ly,lz);
                            }
    }
}
void EarlyRef27Gen::findImage(float ufx,float ufy,float ufz,float lfx,float lfy,float lfz,float * res){
    res[0] = ufx*Ps[0] + lfx*L[0];
    res[1] = ufy*Ps[1] + lfy*L[1];
    res[2] = ufz*Ps[2] + lfz*L[2];
}
bool EarlyRef27Gen::getargs(Unit *unit,bool force){
    float Pst[3],Prt[3],Lt[3],HWt,Bt;
    int Nt;
    float Ntf;
    int ins = 2;
    Pst[0] = ZIN0(ins++);
    Pst[1] = ZIN0(ins++);
    Pst[2] = ZIN0(ins++);
    Prt[0] = ZIN0(ins++);
    Prt[1] = ZIN0(ins++);
    Prt[2] = ZIN0(ins++);
    Lt[0] = ZIN0(ins++);
    Lt[1] = ZIN0(ins++);
    Lt[2] = ZIN0(ins++);
    HWt = ZIN0(ins++);
    Bt = sc_clip(ZIN0(ins++),-1.0,1.0);
    Ntf = ZIN0(ins++);
    Ntf = sc_clip(Ntf,0.0f,5.0f);//if done above dont works!!
    Nt = Ntf;
    bool changed = (Pst[0] != Ps[0] || Pst[1] != Ps[1] || Pst[2] != Ps[2] || Prt[0] != Pr[0] || Prt[1] != Pr[1] || Prt[2] != Pr[2] || Lt[0] != L[0] || Lt[1] != L[1] || Lt[2] != L[2] || HWt != HW || Bt != B || Nt != N);
    if (changed || force) {
        //printf("Nt %d %f \n",Nt,Ntf);
        Ps[0] = Pst[0]; Ps[1] = Pst[1]; Ps[2] = Pst[2]; Pr[0] = Prt[0]; Pr[1] = Prt[1];Pr[2] = Prt[2]; L[0] = Lt[0]; L[1] = Lt[1] ; L[2] = Lt[2] ; HW = HWt ; B = Bt; N = Nt;
        //printf("%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n",Ps[0],Ps[1],Ps[2],Pr[0],Pr[1],Pr[2],L[0],L[1],L[2],HW,B);
        refsCalculation();
    }
    return changed;
}
void EarlyRef27Gen_next(EarlyRef27Gen* unit,int inwrongNumSamples)
{
    int numSamples = unit->mWorld->mFullRate.mBufLength;
    unit->m_pos += numSamples;
    if(unit->m_pos >= unit->framesize){
        unit->let_trig = true;
        unit->m_pos = 0;
    }
    if(unit->let_trig){
        if(unit->getargs(unit)){
            ZOUT0(0) = 1.0;
            //printf("send trig\n");
            unit->let_trig = false;
            unit->m_pos = 0;
        }else{
            ZOUT0(0) = 0.0;
        }
    }else{
            ZOUT0(0) = 0.0;
    }
}
////////////ambisonics bformat
struct EarlyRefAtkGen:public Unit
{
	EarlyRefAtkGen(Unit* unit);
    void CalcOne(int n,float exp,float ux,float uy,float uz,float lx,float ly,float lz);
    void refsCalculation();
    void findImage(float ufx,float ufy,float ufz,float lfx,float lfy,float lfz,float * res);
    bool getargs(Unit * unit, bool force=false);
	float L[3];
    float Ps[3];
    float Pr[3];
    float B,HW,d0;
    int N;
    float samprate;
    float bufnumW,bufnumX,bufnumY,bufnumZ;
    bool let_trig;
    SndBuf *sndbufW, *sndbufX,*sndbufY, *sndbufZ;
    unsigned int framesize,framesize_1, m_pos;
};
SCWrapClass(EarlyRefAtkGen);
EarlyRefAtkGen::EarlyRefAtkGen(Unit *unit)
{
    int ins = 0;
    bufnumW = IN0(ins++);
    bufnumX = IN0(ins++);
    bufnumY = IN0(ins++);
    bufnumZ = IN0(ins++);
    Ps[0] = IN0(ins++) - 1.0;//hack to force getargs
    samprate = FULLRATE;
    sndbufW = GetBuffer(unit,bufnumW);
    sndbufX = GetBuffer(unit,bufnumX);
    sndbufY = GetBuffer(unit,bufnumY);
    sndbufZ = GetBuffer(unit,bufnumZ);
    if (!sndbufW || !sndbufX || !sndbufY || !sndbufZ){
        SETCALC(*ClearUnitOutputs);
        unit->mDone = true;
    }
    Clear(sndbufW->samples, sndbufW->data);
    Clear(sndbufX->samples, sndbufX->data);
    Clear(sndbufY->samples, sndbufY->data);
    Clear(sndbufZ->samples, sndbufZ->data);
    framesize = sc_min(sc_min(sc_min(sndbufW->frames,sndbufX->frames),sndbufY->frames),sndbufZ->frames);
    framesize_1 = framesize -1;
    let_trig = true;
    m_pos = 0;
    //printf("%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n",Ps[0],Ps[1],Ps[2],Pr[0],Pr[1],Pr[2],L[0],L[1],L[2],HW,B);
    //refsCalculation();
    ZOUT0(0) = 0.0;
    SETCALC(EarlyRefAtkGen_next);
}
inline void setfracdel(int pos,float frac,float val,float *buf){
	buf[pos] += val*(1.0-frac);
	buf[pos + 1] += val*frac;
}
static double invsqrt2 = 1.0/sqrt(2.0);
void EarlyRefAtkGen::CalcOne(int n,float exp,float ux,float uy,float uz,float lx,float ly,float lz)
{
    float image[3];
    float del,amp;
    findImage(ux,uy,uz,lx,ly,lz,image);
    /*
    double vec[3];
    vec[0] = image[0] - Pr[0];
    vec[1] = image[1] - Pr[1];
    vec[2] = image[2] - Pr[2];
    
    double planeDist = hypot(vec[1], vec[0]);
	double azim = atan2(vec[1], vec[0]) - M_PI*0.5; 
	double elev = atan2(vec[2], planeDist);
	double dista = hypot(planeDist, vec[2]);
    //float dista = sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
    float preA = pow(B,exp);
    amp = preA/(dista + 0.1);//d0*preA/(distL);
    del = samprate*dista/340.0;
    double cosa = cos(azim);
    double sina = sin(azim);
    double cosb = cos(elev);
    double sinb = sin(elev);

    fracdel2(del,amp*invsqrt2,sndbufW->data,sndbufW->frames);
    fracdel2(del,amp*cosa*cosb,sndbufX->data,sndbufX->frames);
    fracdel2(del,amp*sina*cosb,sndbufY->data,sndbufY->frames);
    fracdel2(del,amp*sinb,sndbufZ->data,sndbufZ->frames);
    return dista;
    */
    //can be done  without trig
    //rotation(x,y,z) -> (y,-x,z)
    double x = image[1] - Pr[1];
    double y = Pr[0] - image[0];
    double z = image[2] - Pr[2];
    
    double dista_2 = x*x + y*y +z*z; 
    double dista = sqrt(dista_2);
        
    del = samprate*dista/340.0;
    if (dista < 1e-9 || (unsigned int)del > framesize_1)
        return;
    double preA = pow(B,exp);
    amp = preA/dista;
    float amp_d2 = preA/dista_2;
    int pos = (int)floor(del);
	float frac = del - (float)pos;
    setfracdel(pos,frac,amp*invsqrt2,sndbufW->data);
    setfracdel(pos,frac,amp_d2*x,sndbufX->data);
    setfracdel(pos,frac,amp_d2*y,sndbufY->data);
    setfracdel(pos,frac,amp_d2*z,sndbufZ->data);
}
/*
void EarlyRefAtkGen::predist(u,v,w,l,m,n,L,Ps){
	local x = (2*u-1)*Ps[1] - 2*l*L[1] 
	local y = (2*v-1)*Ps[2] - 2*m*L[2] 
	local z = (2*w-1)*Ps[3] - 2*n*L[3]
	return {x,y,z},-TA{2*u-1,2*v-1,2*w-1},-TA{-2*l,-2*m,-2*n}
}
*/
void EarlyRefAtkGen::refsCalculation(){
    
    Clear(sndbufW->samples,sndbufW->data);
    Clear(sndbufX->samples,sndbufX->data);
    Clear(sndbufY->samples,sndbufY->data);
    Clear(sndbufZ->samples,sndbufZ->data);
    d0 = 1.0;
    if(N == 0){
        //printf("EarlyRefAtkGen::refsCalculation\n");
        CalcOne(0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0);	
        CalcOne(1, 1.0, 1.0, 1.0,-1.0, 0.0, 0.0, 0.0);	
        CalcOne(2, 1.0, 1.0,-1.0, 1.0, 0.0, 0.0, 0.0);	
        CalcOne(3, 2.0, 1.0,-1.0,-1.0, 0.0, 0.0, 0.0);	
        CalcOne(4, 1.0,-1.0, 1.0, 1.0, 0.0, 0.0, 0.0);	
        CalcOne(5, 2.0,-1.0, 1.0,-1.0, 0.0, 0.0, 0.0);	
        CalcOne(6, 2.0,-1.0,-1.0, 1.0, 0.0, 0.0, 0.0);	
        CalcOne(7, 3.0,-1.0,-1.0,-1.0, 0.0, 0.0, 0.0);	
        CalcOne(8, 1.0, 1.0, 1.0,-1.0, 0.0, 0.0, 2.0);	
        CalcOne(9, 2.0, 1.0,-1.0,-1.0, 0.0, 0.0, 2.0);	
        CalcOne(10, 2.0,-1.0, 1.0,-1.0, 0.0, 0.0, 2.0);	
        CalcOne(11, 3.0,-1.0,-1.0,-1.0, 0.0, 0.0, 2.0);	
        CalcOne(12, 1.0, 1.0,-1.0, 1.0, 0.0, 2.0, 0.0);	
        CalcOne(13, 2.0, 1.0,-1.0,-1.0, 0.0, 2.0, 0.0);	
        CalcOne(14, 2.0,-1.0,-1.0, 1.0, 0.0, 2.0, 0.0);	
        CalcOne(15, 3.0,-1.0,-1.0,-1.0, 0.0, 2.0, 0.0);	
        CalcOne(16, 2.0, 1.0,-1.0,-1.0, 0.0, 2.0, 2.0);	
        CalcOne(17, 3.0,-1.0,-1.0,-1.0, 0.0, 2.0, 2.0);	
        CalcOne(18, 1.0,-1.0, 1.0, 1.0, 2.0, 0.0, 0.0);	
        CalcOne(19, 2.0,-1.0, 1.0,-1.0, 2.0, 0.0, 0.0);	
        CalcOne(20, 2.0,-1.0,-1.0, 1.0, 2.0, 0.0, 0.0);	
        CalcOne(21, 3.0,-1.0,-1.0,-1.0, 2.0, 0.0, 0.0);	
        CalcOne(22, 2.0,-1.0, 1.0,-1.0, 2.0, 0.0, 2.0);	
        CalcOne(23, 3.0,-1.0,-1.0,-1.0, 2.0, 0.0, 2.0);	
        CalcOne(24, 2.0,-1.0,-1.0, 1.0, 2.0, 2.0, 0.0);	
        CalcOne(25, 3.0,-1.0,-1.0,-1.0, 2.0, 2.0, 0.0);	
        CalcOne(26, 3.0,-1.0,-1.0,-1.0, 2.0, 2.0, 2.0);	
    }else{
        //printf("EarlyRefAtkGen::refsCalculation N:%d\n",N);
        for(int l = -N;l <=N;l++)
            for(int m = -N;m <=N;m++)
                for(int n = -N;n <=N;n++)
                    for(int u = 0;u <=1;u++)
                        for(int v = 0;v <=1;v++)
                            for(int w = 0;w <=1;w++){
                                float exp = abs(l-u)+abs(l)+abs(m-v)+abs(m)+abs(n-w)+abs(n);
                                float ux = 2*u-1;
                                float uy = 2*v-1;
                                float uz = 2*w-1;
                                float lx = 2*l;
                                float ly = 2*m;
                                float lz = 2*n;
                                CalcOne(0,exp,-ux,-uy,-uz,lx,ly,lz);
                            }
    }
}
void EarlyRefAtkGen::findImage(float ufx,float ufy,float ufz,float lfx,float lfy,float lfz,float * res){
    res[0] = ufx*Ps[0] + lfx*L[0];
    res[1] = ufy*Ps[1] + lfy*L[1];
    res[2] = ufz*Ps[2] + lfz*L[2];
}
bool EarlyRefAtkGen::getargs(Unit *unit,bool force){
    float Pst[3],Prt[3],Lt[3],HWt,Bt;
    int Nt;
    float Ntf;
    int ins = 4;
    Pst[0] = ZIN0(ins++);
    Pst[1] = ZIN0(ins++);
    Pst[2] = ZIN0(ins++);
    Prt[0] = ZIN0(ins++);
    Prt[1] = ZIN0(ins++);
    Prt[2] = ZIN0(ins++);
    Lt[0] = ZIN0(ins++);
    Lt[1] = ZIN0(ins++);
    Lt[2] = ZIN0(ins++);
    HWt = ZIN0(ins++);
    Bt = sc_clip(ZIN0(ins++),-1.0,1.0);
    Ntf = ZIN0(ins++);
    Ntf = sc_clip(Ntf,0.0f,5.0f);//if done above dont works!!
    Nt = Ntf;
    bool changed = (Pst[0] != Ps[0] || Pst[1] != Ps[1] || Pst[2] != Ps[2] || Prt[0] != Pr[0] || Prt[1] != Pr[1] || Prt[2] != Pr[2] || Lt[0] != L[0] || Lt[1] != L[1] || Lt[2] != L[2] || HWt != HW || Bt != B || Nt != N);
    if (changed || force) {
        //printf("Nt %d %f \n",Nt,Ntf);
        Ps[0] = Pst[0]; Ps[1] = Pst[1]; Ps[2] = Pst[2]; Pr[0] = Prt[0]; Pr[1] = Prt[1];Pr[2] = Prt[2]; L[0] = Lt[0]; L[1] = Lt[1] ; L[2] = Lt[2] ; HW = HWt ; B = Bt; N = Nt;
        //printf("%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n",Ps[0],Ps[1],Ps[2],Pr[0],Pr[1],Pr[2],L[0],L[1],L[2],HW,B);
        refsCalculation();
    }
    return changed;
}
void EarlyRefAtkGen_next(EarlyRefAtkGen* unit,int inwrongNumSamples)
{
    int numSamples = unit->mWorld->mFullRate.mBufLength;
    unit->m_pos += numSamples;
    if(unit->m_pos >= unit->framesize){
        unit->let_trig = true;
        unit->m_pos = 0;
    }
    if(unit->let_trig){
        if(unit->getargs(unit)){
            ZOUT0(0) = 1.0;
            //printf("send trig\n");
            unit->let_trig = false;
            unit->m_pos = 0;
        }else{
            ZOUT0(0) = 0.0;
        }
    }else{
            ZOUT0(0) = 0.0;
    }
}
PluginLoad(DWGReverb){
	ft = inTable;
    DefineDtorUnit(EarlyRef);
	DefineDtorUnit(EarlyRef27);
    DefineDtorUnit(EarlyRef27Gen);
    DefineDtorUnit(EarlyRefAtkGen);
    DefineDtorUnit(DWGReverb);
}

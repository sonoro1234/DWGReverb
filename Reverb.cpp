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


///////////////////DWGAllpass
struct DWGAllpass : public Unit
{
    AllPass2powT<1024*2> allpass;
    DWGAllpass(Unit *unit);
};

SCWrapClassT(DWGAllpass);

DWGAllpass::DWGAllpass(Unit *unit){
    //SETCALC(DWGAllpass_next);
}

void DWGAllpass_next(DWGAllpass *unit,int numsamples)
{
    float* in = ZIN(0);
    float* out = ZOUT(0);
    int M = ZIN0(1);
    float a = sc_clip(ZIN0(2),-1,1);
    unit->allpass.set_coeffs(M,a);
    for(int i=0;i<numsamples;i++){
        out[i] = unit->allpass.filter(in[i]);
    }
}


//////////////////////////////Reverb FDN
template<template<int,int>class Tfdn,int size,int sizedel>
struct DWGReverbBase : public Unit
{
    Tfdn<size,sizedel> fdn;

	void setlengths();
	//void setlengths2();
    float len,mix;
    float SR;
	float factors[size];
	float doprime;
	DWGReverbBase(DWGReverbBase *unit);
};
template<template<int,int>class Tfdn,int size,int sizedel>
DWGReverbBase<Tfdn,size,sizedel>::DWGReverbBase(DWGReverbBase *unit){
    SR = SAMPLERATE;
}

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
template<template<int,int>class Tfdn,int size,int sizedel>
void DWGReverbBase<Tfdn,size,sizedel>::setlengths(){
    //Print("setlengths\n");
	float lent;
    float lengths[size];
	for(int i = 0;i < size ;i++){
		if(doprime==1){
			lent = nearestprime(len * factors[i],0.5);
            if(lent == -1){
                Print("DWGReverb error nearestprime\n");
                lent = len * factors[i];
            }else{
                //Print("nearest %g",lent);
            }
            //Print("\n");
		}else
			lent = len * factors[i];
		//Print("len %d: %g factor:%g \n",i,lent,factors[i]);
		lengths[i] = lent;
	}
    fdn.setlengths(lengths);
	//Print("\n");
}
/*
template<typename Tfdn>
void DWGReverbBase<Tfdn>::setlengths2(){
	int primes[8] = {19,17,13,11,7,5,3,2};
	//float factors[8] = {1,0.9464,0.87352,0.83,0.8123,0.7398,0.69346,0.6349};

	float lent;
	int lent2;
	for(int i = 0;i < 8 ;i++){
		lent = len * factors[i];
		int mi = floor(0.5 + (log(lent)/log(primes[i])));
		lent2 = pow(primes[i],mi);
		Print("nearest %g %d %d\n",lent,lent2,mi);
		fdn.lengths[i] = lent2;
	}
	Print("\n");
}
*/
///C1C3 derivation
template<int size,int sizedel>
struct DWGReverbC1C3T : public DWGReverbBase<FDNC1C3,size,sizedel>
{
    float c1,c3;//for c1,c3 implementation
	void setcoeffs();
    void get_args(Unit* unit,bool force = false);
	DWGReverbC1C3T(DWGReverbC1C3T *unit):DWGReverbBase<FDNC1C3,size,sizedel>(unit){get_args(unit,true);};
};
template<int size,int sizedel>
void DWGReverbC1C3T<size,sizedel>::get_args(Unit*unit,bool force)
{

    float lent = sc_clip(ZIN0(1),2.0,1024.0*8.0);
	float c1_t = sc_clip(ZIN0(2),0.1,1000);
	float c3_t = sc_clip(ZIN0(3),0.1,1000);
    bool docoefs = false;
	float mix = this->fdn.mix = sc_clip(ZIN0(4),0.0,1.0);
    this->doprime = ZIN0(5);
	for(int i=0;i<size;i++)
		this->factors[i] = ZIN0(6+i);
    for(int i=0;i<size;i++)
		this->fdn.o_perm[i] = ZIN0(6+size+i);
    if(lent != this->len || force){
        this->len = lent;
        this->setlengths();
        docoefs = true;
    }
    if(docoefs || c1_t != c1 || c3_t != c3 || force){
        c1 = c1_t;
        c3 = c3_t;
        this->setcoeffs();
    }
}
template<int size,int sizedel>
void DWGReverbC1C3T<size,sizedel>::setcoeffs(){
		this->fdn.setcoeffs(c1,c3,1.0,this->SR);
}

////////////instantiation DWGReverbC1C3
typedef  DWGReverbC1C3T<8,1024*8> DWGReverbC1C3;
SCWrapClassT(DWGReverbC1C3);
void DWGReverbC1C3_next(DWGReverbC1C3  *unit, int inNumSamples)
{
	float *out[2];
	out[0] = OUT(0);
	out[1] = OUT(1);
	float *in = IN(0);
    unit->get_args(unit);

	unit->fdn.go_st(in,out,inNumSamples);
}
////////////instantiation DWGReverbC1C3_16
typedef  DWGReverbC1C3T<16,1024*8> DWGReverbC1C3_16;
SCWrapClassT(DWGReverbC1C3_16);
void DWGReverbC1C3_16_next(DWGReverbC1C3_16  *unit, int inNumSamples)
{
	float *out[2];
	out[0] = OUT(0);
	out[1] = OUT(1);
	float *in = IN(0);
    unit->get_args(unit);

	unit->fdn.go_st(in,out,inNumSamples);
}

//////////////////////////implementation Zitarev1 filter
class Filt1
{
public:
    
    Filt1 (void) : _slo (0), _shi (0) {}
    ~Filt1 (void) {}
    void calc_chi(float _fsamp,float _xover,float _fdamp,float &wlo,float & chi){
        wlo = 6.2832f * _xover / _fsamp;
        if (_fdamp > 0.49f * _fsamp)
            chi = 2.0f;
        else 
            chi = 1.0f - cosf (6.2832f * _fdamp / _fsamp);
    }
    //void  set_params (float del, float tmf, float tlo, float wlo, float thi, float chi);
    void set_params (float del, float tmf, float tlo, float wlo, float thi, float chi)
    {
        float g, t;

        _gmf = powf (0.001f, del / tmf);
        _glo = powf (0.001f, del / tlo) / _gmf - 1.0f;
        _wlo = wlo;    
        g = powf (0.001f, del / thi) / _gmf;
        t = (1 - g * g) / (2 * g * g * chi);
        _whi = (sqrtf (1 + 4 * t) - 1) / (2 * t); 
    } 
    float filter (float x)
    {
        _slo += _wlo * (x - _slo); //+ 1e-10f;
        x += _glo * _slo;
        _shi += _whi * (x - _shi);
        return _gmf * _shi;
    }
    float   _gmf;
    float   _glo;
    float   _wlo;
    float   _whi;
    float   _slo;
    float   _shi;    
};

template<int size,int sizedel>
struct FDNZita: public FDN_HH_Base<Filt1,size,sizedel>{
    //float xover,rtlow,rtmid,fdamp;
    void setcoeffs(float xover,float rtlow, float rtmid,float fdamp,float SR){
        float wlo,chi;
        this->decay[0].calc_chi(SR, xover, fdamp, wlo, chi);
        for (int i = 0; i < size; i++){
             this->decay[i].set_params(this->lengths[i]/SR, rtmid, rtlow, wlo, 0.5f * rtmid, chi);
        }
    }        
};

//Zita derivation
template<int size,int sizedel>
struct DWGReverbZitaT : public DWGReverbBase<FDNZita,size,sizedel>
{
    float xover,rtlow,rtmid,fdamp;
	void setcoeffs();
    void get_args(Unit* unit,bool force = false);
	DWGReverbZitaT(DWGReverbZitaT *unit):DWGReverbBase<FDNZita,size,sizedel>(unit){get_args(unit,true);};
};
template<int size,int sizedel>
void DWGReverbZitaT<size,sizedel>::setcoeffs(){
        this->fdn.setcoeffs(xover,rtlow,rtmid,fdamp,this->SR);
}
//len,xover,rtlow,rtmid,fdamp,doprime
template<int size,int sizedel>
void DWGReverbZitaT<size,sizedel>::get_args(Unit *unit,bool force)
{
    float lent = sc_clip(ZIN0(1),2.0,1024.0*8.0);
	float xover_t = sc_clip(ZIN0(2),50,1500);
	float rtlow_t = sc_clip(ZIN0(3),0.25,12);
	float rtmid_t = sc_clip(ZIN0(4),0.25,12);
	float fdamp_t = sc_clip(ZIN0(5),1500,24000); 
    bool docoefs = false;
	this->mix = this->fdn.mix= 1.0;//sc_clip(ZIN0(4),0.0,1.0);
    this->doprime = ZIN0(6);
	for(int i=0;i<size;i++)
		this->factors[i] = ZIN0(7+i);
    for(int i=0;i<size;i++)
		this->fdn.o_perm[i] = ZIN0(7 + size + i);
    if(lent != this->len || force){
        this->len = lent;
        this->setlengths();
        docoefs = true;
    }
    if(docoefs || xover_t != xover || rtlow_t != rtlow || rtmid_t != rtmid || fdamp_t != fdamp || force){
        xover = xover_t;
        rtlow = rtlow_t;
        rtmid = rtmid_t;
        fdamp = fdamp_t;
        this->setcoeffs();
    }
	/*
    if(force){
        for(int i=1;i<7+2*size; i++)
            printf("ZIN(%d) : %g\n",i,ZIN0(i)); 
    }*/
}
//instantiations Zita8
typedef  DWGReverbZitaT<8,1024*8> DWGReverb3Band;
SCWrapClassT(DWGReverb3Band);
void DWGReverb3Band_next(DWGReverb3Band  *unit, int inNumSamples)
{
	float *out[2];
	out[0] = OUT(0);
	out[1] = OUT(1);
	float *in = IN(0);
    unit->get_args(unit);

	unit->fdn.go_st(in,out,inNumSamples);
}

////////////instantiation Zita_16
typedef  DWGReverbZitaT<16,1024*8> DWGReverb3Band_16;
SCWrapClassT(DWGReverb3Band_16);
void DWGReverb3Band_16_next(DWGReverb3Band_16  *unit, int inNumSamples)
{
	float *out[2];
	out[0] = OUT(0);
	out[1] = OUT(1);
	float *in = IN(0);
    unit->get_args(unit);

	unit->fdn.go_st(in,out,inNumSamples);
}
////////////////////////////////////////////EarlyRef Ugens/////////////////////////////////

float dist(float im[3],float r[3]){
    float x = im[0] - r[0];
    float y = im[1] - r[1];
    float z = im[2] - r[2];
    return sqrt(x*x+y*y+z*z);
}
float dist(float im[3],float r0,float r1,float r2){
    float x = im[0] - r0;
    float y = im[1] - r1;
    float z = im[2] - r2;
    return sqrt(x*x+y*y+z*z);
}
////////////Kendal 1986
//double feedback filtered comb
template<int size>
struct R2{
	int M[2];
	float p,f[2];
	float lastout;
	CircularBuffer2POWSizedT<size> delays[2];
	LTITv<1,1> filters[2];
	R2():lastout(0.0){};
	void set_coeffs(int M1,int M2,float p_t, float f_t1,float f_t2){
		//printf("%d, %d, %f,%f,%f\n",M1,M2,p_t,f_t1,f_t2);
		M[0] = std::min(size,std::max(1,M1)) - 1;
		M[1] = std::min(size,std::max(1,M2)) - 1;
        p = std::min(1.0f,std::max(0.0f,p_t));
		f[0] = std::min(1.0f,std::max(0.0f,f_t1));
		f[1] = std::min(1.0f,std::max(0.0f,f_t2));
		for(int i=0; i< 2; i++){
			filters[i].KernelA = -p;
			filters[i].KernelB = 1.0 - p;
		}
	}  
	inline float filter(float in){
		float odel1 = delays[0].delay(M[0])*f[0];
		odel1 = filters[0].filter(odel1);
		float odel2 = delays[1].delay(M[1])*f[1];
		odel2 = filters[1].filter(odel2);
		delays[0].push(in + odel2);
		delays[1].push(odel1);
		lastout = in + odel1 + odel2;
		return lastout;
	}
};

struct Kendall:public Unit{
	R2<4096*2> R2dels[18];
	R2<4096*2> R2delsR[18];
	float rooms[18][3] = {{0,0,1},{-1,0,0},{0,1,0},{1,0,0},{0,-1,0},{0,0,-1},
		{-1,0,1},{0,1,1},{1,0,1},{0,-1,1},{-1,-1,0},{-1,1,0},{1,1,0},{1,-1,0},
			{-1,0,-1},{0,1,-1},{1,0,-1},{0,-1,-1}};
	int inputs[6][4] = {{6,8,7,9},{10,11,6,14},{11,12,7,15},{12,13,8,16},{13,10,9,17},{14,16,15,17}};//,
		//{0,1},{0,2},{0,3},{0,4}};
	CircularBuffer2POWSizedT<4096*8> MtapDel;
	Kendall(Unit *unit);
	void findDist(float room[3],float rec[2]);
	void CalcOne(int n,float room[3]);
	void refsCalculation();
	void get_args(Unit *unit,bool force = false);
	void go(float *in,float *outL, float *outR,int numsamp);
	float L[3];
    float Ps[3];
    float Pr[3];
    float Ps_[3];//center room is 0,0
    float Pr_[3];

    float B,HW,d0,p;

    float SR;
	float deldirL,deldirR,ampdirL,ampdirR;
    float delL[18],delR[18],ampL[18],ampR[18];
    int RefN[18];
};
SCWrapClassT(Kendall);
Kendall::Kendall(Unit *unit){

	SR = SAMPLERATE;
	get_args(unit,true);
	//SETCALC(Kendall_next);
}

void Kendall::findDist(float room[3],float dists[2]){
    float u[3],res[3];
    for(int i=0; i<3; i++){
        u[i] = 1.0 - 2.0*float((int)room[i] & 1);
        res[i] = u[i]*Ps_[i] + room[i]*L[i];
    }
	dists[0] = dist(res,Pr_[0] - HW,Pr_[1],Pr_[2]);
	dists[1] = dist(res,Pr_[0] + HW,Pr_[1],Pr_[2]);
}
void Kendall::CalcOne(int n,float room[3]){
    float r2[3],r3[3];
    for(int i=0; i<3; i++){
        r2[i] = room[i]*2.0;
        r3[i] = room[i]*3.0;
    }
	float exp = 0.0;
	for(int i=0; i<3; i++)
		exp += abs(room[i]);
	
	float dist1[2],dist2[2],dist3[2];
    findDist(room,dist1);
    findDist(r2,dist2);
    findDist(r3,dist3);

	float d1 = dist2[0] - dist1[0];
	float d2 = dist3[0] - dist2[0];
	//printf("coeffsL %d\n",n);
	R2dels[n].set_coeffs(SR*d1/340.0,SR*d2/340.0,p,B*dist1[0]/dist2[0],B*dist2[0]/dist3[0]);
	d1 = dist2[1] - dist1[1];
	d2 = dist3[1] - dist2[1];
	//printf("coeffsR %d\n",n);
	R2delsR[n].set_coeffs(SR*d1/340.0,SR*d2/340.0,p,B*dist1[1]/dist2[1],B*dist2[1]/dist3[1]);
	
	delL[n] = SR*dist1[0]/340.0;
	ampL[n] = pow(B,exp)/dist1[0];
	delR[n] = SR*dist1[1]/340.0;
	ampR[n] = pow(B,exp)/dist1[1];
}



void Kendall::refsCalculation(){
	Ps_[0] = Ps[0] - L[0]*0.5;
    Ps_[1] = Ps[1] - L[1]*0.5;
    Ps_[2] = Ps[2] - L[2]*0.5;
    Pr_[0] = Pr[0] - L[0]*0.5;
    Pr_[1] = Pr[1] - L[1]*0.5;
    Pr_[2] = Pr[2] - L[2]*0.5;
	deldirL = dist(Ps_,Pr_[0] - HW,Pr_[1],Pr_[2]);
	deldirR = dist(Ps_,Pr_[0] + HW,Pr_[1],Pr_[2]);
	ampdirL = 1.0/deldirL;
	ampdirR = 1.0/deldirR;
	deldirL *=SR/340.0;
	deldirR *=SR/340.0;
	for(int i=0; i <18; i++)
		CalcOne(i,rooms[i]);
	
}
void Kendall::go(float *in,float *outL, float *outR,int num){
	float oo[18],ooR[18];
	for(int i=0; i<num; i++){
		MtapDel.push(in[i]);
		//prepare inputs
		for(int j=0; j <18; j++){
			oo[j] = MtapDel.delay(delL[j])*ampL[j];
			ooR[j] = MtapDel.delay(delR[j])*ampR[j];
			
			if(j < 6)
			for(int k=0; k <4; k++){
				oo[j] += R2dels[inputs[j][k]].lastout;
				ooR[j] += R2delsR[inputs[j][k]].lastout;
			}
		}
		//getouts
		outL[i] = MtapDel.delay(deldirL)*ampdirL;
		outR[i] = MtapDel.delay(deldirR)*ampdirR;
		for(int j=0; j <18; j++){
			outL[i] += R2dels[j].filter(oo[j]);
			outR[i] += R2delsR[j].filter(ooR[j]);
		}
	}
}
void Kendall::get_args(Unit *unit,bool force){
	float Pst[3],Prt[3],Lt[3],HWt,Bt,p_t;
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
    Bt = sc_clip(ZIN0(11),-1.0,1.0);
    p_t = sc_clip(ZIN0(12),0.0,1.0);

    bool changed = (Pst[0] != Ps[0] || Pst[1] != Ps[1] || Pst[2] != Ps[2] || Prt[0] != Pr[0] || Prt[1] != Pr[1] || Prt[2] != Pr[2] || Lt[0] != L[0] || Lt[1] != L[1] || Lt[2] != L[2] || HWt != HW || Bt != B || p_t != p);
    if (changed || force) {
        Ps[0] = Pst[0]; Ps[1] = Pst[1]; Ps[2] = Pst[2]; Pr[0] = Prt[0]; Pr[1] = Prt[1];Pr[2] = Prt[2]; L[0] = Lt[0]; L[1] = Lt[1] ; L[2] = Lt[2] ; HW = HWt ; B = Bt;p = p_t;
        //printf("%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n",Ps[0],Ps[1],Ps[2],Pr[0],Pr[1],Pr[2],L[0],L[1],L[2],HW,B);
        refsCalculation();
    }
}

void Kendall_next(Kendall *unit,int numsamples){
	unit->get_args(unit);
    float * outL = OUT(0);
    float * outR = OUT(1);
    float * in = IN(0);
	unit->go(in,outL,outR,numsamples);
}
/////////////EarlyRef
const int MaxNits = 6;
//const int Nits = MaxNits;
//const int Nrefs = pow((2*Nits + 1),3)*8;
const int MaxNrefs = pow((2*MaxNits + 1),3);
//const int Nrefs = MaxNrefs;
struct EarlyRef:public Unit
{
	Unit * unit;
    EarlyRef(Unit* unit);
    float CalcOne(int n,float exp,float ux,float uy,float uz,float lx,float ly,float lz);
    void refsCalculation();
    void filters_init();
    void allpass_init();
    void filters_tick(float in);
    void filters_tickLR(float *inL,float *inR,float &outL,float &outR);
    void findImage(float ufx,float ufy,float ufz,float lfx,float lfy,float lfz,float * res);
    void getargs(bool force=false);
    CircularBuffer2POWSizedT<4096*16> MtapDel;
    //CircularBuffer2POWSizedT<4096*16> DelR;
    //CircularBuffer2POWSizedT<4096*16> DelL;
    LTITv<1,1> filters[4];//for version 2
    LTITv<1,1> filtersL[4];
    LTITv<1,1> filtersR[4];
    float filters_out[5];
	float L[3];
    float Ps[3];
    float Pr[3];
    float Ps_[3];//center room is 0,0
    float Pr_[3];
    AllPass2powT<1024*2> allpassL[3];
    AllPass2powT<1024*2> allpassR[3];
    int allp_lens[3];
    float allp_coef;
    float B,HW,d0,p;
    int N,Nrefs;
    float samprate;
    float delL[MaxNrefs],delR[MaxNrefs],ampL[MaxNrefs],ampR[MaxNrefs];
    int RefN[MaxNrefs];
};
SCWrapClassT(EarlyRef);
EarlyRef::EarlyRef(Unit *unit)
{
    this->unit = unit;
    samprate = SAMPLERATE;
    getargs(true);
   // SETCALC(EarlyRef_next);
}
void EarlyRef::filters_init(){
    for(int i= 0;i<4;i++){
        filtersL[i].KernelA = filtersR[i].KernelA = filters[i].KernelA = -p;
        filtersL[i].KernelB = filtersR[i].KernelB = filters[i].KernelB = 1.0 - p;
    }
}
void EarlyRef::allpass_init(){
    for(int i=0;i < 3;i++){
        allpassL[i].set_coeffs(allp_lens[i],allp_coef);
        allpassR[i].set_coeffs(allp_lens[i],allp_coef);
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
    if(allp_coef!=0.0)
    for(int i =0;i<3;i++){
        outL = allpassL[i].filter(outL);
        outR = allpassR[i].filter(outR);
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
    //printf("num is %d MaNrefs is %d Nrefs is %d\n",num,MaxNrefs,Nrefs);
}
void EarlyRef::findImage(float ufx,float ufy,float ufz,float lfx,float lfy,float lfz,float * res){
    res[0] = ufx*Ps_[0] + lfx*L[0];
    res[1] = ufy*Ps_[1] + lfy*L[1];
    res[2] = ufz*Ps_[2] + lfz*L[2];
}
void EarlyRef::getargs(bool force){
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
    Bt = sc_clip(ZIN0(11),-1.0,1.0);
    Nt = sc_clip(ZIN0(12),0.0,(float)MaxNits);

    bool changed = (Pst[0] != Ps[0] || Pst[1] != Ps[1] || Pst[2] != Ps[2] || Prt[0] != Pr[0] || Prt[1] != Pr[1] || Prt[2] != Pr[2] || Lt[0] != L[0] || Lt[1] != L[1] || Lt[2] != L[2] || HWt != HW || Bt != B || Nt != N);
    if (changed || force) {
        Ps[0] = Pst[0]; Ps[1] = Pst[1]; Ps[2] = Pst[2]; Pr[0] = Prt[0]; Pr[1] = Prt[1];Pr[2] = Prt[2]; L[0] = Lt[0]; L[1] = Lt[1] ; L[2] = Lt[2] ; HW = HWt ; B = Bt;N = Nt;
        //printf("%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n",Ps[0],Ps[1],Ps[2],Pr[0],Pr[1],Pr[2],L[0],L[1],L[2],HW,B);
        Nrefs = pow((2*N + 1),3);
        refsCalculation();
    }
    
    pt = sc_clip(ZIN0(13),0.0,0.9);
    if(p != pt || force){
        p = pt;
        filters_init();
    }
    //dont check for changes its too cheap
    allp_lens[0] = ZIN0(14);
    allp_lens[1] = ZIN0(15);
    allp_lens[2] = ZIN0(16);
    allp_coef = sc_clip(ZIN0(17),0.0,1.0);
    allpass_init();

}
void EarlyRef_next(EarlyRef* unit,int inNumSamples)
{
    unit->getargs();
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
SCWrapClassT(EarlyRef27);
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
    //SETCALC(EarlyRef27_next);
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
struct EarlyRefGen:public Unit
{
	EarlyRefGen(Unit* unit);
	void CalcFirst();
	float mindist;
    void CalcOne(int n,float exp,float ux,float uy,float uz,float lx,float ly,float lz);
    void refsCalculation();
    void findImage(float ufx,float ufy,float ufz,float lfx,float lfy,float lfz,float * res);
    bool getargs(Unit * unit, bool force=false);
	float L[3];
    float Ps[3];
    float Pr[3];
	float EarL[3],EarR[3];
    float B,HW,d0,Hangle;
    int N;
    float samprate;
    float bufnumL,bufnumR;
    bool let_trig;
    SndBuf *sndbufL, *sndbufR;
    unsigned int framesize,framesize_1, m_pos;
};
SCWrapClass(EarlyRefGen);
EarlyRefGen::EarlyRefGen(Unit *unit)
{
    int ins = 0;
    bufnumL = IN0(ins++);
    bufnumR = IN0(ins++);
    Ps[0] = IN0(ins++) - 1.0;//hack to force getargs
    samprate = FULLRATE;
    sndbufL = GetBuffer(unit,bufnumL);
    sndbufR = GetBuffer(unit,bufnumR);
    if (!sndbufL || !sndbufR){
		Print("EarlyRefGen cant allocate buffers.\n");
        SETCALC(*ClearUnitOutputs);
        unit->mDone = true;
		return;
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
    SETCALC(EarlyRefGen_next);
}
inline void fracdel2(float fsample,float val,float *buf, unsigned int size){
    if(fsample >=size)
        return;
	int pos1 = (int)floor(fsample);
	float frac = fsample - (float)pos1;
	buf[pos1] += val*(1.0-frac);
	buf[pos1 + 1] += val*frac;
}
void EarlyRefGen::CalcFirst()
{
	float sinHW = sin(Hangle)*HW;
	float cosHW = cos(Hangle)*HW;
    EarL[0] = Pr[0] - cosHW;
    EarL[1] = Pr[1] + sinHW;
    EarL[2] = Pr[2];
	float distL = dist(Ps,EarL);
	EarR[0] = Pr[0] + cosHW;
    EarR[1] = Pr[1] - sinHW;
    EarR[2] = Pr[2];
    float distR = dist(Ps,EarR);
	mindist =  std::min(distR,distL);
}
void EarlyRefGen::CalcOne(int n,float exp,float ux,float uy,float uz,float lx,float ly,float lz)
{
    float image[3];
    float dell,delr,ampl,ampr;
    findImage(ux,uy,uz,lx,ly,lz,image);

    float distL = dist(image,EarL);
    float distR = dist(image,EarR);
    float preA = pow(B,exp);
    ampl = preA/(distL + 0.001);//d0*preA/(distL);
    ampr = preA/(distR + 0.001);//d0*preA/(distR);
    dell = samprate*(distL-mindist)/340.0;
    delr = samprate*(distR-mindist)/340.0;
    //printf("%g, %g, %g, %g, %g\n",d0,delL[n],delR[n],ampL[n],ampR[n]);
    /*
    if (dell < sndbufL->samples)
        sndbufL->data[(int)dell] += ampl;
    if (delr < sndbufR->samples)
        sndbufR->data[(int)delr] += ampr;
    */
    fracdel2(dell,ampl,sndbufL->data,sndbufL->frames);
    fracdel2(delr,ampr,sndbufR->data,sndbufR->frames);

}
/*
void EarlyRefGen::predist(u,v,w,l,m,n,L,Ps){
	local x = (2*u-1)*Ps[1] - 2*l*L[1] 
	local y = (2*v-1)*Ps[2] - 2*m*L[2] 
	local z = (2*w-1)*Ps[3] - 2*n*L[3]
	return {x,y,z},-TA{2*u-1,2*v-1,2*w-1},-TA{-2*l,-2*m,-2*n}
}
*/
void EarlyRefGen::refsCalculation(){
    
    Clear(sndbufL->samples,sndbufL->data);
    Clear(sndbufR->samples,sndbufR->data);
    d0 = 1.0;
	CalcFirst();
    if(N == 0){
        //printf("EarlyRefGen::refsCalculation\n");
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
        //printf("EarlyRefGen::refsCalculation N:%d\n",N);
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
void EarlyRefGen::findImage(float ufx,float ufy,float ufz,float lfx,float lfy,float lfz,float * res){
    res[0] = ufx*Ps[0] + lfx*L[0];
    res[1] = ufy*Ps[1] + lfy*L[1];
    res[2] = ufz*Ps[2] + lfz*L[2];
}
bool EarlyRefGen::getargs(Unit *unit,bool force){
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
	float Hangle_t = ZIN0(ins++);
    bool changed = (Pst[0] != Ps[0] || Pst[1] != Ps[1] || Pst[2] != Ps[2] || Prt[0] != Pr[0] || Prt[1] != Pr[1] || Prt[2] != Pr[2] || Lt[0] != L[0] || Lt[1] != L[1] || Lt[2] != L[2] || HWt != HW || Bt != B || Nt != N || Hangle_t!= Hangle);
    if (changed || force) {
        //printf("Nt %d %f \n",Nt,Ntf);
        Ps[0] = Pst[0]; Ps[1] = Pst[1]; Ps[2] = Pst[2]; Pr[0] = Prt[0]; Pr[1] = Prt[1];Pr[2] = Prt[2]; L[0] = Lt[0]; L[1] = Lt[1] ; L[2] = Lt[2] ; HW = HWt ; B = Bt; N = Nt; Hangle = Hangle_t;
        //printf("%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n",Ps[0],Ps[1],Ps[2],Pr[0],Pr[1],Pr[2],L[0],L[1],L[2],HW,B);
        refsCalculation();
    }
    return changed;
}
void EarlyRefGen_next(EarlyRefGen* unit,int inwrongNumSamples)
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
		Print("EarlyRefAtkGen cant allocate buffers.\n");
        SETCALC(*ClearUnitOutputs);
        unit->mDone = true;
		return;
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

extern void initPartConvT(InterfaceTable *);
PluginLoad(DWGReverb){
    /*
    Base<int> der;
    der.inc();
    der.show();
    Base<char> derc;
    derc.inc();
    derc.show();
    */
	ft = inTable;
	
	initPartConvT(inTable);
    DefineDtorUnit(DWGAllpass);
	DefineDtorUnit(Kendall);
    DefineDtorUnit(EarlyRef);
	DefineDtorUnit(EarlyRef27);
    DefineDtorUnit(EarlyRefGen);
    DefineDtorUnit(EarlyRefAtkGen);
    DefineDtorUnit(DWGReverbC1C3);
    DefineDtorUnit(DWGReverbC1C3_16);
    DefineDtorUnit(DWGReverb3Band);
    DefineDtorUnit(DWGReverb3Band_16);
}

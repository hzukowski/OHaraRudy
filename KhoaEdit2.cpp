

#include <math.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <time.h>
#include <string>
#include <sstream>
using std::ifstream;

using namespace std;

// %extracellular ionic concentrations
const double nao = 140.0;
const double cao = 1.8;
const double ko = 5.4;

// %physical constants
const double R = 8314.0;
const double T = 310.0;
const double F = 96485.0;

// %cell geometry
const double L = 0.01;
const double rad = 0.0011;
const double vcell = 1000*3.14*rad*rad*L;
const double Ageo = 2*3.14*rad*rad+2*3.14*rad*L;
const double Acap = 2*Ageo;
const double vmyo = 0.68*vcell;
const double vnsr = 0.0552*vcell;
const double vjsr = 0.0048*vcell;
const double vss = 0.02*vcell;


//// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

const int celltype = 0; //Cell type = (0) for endo; (1) epi (2) for M,

const int num_simulations = 3;
double BCL_values[num_simulations] = {200, 1000, 2000};


const int stimuli = 1000;

//const double Drug = 2.72;  //2.72 (ng/ml) DOF
//const int drugname = Drug;

const double base_dt = 0.005;
const int fold = 1.0/base_dt;


typedef struct cell {
    double CaMKa, CaMKt;
    double m, j, hsp, hp;
    double hf, hs, h, jp;
    double mL, hL, hLp;
    double d, ff, fs, f, fcaf, fcas, fca, jca;
    double ffp, fp, fcafp, fcap, nca, ap;
    double xs1, xs2, xr, xrf, xrs, xk1;
    double a, iF, iS, i, iFp, iSp, ip;
    double nai, nass, ki, kss;
    double cai, cass, cansr, cajsr;
    double Jrelp, Jrelnp;
    double v, v_new, dvdt;
    double I_Kr, C3_Kr, C2_Kr, O_Kr, C1_Kr;
    double ICaL;
    //double In_drug, C3n_drug, C2n_drug, On_drug, C1n_drug;
    //double Ic_drug, C3c_drug, C2c_drug, Oc_drug, C1c_drug;
} Cell;

typedef struct simState {
    double t, tt;
    int tstep, counter, beat;
    Cell cellData;
} SimState;

SimState theState;
SimState *S = &theState;


double Calcu_I_Total(Cell *theCell, double dt, double I_inj);


    int main () {
        
        for (int sim_index = 0; sim_index < num_simulations; ++sim_index) {

        const double BCL = BCL_values[sim_index];
        
        char name[30];
        
      //output files
            
            // Construct output file name using BCL
            std::ostringstream output_name_stream;
            output_name_stream << "Output_CL" << BCL << ".txt";
            std::string output_name = output_name_stream.str();

            ofstream myfile(output_name); // Creating the output file stream with the constructed file name
        
        
        
        double V;
        Cell *theCell;
        double dt = base_dt;
        double I_Total;
        double dv, dvdt;
        double cycle_length;
        double csc, randomN;
        double I_inj;
        
        
        double vmax = 0;
        double tt1, tt2, tt3, tt0, highestV, lowestV;
        double Vm = -86;
        
        
        time_t startTime;
        time_t previousTime;
        
        const time_t timeSave = 2*60*60;
        startTime = time(NULL);
        previousTime = startTime;
        
        //str = "CL =";
        //std::string str = "CL" + std::to_string(BCL);
        
        
        
        
        if ( celltype == 1) {
            cout << " Cell type = EPI " << endl;
        } else if ( celltype == 2) {
            cout << " Cell type = M-cell " << endl;
        }else if ( celltype == 0) {
            cout << " Cell type = ENDO " << endl;
        }
        
        
        
        theCell = &(S->cellData);
        
        
        //// Initial conditions
        theCell->v_new=-87.;
        theCell->nai=7.;
        theCell->nass=theCell->nai;
        theCell->ki=145.;
        theCell->kss=theCell->ki;
        theCell->cai=1.0e-4;
        theCell->cass=theCell->cai;
        theCell->cansr=1.2;
        theCell->cajsr=theCell->cansr;
        theCell->m=0.;
        theCell->hf=1.;
        theCell->hs=1.;
        theCell->j=1.;
        theCell->hsp=1.;
        theCell->jp=1.;
        theCell->mL=0.;
        theCell->hL=1.;
        theCell->hLp=1.;
        theCell->a=0.;
        theCell->iF=1.;
        theCell->iS=1.;
        theCell->ap=0.;
        theCell->iFp=1.;
        theCell->iSp=1.;
        theCell->d=0.;
        theCell->ff=1.;
        theCell->fs=1.;
        theCell->fcaf=1.;
        theCell->fcas=1.;
        theCell->jca=1.;
        theCell->nca=0.;
        theCell->ffp=1.;
        theCell->fcafp=1.;
        theCell->xrf=0.;
        theCell->xrs=0.;
        theCell->xs1=0.;
        theCell->xs2=0.;
        theCell->xk1=1.;
        theCell->Jrelnp=0.;
        theCell->Jrelp=0.;
        theCell->CaMKt=0.;
        theCell->CaMKa=0.;
        
        
        S->counter=0;
        
        S->beat=1;
        S->t = 0.0;
        S->tt = 0.0;
        S->tstep = 0;
        
        tt0 = 6E20;
        tt1 = 6E20;
        tt2 = 6E20;
        tt3 = 6E20;
        
        highestV = -86;
        lowestV = 0;
        dvdt = 0;
        
        
        srand48(startTime);
        
        int done = 0;
        
        while (!done) {
            
            
            cycle_length = BCL;
            
            if (S->t > 20 ) {dt=base_dt;}
            else {dt=base_dt/100;}
            
            
            theCell->v = theCell->v_new;
            
            theCell = &(S->cellData);
            
            
            
            if ( S->t <= 0.5  ) {
                I_inj = -80.0;
            } else {
                I_inj=0.0;
            }
            
            V=theCell->v;
            
            I_Total=Calcu_I_Total(theCell, dt, I_inj);
            
            dv=dt*(-I_Total) ;
            
            theCell->dvdt=fabs(dv/dt);
            theCell->v_new = theCell->v+dv;
            
            
            
            if(  S->beat>900 && ( S->counter % (fold*1) == 0 ) ) {
                myfile << S->tt << "\t" << S->cellData.v << endl;
            }
            
            
            
            S->t += dt;
            S->tt += dt;
            S->counter += 1;
            
            
            
            if(S->t >= cycle_length) {
                
                cout << "now the beat is: " << S->beat << endl;
                S->t = 0;
                S->beat++;
                
                
            }
            if (S->beat > stimuli) { done = 1;}
            
            
        } // end while
        
        myfile.close(); // Close the file
        //fclose ( output );
    }//end for
    
    return 0;
}


// total current
double Calcu_I_Total(Cell *theCell, double dt, double I_inj){
    double v = theCell->v;
    
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    // %CaMK constants
    const double KmCaMK=0.15;
    
    const double aCaMK=0.05;
    const double bCaMK=0.00068;
    const double CaMKo=0.05;
    const double KmCaM=0.0015;
    // %update CaMK
    const double CaMKb=CaMKo*(1.0-theCell->CaMKt)/(1.0+KmCaM/theCell->cass);
    //const double CaMKa=CaMKb+CaMKt;
    theCell->CaMKa=CaMKb+theCell->CaMKt;
    const double dCaMKt=aCaMK*CaMKb*(CaMKb+theCell->CaMKt)-bCaMK*theCell->CaMKt;
    theCell->CaMKt = theCell->CaMKt + dt*dCaMKt;
    
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    // %reversal potentials
    double ENa=(R*T/F)*log(nao/theCell->nai);
    double EK=(R*T/F)*log(ko/theCell->ki);
    double PKNa=0.01833;
    double EKs=(R*T/F)*log((ko+PKNa*nao)/(theCell->ki+PKNa*theCell->nai));
    
    // %convenient shorthand calculations
    double vffrt=v*F*F/(R*T);
    double vfrt=v*F/(R*T);
    
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    // %calculate INa
    
    double mss=1.0/(1.0+exp((-(v+39.57))/9.871));
    double tm=1.0/(6.765*exp((v+11.64)/34.77)+8.552*exp(-(v+77.42)/5.955));
    //double dm=(mss-m)/tm;
    theCell->m = mss-(mss-theCell->m)*exp(-dt/tm);
    
    double hss=1.0/(1+exp((v+82.90)/6.086));
    double thf=1.0/(1.432e-5*exp(-(v+1.196)/6.285)+6.149*exp((v+0.5096)/20.27));
    double ths=1.0/(0.009794*exp(-(v+17.95)/28.05)+0.3343*exp((v+5.730)/56.66));
    double Ahf=0.99;
    double Ahs=1.0-Ahf;
    //double dhf=(hss-hf)/thf;
    //double dhs=(hss-hs)/ths;
    theCell->hf = hss-(hss-theCell->hf)*exp(-dt/thf);
    theCell->hs = hss-(hss-theCell->hs)*exp(-dt/ths);
    theCell->h = Ahf * theCell->hf + Ahs * theCell->hs;
    
    double jss=hss;
    double tj=2.038+1.0/(0.02136*exp(-(v+100.6)/8.281)+0.3052*exp((v+0.9941)/38.45));
    //double dj=(jss-j)/tj;
    theCell->j = jss-(jss-theCell->j)*exp(-dt/tj);
    
    double hssp=1.0/(1+exp((v+89.1)/6.086));
    double thsp=3.0*ths;
    //double dhsp=(hssp-hsp)/thsp;
    theCell->hsp = hssp-(hssp-theCell->hsp)*exp(-dt/thsp);
    theCell->hp=Ahf*theCell->hf+Ahs*theCell->hsp;
    
    double tjp=1.46*tj;
    //double djp=(jss-jp)/tjp;
    theCell->jp = jss-(jss-theCell->jp)*exp(-dt/tjp);
    
    double GNa=75.;
    double fINap=(1.0/(1.0+KmCaMK/theCell->CaMKa));
    
    double INa=GNa*(v-ENa)*theCell->m*theCell->m*theCell->m*((1.0-fINap)*theCell->h*theCell->j+fINap*theCell->hp*theCell->jp);
    
    // %calculate INaL
    double mLss=1.0/(1.0+exp((-(v+42.85))/5.264));
    double tmL=1.0/(6.765*exp((v+11.64)/34.77)+8.552*exp(-(v+77.42)/5.955));
    //double dmL=(mLss-mL)/tmL;
    theCell->mL = mLss-(mLss-theCell->mL)*exp(-dt/tmL);
    
    double hLss=1.0/(1.0+exp((v+87.61)/7.488));
    double thL=200.0;
    //double dhL=(hLss-hL)/thL;
    theCell->hL = hLss-(hLss-theCell->hL)*exp(-dt/thL);
    
    double hLssp=1.0/(1.0+exp((v+93.81)/7.488));
    double thLp=3.0*thL;
    //double dhLp=(hLssp-hLp)/thLp;
    theCell->hLp = hLssp-(hLssp-theCell->hLp)*exp(-dt/thLp);
    
    double GNaL=0.0075;
    if( celltype==1 ) {
        GNaL=GNaL*0.6;
    }
    double fINaLp=(1.0/(1.0+KmCaMK/theCell->CaMKa));
    
    double INaL=GNaL*(v-ENa)*theCell->mL*((1.0-fINaLp)*theCell->hL+fINaLp*theCell->hLp);
    
    // %calculate Ito
    double ass=1.0/(1.0+exp((-(v-14.34))/14.82));
    double ta=1.0515/(1.0/(1.2089*(1.0+exp(-(v-18.4099)/29.3814)))+3.5/(1.0+exp((v+100.0)/29.3814)));
    //double da=(ass-a)/ta;
    theCell->a = ass - (ass - theCell->a)*exp(-dt/ta);
    
    double iss=1.0/(1.0+exp((v+43.94)/5.711));
    double delta_epi;
    if( celltype==1 ) {
        delta_epi=1.0-(0.95/(1.0+exp((v+70.0)/5.0)));
    } else {
        delta_epi=1.0;
    }
    double tiF=4.562+1/(0.3933*exp((-(v+100.0))/100.0)+0.08004*exp((v+50.0)/16.59));
    double tiS=23.62+1/(0.001416*exp((-(v+96.52))/59.05)+1.780e-8*exp((v+114.1)/8.079));
    tiF=tiF*delta_epi;
    tiS=tiS*delta_epi;
    double AiF=1.0/(1.0+exp((v-213.6)/151.2));
    double AiS=1.0-AiF;
    //double diF=(iss-iF)/tiF;
    //double diS=(iss-iS)/tiS;
    theCell->iF = iss - (iss - theCell->iF)*exp(-dt/tiF);
    theCell->iS = iss - (iss - theCell->iS)*exp(-dt/tiS);
    theCell->i = AiF * theCell->iF + AiS * theCell->iS;
    
    double assp=1.0/(1.0+exp((-(v-24.34))/14.82));
    //double dap=(assp-ap)/ta;
    theCell->ap = assp - (assp - theCell->ap)*exp(-dt/ta);
    
    double dti_develop=1.354+1.0e-4/(exp((v-167.4)/15.89)+exp(-(v-12.23)/0.2154));
    double dti_recover=1.0-0.5/(1.0+exp((v+70.0)/20.0));
    double tiFp=dti_develop*dti_recover*tiF;
    double tiSp=dti_develop*dti_recover*tiS;
    //double diFp=(iss-iFp)/tiFp;
    //double diSp=(iss-iSp)/tiSp;
    theCell->iFp = iss - (iss - theCell->iFp)*exp(-dt/tiFp);
    theCell->iSp = iss - (iss - theCell->iSp)*exp(-dt/tiSp);
    theCell->ip = AiF * theCell->iFp + AiS * theCell->iSp;
    
    double Gto=0.02;
    if( celltype==1 ){
        Gto=Gto*4.0;
    } else if( celltype==2 ){
        Gto=Gto*4.0;
    }
    double fItop=(1.0/(1.0+KmCaMK/theCell->CaMKa));
    
    double Ito=Gto*(v-EK)*((1.0-fItop)*theCell->a*theCell->i+fItop*theCell->ap*theCell->ip);
    
    // %calculate ICaL, ICaNa, ICaK
    double dss=1.0/(1.0+exp((-(v+3.940))/4.230));
    double td=0.6+1.0/(exp(-0.05*(v+6.0))+exp(0.09*(v+14.0)));
    //double dd=(dss-d)/td;
    theCell->d = dss-(dss-theCell->d)*exp(-dt/td);
    
    double fss=1.0/(1.0+exp((v+19.58)/3.696));
    double tff=7.0+1.0/(0.0045*exp(-(v+20.0)/10.0)+0.0045*exp((v+20.0)/10.0));
    double tfs=1000.0+1.0/(0.000035*exp(-(v+5.0)/4.0)+0.000035*exp((v+5.0)/6.0));
    double Aff=0.6;
    double Afs=1.0-Aff;
    //double dff=(fss-ff)/tff;
    //double dfs=(fss-fs)/tfs;
    theCell->ff = fss-(fss-theCell->ff)*exp(-dt/tff);
    theCell->fs = fss-(fss-theCell->fs)*exp(-dt/tfs);
    theCell->f = Aff * theCell->ff + Afs * theCell->fs;
    
    double fcass=fss;
    double tfcaf=7.0+1.0/(0.04*exp(-(v-4.0)/7.0)+0.04*exp((v-4.0)/7.0));
    double tfcas=100.0+1.0/(0.00012*exp(-v/3.0)+0.00012*exp(v/7.0));
    double Afcaf=0.3+0.6/(1.0+exp((v-10.0)/10.0));
    double Afcas=1.0-Afcaf;
    //double dfcaf=(fcass-fcaf)/tfcaf;
    //double dfcas=(fcass-fcas)/tfcas;
    theCell->fcaf = fcass - (fcass - theCell->fcaf)*exp(-dt/tfcaf);
    theCell->fcas = fcass - (fcass - theCell->fcas)*exp(-dt/tfcas);
    theCell->fca = Afcaf * theCell->fcaf + Afcas * theCell->fcas;
    
    double tjca=75.0;
    //double djca=(fcass-jca)/tjca;
    theCell->jca = fcass - (fcass - theCell->jca)*exp(-dt/tjca);
    
    double tffp=2.5*tff;
    //double dffp=(fss-ffp)/tffp;
    theCell->ffp = fss - (fss - theCell->ffp)*exp(-dt/tffp);
    theCell->fp = Aff * theCell->ffp + Afs * theCell->fs;
    
    double tfcafp=2.5*tfcaf;
    //double dfcafp=(fcass-fcafp)/tfcafp;
    theCell->fcafp = fcass - (fcass - theCell->fcafp)*exp(-dt/tfcafp);
    theCell->fcap = Afcaf * theCell->fcafp + Afcas * theCell->fcas;
    
    double Kmn=0.002;
    double k2n=1000.0;
    double km2n=theCell->jca*1.0;
    double anca=1.0 / ( k2n/km2n+ pow((1.0+Kmn/theCell->cass),4) );
    //double dnca=anca*k2n-nca*km2n;
    theCell->nca = anca*(k2n/km2n)-(anca*(k2n/km2n)-theCell->nca)*exp(-km2n*dt);
    
    double PhiCaL=4.0*vffrt*(theCell->cass*exp(2.0*vfrt)-0.341*cao)/(exp(2.0*vfrt)-1.0);
    double PhiCaNa=1.0*vffrt*(0.75*theCell->nass*exp(1.0*vfrt)-0.75*nao)/(exp(1.0*vfrt)-1.0);
    double PhiCaK=1.0*vffrt*(0.75*theCell->kss*exp(1.0*vfrt)-0.75*ko)/(exp(1.0*vfrt)-1.0);
    double zca=2.0;
    double PCa=0.0001;
    if( celltype==1 ) {
        PCa=PCa*1.2;
    } else if( celltype==2 ) {
        PCa=PCa*2.5;
    }
    
    double PCap=1.1*PCa;
    double PCaNa=0.00125*PCa;
    double PCaK=3.574e-4*PCa;
    double PCaNap=0.00125*PCap;
    double PCaKp=3.574e-4*PCap;
    double fICaLp=(1.0/(1.0+KmCaMK/theCell->CaMKa));
    double ICaL=(1.0-fICaLp)*PCa*PhiCaL*theCell->d*(theCell->f*(1.0-theCell->nca)+theCell->jca*theCell->fca*theCell->nca)+fICaLp*PCap*PhiCaL*theCell->d*(theCell->fp*(1.0-theCell->nca)+theCell->jca*theCell->fcap*theCell->nca);
    theCell->ICaL = ICaL;
    double ICaNa=(1.0-fICaLp)*PCaNa*PhiCaNa*theCell->d*(theCell->f*(1.0-theCell->nca)+theCell->jca*theCell->fca*theCell->nca)+fICaLp*PCaNap*PhiCaNa*theCell->d*(theCell->fp*(1.0-theCell->nca)+theCell->jca*theCell->fcap*theCell->nca);
    double ICaK=(1.0-fICaLp)*PCaK*PhiCaK*theCell->d*(theCell->f*(1.0-theCell->nca)+theCell->jca*theCell->fca*theCell->nca)+fICaLp*PCaKp*PhiCaK*theCell->d*(theCell->fp*(1.0-theCell->nca)+theCell->jca*theCell->fcap*theCell->nca);
    
    
      
    
    //added section below
    double xrss=1.0/(1.0+exp((-(v+8.337))/6.789));
    double txrf=12.98+1.0/(0.3652*exp((v-31.66)/3.869)+4.123e-5*exp((-(v-47.78))/20.38));
    double txrs=1.865+1.0/(0.06629*exp((v-34.70)/7.355)+1.128e-5*exp((-(v-29.74))/25.94));
    double Axrf=1.0/(1.0+exp((v+54.81)/38.21));
    double Axrs=1.0-Axrf;
    double xrf=xrss-(xrss-xrf)*exp(-dt/txrf);
    double xrs=xrss-(xrss-xrs)*exp(-dt/txrs);
    double xr=Axrf*xrf+Axrs*xrs;
    double rkr=1.0/(1.0+exp((v+55.0)/75.0))*1.0/(1.0+exp((v-10.0)/30.0));
    //
    
    double GKr = 0.024 * 1.7575;
    
    if (celltype == 1)
    {
        GKr *= 1.3;
    }
    if (celltype == 2)
    {
        GKr *= 0.8;
    }
    
    //    if (fabs (sum2 -1.0) > 0.001) {cout << "Error in Sum at t = " << S->t << endl;}
    
    //double IKr = GKr*sqrt(ko/5.4)*theCell->O_Kr*(v-EK);
    double IKr = GKr*sqrt(ko/5.4)*rkr*(v-EK);
    theCell->xr = IKr;
    
    
    //%calculate IKs
    double xs1ss=1.0/(1.0+exp((-(v+11.60))/8.932));
    double txs1=817.3+1.0/(2.326e-4*exp((v+48.28)/17.80)+0.001292*exp((-(v+210.0))/230.0));
    //double dxs1=(xs1ss-xs1)/txs1;
    theCell->xs1 = xs1ss - (xs1ss - theCell->xs1)*exp(-dt/txs1);
    
    double xs2ss=xs1ss;
    double txs2=1.0/(0.01*exp((v-50.0)/20.0)+0.0193*exp((-(v+66.54))/31.0));
    //double dxs2=(xs2ss-xs2)/txs2;
    theCell->xs2 = xs2ss - (xs2ss - theCell->xs2)*exp(-dt/txs2);
    
    double KsCa=1.0+0.6/(1.0+ pow((3.8e-5/theCell->cai),1.4) );
    double GKs=0.0034;
    if( celltype==1 ) {
        GKs=GKs*1.4;
    }
    
    double IKs=GKs*KsCa*theCell->xs1*theCell->xs2*(v-EKs);
    
    // %calculate IK1
    double xk1ss=1.0/(1.0+exp(-(v+2.5538*ko+144.59)/(1.5692*ko+3.8115)));
    double txk1=122.2/(exp((-(v+127.2))/20.36)+exp((v+236.8)/69.33));
    //double dxk1=(xk1ss-xk1)/txk1;
    theCell->xk1 = xk1ss - (xk1ss - theCell->xk1)*exp(-dt/txk1);
    
    double rk1=1.0/(1.0+exp((v+105.8-2.6*ko)/9.493));
    double GK1=0.1908;
    if( celltype==1 ) {
        GK1=GK1*1.2;
    } else if( celltype==2 ) {
        GK1=GK1*1.3;
    }
    
    double IK1=GK1*sqrt(ko)*rk1*theCell->xk1*(v-EK);
    
    // %calculate INaCa_i
    double kna1=15.0;
    double kna2=5.0;
    double kna3=88.12;
    double kasymm=12.5;
    double wna=6.0e4;
    double wca=6.0e4;
    double wnaca=5.0e3;
    double kcaon=1.5e6;
    double kcaoff=5.0e3;
    double qna=0.5224;
    double qca=0.1670;
    double hca=exp((qca*v*F)/(R*T));
    double hna=exp((qna*v*F)/(R*T));
    double h1=1+theCell->nai/kna3*(1+hna);
    double h2=(theCell->nai*hna)/(kna3*h1);
    double h3=1.0/h1;
    double h4=1.0+theCell->nai/kna1*(1+theCell->nai/kna2);
    double h5=theCell->nai*theCell->nai/(h4*kna1*kna2);
    double h6=1.0/h4;
    double h7=1.0+nao/kna3*(1.0+1.0/hna);
    double h8=nao/(kna3*hna*h7);
    double h9=1.0/h7;
    double h10=kasymm+1.0+nao/kna1*(1.0+nao/kna2);
    double h11=nao*nao/(h10*kna1*kna2);
    double h12=1.0/h10;
    double k1=h12*cao*kcaon;
    double k2=kcaoff;
    double k3p=h9*wca;
    double k3pp=h8*wnaca;
    double k3=k3p+k3pp;
    double k4p=h3*wca/hca;
    double k4pp=h2*wnaca;
    double k4=k4p+k4pp;
    double k5=kcaoff;
    double k6=h6*theCell->cai*kcaon;
    double k7=h5*h2*wna;
    double k8=h8*h11*wna;
    double x1=k2*k4*(k7+k6)+k5*k7*(k2+k3);
    double x2=k1*k7*(k4+k5)+k4*k6*(k1+k8);
    double x3=k1*k3*(k7+k6)+k8*k6*(k2+k3);
    double x4=k2*k8*(k4+k5)+k3*k5*(k1+k8);
    double E1=x1/(x1+x2+x3+x4);
    double E2=x2/(x1+x2+x3+x4);
    double E3=x3/(x1+x2+x3+x4);
    double E4=x4/(x1+x2+x3+x4);
    double KmCaAct=150.0e-6;
    double allo=1.0/(1.0+ pow((KmCaAct/theCell->cai),2));
    double zna=1.0;
    double JncxNa=3.0*(E4*k7-E1*k8)+E3*k4pp-E2*k3pp;
    double JncxCa=E2*k2-E1*k1;
    double Gncx=0.0008;
    if( celltype==1 ) {
        Gncx=Gncx*1.1;
    } else if( celltype==2 ) {
        Gncx=Gncx*1.4;
    }
    double INaCa_i=0.8*Gncx*allo*(zna*JncxNa+zca*JncxCa);
    
    // %calculate INaCa_ss
    h1=1+theCell->nass/kna3*(1+hna);
    h2=(theCell->nass*hna)/(kna3*h1);
    h3=1.0/h1;
    h4=1.0+theCell->nass/kna1*(1+theCell->nass/kna2);
    h5=theCell->nass*theCell->nass/(h4*kna1*kna2);
    h6=1.0/h4;
    h7=1.0+nao/kna3*(1.0+1.0/hna);
    h8=nao/(kna3*hna*h7);
    h9=1.0/h7;
    h10=kasymm+1.0+nao/kna1*(1+nao/kna2);
    h11=nao*nao/(h10*kna1*kna2);
    h12=1.0/h10;
    k1=h12*cao*kcaon;
    k2=kcaoff;
    k3p=h9*wca;
    k3pp=h8*wnaca;
    k3=k3p+k3pp;
    k4p=h3*wca/hca;
    k4pp=h2*wnaca;
    k4=k4p+k4pp;
    k5=kcaoff;
    k6=h6*theCell->cass*kcaon;
    k7=h5*h2*wna;
    k8=h8*h11*wna;
    x1=k2*k4*(k7+k6)+k5*k7*(k2+k3);
    x2=k1*k7*(k4+k5)+k4*k6*(k1+k8);
    x3=k1*k3*(k7+k6)+k8*k6*(k2+k3);
    x4=k2*k8*(k4+k5)+k3*k5*(k1+k8);
    E1=x1/(x1+x2+x3+x4);
    E2=x2/(x1+x2+x3+x4);
    E3=x3/(x1+x2+x3+x4);
    E4=x4/(x1+x2+x3+x4);
    KmCaAct=150.0e-6;
    allo=1.0/(1.0+ pow((KmCaAct/theCell->cass),2));
    JncxNa=3.0*(E4*k7-E1*k8)+E3*k4pp-E2*k3pp;
    JncxCa=E2*k2-E1*k1;
    double INaCa_ss=0.2*Gncx*allo*(zna*JncxNa+zca*JncxCa);
    
    double k1p=949.5;
    double k1m=182.4;
    double k2p=687.2;
    double k2m=39.4;
    k3p=1899.0;
    double k3m=79300.0;
    k4p=639.0;
    double k4m=40.0;
    double Knai0=9.073;
    double Knao0=27.78;
    double delta=-0.1550;
    double Knai=Knai0*exp((delta*v*F)/(3.0*R*T));
    double Knao=Knao0*exp(((1.0-delta)*v*F)/(3.0*R*T));
    double Kki=0.5;
    double Kko=0.3582;
    double MgADP=0.05;
    double MgATP=9.8;
    double Kmgatp=1.698e-7;
    double H=1.0e-7;
    double eP=4.2;
    double Khp=1.698e-7;
    double Knap=224.0;
    double Kxkur=292.0;
    double P=eP/(1.0+H/Khp+theCell->nai/Knap+theCell->ki/Kxkur);
    double a1=(k1p * pow((theCell->nai/Knai),3) )/( pow((1.0+theCell->nai/Knai),3) + pow((1.0+theCell->ki/Kki),2) - 1.0);
    double b1=k1m*MgADP;
    double a2=k2p;
    double b2=(k2m * pow((nao/Knao),3) )/( pow((1.0+nao/Knao),3) + pow((1.0+ko/Kko),2) - 1.0);
    double a3=(k3p * pow((ko/Kko),2) )/( pow((1.0+nao/Knao),3) + pow((1.0+ko/Kko),2) - 1.0);
    double b3=(k3m*P*H)/(1.0+MgATP/Kmgatp);
    double a4=(k4p*MgATP/Kmgatp)/(1.0+MgATP/Kmgatp);
    double b4=(k4m * pow((theCell->ki/Kki),2) )/( pow((1.0+theCell->nai/Knai),3) + pow((1.0+theCell->ki/Kki),2) - 1.0);
    x1=a4*a1*a2+b2*b4*b3+a2*b4*b3+b3*a1*a2;
    x2=b2*b1*b4+a1*a2*a3+a3*b1*b4+a2*a3*b4;
    x3=a2*a3*a4+b3*b2*b1+b2*b1*a4+a3*a4*b1;
    x4=b4*b3*b2+a3*a4*a1+b2*a4*a1+b3*b2*a1;
    E1=x1/(x1+x2+x3+x4);
    E2=x2/(x1+x2+x3+x4);
    E3=x3/(x1+x2+x3+x4);
    E4=x4/(x1+x2+x3+x4);
    double zk=1.0;
    double JnakNa=3.0*(E1*a3-E2*b3);
    double JnakK=2.0*(E4*b1-E3*a1);
    double Pnak=30;
    if( celltype==1 ) {
        Pnak=Pnak*0.9;
    } else if( celltype==2 ) {
        Pnak=Pnak*0.7;
    }
    double INaK=Pnak*(zna*JnakNa+zk*JnakK);
    
    // %calculate IKb
    double xkb=1.0/(1.0+exp(-(v-14.48)/18.34));
    double GKb=0.003;
    if( celltype==1 ) {
        GKb=GKb*0.6;
    }
    double IKb=GKb*xkb*(v-EK);
    
    // %calculate INab
    double PNab=3.75e-10;
    double INab=PNab*vffrt*(theCell->nai*exp(vfrt)-nao)/(exp(vfrt)-1.0);
    
    // %calculate ICab
    double PCab=2.5e-8;
    double ICab=PCab*4.0*vffrt*(theCell->cai*exp(2.0*vfrt)-0.341*cao)/(exp(2.0*vfrt)-1.0);
    
    // %calculate IpCa
    double GpCa=0.0005;
    double IpCa=GpCa*theCell->cai/(0.0005+theCell->cai);
    
    //%%%%%%%%%%%%%%%%%%%%
    double I_total = INa+INaL+Ito+ICaL+ICaNa+ICaK+IKr+IKs+IK1+INaCa_i+INaCa_ss+INaK+INab+IKb+IpCa+ICab+I_inj;
    //%%%%%%%%%%%%%%%%%%%%%
    
    // %calculate diffusion fluxes
    double JdiffNa=(theCell->nass-theCell->nai)/2.0;
    double JdiffK=(theCell->kss-theCell->ki)/2.0;
    double Jdiff=(theCell->cass-theCell->cai)/0.2;
    
    // %calculate ryanodione receptor calcium induced calcium release from the jsr
    double bt=4.75;
    double a_rel=0.5*bt;
    double Jrel_inf=a_rel*(-ICaL)/(1.0 + pow((1.5/theCell->cajsr),8) );
    if( celltype==2 ){
        Jrel_inf=Jrel_inf*1.7;
    }
    double tau_rel=bt/(1.0+0.0123/theCell->cajsr);
    
    if( tau_rel<0.001 ) {
        tau_rel=0.001;
    }
    
    //double dJrelnp=(Jrel_inf-Jrelnp)/tau_rel;
    theCell->Jrelnp = Jrel_inf - (Jrel_inf - theCell->Jrelnp)*exp(-dt/tau_rel);
    
    double btp=1.25*bt;
    double a_relp=0.5*btp;
    double Jrel_infp=a_relp*(-ICaL)/(1.0 + pow((1.5/theCell->cajsr),8) );
    if( celltype==2 ) {
        Jrel_infp=Jrel_infp*1.7;
    }
    double tau_relp=btp/(1.0+0.0123/theCell->cajsr);
    
    if( tau_relp<0.001 ) {
        tau_relp=0.001;
    }
    
    //double dJrelp=(Jrel_infp-Jrelp)/tau_relp;
    theCell->Jrelp = Jrel_infp - (Jrel_infp - theCell->Jrelp)*exp(-dt/tau_relp);
    
    double fJrelp=(1.0/(1.0+KmCaMK/theCell->CaMKa));
    double Jrel = (1.0-fJrelp)*theCell->Jrelnp+fJrelp*theCell->Jrelp;
    // %calculate serca pump, ca uptake flux
    double Jupnp=0.004375*theCell->cai/(theCell->cai+0.00092);
    double Jupp=2.75*0.004375*theCell->cai/(theCell->cai+0.00092-0.00017);
    if( celltype==1 ) {
        Jupnp=Jupnp*1.3;
        Jupp=Jupp*1.3;
    }
    double fJupp=(1.0/(1.0+KmCaMK/theCell->CaMKa));
    double Jleak=0.0039375*theCell->cansr/15.0;
    double Jup=(1.0-fJupp)*Jupnp+fJupp*Jupp-Jleak;
    
    // %calculate tranlocation flux
    double Jtr=(theCell->cansr-theCell->cajsr)/100.0;
    
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    // %calcium buffer constants
    double cmdnmax=0.05;
    if( celltype==1 ) {
        cmdnmax=cmdnmax*1.3;
    }
    double kmcmdn=0.00238;
    double trpnmax=0.07;
    double kmtrpn=0.0005;
    double BSRmax=0.047;
    double KmBSR=0.00087;
    double BSLmax=1.124;
    double KmBSL=0.0087;
    double csqnmax=10.0;
    double kmcsqn=0.8;
    
    // %update intracellular concentrations, using buffers for cai, cass, cajsr
    double dnai=-(INa+INaL+3.0*INaCa_i+3.0*INaK+INab)*Acap/(F*vmyo)+JdiffNa*vss/vmyo;
    double dnass=-(ICaNa+3.0*INaCa_ss)*Acap/(F*vss)-JdiffNa;
    theCell->nai = theCell->nai + dt*dnai;
    theCell->nass = theCell->nass + dt*dnass;
    
    double dki=-(Ito+IKr+IKs+IK1+IKb+I_inj-2.0*INaK)*Acap/(F*vmyo)+JdiffK*vss/vmyo;
    double dkss=-(ICaK)*Acap/(F*vss)-JdiffK;
    theCell->ki = theCell->ki + dt*dki;
    theCell->kss = theCell->kss + dt*dkss;
    
    double Bcai=1.0/(1.0 + cmdnmax * kmcmdn / pow((kmcmdn+theCell->cai),2) + trpnmax * kmtrpn / pow((kmtrpn+theCell->cai),2) );
    double dcai=Bcai*(-(IpCa+ICab-2.0*INaCa_i)*Acap/(2.0*F*vmyo)-Jup*vnsr/vmyo+Jdiff*vss/vmyo);
    theCell->cai = theCell->cai + dt*dcai;
    double Bcass=1.0/(1.0 + BSRmax * KmBSR / pow((KmBSR+theCell->cass),2) + BSLmax * KmBSL / pow((KmBSL+theCell->cass),2) );
    double dcass=Bcass*(-(ICaL-2.0*INaCa_ss)*Acap/(2.0*F*vss)+Jrel*vjsr/vss-Jdiff);
    theCell->cass = theCell->cass + dt*dcass;
    
    double dcansr=Jup-Jtr*vjsr/vnsr;
    theCell->cansr = theCell->cansr + dt*dcansr;
    
    double Bcajsr=1.0/(1.0 + csqnmax * kmcsqn / pow((kmcsqn+theCell->cajsr),2) );
    double dcajsr=Bcajsr*(Jtr-Jrel);
    theCell->cajsr = theCell->cajsr + dt*dcajsr;
    
    
    
    return I_total;
}



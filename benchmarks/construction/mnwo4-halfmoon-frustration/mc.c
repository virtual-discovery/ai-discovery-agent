/* Classical MC for H = -sum_{i!=j} J_ij e_i.e_j + D sum_i (e_i.n)^2  (unit vectors, ordered pairs).
   usage: mc bonds.txt L1 L2 L3 T_meV D_meV nx ny nz neq nmeas every seed qfile out_sq [restart_in|-]
   bonds.txt: lines "i j d1 d2 d3 J_meV" (sublattice 0/1, cell offsets)  (both directions listed)
   qfile: nq lines "h k l"  (r.l.u.), positions from TAU below. Output: per q accumulated <Sxx>,<Syy>,<Szz> (1/N^2 |m(q)|^2)
   also prints mean energy per spin and specific heat. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <complex.h>
static unsigned long long s[2];
static inline unsigned long long rotl(const unsigned long long x,int k){return (x<<k)|(x>>(64-k));}
static inline unsigned long long nxt(void){unsigned long long s0=s[0],s1=s[1],r=s0+s1;s1^=s0;s[0]=rotl(s0,55)^s1^(s1<<14);s[1]=rotl(s1,36);return r;}
static inline double U(void){return (nxt()>>11)*(1.0/9007199254740992.0);}
int L1,L2,L3,N,nb; int *bi,*bj,*b1,*b2,*b3; double *bJ;
int *nbr; double *nJ; int *nn;  /* neighbour lists */
double *e; double D,nvec[3],T;
const double TAU[2][3]={{0.5,0.68490,0.25},{0.5,0.31510,0.75}};
int idx(int s_,int n1,int n2,int n3){n1=((n1%L1)+L1)%L1;n2=((n2%L2)+L2)%L2;n3=((n3%L3)+L3)%L3;return s_+2*(n1+L1*(n2+L2*n3));}
static void field(int i,double h[3]){h[0]=h[1]=h[2]=0;for(int k=0;k<nn[i];k++){int j=nbr[i*64+k];double J=nJ[i*64+k];h[0]+=2*J*e[3*j];h[1]+=2*J*e[3*j+1];h[2]+=2*J*e[3*j+2];}}
/* local energy part: E_i = -e.h + D (e.n)^2 , h = 2 sum_j J e_j */
static void randvec(double v[3]){double z=2*U()-1,p=2*M_PI*U(),r=sqrt(1-z*z);v[0]=r*cos(p);v[1]=r*sin(p);v[2]=z;}
double cone=0.5; long acc=0,tri=0;
static void metro(int i){double h[3];field(i,h);double *ei=e+3*i;double v[3];randvec(v);
  double n_[3]={ei[0]+cone*v[0],ei[1]+cone*v[1],ei[2]+cone*v[2]};double nr=sqrt(n_[0]*n_[0]+n_[1]*n_[1]+n_[2]*n_[2]);n_[0]/=nr;n_[1]/=nr;n_[2]/=nr;
  double po=ei[0]*nvec[0]+ei[1]*nvec[1]+ei[2]*nvec[2],pn=n_[0]*nvec[0]+n_[1]*nvec[1]+n_[2]*nvec[2];
  double dE=-( (n_[0]-ei[0])*h[0]+(n_[1]-ei[1])*h[1]+(n_[2]-ei[2])*h[2]) + D*(pn*pn-po*po);
  tri++; if(dE<=0||U()<exp(-dE/T)){ei[0]=n_[0];ei[1]=n_[1];ei[2]=n_[2];acc++;}}
static void overrelax(int i){double h[3];field(i,h);double hh=h[0]*h[0]+h[1]*h[1]+h[2]*h[2];if(hh<1e-20)return;double *ei=e+3*i;
  double d=(ei[0]*h[0]+ei[1]*h[1]+ei[2]*h[2])/hh;double n_[3]={2*d*h[0]-ei[0],2*d*h[1]-ei[1],2*d*h[2]-ei[2]};
  double po=ei[0]*nvec[0]+ei[1]*nvec[1]+ei[2]*nvec[2],pn=n_[0]*nvec[0]+n_[1]*nvec[1]+n_[2]*nvec[2];
  double dE=D*(pn*pn-po*po); if(dE<=0||U()<exp(-dE/T)){ei[0]=n_[0];ei[1]=n_[1];ei[2]=n_[2];}}
double energy(void){double E=0;for(int i=0;i<N;i++){double h[3];field(i,h);double *ei=e+3*i;double p=ei[0]*nvec[0]+ei[1]*nvec[1]+ei[2]*nvec[2];E+=-0.5*(ei[0]*h[0]+ei[1]*h[1]+ei[2]*h[2])+D*p*p;}return E/N;}
int main(int argc,char**argv){
  if(argc<17){fprintf(stderr,"usage\n");return 1;}
  FILE*f=fopen(argv[1],"r");L1=atoi(argv[2]);L2=atoi(argv[3]);L3=atoi(argv[4]);T=atof(argv[5]);D=atof(argv[6]);
  nvec[0]=atof(argv[7]);nvec[1]=atof(argv[8]);nvec[2]=atof(argv[9]);{double r=sqrt(nvec[0]*nvec[0]+nvec[1]*nvec[1]+nvec[2]*nvec[2]);for(int c=0;c<3;c++)nvec[c]/=r;}
  int neq=atoi(argv[10]),nmeas=atoi(argv[11]),every=atoi(argv[12]);unsigned long long seed=strtoull(argv[13],0,10);s[0]=seed*0x9E3779B97F4A7C15ULL+1;s[1]=seed^0xD1B54A32D192ED03ULL;for(int k=0;k<20;k++)nxt();
  char*qfile=argv[14],*out=argv[15],*rin=argv[16];
  int cap=1000;bi=malloc(cap*sizeof(int));bj=malloc(cap*sizeof(int));b1=malloc(cap*sizeof(int));b2=malloc(cap*sizeof(int));b3=malloc(cap*sizeof(int));bJ=malloc(cap*sizeof(double));nb=0;
  while(fscanf(f,"%d %d %d %d %d %lf",&bi[nb],&bj[nb],&b1[nb],&b2[nb],&b3[nb],&bJ[nb])==6)nb++;fclose(f);
  N=2*L1*L2*L3;e=malloc(3*N*sizeof(double));nbr=malloc(64*N*sizeof(int));nJ=malloc(64*N*sizeof(double));nn=calloc(N,sizeof(int));
  for(int n3=0;n3<L3;n3++)for(int n2=0;n2<L2;n2++)for(int n1=0;n1<L1;n1++)for(int k=0;k<nb;k++){int i=idx(bi[k],n1,n2,n3),j=idx(bj[k],n1+b1[k],n2+b2[k],n3+b3[k]);if(bJ[k]==0)continue;nbr[i*64+nn[i]]=j;nJ[i*64+nn[i]]=bJ[k];nn[i]++;}
  if(strcmp(rin,"-")!=0){FILE*g=fopen(rin,"r");for(int i=0;i<N;i++)if(fscanf(g,"%lf %lf %lf",&e[3*i],&e[3*i+1],&e[3*i+2])!=3){fprintf(stderr,"restart short\n");return 1;}fclose(g);}
  else for(int i=0;i<N;i++)randvec(e+3*i);
  /* q points */
  FILE*qf=fopen(qfile,"r");int nq=0,qc=1024;double*qs=malloc(3*qc*sizeof(double));
  while(fscanf(qf,"%lf %lf %lf",&qs[3*nq],&qs[3*nq+1],&qs[3*nq+2])==3){nq++;if(nq==qc){qc*=2;qs=realloc(qs,3*qc*sizeof(double));}}fclose(qf);
  double*Sq=calloc(3*nq,sizeof(double));
  double complex *p1=malloc(L1*sizeof(double complex)),*p2=malloc(L2*sizeof(double complex)),*p3=malloc(L3*sizeof(double complex));
  int *order=malloc(N*sizeof(int));for(int i=0;i<N;i++)order[i]=i;
  double Es=0,E2s=0;int nsamp=0;
  for(int sweep=0;sweep<neq+nmeas;sweep++){
    for(int i=0;i<N;i++){int j=(int)(U()*N);metro(j);}
    for(int r=0;r<2;r++)for(int i=0;i<N;i++)overrelax(i);
    if(sweep<neq && sweep%50==49){double a=(double)acc/tri; if(a<0.4)cone*=0.8; else if(a>0.6) cone*=1.2; if(cone>2)cone=2; acc=tri=0;}
    if(sweep>=neq && (sweep-neq)%every==0){
      double E=energy();Es+=E;E2s+=E*E;nsamp++;
      for(int iq=0;iq<nq;iq++){double h=qs[3*iq],k=qs[3*iq+1],l=qs[3*iq+2];
        for(int a=0;a<L1;a++)p1[a]=cexp(2*M_PI*I*h*a);for(int a=0;a<L2;a++)p2[a]=cexp(2*M_PI*I*k*a);for(int a=0;a<L3;a++)p3[a]=cexp(2*M_PI*I*l*a);
        double complex m[3]={0,0,0};
        for(int s_=0;s_<2;s_++){double complex ps=cexp(2*M_PI*I*(h*TAU[s_][0]+k*TAU[s_][1]+l*TAU[s_][2]));
          double complex ms[3]={0,0,0};
          for(int n3=0;n3<L3;n3++){double complex m3[3]={0,0,0};
            for(int n2=0;n2<L2;n2++){double complex m2[3]={0,0,0};int base=2*L1*(n2+L2*n3)+s_;
              for(int n1=0;n1<L1;n1++){int i=base+2*n1;double complex ph=p1[n1];m2[0]+=ph*e[3*i];m2[1]+=ph*e[3*i+1];m2[2]+=ph*e[3*i+2];}
              for(int c=0;c<3;c++)m3[c]+=p2[n2]*m2[c];}
            for(int c=0;c<3;c++)ms[c]+=p3[n3]*m3[c];}
          for(int c=0;c<3;c++)m[c]+=ps*ms[c];}
        for(int c=0;c<3;c++)Sq[3*iq+c]+=creal(m[c]*conj(m[c]))/((double)N*N);
      }
    }
  }
  FILE*o=fopen(out,"w");fprintf(o,"# T_meV %g D %g L %d %d %d nsamp %d E %.8f C %.8f\n",T,D,L1,L2,L3,nsamp,Es/nsamp,(E2s/nsamp-(Es/nsamp)*(Es/nsamp))*N/(T*T));
  for(int iq=0;iq<nq;iq++)fprintf(o,"%.6e %.6e %.6e\n",Sq[3*iq]/nsamp,Sq[3*iq+1]/nsamp,Sq[3*iq+2]/nsamp);fclose(o);
  char rp[512];snprintf(rp,512,"%s.spins",out);FILE*rs=fopen(rp,"w");for(int i=0;i<N;i++)fprintf(rs,"%.8f %.8f %.8f\n",e[3*i],e[3*i+1],e[3*i+2]);fclose(rs);
  printf("E %.8f C %.8f acc_cone %.3f\n",Es/nsamp,(E2s/nsamp-(Es/nsamp)*(Es/nsamp))*N/(T*T),cone);
  return 0;}

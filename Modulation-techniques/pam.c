#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>

#define NoE 10

#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])

gsl_rng *r,*rg;
const gsl_rng_type *T,*TG;

//genrates grey codes for constellation diagram
void generateGreyCodes(char A[][5], int n, int len){
    for(int i=0;i<n;i++){
        int num=i;
        for(int j=len-1;j>=0;j--){
            A[i][j]=num%2+'0';
            num=num/2;
        }
        A[i][len]='\0';
    }
}

//integration function
float trapz_integral(float f[], float sm[]){
	int i,n;
	float x0 = 0.01,xn = 1.0,h = 0.01 ,so,ans=0.0;
	n=(xn-x0)/h;
	float y[1001];
	if(n%2==1)
		n=n+1;
	h=(xn-x0)/n;
	for(i=0; i<=n; i++){
		y[i]=(sm[i]*f[i]);
	}
	so=0;
	for(i=1; i<n; i++){
			so=so+y[i];
	}
	ans=(y[0]+y[n]+so)/4.8;
	return ans;
}

//generation of raised cosine pulse
float generateRaisedCosine(float *rc){
	float T=1.0,a=0.25;
    float t=0.01,Eg,exc=T/a;
    exc=exc/2.0;
    int t2=(int)(exc*100)-1;
    double data[2*700];
    float index[700];
    for(int t1=0;t1<700;t1++){
		   	rc[t1]=sin((3.14*t)/T)*cos(3.14*a*t/T);
		    float temp = (2*a*t);
		    temp=temp/T ;
		    temp=1-pow(temp,2);
		    rc[t1]=rc[t1]/temp;
		    rc[t1]=(rc[t1]*T)/(3.14*t);
		   	if(t1==t2) {
		   		rc[t1]=(3.14/(4*T))*sin(3.14/(2*a));
		   		rc[t1]=rc[t1]/(3.14/(2*a));
		   	}
		   	REAL(data,t1)=rc[t1];
		   	IMAG(data,t1) = 0.0;
		   	float temp1=rc[t1];
		    if(temp1<0){
		        temp1=temp1*-1;
		    }
		    if(t1<128){
		   		index[t1]=t1;
		   		REAL(data,t1)=rc[t1];
			   	IMAG(data,t1) = 0.0;	
		   	}
		   	Eg+=pow(temp1,2);
        t=t+0.01;
    }
    //fft
    gsl_fft_complex_radix2_forward(data,1,128);
    //plotting
    char * commandsForGnuplot[] = {"set title \"Raised Cosine Spectrum\"","plot 'data.temp'"};
			FILE * temp = fopen("data.temp", "w");
			FILE * gnuplotPipe = popen ("gnuplot -persistent", "w");
			for (int i=0; i < 128; i++){
				fprintf(temp, "%f %e \n", index[i], REAL(data,i));
			}
			for (int i=0; i < 2; i++){
				fprintf(gnuplotPipe, "%s \n", commandsForGnuplot[i]);
			}
    return Eg;
}

//return nearest point to the estimated point in PAM
int pamNearestPoint(float ans,int M){
	if(ans>(M-1))
		return (M-1);
	else if(ans<(-1*(M-1)))
		return (-1*(M-1));
	else{
		int points[2];
		points[0]=floor(ans);
		points[1]=ceil(ans);
		if(points[0]%2==0){
			if(ans>0)
				points[0]=points[0]-1;
			else
				points[0]=points[0]+1;
		}
		if(points[1]%2==0){
			if(ans>0)
				points[1]=points[1]+1;
			else
				points[1]=points[1]-1;
		}
		float d1=sqrt(pow((ans-points[0]),2));
		float d2=sqrt(pow((ans-points[1]),2));
		if(d1<d2)
			return points[0];
		else
			return points[1];
	}
}
int main(){

	gsl_rng_env_setup();

	T = gsl_rng_default;
  	TG = gsl_rng_default;

  	r = gsl_rng_alloc(T);
  	rg = gsl_rng_alloc(TG);
  	int M;
  	scanf("%d",&M);
  	int k=log(M)/log(2);
  	char A[M][5];
  	generateGreyCodes(A,M,k);
  	int j,count,message,iter=0,f = 20,snr_arr[10],flag=0;
    float *rc = (float*)malloc(1000*sizeof(float));
    float Eg = generateRaisedCosine(rc);
    float sm[2001],f1[2001],N0,std_noise,prob[10];
    int cons2[M];
    int *cons1;
	generateGreyCodes(A,M,k);
	//varying snr
    for(int snr=1;snr<=10;snr++)
    { 
        count=0;
        N0=(2*Eg)/snr;
        snr_arr[snr-1]=snr;
        std_noise = sqrt(N0/2);
        for(iter=0;;iter++)
        {
            char bits[k+1];
            int x;
            //constellation point generation
            if(M==16){
			int temp[16]={-5,-7,-3,-1,-11,-9,-13,-15,5,7,3,1,11,9,13,15};
			cons1=&temp[0];
		}
		else if(M==8){
			int temp[8]={-7,-5,-1,-3,7,5,1,3};
			cons1=&temp[0];
		}
		else if(M==4){
			int temp[4]={-3,-1,3,1};
			cons1=&temp[0];
		}
		else if(M==2){
			int temp[2]={-1,1};
			cons1=&temp[0];
		}
		for(j=0;j<M;j++){
			cons2[j]=cons1[j];	
		}
		//plotting constellation
		if(flag==0){
					char * commandsForGnuplot1[] = {"set title \"PAM constellation points\"","set xrange [-16:16]","set yrange [-1:1]" ,"plot 'cons.temp'"};
			FILE * temp1 = fopen("cons.temp", "w");
			FILE * gnuplotPipe1 = popen ("gnuplot -persistent", "w");
			int z=0;
			for (int i=0; i < M; i++){
				fprintf(temp1, "%d %d \n", cons2[i],z);
			}
			for (int i=0; i < 4; i++){
				fprintf(gnuplotPipe1, "%s \n", commandsForGnuplot1[i]);
    		}
    		flag=1;
		}
			//bits generation
            for(j=0;j<k;j++)
            {
                message = gsl_rng_uniform_int(r,2);
                bits[j] = message+'0';
            }
            bits[j] = '\0';
            	//mapping bits to contellation
		        for(j = 0;j < M;j++)
			    {
			        if(strcmp(bits,A[j]) == 0)
			        {
			            x = cons1[j];
			            break;
			        }
				} 
				//gaussian noise generation
		        float gn=gsl_ran_gaussian(r,std_noise);
		        float t=0.01;
		        for(int t1=0;t1<700;t1++)
		        {
		        	//basis function generation
		            f1[t1]=sqrt(2/Eg) * rc[t1] * cos(2*3.14*f*t);
		            //modulated signal
	            sm[t1]=x*rc[t1]*cos(2*3.14*f*t);
	            	//adding gaussian noise                    
		            sm[t1]+=gn;
		            t=t+0.01;
		        }
		        //demodulation using correlation demodulator
		        float ans=trapz_integral(f1,sm);
		        //applying ML rule
		        int est = pamNearestPoint(ans,M);
		        //counting errors		        
		        if(x!=est){
		            count++;
		        }
		        if(count == NoE)
		            break;
        }
        printf("**********count = %d in iterations = %d at SNR = %d dB\n",NoE,iter,snr);
        printf("**********probability of error is %0.14f\n", (float)NoE/iter);
        prob[snr-1]=(float)NoE/iter;
    }
	//plotting BER vs Eb/N0
	char * commandsForGnuplot[] = {"set title \"BER vs Eb/N0(in dB)\"", "plot 'prob.temp'"};
    FILE * temp = fopen("prob.temp", "w");
    FILE * gnuplotPipe = popen ("gnuplot -persistent", "w");
    for (int i=0; i < 10; i++){
    	fprintf(temp, "%d %0.14f \n", snr_arr[i], prob[i]);
    }
    for (int i=0; i < 2; i++){
    	fprintf(gnuplotPipe, "%s \n", commandsForGnuplot[i]);
    }
    return 0;
}

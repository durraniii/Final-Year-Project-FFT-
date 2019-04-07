
#include <iostream>
#include <cmath>
#include <ctime>
#include <fstream>
#include<iomanip>
#include <cstdlib>
#include<thread>
using namespace std;
const int N=256;
thread mythreads[N];

bool ifft=false;
bool i_im=false;
float xr[N], xi[N] ,sum[N], xr2[N][N],xi2[N][N],orignal[N][N],arr[N][N],real[N][N],imag[N][N],re[N][N],im[N][N],outreal[N][N],outimag[N][N], arr1[N][N];
void fft(float *a,char ch){

int j, k;


float tempi[N], tempr[N], i2, temp, z, complexR[N], complexI[N], z1, z2;


	//bit reversal
	i2 = N >> 1;
	j = 0;
	temp = 0;
	for (int i = 0; i < N - 1; i++) {
		if (i < j) {
			temp = a[i];
			a[i] = a[j];
			a[j] = temp;
		}
		k = i2;
		while (k <= j) {
			j -= k;
			k >>= 1;
		}
		j += k;
	}

	k = 0;
	int n, s = 0;
	int e = 1;
	int clip1, clip2;

	//FFT
	for (int stage = 1; stage <= log2(N); stage++)
	{
		for (int b = 0; b < N / 2; b++)	//b=butterfly
		{
			n = (int)(b / (int)pow(2, stage - 1));
			k = (int)(b % (int)pow(2, stage - 1));
			clip1 = k + n*pow(2, stage);
			clip2 = k + n*pow(2, stage) + pow(2,stage-1);

			if (stage == 1)
			{
			    if(ch == 'i'){
                       
				tempi[clip1] = a[clip1];
				tempr[clip1] = 0;
				tempi[clip2] = a[clip2];
				tempr[clip2] = 0;
				}
				else{

                                tempr[clip1] = a[clip1];
				tempi[clip1] = 0;
				tempr[clip2] = a[clip2];
				tempi[clip2] = 0;
				}
			}
			if (stage > 1)
			{
				tempr[clip1]  = xr[clip1];
				tempi[clip1]  = xi[clip1];
				tempr[clip2] = xr[clip2];
				tempi[clip2] = xi[clip2];
			}
			if(ifft==true)
			{

			z1 = cos(2. *acos(-1) * k / pow(2, stage));
			z2 = sin(2. *acos(-1) * k / pow(2, stage));

			}

			else
			{
			z1 =  cos(2. *acos(-1) * k / pow(2, stage));
			z2 =  -sin(2. *acos(-1) * k / pow(2, stage));
			}
			complexR[clip2] = z1 * tempr[clip2] - z2 * tempi[clip2];
			complexI[clip2] = z2 * tempr[clip2] + z1 * tempi[clip2];

			xr[clip1] = tempr[clip1] + complexR[clip2];
			xi[clip1] = tempi[clip1] + complexI[clip2];

			xr[clip2] = tempr[clip1] - complexR[clip2];
			xi[clip2] = tempi[clip1] - complexI[clip2];

		}

	}


	}
void fftt(float arr[N][N],int c,char ch,char datatype){

float a[N];

if(ch == 'c'){
    for(int i=0;i<N;i++){
        a[i]=arr[i][c];

    }


}
else{
   for(int i=0;i<N;i++){
        a[i]=arr[c][i];

        }


}
if(datatype == 'I')
    fft(a,'i');
else
fft(a,'r');
if(ch == 'c'){

  for(int i=0;i<N;i++){



    xr2[i][c]=xr[i];
    xi2[i][c]=xi[i];


}

}
else{

    for(int i=0;i<N;i++){

    xr2[c][i]=xr[i];
    xi2[c][i]=xi[i];
}}

	}
void ifftt(float r[N][N],float im[N][N],int c,char ch){
float real[N],imag[N];
thread t,t1;
if(ch == 'c'){

      for(int i=0;i<N;i++){
        real[i]=r[i][c];
        imag[i]=im[i][c];
       	}
}
else{
      for(int i=0;i<N;i++){
        real[i]=r[c][i];
        imag[i]=0;
        }

}


fft(real,'r');
for (int i = 0; i < N; i++)
	{
	    xr[i]=xr[i]/N;
             xi[i]=xi[i]/N;
            sum[i]=xr[i] + xi[i];
	}
fft(imag,'i');
	
	

    for (int i = 0; i < N; i++){
	    xr[i]=xr[i]/N;
        xi[i]=xi[i]/N;
        sum[i]+=xr[i] + xi[i];
	}

	if(ch == 'c'){
	for(int i=0;i<N;i++)
	{
        orignal[i][c]=sum[i];
		}
		}
        else{
            	for(int i=0;i<N;i++)
                    orignal[c][i]=sum[i];
            }
}
void inline fun1(float re[N][N],float im[N][N],int i){
 

       ifft=true;
   	ifftt(re,im,i,'c');
       ifft=false;
     
   
}
void  ifft2(float re[N][N],float im[N][N]){
    //inverse fft2 col wise maybe
    
    for(int i=0; i < N ;i++)
        mythreads[i]=thread{fun1,re,im,i};
    for(int i=0;i<N;i++)
           mythreads[i].join();

		
       
//inverse fft2 row wise
    for(int i=0;i<N;i++){

        mythreads[i]=thread{ifftt,orignal,im,i,'r'};
    }
    for(int i=0;i<N;i++)
       mythreads[i].join();

}
void  fft2_image(float arr[N][N]){
   
    for(int col=0;col<N;col++){
        mythreads[col]=thread {fftt,arr,col,'c','R'};

		}

 for(int i=0;i<N;i++)
        mythreads[i].join();
   
for(int i=0;i<N;i++){
 
    for(int j=0;j<N;j++){

            real[i][j]=xr2[i][j];
            imag[i][j]=xi2[i][j];}

    }

 
       //fft row wise real data
        for(int row=0;row < N;row++)
            mythreads[row]=thread{fftt,real,row,'r','R'};
    for(int i=0;i<N;i++)
        mythreads[i].join();
        
        for(int i=0;i<N;i++){
               for(int j=0;j<N;j++){

            re[i][j]=xr2[i][j];
            im[i][j]=xi2[i][j];
			                   }
}


///    fft row wise imag data
        for(int row=0;row< N ;row++)
            mythreads[row]=thread{fftt,imag,row,'r','I'};

    for(int i=0;i<N;i++)
        mythreads[i].join();
        
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){

            re[i][j]+=xr2[i][j];
            im[i][j]+=xi2[i][j];}
            }

 
  }
void  fft2_filter(float arr[N][N],int size){

 
    for(int col=0;col<size;col++){
        mythreads[col]=thread{fftt,arr,col,'c','R'};

		}
    for(int i=0; i < size; i++)
            mythreads[i].join();
    
    
for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){

            real[i][j]=xr2[i][j];
            imag[i][j]=xi2[i][j];}

    }

    
       //fft row wise real data
        for(int row=0;row < N;row++)
           mythreads[row]=thread{fftt,real,row,'r','R'};
     for(int i=0;i<N;i++)
        mythreads[i].join();
       
           
        for(int i=0;i<N;i++){
                for(int j=0;j<N;j++){

            re[i][j]=xr2[i][j];
            im[i][j]=xi2[i][j];
			                   }
}


///    fft row wise imag data
        for(int row=0;row< N ;row++)
            mythreads[row]=thread{fftt,imag,row,'r','I'};
 for(int i=0;i<N;i++)
        mythreads[i].join();
   
    for(int i=0;i<N;i++){
    
    
        for(int j=0;j<N;j++){

            re[i][j]+=xr2[i][j];
            im[i][j]+=xi2[i][j];}
            
    }

 
  }
void complex(float r[N][N],float ii[N][N],float r1[N][N],float i1[N][N]){
for(int i=0;i<N;i++){
	for(int j=0;j<N;j++)
		{
		real[i][j]=((r[i][j]*r1[i][j])-(ii[i][j]*i1[i][j]));
		imag[i][j]=((r[i][j]*i1[i][j])+(ii[i][j]*r1[i][j]));}}
}
void FFT(float image[N][N],float filter[N][N]){

fft2_image(image);
fft2_filter(filter,N/10);


complex(re,im,re,im);
ifft2(real,imag);
}

int main(){
   
    
for(int i=0;i< N;i++){

    for(int j=0;j< N;j++){
       
            arr[i][j]=i*j;
       
            }


}
    float timer=0;
clock_t t=clock();
    for(int i=0; i < 1 ; i++){


 FFT(arr,arr);
    }
 t = clock() - t; 
 timer += ((float)t)/CLOCKS_PER_SEC; // in seconds 
 cout << timer << endl;

}

#include <bits/stdc++.h>
#include <stdlib.h>
#include <omp.h>
using namespace std ;

void print(double** matrix, int n){
	for(int i = 0; i < n; i++){
	   for(int j = 0; j < n; j++){
		   cout << matrix[i][j] << "  ";
	   }
	   cout<< endl;
   }
   cout<<endl;
	
}

double verify(int n,double **a, double** l, double** u, vector<int> pvector){
	double **p; double **result;
	p = new double* [n]; result = new double* [n];
	for (int i=0; i < n; i++){
		p[i] = new double[n];
		result[i] = new double[n];
	}
	for(int i =0 ; i < n; i++){
		for(int j = 0; j < n; j++){
			p[i][j] = 0.0;
			result[i][j] = 0.0;
		}
	}
	for(int i = 0; i < n; i++){
		p[i][pvector[i]] = 1.0;
	}
	
	double tmp = 0.0;
	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			tmp = 0.0;
			for(int k = 0; k < n; k++){
				tmp = tmp + p[i][k]*a[k][j] - l[i][k]*u[k][j]; 
			}
			result[i][j] = tmp;
		}
	}
	double norm=0.0,colnorm = 0.0;
	for(int i = 0; i < n; i++){
		colnorm = 0.0;
		for(int j = 0; j < n; j++){
			colnorm += result[i][j]*result[i][j];
		}
		norm += sqrt(colnorm);
	}
	
	return norm;
	
}

int main(int argc, char *argv[]) {
   // printf() displays the string inside quotation
    string arg1 = argv[1] ;
	int n = stoi(arg1);
	int kdash;
	//printf("n = %d", n);
	double **a;
	a = new double* [n];
	double **a_orig;
	a_orig = new double* [n];
	double **u;
	u = new double* [n];
	double **l;
	l = new double* [n];

	auto start = std::chrono::high_resolution_clock::now();	
#pragma omp parallel for	
	for (int i=0; i < n; i++){
		a[i] = new double[n];
		u[i] = new double[n];
		l[i] = new double[n];
		a_orig[i] = new double[n];
	}
   
   
	vector<int> p(n,0);
	double drand48();
	double maxelem;

	for(int i =0; i < n; i++){
		p[i] = i;		
		for(int j = 0; j < n; j++){
			if( j < i ){
				u[i][j] = 0.0;
				l[i][j] = drand48();
			}
			else if(j == i){
				l[i][j] = 1.0;
				u[i][j] = drand48();
			}
			else{
				u[i][j] = drand48();
				l[i][j] = 0.0;
			}
			a[i][j] = drand48();
			a_orig[i][j] = a[i][j];		   
		}  
	}
   
 cout<<"Time Taken in Initialisation: "<<chrono::duration_cast<chrono::microseconds>(std::chrono::high_resolution_clock::now()-start).count()/1000000.0<<endl ;  
	double sum=0.0 ;
	for(int k = 0; k < n; k++){
		maxelem = 0.0;
	auto startfor = std::chrono::high_resolution_clock::now();	
	
#pragma omp parallel reduction(max:maxelem)
		for(int i = k; i < n; i++){
			if(fabs(a[i][k])> maxelem){
				maxelem = fabs(a[i][k]);
				kdash = i;
			}
		}
		
		if(maxelem==0.0){
			cerr<<"Singular matrix."<<endl;
		}
		
		//swap in p 
		int temp = p[k];
		p[k] = p[kdash];
		p[kdash] = temp;
		
		//swap rows in A
		double* tmp = a[k];
		a[k] = a[kdash];
		a[kdash] = tmp;
		   
		double tmpdbl=0.0;

#pragma omp parallel for private(tmpdbl) if(k>500)
		for(int i = 0; i < k; i++){
			tmpdbl = l[k][i];
			l[k][i] = l[kdash][i];
			l[kdash][i] = tmpdbl;
		}

		u[k][k] = a[k][k];
		
#pragma omp parallel for if(k<n-500)
		for(int i = k+1; i < n; i++){
		   l[i][k] = a[i][k]*1.0/u[k][k];
		   u[k][i] = a[k][i];
		}
		
	int j=0;
	double utemp[n-k-1] ;
	
#pragma omp parallel for
	for(int i=0;i<n-k-1;i++){
		utemp[i]=u[k][i+k+1] ;
	}

	double *tempdub ;
	double *tempdub2;
	double val,val2 ;
	int ei=0 ;
	
	auto endfor = std::chrono::high_resolution_clock::now();	
	sum+=chrono::duration_cast<chrono::microseconds>(endfor-startfor).count() ;
	

// #pragma omp parallel for private(j,val,tempdub,tempdub2,val2) lastprivate(ei)
	// for(int i = k+1; i < n-1; i+=2){
		// ei=i ;
		// val=l[i][k] ;
		// val2=l[i+1][k] ;
		// tempdub=a[i] ;
		// tempdub2=a[i+1] ;
	   // for(j = k+1; j<n; j++){
		   // tempdub[j] -= val*utemp[j-k-1];	
		   // tempdub2[j] -= val2*utemp[j-k-1];	
	   // }
	// }
	
	// if(ei==n-3||k+1==n-1){
		// ei=n-1 ;	
		// tempdub=a[ei] ;
		// val=l[ei][k] ;
		// for(j = k+1; j<n; j++){
		   // tempdub[j] -= val*utemp[j-k-1];	
		// }
	// }
	
#pragma omp parallel for private(j,val,tempdub)
	for(int i = k+1; i < n; i++){
		val=l[i][k] ;
		tempdub=a[i] ;
	   for(j = k+1; j<n; j++){
		   tempdub[j] -= val*utemp[j-k-1];	
	   }
	}
}

	auto end = std::chrono::high_resolution_clock::now();
	cout<<"Time Taken: "<<chrono::duration_cast<chrono::microseconds>(end-start).count()/1000000.0<<endl ;
	cout<<"sum: "<<sum/1000000.0<<endl ;
   //print(a,n);
   //print(l,n);
   //print(u,n);
    //double result = verify(n,a_orig,l,u,p);
   // cout<<result<<endl;
   
   return 0;
}

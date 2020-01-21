#include <bits/stdc++.h>
#include <pthread.h>
using namespace std ;

double **a;
double **a_orig;
double **u;
double **l;
int n;
int cores;

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
				//cout<<p[i][k]*a[k][j] - l[i][k]*u[k][j]<<" " ; 
			}
			result[i][j] = tmp;
			// cout<<tmp<<" " ;
		}
		// cout<<endl ;
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

void* assignParallel(void* pair_id_k){
	int id = (*(pair<int,int> *)pair_id_k).first ;
	int k = (*(pair<int,int> *)pair_id_k).second ;
	int total=n-1-k ;
	int start = (k+1)+id*total/cores ;
	int end = min((k+1)+(id+1)*total/cores,n) ;
	for(int i=start;i<end;i++){		
		for(int j=k+1;j<n;j++){
			a[i][j] -= l[i][k]*u[k][j];
		}
	}
}

int main(int argc, char *argv[]) {
   // printf() displays the string inside quotation
    string arg1= argv[1] ;
    string arg2= argv[2] ;
	cores= stoi(arg1);
	n = stoi(arg2);
	int kdash;
	//printf("n = %d", n);
	a = new double* [n];
	a_orig = new double* [n];
	u = new double* [n];
	l = new double* [n];

	auto start = std::chrono::high_resolution_clock::now();	
	for (int i=0; i < n; i++){
		a[i] = new double[n];
		u[i] = new double[n];
		l[i] = new double[n];
		a_orig[i] = new double[n];
	}
   
   
	vector<int> p(n,0);
	double drand48();
	double max;

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
   
   
for(int k = 0; k < n; k++){
	max = 0.0;
		
	for(int i = k; i < n; i++){
		if(fabs(a[i][k])> max){
			max = fabs(a[i][k]);
			kdash = i;
		}
	}
	
	if(max==0.0){
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

	for(int i = 0; i < k; i++){
		tmpdbl = l[k][i];
		l[k][i] = l[kdash][i];
		l[kdash][i] = tmpdbl;
	}

	u[k][k] = a[k][k];
	
	for(int i = k+1; i < n; i++){
	   l[i][k] = a[i][k]*1.0/u[k][k];
	   u[k][i] = a[k][i];
	}

	pthread_t thread[cores] ;
	pair<int,int> arm[cores] ;
	for(int i=0;i<cores;i++){
		arm[i]= make_pair(i,k) ;
	}
	
	for(int i=0;i<cores;i++){
		pthread_create(&thread[i], NULL, &assignParallel, (void*)&arm[i]); 
	}
	for(int i=0;i<cores;i++)
		pthread_join(thread[i],NULL) ;

}

	auto end = std::chrono::high_resolution_clock::now();
	cout<<"Time Taken: "<<chrono::duration_cast<chrono::microseconds>(end-start).count()/1000000.0<<endl ;
   //print(a,n);
   //print(l,n);
   //print(u,n);
   //double result = verify(n,a_orig,l,u,p);
   //cout<<result<<endl;

	for (int i=0; i < n; i++){
		delete a[i] ;
		delete u[i] ;
		delete l[i] ;
		delete a_orig[i] ;
	}
	delete a,u,l,a_orig ;
   return 0;
}

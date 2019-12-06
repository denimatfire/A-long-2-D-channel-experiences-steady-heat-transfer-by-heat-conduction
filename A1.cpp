//Program to Solve Laplace Equation with Given Boundary Conditions
#include<iostream>
#include<fstream>
#include<iomanip>
#include<math.h>
#include<time.h>
#define L 1.0
#define M 200
#define N 200

using namespace std;
	const float dx=1.0/50.0;
	const float dy=1.0/50.0;
	const int m=((2*L)/dx)+1, n=((3*L)/(2*dy))+1; 
	const int	i2=(((L/2.0)-(dx/2.0))/dx)+2, i3=(((3*L/2.0)-(dx/2.0))/dx)+0.5, j2=((L-(dy/2.0))/dy)+1.5;

	double x[M][N], y[M][N], T[M][N], Tnew[M][N], Told[M][N], res_l[M][N], res_u[M][N], res_l_new[M][N], res_u_new[M][N];
	double ap[M][N], ae[M][N], aw[M][N], an[M][N], as[M][N], b[M][N], Ls[M][N], Lw[M][N], Lp[M][N], Ue[M][N], Un[M][N], S[M][N], delta[M][N];
	double r[M][N], r_star[M][N], r_tilda[M][N], P[M][N], P_tilda[M][N], t[M][N], t_tilda[M][N];
	
	
	double sum, maxerr, maxerrt=0.0, w, err_l, err_u, convergence=1.0e-8, sum_l, sum_u, alpha, alpha_bicg, samay, var;
	double sum_nr, sum_dr, delta1, t_nr, t_dr, rho, bet, dt=1.0e-3, tx=dt/(dx*dx), ty=dt/(dy*dy);
	
	float beta=pow((dy/dx),2), Tb=1.0;
	int nl=0, ml=0, nu=0, mu=0, itr=1;	

//..................................................Grid Generation.......................................................	
void grid()
{
	for(int j=1; j<n; j++)
	for(int i=1; i<m; i++)
		{
			x[i][j]=((i-1.0)*dx)+(dx/2.0);
			y[i][j]=((j-1.0)*dy)+(dy/2.0);
		}
		
	ofstream gridout;
	gridout.open("Computational_Domain.plt");
	gridout<<"ZONE "<<"I="<<m-1<<" J="<<n-1<<endl;
	for(int j=1; j<n; j++)
	for(int i=1; i<m; i++)
		gridout<<x[i][j]<<" "<<y[i][j]<<" "<<T[i][j]<<endl;
	gridout.close();
}

//..................................................Print Solution..............................................................
void print_solution()
{
	ofstream filout;
	filout.open("BiCG_Pre.plt");

	filout<<"ZONE "<<"I="<<ml<<" J="<<nl<<endl;
	for(int j=1; j<j2; j++)			// Lower Block
	for(int i=i2; i<=i3; i++)
		filout<<x[i][j]<<" "<<y[i][j]<<" "<<T[i][j]<<endl;
	
	filout<<"ZONE "<<"I="<<mu<<" J="<<nu<<endl;	
	for(int j=j2; j<n; j++)		// Upper Block
	for(int i=1; i<m; i++)
		filout<<x[i][j]<<" "<<y[i][j]<<" "<<T[i][j]<<endl;
	filout.close();
}

//..................................................Assign Default Coefficients.........................................................
void default_coefficients()
{
	for(int j=1; j<j2; j++)
	for(int i=i2; i<=i3; i++)
	{
		ap[i][j]=-2.0*(1+beta);
		ae[i][j]=beta;
		aw[i][j]=beta;
		an[i][j]=1.0;
		as[i][j]=1.0;
		b[i][j]=0.0;
	}
	for(int j=j2; j<n; j++)
	for(int i=1; i<m; i++)
	{
		ap[i][j]=-2.0*(1+beta);
		ae[i][j]=beta;
		aw[i][j]=beta;
		an[i][j]=1.0;
		as[i][j]=1.0;
		b[i][j]=0.0;
	}
	cout<<"\nDefault Coefficients Assigned !!";
}
//..................................................Assign Default Coefficients-Approach 1.........................................................
void default_coefficients_app1()
{
	for(int j=1; j<j2; j++)
	for(int i=i2; i<=i3; i++)
	{
		ap[i][j] = 1+tx+ty;
		ae[i][j] = -tx/2.0;
		aw[i][j] = -tx/2.0;
		an[i][j] = -ty/2.0;
		as[i][j] = -ty/2.0;
		b[i][j] = ((1-tx-ty)*T[i][j])+(T[i+1][j]*tx/2.0)+(T[i-1][j]*tx/2.0)+(T[i][j+1]*ty/2.0)+(T[i][j-1]*ty/2.0);
	}
	for(int j=j2; j<n; j++)
	for(int i=1; i<m; i++)
	{
		ap[i][j] = 1+tx+ty;
		ae[i][j] = -tx/2.0;
		aw[i][j] = -tx/2.0;
		an[i][j] = -ty/2.0;
		as[i][j] = -ty/2.0;
		b[i][j] = ((1-tx-ty)*T[i][j])+(T[i+1][j]*tx/2.0)+(T[i-1][j]*tx/2.0)+(T[i][j+1]*ty/2.0)+(T[i][j-1]*ty/2.0);
	}
	cout<<"\nDefault Coefficients Assigned !!";
}

//..................................................Assign Boundary Coefficients.........................................................
void boundary_coefficients()
{
	for(int j=1; j<j2; j++)			// Lower Block
	for(int i=i2; i<=i3; i++)
	{
		if(j!=1)
		{
			if(i==i2)
			{
				ap[i][j]=ap[i][j]-aw[i][j];
				b[i][j]=b[i][j]-(2*aw[i][j]*Tb);
				aw[i][j]=0.0;
			}
			else if(i==i3)
			{
				ap[i][j]=ap[i][j]-ae[i][j];
				b[i][j]=b[i][j]-(2*ae[i][j]*Tb);
				ae[i][j]=0.0;
			}
		}
		else if(j==1)
		{
			if((i>i2)&(i<i3))
			{
				ap[i][j]=ap[i][j]+as[i][j];
				as[i][j]=0.0;
			}
			else if(i==i2)
			{
				ap[i][j]=ap[i][j]-aw[i][j]+as[i][j];
				b[i][j]=b[i][j]-(2*aw[i][j]*Tb);
				as[i][j]=0.0;
				aw[i][j]=0.0;
			}
			else if(i==i3)
			{
				ap[i][j]=ap[i][j]-ae[i][j]+as[i][j];
				b[i][j]=b[i][j]-(2*ae[i][j]*Tb);
				as[i][j]=0.0;
				ae[i][j]=0.0;
			}
		}
	}
	
	for(int j=j2; j<n; j++)		// Upper Block
	for(int i=1; i<m; i++)
	{
		if((j!=j2)&(j!=n-1))
		{
			if(i==1)
			{
				ap[i][j]=ap[i][j]+aw[i][j];
				aw[i][j]=0.0;
			}
			else if(i==m-1)
			{
				ap[i][j]=ap[i][j]+ae[i][j];
				ae[i][j]=0.0;
			}
		}
		else if(j==j2)
		{
			if((i<i2)&(i>1))
			{
				ap[i][j]=ap[i][j]-as[i][j];
				b[i][j]=b[i][j]-(2*as[i][j]*Tb);
				as[i][j]=0.0;
			}
			else if((i<m-1)&(i>i3))
			{
				ap[i][j]=ap[i][j]-as[i][j];
				b[i][j]=b[i][j]-(2*as[i][j]*Tb);
				as[i][j]=0.0;
			}
			else if(i==1)
			{
				ap[i][j]=ap[i][j]-as[i][j]+aw[i][j];
				b[i][j]=b[i][j]-(2*as[i][j]*Tb);
				as[i][j]=0.0;
				aw[i][j]=0.0;
			}
			else if(i==m-1)
			{
				ap[i][j]=ap[i][j]-as[i][j]+ae[i][j];
				b[i][j]=b[i][j]-(2*as[i][j]*Tb);
				as[i][j]=0.0;
				ae[i][j]=0.0;
			}
		}
		else if(j==n-1)
		{
			if((i>1)&(i<m-1))
			{
				ap[i][j]=ap[i][j]+(an[i][j]*(1-dy));
				
				an[i][j]=0.0;
			}
			else if(i==1)
			{
				ap[i][j]=ap[i][j]+aw[i][j]+(an[i][j]*(1-dy));
				
				an[i][j]=0.0;
				aw[i][j]=0.0;
			}
			else if(i==m-1)
			{
				ap[i][j]=ap[i][j]+ae[i][j]+(an[i][j]*(1-dy));
				
				an[i][j]=0.0;
				ae[i][j]=0.0;
			}
		}
	}
	cout<<"\nBoundary Coefficients Assigned !!";
}

//....................................................Solution Initialization.....................................................
void initialize()
{
	for(int j=1; j<j2; j++)
	for(int i=i2; i<=i3; i++)
	{
		T[i][j] = 0.1;
		Told[i][j] = 0.1;
	}
	for(int j=j2; j<n; j++)
	for(int i=1; i<m; i++)
	{
		T[i][j] = 0.1;
		Told[i][j] = 0.1;
	}
	for(int j=0; j<M; j++)
	for(int i=0; i<N; i++)
	{
		Ls[i][j] = 0.0;
		Lw[i][j] = 0.0;
		Lp[i][j] = 0.0;
		Ue[i][j] = 0.0;
		Un[i][j] = 0.0;
		
		ap[i][j]=0.0;
		ae[i][j]=0.0;
		aw[i][j]=0.0;
		an[i][j]=0.0;
		as[i][j]=0.0;
		b[i][j]=0.0;
		S[i][j]=0.0;
		delta[i][j]=0.0;

		r[i][j] = 0.0;
		r_star[i][j] = 0.0;
		r_tilda[i][j] = 0.0;
		P[i][j] = 0.0;
		P_tilda[i][j] = 0.0;
		t[i][j] = 0.0;
		t_tilda[i][j] = 0.0;
	}	
}

//................................................Cell Count.....................................................
void cell_count()
{
	nl=0, ml=0, nu=0, mu=0;

	for(int j=1; j<j2; j++) nl=nl+1;
	for(int i=i2; i<=i3; i++) ml=ml+1;
	for(int j=j2; j<n; j++)	 nu=nu+1;
	for(int i=1; i<m; i++) mu=mu+1;	
	
	
}


//....................................................Error Calculation.........................................................
double error()
{
 	sum_l=0.0, sum_u=0.0, maxerr=0.0;
	
	// Lower Block
	for(int j=1; j<j2; j++)			
	for(int i=i2; i<=i3; i++)
		sum_l=sum_l+pow(res_l_new[i][j],2);
			sum_l=sum_l/(ml*nl);
			err_l=sqrt(sum_l);
			
	// Upper Block
	for(int j=j2; j<n; j++)		
	for(int i=1; i<m; i++)
		sum_u=sum_u+pow(res_u_new[i][j],2);
			sum_u=sum_u/(mu*nu);
			err_u=sqrt(sum_u);
		
	if(err_l>=err_u)
		maxerr=err_l;
	else
		maxerr=err_u;
			
	return maxerr;
}

//................................................Gauss Siedel Iteration..........................................................
void gauss_siedel()
{
	cout<<"\nw : ";
	cin>>w;
	
	default_coefficients();
	
	boundary_coefficients();
	
	ofstream errout;
	errout.open("convergence_Gauss.txt");
	
	
	do
	{
		// Lower Block
		for(int j=1; j<j2; j++)		
		for(int i=i2; i<=i3; i++)
		{
			res_l_new[i][j]=b[i][j]-((as[i][j]*T[i][j-1])+(aw[i][j]*T[i-1][j])+(ap[i][j]*T[i][j])+(ae[i][j]*T[i+1][j])+(an[i][j]*T[i][j+1]));
			T[i][j]=T[i][j]+((w*res_l_new[i][j])/ap[i][j]);
		}
		
		// Upper Block
		for(int j=j2; j<n; j++)		
		for(int i=1; i<m; i++)
		{
			res_u_new[i][j]=b[i][j]-((as[i][j]*T[i][j-1])+(aw[i][j]*T[i-1][j])+(ap[i][j]*T[i][j])+(ae[i][j]*T[i+1][j])+(an[i][j]*T[i][j+1]));
			T[i][j]=T[i][j]+((w*res_u_new[i][j])/ap[i][j]);
		}
		
		maxerr=error();
		cout<<"\nitr: "<<itr<<" error: "<<maxerr;
		errout<<fixed<<setprecision(16)<<log(itr)<<" "<<log(maxerr)<<endl;
		itr++;	
		
	}while(maxerr>convergence);
	
	errout.close();
	cout<<"\nGauss Siedel operation done !!";
}

//.......................................................Jacobi..................................................
void jacobi()
{
	default_coefficients();
	
	boundary_coefficients();
	
	ofstream errout;
	errout.open("convergence_Jacobi.txt");
	
	do
	{
		// Lower Block
		for(int j=1; j<j2; j++)		
		for(int i=i2; i<=i3; i++)
		{
			res_l_new[i][j]=b[i][j]-((as[i][j]*T[i][j-1])+(aw[i][j]*T[i-1][j])+(ap[i][j]*T[i][j])+(ae[i][j]*T[i+1][j])+(an[i][j]*T[i][j+1]));
			Tnew[i][j]=T[i][j]+((res_l_new[i][j])/ap[i][j]);
		}
			
		// Upper Block
		for(int j=j2; j<n; j++)		
		for(int i=1; i<m; i++)
		{
			res_u_new[i][j]=b[i][j]-((as[i][j]*T[i][j-1])+(aw[i][j]*T[i-1][j])+(ap[i][j]*T[i][j])+(ae[i][j]*T[i+1][j])+(an[i][j]*T[i][j+1]));	
			Tnew[i][j]=T[i][j]+((res_u_new[i][j])/ap[i][j]);
		}
	
		maxerr=error();		
			
		// Lower Block
		for(int j=1; j<j2; j++)		
		for(int i=i2; i<=i3; i++)
			T[i][j]=Tnew[i][j];
			
		// Upper Block
		for(int j=j2; j<n; j++)		
		for(int i=1; i<m; i++)
			T[i][j]=Tnew[i][j];
			
		cout<<"\nitr: "<<itr<<" error: "<<maxerr;
		errout<<fixed<<setprecision(16)<<log(itr)<<" "<<log(maxerr)<<endl;
		itr++;	
		
	}while(maxerr>convergence);
	
	errout.close();
	cout<<"\nJacobi operation done !!";
}

//......................................................LU.......................................................
void LU(float alpha)
{
	//Lower Block
	for(int j=1; j<j2; j++)
	for(int i=i2; i<=i3; i++)
	{
		Ls[i][j] = as[i][j]/(1+alpha*Ue[i][j-1]);
		Lw[i][j] = aw[i][j]/(1+alpha*Un[i-1][j]);
		Lp[i][j] = ap[i][j]-(Lw[i][j]*Ue[i-1][j])-(Ls[i][j]*Un[i][j-1])+(alpha*Lw[i][j]*Un[i-1][j])+(alpha*Ls[i][j]*Ue[i][j-1]);
		Ue[i][j] = (ae[i][j]-(alpha*Ls[i][j]*Ue[i][j-1]))/Lp[i][j];
		Un[i][j] = (an[i][j]-(alpha*Lw[i][j]*Un[i-1][j]))/Lp[i][j];
	}

	//Upper Block
	for(int j=j2; j<n; j++)		
	for(int i=1; i<m; i++)
	{
		Ls[i][j] = as[i][j]/(1+alpha*Ue[i][j-1]);
		Lw[i][j] = aw[i][j]/(1+alpha*Un[i-1][j]);
		Lp[i][j] = ap[i][j]-(Lw[i][j]*Ue[i-1][j])-(Ls[i][j]*Un[i][j-1])+(alpha*Lw[i][j]*Un[i-1][j])+(alpha*Ls[i][j]*Ue[i][j-1]);
		Ue[i][j] = (ae[i][j]-(alpha*Ls[i][j]*Ue[i][j-1]))/Lp[i][j];
		Un[i][j] = (an[i][j]-(alpha*Lw[i][j]*Un[i-1][j]))/Lp[i][j];
	}
}

//......................................LS (Forward Substitution)................................................
void LS()
{
	//Lower Block
	for(int j=1; j<j2; j++)
	for(int i=i2; i<=i3; i++)
		S[i][j] = (res_l[i][j]-(Ls[i][j]*S[i][j-1])-(Lw[i][j]*S[i-1][j]))/Lp[i][j];

	//Upper Block
	for(int j=j2; j<n; j++)		
	for(int i=1; i<m; i++)
		S[i][j] = (res_u[i][j]-(Ls[i][j]*S[i][j-1])-(Lw[i][j]*S[i-1][j]))/Lp[i][j];
}

//......................................UD (Backward Substitution)...............................................
void UD()
{
	//Upper Block
	for(int j=n-1; j>=j2; j--)
	for(int i=m-1; i>=1; i--)
		delta[i][j] = S[i][j]-(Ue[i][j]*delta[i+1][j])-(Un[i][j]*delta[i][j+1]);
			
	//Lower Block
	for(int j=j2-1; j>=1; j--)
	for(int i=i3; i>=i2; i--)
		delta[i][j] = S[i][j]-(Ue[i][j]*delta[i+1][j])-(Un[i][j]*delta[i][j+1]);
}

//.....................................................SIP........................................................
void SIP()
{
	ofstream errout;
	errout.open("convergence_SIP.txt");
	
	default_coefficients();
	
	boundary_coefficients();
	
	cout<<"\nalpha: ";
	cin>>alpha;
	
	do
	{
		maxerr=0.0;
		
		LU(alpha);
		
		//Lower Block
		for(int j=1; j<j2; j++)
		for(int i=i2; i<=i3; i++)
				res_l[i][j]=b[i][j]-((as[i][j]*T[i][j-1])+(aw[i][j]*T[i-1][j])+(ap[i][j]*T[i][j])+(ae[i][j]*T[i+1][j])+(an[i][j]*T[i][j+1]));
		
		//Upper Block
		for(int j=j2; j<n; j++)		
		for(int i=1; i<m; i++)
				res_u[i][j]=b[i][j]-((as[i][j]*T[i][j-1])+(aw[i][j]*T[i-1][j])+(ap[i][j]*T[i][j])+(ae[i][j]*T[i+1][j])+(an[i][j]*T[i][j+1]));
		
		LS();

		UD();	

		//Lower Block
		for(int j=1; j<j2; j++)
		for(int i=i2; i<=i3; i++)
			Tnew[i][j] = T[i][j]+delta[i][j];

		//Upper Block
		for(int j=j2; j<n; j++)		
		for(int i=1; i<m; i++)
			Tnew[i][j] = T[i][j]+delta[i][j];

		//Lower Block
		for(int j=1; j<j2; j++)
		for(int i=i2; i<=i3; i++)
				res_l_new[i][j]=b[i][j]-((as[i][j]*Tnew[i][j-1])+(aw[i][j]*Tnew[i-1][j])+(ap[i][j]*Tnew[i][j])+(ae[i][j]*Tnew[i+1][j])+(an[i][j]*Tnew[i][j+1]));
		
		//Upper Block
		for(int j=j2; j<n; j++)		
		for(int i=1; i<m; i++)
				res_u_new[i][j]=b[i][j]-((as[i][j]*Tnew[i][j-1])+(aw[i][j]*Tnew[i-1][j])+(ap[i][j]*Tnew[i][j])+(ae[i][j]*Tnew[i+1][j])+(an[i][j]*Tnew[i][j+1]));
		
		//Lower Block
		for(int j=1; j<j2; j++)
		for(int i=i2; i<=i3; i++)
			T[i][j] = Tnew[i][j];

		//Upper Block
		for(int j=j2; j<n; j++)		
		for(int i=1; i<m; i++)
			T[i][j] = Tnew[i][j];

		maxerr=error();
		cout<<"\nitr: "<<itr<<" error: "<<maxerr;
		errout<<fixed<<setprecision(16)<<log(itr)<<" "<<log(maxerr)<<endl;
		itr++;

	}while(maxerr>convergence);
	
	errout.close();
	cout<<"\nSIP operation done !!";
}

//....................................................BiCGSTAB...........................................................................
void BiCGSTAB()
{
	ofstream errout;
	errout.open("convergence_BiCGSTAB.txt");

	default_coefficients();
	boundary_coefficients();

	bet=0.0, rho=0.0;

	//Initial Residual
	//Lower Block
	for(int j=1; j<j2; j++)
	for(int i=i2; i<=i3; i++)
		r_tilda[i][j]=b[i][j]-((as[i][j]*T[i][j-1])+(aw[i][j]*T[i-1][j])+(ap[i][j]*T[i][j])+(ae[i][j]*T[i+1][j])+(an[i][j]*T[i][j+1]));
		
	//Upper Block
	for(int j=j2; j<n; j++)		
	for(int i=1; i<m; i++)
		r_tilda[i][j]=b[i][j]-((as[i][j]*T[i][j-1])+(aw[i][j]*T[i-1][j])+(ap[i][j]*T[i][j])+(ae[i][j]*T[i+1][j])+(an[i][j]*T[i][j+1]));
	

	//r,P,r_star
	//Lower Block
	for(int j=1; j<j2; j++)
	for(int i=i2; i<=i3; i++)
	{
		r[i][j] = r_tilda[i][j];
		P[i][j] = r[i][j];
		r_star[i][j] = r[i][j];
	}

	//Upper Block
	for(int j=j2; j<n; j++)		
	for(int i=1; i<m; i++)
	{
		r[i][j] = r_tilda[i][j];
		P[i][j] = r[i][j];
		r_star[i][j] = r[i][j];
	}

	do
	{
		//P_tilda
		//Lower Block
		for(int j=1; j<j2; j++)
		for(int i=i2; i<=i3; i++)
			P_tilda[i][j] = (as[i][j]*P[i][j-1])+(aw[i][j]*P[i-1][j])+(ap[i][j]*P[i][j])+(ae[i][j]*P[i+1][j])+(an[i][j]*P[i][j+1]);
		
		//Upper Block
		for(int j=j2; j<n; j++)		
		for(int i=1; i<m; i++)
			P_tilda[i][j] = (as[i][j]*P[i][j-1])+(aw[i][j]*P[i-1][j])+(ap[i][j]*P[i][j])+(ae[i][j]*P[i+1][j])+(an[i][j]*P[i][j+1]);

		//Finding alpha_bicg
		sum_dr=0.0;
		delta1=0.0;
		//Lower Block
		for(int j=1; j<j2; j++)
		for(int i=i2; i<=i3; i++)
		{
			delta1 = delta1+(r_star[i][j]*r[i][j]);
			sum_dr = sum_dr+(r_star[i][j]*P_tilda[i][j]);
		}
		//Upper Block
		for(int j=j2; j<n; j++)		
		for(int i=1; i<m; i++)
		{
			delta1 = delta1+(r_star[i][j]*r[i][j]);
			sum_dr = sum_dr+(r_star[i][j]*P_tilda[i][j]);
		}
		alpha_bicg=delta1/sum_dr;


		//Finding t
		//Lower Block
		for(int j=1; j<j2; j++)
		for(int i=i2; i<=i3; i++)
			t[i][j] = r[i][j]-(alpha_bicg*P_tilda[i][j]);
		//Upper Block
		for(int j=j2; j<n; j++)		
		for(int i=1; i<m; i++)
			t[i][j] = r[i][j]-(alpha_bicg*P_tilda[i][j]);
		
		//Finding t_tilda
		//Lower Block
		for(int j=1; j<j2; j++)
		for(int i=i2; i<=i3; i++)
			t_tilda[i][j] = (as[i][j]*t[i][j-1])+(aw[i][j]*t[i-1][j])+(ap[i][j]*t[i][j])+(ae[i][j]*t[i+1][j])+(an[i][j]*t[i][j+1]);
		//Upper Block
		for(int j=j2; j<n; j++)		
		for(int i=1; i<m; i++)	
			t_tilda[i][j] = (as[i][j]*t[i][j-1])+(aw[i][j]*t[i-1][j])+(ap[i][j]*t[i][j])+(ae[i][j]*t[i+1][j])+(an[i][j]*t[i][j+1]);

		//Finding rho
		t_nr=0.0, t_dr=0.0;
		//Lower Block
		for(int j=1; j<j2; j++)
		for(int i=i2; i<=i3; i++)
		{
			t_nr = t_nr+(t_tilda[i][j]*t[i][j]);
			t_dr = t_dr+(t_tilda[i][j]*t_tilda[i][j]);
		}
		//Upper Block
		for(int j=j2; j<n; j++)		
		for(int i=1; i<m; i++)
		{
			t_nr = t_nr+(t_tilda[i][j]*t[i][j]);
			t_dr = t_dr+(t_tilda[i][j]*t_tilda[i][j]);
		}
		rho = t_nr/t_dr;

		//Updating Temperature (T)
		//Lower Block
		for(int j=1; j<j2; j++)
		for(int i=i2; i<=i3; i++)
		{
			T[i][j] = T[i][j]+(alpha_bicg*P[i][j])+(rho*t[i][j]);
			r[i][j] = t[i][j]-(rho*t_tilda[i][j]);
		}
		//Upper Block
		for(int j=j2; j<n; j++)		
		for(int i=1; i<m; i++)
		{
			T[i][j] = T[i][j]+(alpha_bicg*P[i][j])+(rho*t[i][j]);
			r[i][j] = t[i][j]-(rho*t_tilda[i][j]);
		}
		
		//Finding bet
		sum_nr=0.0;
		//Lower Block
		for(int j=1; j<j2; j++)
		for(int i=i2; i<=i3; i++)
			sum_nr = sum_nr+(r_star[i][j]*r[i][j]);	
		//Upper Block
		for(int j=j2; j<n; j++)		
		for(int i=1; i<m; i++)
			sum_nr = sum_nr+(r_star[i][j]*r[i][j]);	
		bet = (alpha_bicg*sum_nr)/(rho*delta1);


		//Finding P
		//Lower Block
		for(int j=1; j<j2; j++)
		for(int i=i2; i<=i3; i++)
			P[i][j] = r[i][j]+(bet*(P[i][j]-(rho*P_tilda[i][j])));
		//Upper Block
		for(int j=j2; j<n; j++)		
		for(int i=1; i<m; i++)
			P[i][j] = r[i][j]+(bet*(P[i][j]-(rho*P_tilda[i][j])));
		
		//Error Calculation
		//Lower Block
		for(int j=1; j<j2; j++)
		for(int i=i2; i<=i3; i++)
		res_l_new[i][j]=b[i][j]-((as[i][j]*T[i][j-1])+(aw[i][j]*T[i-1][j])+(ap[i][j]*T[i][j])+(ae[i][j]*T[i+1][j])+(an[i][j]*T[i][j+1]));
		//Upper Block
		for(int j=j2; j<n; j++)		
		for(int i=1; i<m; i++)
		res_u_new[i][j]=b[i][j]-((as[i][j]*T[i][j-1])+(aw[i][j]*T[i-1][j])+(ap[i][j]*T[i][j])+(ae[i][j]*T[i+1][j])+(an[i][j]*T[i][j+1]));
		
		maxerr=error();

		
		cout<<"\nitr: "<<itr<<" error: "<<maxerr;
		errout<<fixed<<setprecision(16)<<log(itr)<<" "<<log(maxerr)<<endl;
		itr++;	

	}while(maxerr>convergence);
	
	errout.close();
	cout<<"\nBiCG operation done !!";
}

//....................................................BiCGSTAB_PRE...........................................................................
void BiCGSTAB_PRE()
{
	ofstream errout;
	errout.open("convergence_BiCGSTAB_Pre.txt");

	default_coefficients();
	boundary_coefficients();

	cout<<"\nEnter alpha:";
	cin>>alpha;
	LU(alpha);

	bet=0.0, rho=0.0;

	//Initial Residual
	//Lower Block
	for(int j=1; j<j2; j++)
	for(int i=i2; i<=i3; i++)
		r_tilda[i][j]=b[i][j]-((as[i][j]*T[i][j-1])+(aw[i][j]*T[i-1][j])+(ap[i][j]*T[i][j])+(ae[i][j]*T[i+1][j])+(an[i][j]*T[i][j+1]));
		
	//Upper Block
	for(int j=j2; j<n; j++)		
	for(int i=1; i<m; i++)
		r_tilda[i][j]=b[i][j]-((as[i][j]*T[i][j-1])+(aw[i][j]*T[i-1][j])+(ap[i][j]*T[i][j])+(ae[i][j]*T[i+1][j])+(an[i][j]*T[i][j+1]));
	
	//LS=r_tilda:Forward Substitution
	//Lower Block
	for(int j=1; j<j2; j++)
	for(int i=i2; i<=i3; i++)
		S[i][j] = (r_tilda[i][j]-(Ls[i][j]*S[i][j-1])-(Lw[i][j]*S[i-1][j]))/Lp[i][j];

	//Upper Block
	for(int j=j2; j<n; j++)		
	for(int i=1; i<m; i++)
		S[i][j] = (r_tilda[i][j]-(Ls[i][j]*S[i][j-1])-(Lw[i][j]*S[i-1][j]))/Lp[i][j];
	
	//Ur_tilda=S:Backward Substitution
	//Upper Block
	for(int j=n-1; j>=j2; j--)
	for(int i=m-1; i>=1; i--)
		r_tilda[i][j] = S[i][j]-(Ue[i][j]*r_tilda[i+1][j])-(Un[i][j]*r_tilda[i][j+1]);
			
	//Lower Block
	for(int j=j2-1; j>=1; j--)
	for(int i=i3; i>=i2; i--)
		r_tilda[i][j] = S[i][j]-(Ue[i][j]*r_tilda[i+1][j])-(Un[i][j]*r_tilda[i][j+1]);


	//r,P,r_star
	//Lower Block
	for(int j=1; j<j2; j++)
	for(int i=i2; i<=i3; i++)
	{
		r[i][j] = r_tilda[i][j];
		P[i][j] = r[i][j];
		r_star[i][j] = r[i][j];
	}

	//Upper Block
	for(int j=j2; j<n; j++)		
	for(int i=1; i<m; i++)
	{
		r[i][j] = r_tilda[i][j];
		P[i][j] = r[i][j];
		r_star[i][j] = r[i][j];
	}

	do
	{
		//P_tilda
		//Lower Block
		for(int j=1; j<j2; j++)
		for(int i=i2; i<=i3; i++)
			P_tilda[i][j] = (as[i][j]*P[i][j-1])+(aw[i][j]*P[i-1][j])+(ap[i][j]*P[i][j])+(ae[i][j]*P[i+1][j])+(an[i][j]*P[i][j+1]);
		
		//Upper Block
		for(int j=j2; j<n; j++)		
		for(int i=1; i<m; i++)
			P_tilda[i][j] = (as[i][j]*P[i][j-1])+(aw[i][j]*P[i-1][j])+(ap[i][j]*P[i][j])+(ae[i][j]*P[i+1][j])+(an[i][j]*P[i][j+1]);

		//LS=P_tilda:Forward Substitution
		//Lower Block
		for(int j=1; j<j2; j++)
		for(int i=i2; i<=i3; i++)
			S[i][j] = (P_tilda[i][j]-(Ls[i][j]*S[i][j-1])-(Lw[i][j]*S[i-1][j]))/Lp[i][j];

		//Upper Block
		for(int j=j2; j<n; j++)		
		for(int i=1; i<m; i++)
			S[i][j] = (P_tilda[i][j]-(Ls[i][j]*S[i][j-1])-(Lw[i][j]*S[i-1][j]))/Lp[i][j];
	
		//UP_tilda=S:Backward Substitution
		//Upper Block
		for(int j=n-1; j>=j2; j--)
		for(int i=m-1; i>=1; i--)
			P_tilda[i][j] = S[i][j]-(Ue[i][j]*P_tilda[i+1][j])-(Un[i][j]*P_tilda[i][j+1]);
			
		//Lower Block
		for(int j=j2-1; j>=1; j--)
		for(int i=i3; i>=i2; i--)
			P_tilda[i][j] = S[i][j]-(Ue[i][j]*P_tilda[i+1][j])-(Un[i][j]*P_tilda[i][j+1]);		

		//Finding alpha_bicg
		sum_dr=0.0;
		delta1=0.0;
		//Lower Block
		for(int j=1; j<j2; j++)
		for(int i=i2; i<=i3; i++)
		{
			delta1 = delta1+(r_star[i][j]*r[i][j]);
			sum_dr = sum_dr+(r_star[i][j]*P_tilda[i][j]);
		}
		//Upper Block
		for(int j=j2; j<n; j++)		
		for(int i=1; i<m; i++)
		{
			delta1 = delta1+(r_star[i][j]*r[i][j]);
			sum_dr = sum_dr+(r_star[i][j]*P_tilda[i][j]);
		}
		alpha_bicg=delta1/sum_dr;


		//Finding t
		//Lower Block
		for(int j=1; j<j2; j++)
		for(int i=i2; i<=i3; i++)
			t[i][j] = r[i][j]-(alpha_bicg*P_tilda[i][j]);
		//Upper Block
		for(int j=j2; j<n; j++)		
		for(int i=1; i<m; i++)
			t[i][j] = r[i][j]-(alpha_bicg*P_tilda[i][j]);
		
		//Finding t_tilda
		//Lower Block
		for(int j=1; j<j2; j++)
		for(int i=i2; i<=i3; i++)
			t_tilda[i][j] = (as[i][j]*t[i][j-1])+(aw[i][j]*t[i-1][j])+(ap[i][j]*t[i][j])+(ae[i][j]*t[i+1][j])+(an[i][j]*t[i][j+1]);
		//Upper Block
		for(int j=j2; j<n; j++)		
		for(int i=1; i<m; i++)	
			t_tilda[i][j] = (as[i][j]*t[i][j-1])+(aw[i][j]*t[i-1][j])+(ap[i][j]*t[i][j])+(ae[i][j]*t[i+1][j])+(an[i][j]*t[i][j+1]);

		//LS=t_tilda:Forward Substitution
		//Lower Block
		for(int j=1; j<j2; j++)
		for(int i=i2; i<=i3; i++)
			S[i][j] = (t_tilda[i][j]-(Ls[i][j]*S[i][j-1])-(Lw[i][j]*S[i-1][j]))/Lp[i][j];

		//Upper Block
		for(int j=j2; j<n; j++)		
		for(int i=1; i<m; i++)
			S[i][j] = (t_tilda[i][j]-(Ls[i][j]*S[i][j-1])-(Lw[i][j]*S[i-1][j]))/Lp[i][j];
	
		//Ut_tilda=S:Backward Substitution
		//Upper Block
		for(int j=n-1; j>=j2; j--)
		for(int i=m-1; i>=1; i--)
			t_tilda[i][j] = S[i][j]-(Ue[i][j]*t_tilda[i+1][j])-(Un[i][j]*t_tilda[i][j+1]);
			
		//Lower Block
		for(int j=j2-1; j>=1; j--)
		for(int i=i3; i>=i2; i--)
			t_tilda[i][j] = S[i][j]-(Ue[i][j]*t_tilda[i+1][j])-(Un[i][j]*t_tilda[i][j+1]);	

		//Finding rho
		t_nr=0.0, t_dr=0.0;
		//Lower Block
		for(int j=1; j<j2; j++)
		for(int i=i2; i<=i3; i++)
		{
			t_nr = t_nr+(t_tilda[i][j]*t[i][j]);
			t_dr = t_dr+(t_tilda[i][j]*t_tilda[i][j]);
		}
		//Upper Block
		for(int j=j2; j<n; j++)		
		for(int i=1; i<m; i++)
		{
			t_nr = t_nr+(t_tilda[i][j]*t[i][j]);
			t_dr = t_dr+(t_tilda[i][j]*t_tilda[i][j]);
		}
		rho = t_nr/t_dr;

		//Updating Temperature (T)
		//Lower Block
		for(int j=1; j<j2; j++)
		for(int i=i2; i<=i3; i++)
		{
			T[i][j] = T[i][j]+(alpha_bicg*P[i][j])+(rho*t[i][j]);
			r[i][j] = t[i][j]-(rho*t_tilda[i][j]);
		}
		//Upper Block
		for(int j=j2; j<n; j++)		
		for(int i=1; i<m; i++)
		{
			T[i][j] = T[i][j]+(alpha_bicg*P[i][j])+(rho*t[i][j]);
			r[i][j] = t[i][j]-(rho*t_tilda[i][j]);
		}
		
		//Finding bet
		sum_nr=0.0;
		//Lower Block
		for(int j=1; j<j2; j++)
		for(int i=i2; i<=i3; i++)
			sum_nr = sum_nr+(r_star[i][j]*r[i][j]);	
		//Upper Block
		for(int j=j2; j<n; j++)		
		for(int i=1; i<m; i++)
			sum_nr = sum_nr+(r_star[i][j]*r[i][j]);	
		bet = (alpha_bicg*sum_nr)/(rho*delta1);


		//Finding P
		//Lower Block
		for(int j=1; j<j2; j++)
		for(int i=i2; i<=i3; i++)
			P[i][j] = r[i][j]+(bet*(P[i][j]-(rho*P_tilda[i][j])));
		//Upper Block
		for(int j=j2; j<n; j++)		
		for(int i=1; i<m; i++)
			P[i][j] = r[i][j]+(bet*(P[i][j]-(rho*P_tilda[i][j])));
		
		//Error Calculation
		//Lower Block
		for(int j=1; j<j2; j++)
		for(int i=i2; i<=i3; i++)
		res_l_new[i][j]=b[i][j]-((as[i][j]*T[i][j-1])+(aw[i][j]*T[i-1][j])+(ap[i][j]*T[i][j])+(ae[i][j]*T[i+1][j])+(an[i][j]*T[i][j+1]));
		//Upper Block
		for(int j=j2; j<n; j++)		
		for(int i=1; i<m; i++)
		res_u_new[i][j]=b[i][j]-((as[i][j]*T[i][j-1])+(aw[i][j]*T[i-1][j])+(ap[i][j]*T[i][j])+(ae[i][j]*T[i+1][j])+(an[i][j]*T[i][j+1]));
		
		maxerr=error();

		
		cout<<"\nitr: "<<itr<<" error: "<<maxerr;
		errout<<fixed<<setprecision(16)<<log(itr)<<" "<<log(maxerr)<<endl;
		itr++;	

	}while(maxerr>convergence);
	
	errout.close();
	cout<<"\nBiCG-Pre operation done !!";
}

//......................................................Approach-1 (Parabolic Equation).............................................................................
void approach_1()
{
	cout<<"\nEnter w: ";
	cin>>w;
	
	ofstream errout;
	errout.open("convergence_BiCGSTAB_Pre.txt");
	
	// Time Loop Starts
	do			
	{
		int k=1;
			
		//Lower Block
		for(int j=1; j<j2; j++)
		for(int i=i2; i<=i3; i++)
			Told[i][j]=T[i][j];
		
		//Upper Block
		for(int j=j2; j<n; j++)		
		for(int i=1; i<m; i++)
			Told[i][j]=T[i][j];
			
		default_coefficients_app1();
		boundary_coefficients();
		
		do
		{
			sum=0.0, maxerr=0.0;
		
			//Lower Block
			for(int j=1; j<j2; j++)
			for(int i=i2; i<=i3; i++)
				{
					res_l_new[i][j]=b[i][j]-((as[i][j]*T[i][j-1])+(aw[i][j]*T[i-1][j])+(ap[i][j]*T[i][j])+(ae[i][j]*T[i+1][j])+(an[i][j]*T[i][j+1]));
					var=T[i][j]+((w*res_l_new[i][j])/ap[i][j]);
					sum = sum+(pow((var-T[i][j]),2));
					T[i][j]=var;
				}
				sum=sum/(ml*nl);
				err_l=sqrt(sum);
				
			sum=0.0;	
			//Upper Block
			for(int j=j2; j<n; j++)		
			for(int i=1; i<m; i++)
			{
				res_u_new[i][j]=b[i][j]-((as[i][j]*T[i][j-1])+(aw[i][j]*T[i-1][j])+(ap[i][j]*T[i][j])+(ae[i][j]*T[i+1][j])+(an[i][j]*T[i][j+1]));
				var=T[i][j]+((w*res_u_new[i][j])/ap[i][j]);
				sum = sum+(pow((var-T[i][j]),2));
				T[i][j]=var;
			}
			sum=sum/(mu*nu);
			err_u=sqrt(sum);
			
			if(err_l>=err_u)
				maxerr=err_l;
			else
				maxerr=err_u;
			
		//	maxerr=error();
			cout<<"\nitr: "<<k<<" error: "<<maxerr;
			errout<<fixed<<setprecision(16)<<log(itr)<<" "<<log(maxerr)<<endl;
			itr++;
			k++;
			
		}while(maxerr>convergence);
		
		sum_l=0.0, sum_u=0.0, maxerrt=0.0;
		
		sum=0.0;
		//Lower Block
		for(int j=1; j<j2; j++)
		for(int i=i2; i<=i3; i++)
		//	sum_l = sum_l+(pow((T[i][j]-Told[i][j]),2));
		//	sum_l = sum_l/(ml*nl);
		sum = sum+ pow((T[i][j]-Told[i][j]), 2);
			err_l = sqrt(sum/(ml*nl));
			
			sum=0.0;
		//Upper Block
		for(int j=j2; j<n; j++)		
		for(int i=1; i<m; i++)
		//	sum_u = sum_u+(pow((T[i][j]-Told[i][j]),2));
		//	sum_u = sum_u/(mu*nu);
		//	err_u = sqrt(sum_u);
			
			sum = sum+ pow((T[i][j]-Told[i][j]), 2);
			err_u = sqrt(sum/(mu*nu));
		
		if(err_l>=err_u)
		maxerrt=err_l;
		else
		maxerrt=err_u;
		
		cout<<"\n\t\tTime: "<<samay<<" error: "<<maxerrt;
		samay=samay+dt;
			
	}while(maxerrt>convergence);
}

int main()
{
	clock_t s;
	
	grid();
	cout<<"\nGrid Generated !!\nGrid printed !!";

	cout<<"\nm="<<m<<"\nn="<<n<<"\ni2="<<i2<<"\ni3="<<i3<<"\nj2="<<j2;
	
	initialize();
	cout<<"\nSolution Initialization done !!";

	cell_count();
	
	s=clock();
	
	//gauss_siedel();	
	//jacobi();
	//SIP();
	//BiCGSTAB();
	//BiCGSTAB_PRE();
	approach_1();
	
	s=clock()-s;
	
	double time_taken=((double)s)/CLOCKS_PER_SEC;
	ofstream ftime;
	ftime.open("Time taken.txt");
	ftime<<"Time taken to converge: "<<time_taken<<" sec";
	
	cout<<"\n SOLUTION COVERGED !!"<<"\nComputational Time: "<<time_taken<<" sec";
	
	print_solution();
		
	return 0;
}

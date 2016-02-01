//><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
//Programer: UBC ADMB May Workshop													 
//Date:	May 12-16, 2013														 
//Purpose:Delay-Difference.											 
//Notes: 				 
//							 
//																 
//><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>

DATA_SECTION
	init_int syr;
	init_int eyr;
	init_number surv;
	init_number rho;
	init_number alpha;
	init_int agek;
	init_int nage;
	init_vector wa(1,nage)
	init_vector yt(syr,eyr);
	init_vector ct(syr,eyr);
	init_vector wbart(syr,eyr);
	init_int eof;
	number wk;
	!!wk=wa(agek);
	int iter;
	!!iter=0;
	LOCAL_CALCS
		if(eof!=999)
		{
			cout<<"Error reading data.\n Fix it."<<endl;
			ad_exit(1);
		}
	END_CALCS


PARAMETER_SECTION
	init_number log_ro;
	!!log_ro=log(4000);
	init_number log_cr;
	!!log_cr=log(10);
	init_number log_ptau;
	!!log_ptau=log(1./(0.3*0.3));
	objective_function_value nll;
	
	number tau;
	number ro;
	number bo;
	number no;
	number cr;
	number q;
	number a;
	number b;
	number bpr;
	number fpen;
	!!fpen=0;
	vector yt_resid(syr,eyr);
	vector Bt(syr,eyr);
	vector Nt(syr,eyr);
	vector Rt(syr-agek+1,eyr);
	vector Ct(syr,eyr);
	vector Wbart(syr,eyr);

PROCEDURE_SECTION
	initialization();
	statedynamics();
	observation_model();
	objective_function();

	if(mceval_phase())
	{ 
		forecast();
		mcmc_output();
	}
	if(last_phase())
	{
		forecast();
	} 

FUNCTION initialization
	ro=mfexp(log_ro);
	cr=mfexp(log_cr)+1.;
	bpr=(surv*alpha/(1.-surv)+wk)/(1.-surv*rho);
	bo=ro*bpr;
	no=ro/(1.-surv);
	a=cr/bpr;
	b=(cr-1.)/(ro*bpr);
	Rt(syr-agek+1,syr)=ro;
	Nt(syr)=no;
	Bt(syr)=bo;

FUNCTION statedynamics
	for(int t=syr+1;t<=eyr;t++)
	{
		dvariable hr;
		hr=.99-posfun(.99-ct(t-1)/Bt(t-1),0.0001,fpen);
		dvariable surv2=surv*(1.-hr);
		Bt(t)=surv2*(alpha*Nt(t-1)+rho*Bt(t-1))+wk*Rt(t-agek);
		Nt(t)=surv2*Nt(t-1)+Rt(t-agek);
		Rt(t)=a*(1.-hr)*Bt(t)/(1.+b*(1.-hr)*Bt(t));
		Ct(t)=hr*Bt(t);
		Wbart(t)=Bt(t)/Nt(t);
	}

FUNCTION observation_model
	yt_resid=log(yt)-log(Bt);
	q=mfexp(mean(yt_resid));
	yt_resid-=mean(yt_resid);


FUNCTION objective_function 
	tau=sqrt(1./mfexp(log_ptau));
	nll=dnorm(yt_resid,tau);
	nll+=1.e5*fpen;
	//if(fpen>0)cout<<fpen<<endl;

FUNCTION mcmc_output
	if(iter==0)
	{
		//ofstream ofs("refpar.mcmc");
		//ofs<<"fmsy\t bmsy\t msy\t b/bmsy\t f/fmsy"<<endl;
		//ofs<<"r\t k\t q\t sig\t"<<endl;
	}
	iter++;
	//double fratio=value(ft(eyr)/fmsy);
	//double bratio=value(Nt(eyr)*wa/bmsy);
	//ofstream ofs("refpar.mcmc",ios::app);
	//ofs<<fmsy<<"\t"<<bmsy<<"\t"<<msy<<"\t"<<bratio<<"\t"<<fratio<<endl;

FUNCTION forecast


REPORT_SECTION
	REPORT(ro);	
	REPORT(bpr);
	REPORT(bo);	
	REPORT(no);	
	REPORT(cr);	
	REPORT(yt);
	REPORT(Bt*q);
	REPORT(tau);
	REPORT(ct);
	REPORT(elem_div(ct,value(Bt)));
	REPORT(yt_resid);
TOP_OF_MAIN_SECTION
	
	time(&start);
	arrmblsize = 50000000;
	gradient_structure::set_GRADSTACK_BUFFER_SIZE(1.e7);
	gradient_structure::set_CMPDIF_BUFFER_SIZE(1.e7);
	gradient_structure::set_MAX_NVAR_OFFSET(5000);
	gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000);


GLOBALS_SECTION
	/**
	\def REPORT(object)
	Prints name and value of \a object on ADMB report %ofstream file.
	*/
	#undef REPORT
	#define REPORT(object) report << #object "\n" << object << endl;

	#include <admodel.h>
	#include <time.h>
	#include <contrib.h>//IF you have ADMB-11
	//#include<stats.cxx>//If you have ADMB-10 and make sure stats.cxx is in your working directory
	time_t start,finish;
	long hour,minute,second;
	double elapsed_time;

FINAL_SECTION
	time(&finish);
	elapsed_time=difftime(finish,start);
	hour=long(elapsed_time)/3600;
	minute=long(elapsed_time)%3600/60;
	second=(long(elapsed_time)%3600)%60;
	cout<<"*******************************************"<<endl;
	cout<<"--Start time: "<<ctime(&start)<<endl;
	cout<<"--Finish time: "<<ctime(&finish)<<endl;
	cout<<"--Runtime: ";
	cout<<hour<<" hours, "<<minute<<" minutes, "<<second<<" seconds"<<endl;
	cout<<"*******************************************"<<endl;



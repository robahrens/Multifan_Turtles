//><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>
//Programer: Robert Ahrens MULTIFAN type Code													 
//Date:	Feb 3, 2014														 
//Purpose:Barazilian Sea Turtle LF Data.											 
//Notes: Alan and Karen Bolten
//Schnute & Fournier 1980 CJFAS 37: 1337-1351.
//Fournier et al. 1990 CJFAS 47 301-317		
//Note that formulas have been modified so that the first age class does not need to be one 							 
//																 
//><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>><>

DATA_SECTION
	init_int nsets;
	init_int nbins;
	init_number minbin;
	init_number binw;
	init_ivector years(1,nsets);
	init_ivector months(1,nsets);
	init_matrix lfdata(1,nsets,1,nbins);
	init_int eof;
	int iter;
	!!iter=0;
	LOCAL_CALCS
		if(eof!=999)
		{
			cout<<"Error reading data.\n Fix it."<<endl;
			ad_exit(1);
		}
	END_CALCS

	!!ad_comm::change_datafile_name("turtleLF.ctl");
	init_number fage;
	init_int nage;
	init_number ivbk;
	init_number ilfirst;
	init_number idl;
	init_number ilam1;
	init_number ilam2;
	init_number itau;
	
	init_number pmeanlone; //prior age 1
	init_number psdmeanlone;//prior age 1
	init_number pmeanlzero;//prior hatchling
	init_number psdmeanlzero;//prior hatchling
	init_number pmeanvbk;//prior vbk
	init_number psdmeanvbk;//prior vbk
	init_number pmeanlinf;//prior linf
	init_number psdmeanlinf;//prior linf
	init_int eof2;
	LOCAL_CALCS
		if(eof2!=999)
		{
			cout<<"Error reading control.\n Fix it."<<endl;
			ad_exit(1);
		}
	END_CALCS
	vector age(1,nage);
	vector lbinmids(1,nbins);
	matrix opl(1,nsets,1,nbins);
	LOCAL_CALCS
		age.fill_seqadd(fage,1.);
		lbinmids.fill_seqadd(minbin,binw);
		for(int i=1;i<=nsets;i++)opl(i)=lfdata(i)/(rowsum(lfdata)(i));
	END_CALCS

PARAMETER_SECTION
	init_bounded_number vbk(0,1.);
	init_number log_lfirst;
	init_number log_dl;
	init_number log_lam1;
	init_number log_lam2;
	init_number log_tau;
	init_bounded_matrix page(1,nsets,1,nage,0,1);
	objective_function_value nll;
	!!vbk=ivbk;
	!!log_lfirst=log(ilfirst);
	!!log_dl=log(idl);
	!!log_lam1=log(ilam1);
	!!log_lam2=log(ilam2);
	!!log_tau=log(itau);
	number npar;
	number LL;
	number AIC;
	number AICc;
	number BIC
	number lfirst;
	number lone;
	number lzero;
	number lmax;
	//number cv;
	number lam1;
	number lam2;
	number tau;
	number tzero;
	number rho;
	number linf;
	number fpen;
	sdreport_number ratio;

	vector la(1,nage);
	vector sdla(1,nage);
	matrix npage(1,nsets,1,nage);
	matrix pla(1,nbins,1,nage);
	matrix ppl(1,nsets,1,nbins);
	matrix pnl(1,nsets,1,nbins)
	matrix meanl(1,nsets,1,nage)
	matrix sdl(1,nsets,1,nage)

	matrix epsilon(1,nsets,1,nbins);

PRELIMINARY_CALCS_SECTION

PROCEDURE_SECTION
	initialization();
	lfset_calculations();
	objective_function();

	if(mceval_phase())
	{
		mcmc_output(); 
	}
	if(last_phase())
	{
	} 

FUNCTION initialization
	fpen=0;
	lmax=mfexp(log_lfirst)+mfexp(log_dl);
	lfirst=mfexp(log_lfirst);
	lam1=mfexp(log_lam1);
	lam2=mfexp(log_lam2);
	//cv=mfexp(log_cv);
	npage=page/sum(page);
	//cout<<npage<<endl;
	tau=mfexp(log_tau);
	rho=mfexp(-vbk);
	linf=(lmax-lfirst*pow(rho,(nage-fage)))/(1-pow(rho,(nage-fage)));
	tzero=-(posfun(-(fage-1./log(rho)*log((lmax-lfirst)/(lmax-lfirst*pow(rho,(nage-fage))))),0.001,fpen));

FUNCTION lfset_calculations
	dvariable z1;
	dvariable z2;	
	for(int i=1;i<=nsets;i++)
	{
		npage(i)=page(i)/sum(page(i));
		la=lfirst+(lmax-lfirst)*(1.-pow(rho,age-fage+(months(i)-1.)/12.))/(1.-pow(rho,nage-fage));
		meanl(i)=la;
		sdla=lam1*mfexp(lam2*(-1.+2.*(1.-pow(rho,age-fage+(months(i)-1.)/12.))/(1.-pow(rho,nage-fage))));
		sdl(i)=sdla;
		for(int j=1;j<=nage;j++) //loop over ages
		{
			 for(int k=1;k<=nbins;k++) //loop over length bins
			{
				z1=((lbinmids(k)-0.5*binw)-la(j))/sdla(j);
				z2=((lbinmids(k)+0.5*binw)-la(j))/sdla(j);
				pla(k,j)=cumd_norm(z2)-cumd_norm(z1);
			}//end nbins
		}//end nage
		dvar_vector pl=pla*npage(i);
		ppl(i)=pl/sum(pl);
		pnl(i)=ppl(i)*(rowsum(lfdata)(i));
		epsilon(i)=elem_prod((1-ppl(i)),ppl(i))+0.1/nbins;
	}
	lzero=linf*(1.-mfexp(-vbk*(0.-tzero)));
	lone=linf*(1.-mfexp(-vbk*(1.-tzero)));

FUNCTION objective_function 
	dvariable comp1;
	dvariable comp2;
	dvariable comp3;
	dvariable comp4;
	dvariable comp5;
	dvariable comp6;
	dvariable comp7;

	//tau=sqrt(1./(nsets*nbins)*sum(elem_div(elem_prod(opl-ppl,opl-ppl),elem_prod(ppl,1-ppl))));
	comp1=-0.5*sum(log(epsilon));
	comp2=-nbins*log(tau);
	comp3=sum(log(mfexp(-elem_div(elem_prod((opl-ppl),(opl-ppl)),(2*epsilon*tau*tau))+0.01)));
	comp4=dlnorm(lone,log(pmeanlone),psdmeanlone);
	comp5=dlnorm(lzero,log(pmeanlzero),psdmeanlzero);
	comp6=dlnorm(vbk,log(pmeanvbk),psdmeanvbk);
	comp7=dlnorm(linf,log(pmeanlinf),psdmeanlinf);
	nll=-(comp1+comp2+comp3)+comp4+comp5+comp6+comp7+fpen;
	npar=nsets*nage+6;
	LL=comp1+comp2+comp3;
	AIC=2*npar-2*(comp1+comp2+comp3);
	AICc=AIC+(2*npar*(npar-1.))/(nsets*nbins-npar-1.);
	BIC=npar*log(nsets*nbins)-2*(comp1+comp2+comp3);

FUNCTION mcmc_output

	
	if(iter==0)
	{
		ofstream ofs("par.mcmc");
		ofs<<"vbk\t linf\t to\t"<<endl;
	}
	iter++;
	ofstream ofs("par.mcmc",ios::app);
	ofs<<vbk<<"\t"<<linf<<"\t"<<tzero<<endl;

REPORT_SECTION

	REPORT(LL);
	REPORT(npar);
	REPORT(AIC);
	REPORT(AICc);
	REPORT(BIC);
	REPORT(fpen);
	REPORT(fage);
	REPORT(nage);
	REPORT(rho);
	REPORT(lfirst);
	REPORT(lmax);
	REPORT(lam1);
	REPORT(lam2);
	REPORT(tau);
	REPORT(vbk);
	REPORT(linf);
	REPORT(tzero);
	REPORT(lzero);
	REPORT(lone);
	REPORT(pmeanlone);
	REPORT(psdmeanlone);
	REPORT(pmeanlzero);//prior hatchling
	REPORT(psdmeanlzero);//prior hatchling
	REPORT(pmeanvbk);//prior vbk
	REPORT(psdmeanvbk);//prior vbk
	REPORT(pmeanlinf);//prior linf
	REPORT(psdmeanlinf);//prior linf
	REPORT(npage);
	REPORT(lbinmids);
	REPORT(opl);
	REPORT(ppl);
	REPORT(lfdata);
	REPORT(pnl);
	REPORT(meanl);
	REPORT(sdl);
	ofstream ofs("short.rep");
	ofs<<vbk<<"\t"<<linf<<"\t"<<tzero<<"\t"<<lzero<<"\t"<<lone<<"\t"<<LL<<"\t"<<npar<<"\t"<<AIC<<"\t"<<AICc<<"\t"<<BIC<<endl;
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



/**All SAS code is based on the code provided in 
Hernán MA, Robins JM (2020). Causal Inference: What If. Boca Raton: Chapman & Hall/CRC.
*/

/*Set up your library*/
libname lab "C:\Users\ejmurray\Dropbox\Teaching\CausalBook\nhefs_sas";

/*Read in the data*/
data nhefs;
set lab.nhefs;
run;

/*create an indicator for missing outcomes*/
data nhefs;
set nhefs;
cens = (wt82 eq .);
run;

/**************************************************************
/*Standardization over age group*/
/**************************************************************/;
/*Create age group and product term*/
data nhefs;
set nhefs;
older = .;
if age > 50 then older = 1;
else if .< age le 50 then older = 0;
qsmkolder = qsmk*older;
run;

/*Fit a model for weight given smoking cessation and age: E[Y|A, L]= ?0 + ?1 A + ?2 L + ?3 AL*/
proc genmod data = nhefs;
	model wt82_71 = qsmk older qsmkolder ;
run;
quit;
/*Calculate strata prevalence*/
proc freq data = nhefs;
where cens = 0;
tables older /norow nocol;
run;

/**************************************************************
/*Standardization over sex and age group*/
/**************************************************************/;
/*Create product terms*/
data nhefs; 
set nhefs; 
	qsmksex = sex*qsmk; 
	sexolder = sex*older; 
	qsmksexolder = qsmk*sex*older; 
run;
/*Fit a model for weight given smoking cessation, age group, sex, and all possible product terms*/
proc genmod data = nhefs;
	model wt82_71 = qsmk older sex qsmkolder qsmksex sexolder qsmksexolder;
run;
quit;

/*Calculate strata prevalence*/
proc freq data = nhefs;
where cens = 0;
tables sex*older /norow nocol;
run;

/**************************************************************
Standardizing the mean outcome to the baseline confounders
Data from NHEFS
***************************************************************/;

/* create a dataset with 3 copies of each subject */  
data onesample ;
  set nhefs ;
  label interv= "Intervention"; 
  interv = -1 ;    /* 1st copy: equal to original one */
  	output ; 
  interv = 0 ;     /* 2nd copy: treatment set to 0, outcome to missing */
  	qsmk = 0 ;
  	wt82_71 = . ;
  	output ;  
  interv = 1 ;     /* 3rd copy: treatment set to 1, outcome to missing*/
  	qsmk = 1 ;
  	wt82_71 = . ;
  	output ;    
run;


* linear model to estimate mean outcome conditional on treatment and confounders;
* parameters are estimated using original observations only (interv= -1) ;
* parameter estimates are used to predict mean outcome for observations with 
  treatment set to 0 (interv=0) and to 1 (innterv=1);
proc genmod data = onesample;
	class exercise active education;
	model wt82_71 = qsmk sex race age age*age education
				smokeintensity smokeintensity*smokeintensity smokeyrs smokeyrs*smokeyrs
				exercise active wt71 wt71*wt71
				qsmk*smokeintensity
				;
    output out = predicted_mean p = meanY ;
run;

* estimate mean outcome in each of the groups interv=0, and interv=1;
* this mean outcome is a weighted average of the mean outcomes in each combination 
	of values of treatment and confounders, that is, the standardized outcome;
proc means data = predicted_mean mean noprint;
  class interv ;
  var meanY ;
  types interv ;
  output out = results (keep = interv mean ) mean = mean ;
run;

proc print data = results noobs label ;
  title "Parametric g-formula";
  var interv mean ;
run;

/*****************************
G-formula for survival 
******************************/

/* some preprocessing of the data */
data nhefs;
    set lab.nhefs;
    if death=0 then survtime=120; 
	else if death=1 then survtime= (yrdth-83)*12 + modth; * yrdth ranges from 83 to 92;
run;

/* creation of person-month data */
data nhefs_surv;
	length seqn 8. time 8. event 8.;
   	set nhefs;
	do time= 0 to (survtime-1);
		event= (death=1 and time=survtime-1);
		timesq= time*time;
		output;
	end;
run;


/* fit a pooled logistic regression model with covariates to estimate the 
observed conditional hazard of mortality over time; here we assume a quadratic
function of time. Include an interaction between time & exposure to relax the 
proportional hazards assumption*/
proc logistic data= nhefs_surv outmodel = gf_model;
	ods exclude ClassLevelInfo Type3 Association FitStatistics GlobalTests;
	class exercise active education;
	model event = qsmk qsmk*time qsmk*timesq time timesq
				sex race age age*age education
				smokeintensity smokeintensity*smokeintensity smkintensity82_71
				smokeyrs smokeyrs*smokeyrs
				exercise active wt71 wt71*wt71;
	output out=est_prob p=p_noevent;
run;

/* make a dataset that has every individual with every time point under every treatment;
i.e. "if no one were censored and if everyone received each treatment"*/
data for_gf ;
	set nhefs_surv ( where = (time = 0));
	do time = 0 to 119;
	    timesq = time*time ;
		qsmk = 0 ;
		output ;
		qsmk = 1 ;
		output ;
	end;
run;

/*use the model results to get predicted hazards for each individual at each time under each treatment*/
proc logistic inmodel=gf_model ;
	score  data=for_gf out=gf_pred  (keep = seqn time qsmk P_0 rename = (P_0=p_noevent)) ;
run;
proc sort data = gf_pred ; by seqn qsmk time ; run;

/*Calculate predicted survival for each individual in our simulated dataset, 
for each exposure level, using a Kaplan-Meier product estimate method*/
data gf_pred ;
	set gf_pred ;
	by seqn qsmk ;
	retain surv ;
	if first.qsmk then surv = 1 ;
	surv = surv * p_noevent ;
run;

/*Summarize by time point and exposure level*/
proc means data = gf_pred noprint  ;
	class qsmk time ;
	var surv ;
	types qsmk*time ;
	output out= mysurv_gf (keep  = qsmk time surv) mean=surv ;
run;

/*Output a graph-friendly version of survival over time*/
data qsmk0_gf qsmk1_gf ;
	set mysurv_gf ;
	if qsmk = 0 then output qsmk0_gf ;
	if qsmk = 1 then output qsmk1_gf ;
	keep time surv;
run;
data gf_graph ;
	merge qsmk0_gf (rename = (surv = surv0)) qsmk1_gf (rename = (surv = surv1));
	survdiff= surv1 - surv0;
run;
proc print data= gf_graph; id time; run;

/*Plot estimated standardized survival over time*/
proc sgplot data = gf_graph ;
   series x = time y = surv0 /legendlabel='A = 0';
   series x = time y = surv1 /legendlabel='A = 1';
   title 'Survival from g-formula';
   yaxis label = 'Survival '  values = (0.6 to 1 by 0.2) ;
   xaxis label = 'Months '  values = (0 to 120 by 12) ;
   keylegend/noborder title='A: ';
run;
title;



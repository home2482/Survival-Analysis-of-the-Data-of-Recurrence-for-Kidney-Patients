data kidrecurr;
input patient $ time1 infect1 time2	infect2	age	gender $	gn $	an $	pkd $;
age10=age-10;
cards;
5	30	1	12	1	10	0	0	0	0
27	132	1	156	1	10	1	1	0	0
6	24	1	245	1	16.5	1	0	0	0
15	536	1	25	0	17	1	0	0	0
35	119	1	8	1	22	1	0	0	0
1	16	1	8	1	28	0	0	0	0
22	24	0	402	1	30	1	0	0	0
4	318	1	447	1	31.5	1	0	0	0
3	22	1	28	1	32	0	0	0	0
12	141	1	8	0	34	1	0	0	0
13	38	1	96	1	35	1	0	1	0
14	70	0	149	0	42	1	0	1	0
36	54	0	16	0	42	1	0	0	0
24	39	1	46	0	42.5	1	0	1	0
25	40	1	12	1	43	0	0	1	0
18	114	1	292	1	43.5	1	0	0	0
11	7	1	33	1	44	1	0	1	0
20	108	0	15	1	44	1	0	0	0
34	190	1	5	0	44.5	1	1	0	0
21	562	1	152	1	46.5	0	0	0	1
2	13	0	23	1	48	1	1	0	0
32	5	0	43	1	50.5	1	0	1	0
7	9	1	7	1	51	0	1	0	0
10	154	1	15	1	51.5	0	1	0	0
28	34	1	30	1	52	1	0	1	0
37	6	0	78	1	52	1	0	0	1
19	159	0	22	0	53	1	1	0	0
29	2	1	25	1	53	0	1	0	0
30	26	1	130	1	54	1	1	0	0
8	30	1	511	1	55.5	1	1	0	0
31	27	1	58	1	56	1	0	1	0
33	152	1	30	1	57	1	0	0	1
26	113	0	201	1	57.5	1	0	1	0
16	4	0	17	1	60	0	0	1	0
17	185	1	177	1	60	1	0	0	0
38	8	0	63	1	60	0	0	0	1
23	66	1	13	1	62.5	1	0	1	0
9	53	1	196	1	69	1	0	1	0
;
run;
proc lifetest data=kidrecurr outsurv=out;   
time time1*infect1(0); 
strata gender;/*same as gn an pkd*/
run;
proc lifetest data=kidrecurr outsurv=out;   
time time1*infect1(0); 
strata  pkd;
run; 
proc lifetest data=kidrecurr outsurv=out;   
time time1*infect1(0); 
strata gn ;
run; 
proc lifetest data=kidrecurr outsurv=out;   
time time1*infect1(0); 
strata an ;
run;  
proc phreg data=kidrecurr;  
model time1*infect1(0)=age10 ; 
run;
proc phreg data = kidrecurr;
class gender gn an pkd; 
model time1*infect1(0)= age10 gender gn an pkd  /ties=efron;
run;
proc phreg data = kidrecurr;
class gender gn an pkd; 
model time1*infect1(0)= age10 gender gn an pkd  /ties=breslow;
run;
proc phreg data = kidrecurr;
class gender gn an pkd; 
model time1*infect1(0)= age10 gender gn an pkd  /ties=exact;
run;
proc phreg data = kidrecurr;
class gender gn an pkd; 
model time1*infect1(0)= age10|gender|gn|an |pkd  /ties=efron;
run;
proc phreg data=kidrecurr;
class gender gn an pkd;
model time1*infect1(0)= gender gn an pkd /ties=efron;
output out=diagnostic xbeta=risk_score resdev=dev resmart=mart;
run; 
/*Residual analysis on the original model (martingale and deviance residuals*/
proc loess data=diagnostic;
model mart=age10;
ods output OutputStatistics=test_age10;
run;
proc sort data=test_age10;
by age10;
run;
proc gplot data=test_age10;
title 'Martingale Residuals with Loess Fit';
symbol1 color=black i=none value=dot;
symbol2 color=red i=join value=circle;
plot depvar*age10 pred*age10/overlay;
run;
/*creating piecewise age10 function*/
data kidney;set kidrecurr;
if age10 le 22 then age32=age10;
else age32=0;
run;
/*full model with all covariates + piecewise age10 function*/
proc phreg data=kidney;
class gender gn an pkd;
model time1*infect1(0)= age10 age32 gender gn an pkd/ties=efron;
output out=diagnostic1 resdev=dev resmart=mart;
run;
proc gplot data=diagnostic1;  
title 'Deviance Residual Plot';  
symbol value=dot i=none;  
plot dev*age10;
run; 
/* checking proportional hazard assumption */
/* wtressch gives scaled Sch. residual while ressch gives un-scaled one */ 
/* If use “class” statement, be very careful to name the wtressch */ 
/* The order of wtressch residuals is the same as the order of parameters */ 
/* If a categorical covariate A has 3 d.f., it has 3 wtressch residuals */  
proc phreg data=kidney;
class gender gn an pkd; 
model time1*infect1(0)=age10 age32 gender gn an pkd/ ties=efron;  
assess ph/resample; 
output out=ph_diag 
wtressch=wtschage10 wtschage32 a b c d;
run;
proc gplot data=ph_diag; 
title'Schoenfeld residuals of age10 versus survial time';
plot wtschage10*time1; 
symbol value=dot i=none; 
run;  
proc gplot data=ph_diag;
plot a*time1 b*time1 c*time1 d*time1;
symbol1 value=dot i=none;
run;
proc reg data=ph_diag;
model wtschage10=time1;
model a=time1;
model b=time1;
model c=time1;
model d=time1;
run;
proc phreg data = kidney; 
class gender gn an pkd ; 
model time1*infect1(0)= age10 gender gn an pkd /ties=efron selection=backward;
run;
proc phreg data=kidney;
class gender pkd;
model time1*infect1(0)=gender pkd/ties=efron;
run;
proc phreg data=kidney;
class gender pkd;
model time1*infect1(0)= gender pkd/ties=efron;
assess PH/resample;
run;

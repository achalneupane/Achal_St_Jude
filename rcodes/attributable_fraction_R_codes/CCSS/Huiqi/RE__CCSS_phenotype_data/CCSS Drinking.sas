%let ccssdata=%str(Z:\SJShare\SJCOMMON\ECC\Ness Research Team\CCSS\CCSS 20200205);

libname cfu5 "&ccssdata.\combined\followup5" access=readonly;
libname cfu6 "&ccssdata.\combined\followup6" access=readonly;
libname cothers "&ccssdata.\combined\others" access=readonly;
libname cpsycho "&ccssdata.\combined\psycho" access=readonly;

libname ebase "&ccssdata.\expansion\baseline" access=readonly;
libname emraf "&ccssdata.\expansion\mraf" access=readonly;
libname eothers "&ccssdata.\expansion\others" access=readonly;

libname obase "&ccssdata.\original\baseline" access=readonly;
libname ofu1 "&ccssdata.\original\followup1" access=readonly;
libname ofu2 "&ccssdata.\original\followup2" access=readonly;
libname ofu3 "&ccssdata.\original\followup3" access=readonly;
libname ofu2007 "&ccssdata.\original\fu2007" access=readonly;
libname omraf "&ccssdata.\original\mraf" access=readonly;
libname oothers "&ccssdata.\original\others" access=readonly;

libname datafmt "&ccssdata." access=readonly;
libname icd "&ccssdata.\icdformat" access=readonly;
options fmtsearch=(datafmt icd);

%let mydata=%str(C:\Users\hwang3\OneDrive - St. Jude Children%'s Research Hospital\Project);
libname qi "&mydata.\0. Files\CCSS" access=readonly;
libname whq "&mydata.\13. Requet from Achal and Yadav";

%macro add_dummy_table();
	ods excel options(sheet_interval="table");
	ods exclude all;
	data _null_;
	file print;
	put _all_;
	run;
	ods select all;
%mend add_dummy_table;

/* 
For drinking, baseline items were very limited. When need to compare the baseline vs FU, needed to ignore some information from FU to 
make the drinking variable comparable. 
Now, with this project, no need to compare with the baseline. If additional items from FU said it is heavy drinking, then make it as heavy.
*/

/** Note: 
In my project, <18 were not used for drikning, so deleted in these codes. You may need to use those people depending on your project.
**/

/*Will use
Risky drinking definition: At least >4 drinks a day or 14 per week (male), at least >3 drinks per day or 7 drinks per week (women).
Heavy drinking is: as =6 drinks per day for men and =5 drinks per day for women, at least once per month.*/

/*#################################################################################################################*/
/*#################################################################################################################*/
/*#################################################################################################################*/
proc freq data=obase.baseN; table EVDRINK*DRINKLYR /norow nocol nopercent missing; run;
/*DRINKDAY: N7-A On days when drink, # of drinks
   1= no drinks in the past 2 years
   2= 1 drink/day
   3= 2 drinks/day
   4= 3 drinks/day
   5= 4 drinks/day
   6= 5 drinks/day
   7= 6+ drinks/day*/
data drkbase; 
merge obase.basea(keep=ccssid d_compq in=aa)
	  obase.baseN(rename=(EVDRINK=everdrink) keep=ccssid EVDRINK--DRINKLYR)
	  oothers.casestat(keep=ccssid sex d_birth);
by ccssid;
if aa;
/*everdrink: N3 Had at least one alcoholic drink  -- entire life. 
DRINKLYR: N8 At least one drink in past year*/
if DRINKLYR=1 then everdrink=1;  	/*Qi added: some yes to drikin in past year but missing in driking in entire life, should be yes*/
if everdrink=2 then DRINKLYR=2; /*No to drinking in entire life should be No to last year drinking.*/
if everdrink=2 or DRINKLYR=2 or DRINKDAY=1 then do; /*No to past year, or to the entire life means to No to heavry/risky.DRINKDAY=1 for no drinks.*/
	heavy=2; risky=2;
end;
else /*if everdrink in (.,1,3) then*/ do;
	weekdr=(sum(of n_wine n_beer n_mixdr))*12/52; /**from per months to per week**/
	if sex=1 then do;
		/*heavy drinking:  =6 drinks per day for men and =5 drinks per day for women, at least once per month*/
		if drinkday ge 7 then heavy=1;
		else if /*drinkday eq .*/ missing(drinkday) then heavy=.; /*Qi modified to handle - as missing*/
		else heavy=2;
		/**risky drinking: A day at least 4 drinks or weekdr least 14 for male**/
		if (weekdr > 14 or drinkday > 5) then risky=1;  /*drinkday > 5 means >4 drinks*/
		else if weekdr=. or drinkday=. then risky=.; 
		else risky=2;
	end;

	if sex=2 then do;
		/*heavy drinking:  =6 drinks per day for men and =5 drinks per day for women, at least once per month*/
		if drinkday ge 6 then heavy=1;
		else if drinkday eq . then heavy=.;
		else heavy=2;
		/**A day at least 3 drinks or weekdr least 7 for female.**/
		if (weekdr > 7 or drinkday > 4) then risky=1;  /*drinkday > 4 means >3 drinks*/
		else if weekdr=. or drinkday=. then risky=.; 
		else risky=2;
	end;
end;

/*Heavy drinking is a subset of risky drinking, in terms of the defintion. So for heavy=1 and risky=. should make risky=1*/
if heavy=1 then risky=1;
format risky yesnof.;
/*Can add No to those age<18, or delete <18*/
ageb=(d_compq-d_birth)/365.25;
if ageb<18 then delete;
run;
proc means data=drkbase; var ageb;; run;
proc freq data=drkbase; table sex*(heavy risky) heavy*risky DRINKDAY*(sex heavy risky)/nocol nopercent missing; run;
data miss; set drkbase(where=(heavy=. or risky=.)); run;


/*##########################################################################################*/
/*######code for FU2007, FU5, and expansion cohort baseline are the same ###################*/
/*##########################################################################################*/
/*MultDrinkFreq:N6 Number of times had 5+ (m) or 4+ (f) alcoholic drinks in last 12 months.
AnyDrinkFreq: N5 Number of times had any alcoholic drink in last 12 months
*/
proc freq data=ofu2007.fu2007N; table AnyDrinkFreq; run;
/*DFRQF
1=Every day                                        
2=5 to 6 times a week                               
3=3 to 4 times a week                              
4=twice a week                                      
5=once a week                                        
6=2 to 3 times a month                               
7=once a month                                       
8=3 to 11 times in the past year                     
9=1 or 2 times in the past year                     
10=Never in the past year*/

data drkf07/*(keep=ccssid alcohol heavy everdrink)*/;
merge ofu2007.fu2007a(in=aa keep=ccssid d_fu2007) 
	  ofu2007.fu2007N(keep=ccssid--MultDrinkFreq) 
	  oothers.casestat(keep=ccssid sex d_birth);
by ccssid;
if aa;
daydrink=sum(n_wine, n_beer, n_mixdr); /*This is per day, while baseline was per month. This "daydrink" is numeric drink amount, no format.*/

if anydrinkfreq eq 1 then drinkweek=daydrink*7;
else if anydrinkfreq eq 2 then drinkweek=daydrink*5.5;
else if anydrinkfreq eq 3 then drinkweek=daydrink*3.5;
else if anydrinkfreq eq 4 then drinkweek=daydrink*2;
else if anydrinkfreq eq 5 then drinkweek=daydrink*1;
else if anydrinkfreq eq 6 then drinkweek=daydrink*2.5/(52/12);
else if anydrinkfreq eq 7 then drinkweek=daydrink*1/(52/12);
else if anydrinkfreq eq 8 then drinkweek=daydrink*7/52;
else if anydrinkfreq eq 9 then drinkweek=daydrink*1.5/52;
else if anydrinkfreq eq 10 then drinkweek=daydrink*0;

/*Kyla and Zhenhong, St Jude. slightly different coefficients. don't know how it came up. Not far from the above when median was taken, but generally lowere here.
if AnyDrinkFreq=1 then wkdrink=(daydrink*1)*7;
if AnyDrinkFreq=2 then wkdrink=(daydrink*0.71)*7;
if AnyDrinkFreq=3 then wkdrink=(daydrink*0.43)*7;
if AnyDrinkFreq=4 then wkdrink=(daydrink*0.29)*7;
if AnyDrinkFreq=5 then wkdrink=(daydrink*0.14)*7;
if AnyDrinkFreq=6 then wkdrink=(daydrink*0.07)*7;
if AnyDrinkFreq=7 then wkdrink=(daydrink*0.03)*7;
if AnyDrinkFreq=8 then wkdrink=(daydrink*0.008)*7;
if AnyDrinkFreq=9 then wkdrink=(daydrink*0.003)*7;
if AnyDrinkFreq=10 then wkdrink=(daydrink*0)*7;*/

/*EverDrink: N1 Ever had at least 2 alcoholic drinks. This can be used to define No. 
Qi added the last 3 conditions: These are Never in the past year, which can define No to heavy/risky was for the last year, too, not for the entire life. */
if EverDrink=2 or MaxDrinkDay=0 or MultDrinkFreq=10 or AnyDrinkFreq=10 then do;
	heavy=2; risky=2;
end; 
else do;
	*heavy drinking;
	if sex=1 then do;
		if daydrink ge 6 then heavy=1;
		else if daydrink eq . then heavy=.;
		else heavy=2;
	end;
	if sex=2 then do;
		if daydrink ge 5 then heavy=1;
		else if daydrink eq . then heavy=.;
		else heavy=2;
	end;
	/*Qi questioned: In the definiton of the paper, it says "at least once per month" for heavy drinking. 
	"risky drinking was defined as >4 drinks per day or 14 drinks per week for men and >3 drinks per day or 7 drinks per week for women. 
	Heavy drinking was defined as =6 drinks per day for men and =5 drinks per day for women, at least once per month. "
	===== The code above used "daydrink" which is the drinks amount on days when drank, but did not use at least once per month. 
	I would add:*/
	if heavy=1 and AnyDrinkFreq in (8,9,10) then heavy=2; /*if < once per month, then heavy is No. -- no such a once per month for risky definiton.
	Qi added this based on the "per month" in the definition. */

	*risky drinking;
	if sex=1 then do;
		if daydrink gt 4 or drinkweek gt 14 then risky=1;
		else if daydrink eq . and drinkweek eq . then risky=.;
		else risky=2;
	end;
	if sex=2 then do;
		if daydrink gt 3 or drinkweek gt 7 then risky=1;
		else if daydrink eq . and drinkweek eq . then risky=.;
		else risky=2;
	end;
end;
agef7=(d_fu2007-d_birth)/365.25;
run;
proc freq data=drkf07; table heavy*risky /norow nocol nopercent missing; run;
data miss; set drkf07; if heavy=. or risky=.; run;
proc means data=drkf07; var agef7; run; /*all are >18 at FU2007.*/


/* These was not asked in the baseline. This can be used to define risky, too. The heavy below is heavy in Kayla's code. 
MaxDrinkDay: N4 Largest number of drinks on any day in last 12 months. 
multdrinkFreq: During the last 12 months, how often did you have 5 or more (males) or 4 or more (females) drinks containing any kind of alcohol in a single day?
 1 every day, …, 3 means 3-4 days a week, 4 means 2 days a week, 5 means 1 day a week. 
*/
data drkf07;
set drkf07;
if 0 LE multdrinkFreq LE 4 then alcohol=3; /***heavy drinking***/
else if (sex=1 and (MaxDrinkDay GT 4 or drinkweek GT 14)) or (sex=2 and (MaxDrinkDay GT 3 or drinkweek GT 7)) then alcohol=2; /*** risky drinking**/
else if MaxDrinkDay gt 0 then alcohol=1;
else if MaxDrinkDay=0 or everdrink=2 /*Qi added these 2 conditons:*/or MultDrinkFreq=10 or AnyDrinkFreq=10 then alcohol=0;
/*many people had missing on max day drikining and Yes to ever drink (entire life), but never in the past year, which should be No instead of missing,
since risky/heavy are for the past 12 months.*/
run;
proc freq data=drkf07; table alcohol*(heavy risky) heavy*risky/norow nocol missing; run;

/*These are the people who are heavy/risky based on N4 Largest number of drinks on any day in last 12 months.; or; frequency of 5 or more (males) or 4 or more (females) drinks
in a single day; but they are No to risky drinking, if using the beer/wine/mix drink.
===== don't know what to do!!!!!!!!!!!!!!!!!!! ==to be consistent with baseline, not use. ====if no need to be consistent with the baselin,
then we can ue them for the timepoint.*/
data chk;
set drkf07(where=(alcohol in (2,3) and risky^=1));
run;
proc freq data=drkf07; table (heavy risky)*alcohol /norow nocol nopercent missing; run;
data drkf07;
set drkf07; /*This made it not consistent with the baseline. == use items not asked in the baseline to further define heavy/risky.*/
if alcohol=3 then heavy=1; 
if alcohol in (2,3) then risky=1;
run;

/*########################################*/
/*############### FU5 #####################*/
data drkfu5/*(keep=ccssid alcohol heavy everdrink)*/;
merge cfu5.fu5a(in=aa keep=ccssid d_fu5) 
	  cfu5.fu5N(keep=ccssid--MultDrinkFreq)
	  oothers.casestat(keep=ccssid sex d_birth in=oo) 
	  eothers.Ecasestat(keep=ccssid sex d_birth in=ee) ;
by ccssid;
if aa;
if oo then group=1; else if ee then group=2;
daydrink=sum(n_wine, n_beer, n_mixdr); /*This is per day, while baseline was per month. This "daydrink" is numeric drink amount, no format.*/

if anydrinkfreq eq 1 then drinkweek=daydrink*7;
else if anydrinkfreq eq 2 then drinkweek=daydrink*5.5;
else if anydrinkfreq eq 3 then drinkweek=daydrink*3.5;
else if anydrinkfreq eq 4 then drinkweek=daydrink*2;
else if anydrinkfreq eq 5 then drinkweek=daydrink*1;
else if anydrinkfreq eq 6 then drinkweek=daydrink*2.5/(52/12);
else if anydrinkfreq eq 7 then drinkweek=daydrink*1/(52/12);
else if anydrinkfreq eq 8 then drinkweek=daydrink*7/52;
else if anydrinkfreq eq 9 then drinkweek=daydrink*1.5/52;
else if anydrinkfreq eq 10 then drinkweek=daydrink*0;

/*EverDrink: N1 Ever had at least 2 alcoholic drinks. This can be used to define No. 
Qi added the last 3 conditons: These are Never in the past year, which can define No to heavy/risky was for the last year, too, not for the entire life. */
if EverDrink=2 or MaxDrinkDay=0 or MultDrinkFreq=10 or AnyDrinkFreq=10 then do;
	heavy=2; risky=2;
end; 
else do;
	*heavy drinking;
	if sex=1 then do;
		if daydrink ge 6 then heavy=1;
		else if daydrink eq . then heavy=.;
		else heavy=2;
	end;
	if sex=2 then do;
		if daydrink ge 5 then heavy=1;
		else if daydrink eq . then heavy=.;
		else heavy=2;
	end;
	/*Qi questioned: In the definiton of the paper, it says "at least once per month" for heavy drinking. 
	" risky drinking was defined as >4 drinks per day or 14 drinks per week for men and >3 drinks per day or 7 drinks per week for women. 
	Heavy drinking was defined as =6 drinks per day for men and =5 drinks per day for women, at least once per month. "
	===== The code above used "daydrink" which is the drinks amount on days when drank, but did not use at least once per month. 
	I would add:*/
	if heavy=1 and AnyDrinkFreq in (8,9,10) then heavy=2; /*if < once per month, then heavy is No. -- no such a once per month for risky definiton.
	Qi added this based on the "per month" in the definition. */

	*risky drinking;
	if sex=1 then do;
		if daydrink gt 4 or drinkweek gt 14 then risky=1;
		else if daydrink eq . and drinkweek eq . then risky=.;
		else risky=2;
	end;
	if sex=2 then do;
		if daydrink gt 3 or drinkweek gt 7 then risky=1;
		else if daydrink eq . and drinkweek eq . then risky=.;
		else risky=2;
	end;
end;
agefu5=(d_fu5-d_birth)/365.25;
run;
proc freq data=drkfu5; table group*(heavy risky) /nocol nopercent; run;
proc means data=drkfu5; class group; var agefu5; run;
data drkfu5; set drkfu5; if agefu5<18 then delete; run; /*	11272/11337 remained.	*/

/*This made it not consistent with the baseline. == use items not asked in the baseline to further define heavy/risky.*/
data drkfu5;
set drkfu5;
if 0 LE multdrinkFreq LE 4 then alcohol=3; /***heavy drinking***/
else if (sex=1 and (MaxDrinkDay GT 4 or drinkweek GT 14)) or (sex=2 and (MaxDrinkDay GT 3 or drinkweek GT 7)) then alcohol=2; /*** risky drinking**/
else if MaxDrinkDay gt 0 then alcohol=1;
else if MaxDrinkDay=0 or everdrink=2 /*Qi added these 2 conditons:*/or MultDrinkFreq=10 or AnyDrinkFreq=10 then alcohol=0;
/*many people had missing on max day drikining and Yes to ever drink (entire life), but never in the past year, which should be No instead of missing,
since risky/heavy are for the past 12 months.*/
run;
proc freq data=drkfu5; table alcohol*(heavy risky) heavy*risky/norow nocol missing; run;
data drkfu5;
set drkfu5; 
if alcohol=3 then heavy=1; 
if alcohol in (2,3) then risky=1;
run;

/*#########################################*/
/*########### Expansion baseline ##########*/

data exbase/*(keep=ccssid alcohol heavy everdrink)*/;
merge ebase.Ebasea(keep=ccssid d_compq) 
	  ebase.ebaseO(in=aa keep=ccssid EverDrink--MultDrinkFreq in=aa)
	  eothers.Ecasestat(keep=ccssid sex d_birth) ;
by ccssid;
if aa;
group=2;
daydrink=sum(n_wine, n_beer, n_mixdr); /*This is per day, while baseline was per month. This "daydrink" is numeric drink amount, no format.*/

if anydrinkfreq eq 1 then drinkweek=daydrink*7;
else if anydrinkfreq eq 2 then drinkweek=daydrink*5.5;
else if anydrinkfreq eq 3 then drinkweek=daydrink*3.5;
else if anydrinkfreq eq 4 then drinkweek=daydrink*2;
else if anydrinkfreq eq 5 then drinkweek=daydrink*1;
else if anydrinkfreq eq 6 then drinkweek=daydrink*2.5/(52/12);
else if anydrinkfreq eq 7 then drinkweek=daydrink*1/(52/12);
else if anydrinkfreq eq 8 then drinkweek=daydrink*7/52;
else if anydrinkfreq eq 9 then drinkweek=daydrink*1.5/52;
else if anydrinkfreq eq 10 then drinkweek=daydrink*0;

/*EverDrink: N1 Ever had at least 2 alcoholic drinks. This can be used to define No. 
Qi added the last 3 conditons: These are Never in the past year, which can define No to heavy/risky was for the last year, too, not for the entire life. */
if EverDrink=2 or MaxDrinkDay=0 or MultDrinkFreq=10 or AnyDrinkFreq=10 then do;
	heavy=2; risky=2;
end; 
else do;
	*heavy drinking;
	if sex=1 then do;
		if daydrink ge 6 then heavy=1;
		else if daydrink eq . then heavy=.;
		else heavy=2;
	end;
	if sex=2 then do;
		if daydrink ge 5 then heavy=1;
		else if daydrink eq . then heavy=.;
		else heavy=2;
	end;
	/*Qi questioned: In the definiton of the paper, it says "at least once per month" for heavy drinking. 
	" risky drinking was defined as >4 drinks per day or 14 drinks per week for men and >3 drinks per day or 7 drinks per week for women. 
	Heavy drinking was defined as =6 drinks per day for men and =5 drinks per day for women, at least once per month. "
	===== The code above used "daydrink" which is the drinks amount on days when drank, but did not use at least once per month. 
	I would add:*/
	if heavy=1 and AnyDrinkFreq in (8,9,10) then heavy=2; /*if < once per month, then heavy is No. -- no such a once per month for risky definiton.
	Qi added this based on the "per month" in the definition. */

	*risky drinking;
	if sex=1 then do;
		if daydrink gt 4 or drinkweek gt 14 then risky=1;
		else if daydrink eq . and drinkweek eq . then risky=.;
		else risky=2;
	end;
	if sex=2 then do;
		if daydrink gt 3 or drinkweek gt 7 then risky=1;
		else if daydrink eq . and drinkweek eq . then risky=.;
		else risky=2;
	end;
end;
ageb=(d_compq-d_birth)/365.25;
run;
proc freq data=exbase; table group*(heavy risky) /nocol nopercent; run;
proc means data=exbase; class group; var ageb; run;
data exbase; set exbase; if ageb<18 then delete; run; /*9270/10756 remained*/
/*This made it not consistent with the baseline. == use items not asked in the baseline to further define heavy/risky.*/
data exbase;
set exbase;
if 0 LE multdrinkFreq LE 4 then alcohol=3; /***heavy drinking***/
else if (sex=1 and (MaxDrinkDay GT 4 or drinkweek GT 14)) or (sex=2 and (MaxDrinkDay GT 3 or drinkweek GT 7)) then alcohol=2; /*** risky drinking**/
else if MaxDrinkDay gt 0 then alcohol=1;
else if MaxDrinkDay=0 or everdrink=2 /*Qi added these 2 conditons:*/or MultDrinkFreq=10 or AnyDrinkFreq=10 then alcohol=0;
/*many people had missing on max day drikining and Yes to ever drink (entire life), but never in the past year, which should be No instead of missing,
since risky/heavy are for the past 12 months.*/
run;
proc freq data=exbase; table alcohol*(heavy risky) heavy*risky/norow nocol missing; run;
data exbase;
set exbase; 
if alcohol=3 then heavy=1; 
if alcohol in (2,3) then risky=1;
run;

data drink;
merge drkbase(keep=ccssid heavy risky ageb rename=(heavy=heavyb risky=riskyb) in=aa)
	  exbase(keep=ccssid heavy risky ageb rename=(heavy=heavyb risky=riskyb) in=bb)
	  drkf07(keep=ccssid heavy risky agef7 rename=(heavy=heavyf7 risky=riskyf7) in=cc)
	  drkfu5(keep=ccssid heavy risky agefu5 rename=(heavy=heavyfu5 risky=riskyfu5) in=dd)
	  /*also need who completed the survey   === same format on person competing survy across the 4 Q's*/
	  obase.basea(keep=ccssid RELTOPT/*format RELTOPF*/ rename=(RELTOPT=selfb))
	  ofu2007.fu2007a(keep=ccssid PersonCompletingRelationship/*format RELTOPF*/ rename=(PersonCompletingRelationship=selff7))
	  ebase.ebasea(keep=ccssid reltop/*format RELTOPF*/ rename=(reltop=selfb) )
	  cfu5.fu5a(keep=ccssid reltop/*format RELTOPF*/  rename=(reltop=selffu5))
	  oothers.casestat(keep=ccssid sex d_birth in=oo) 
	  eothers.Ecasestat(keep=ccssid sex d_birth in=ee);
by ccssid;
if aa or bb or cc or dd; /*18+ at each Q*/
if oo then group=1; else if ee then group=2;
if selfb=8 then selfb=1; else selfb=2;
if selff7=8 then selff7=1; else selff7=2;
if selffu5=8 then selffu5=1; else selffu5=2;
run;
proc means data=drink; var age:; run;
proc means data=drink; class group; var age:; run;

proc freq data=drink; tables heavyb*riskyb heavyf7*riskyf7 heavyfu5*riskyfu5; run;

data whq.drink;
set drink;
risky_mostrecent=coalesce(riskyfu5, riskyf7, riskyb);
label risky_mostrecent="Heavy/Risky drinking from most recent data";
format risky_mostrecent riskyfu5 riskyf7 heavyb heavyf7 heavyfu5 YESNOF.;
run;
proc freq data=whq.drink; tables risky_mostrecent; run;


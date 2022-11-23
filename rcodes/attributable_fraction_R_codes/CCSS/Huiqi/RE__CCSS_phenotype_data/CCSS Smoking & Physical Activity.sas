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

/*----------------------|
|         Smoke         |
-----------------------*/
/*#############Smoking status#############*/
/*In baseline: "N1 Smoked at least 100 cigarettes"*/
proc freq data=obase.baseN; table EVSMOKE*smokenow/norow nocol nopercent missing; run;
data smkbase; 
merge obase.basea(keep=ccssid d_compq)
	  obase.baseN(keep=ccssid EVSMOKE SMOKENOW A_SMOKE N_CIGDAY N_YSMOKE A_SMKEND);
by ccssid;
if EVSMOKE=3 then EVSMOKE=2; /*make not sure to No*/
if EVSMOKE=2 then SMOKENOW=2; /*Eversmoke=No then smokenow=No*/
if SMOKENOW=1 and EVSMOKE=. then do;  EVSMOKE=1; d_smokeb=d_compq; end; /*some answered Now smoking but missing in ever smoked. Use the date now as the date for ever smoking*/
run;
proc freq data=smkbase; table EVSMOKE*smokenow/norow nocol nopercent missing; run;

/*In fu2003, it is asked "L1 Ever smoked at least 100 cigarettes"*/
proc freq data=ofu2.f2main; table EVSMOKE*SMOKENOW/norow nocol nopercent missing; run;
data smkf2; 
set ofu2.f2main(keep=ccssid d_fu2 EVSMOKE SMOKENOW A_SMOKE N_CIGDAY N_YSMOKE QUIT12MO);
if EVSMOKE=2 then SMOKENOW=2; /*Eversmoke=No then smokenow=No*/
run; 
proc freq data=smkf2; table EVSMOKE*SMOKENOW/norow nocol nopercent missing; run;

/*In Fu2007, "N7 Smoked at least 100 cigarettes in previous two years"*/
proc freq data=ofu2007.fu2007N; table EvSmoke2yr*SmokeNow/norow nocol nopercent missing; run;
data smkf07; 
merge ofu2007.fu2007a(keep=ccssid d_fu2007)
	  ofu2007.fu2007N(keep=ccssid EvSmoke2yr SmokeNow A_Smoke N_CIGDAY N_YSMOKE N_QUITSMOKE); 
by ccssid;
if EvSmoke2yr=2 then smokenow=2; /*Eversmoke in 2 years=No then smokenow=No*/
if smokenow=1 and EvSmoke2yr=. then do; EvSmoke2yr=1; d_smokef07=d_fu2007; end; /*some answered Now smoking but missing in ever smoked.*/
run;
proc freq data=smkf07; table EvSmoke2yr*SmokeNow/norow nocol nopercent missing; run;

/*In Fu5, "N7 Smoked at least 100 cigarettes in previous two years"*/
proc freq data=cfu5.fu5N; table SMOKESINCE*SmokeNow/norow nocol nopercent missing; run;
data smkfu5; 
merge cfu5.fu5a(keep=ccssid d_fu5)
	  cfu5.fu5N(keep=ccssid SMOKESINCE SmokeNow A_Smoke N_CIGDAY N_YSMOKE N_QUITSMOKE); 
by ccssid;
if SMOKESINCE=2 then smokenow=2; /*Eversmoke in 2 years=No then smokenow=No*/
if smokenow=1 and SMOKESINCE=. then do; SMOKESINCE=1; d_smokefu5=d_fu5; end; /*some answered Now smoking but missing in ever smoked.*/
run;
proc freq data=smkfu5; table SMOKESINCE*SmokeNow/norow nocol nopercent missing; run;

data smk;
merge smkbase(rename=(EVSMOKE=evsmokeb SMOKENOW=smokenowb a_smoke=a_smokeb) in=aa)
	  smkf2(rename=(EVSMOKE=evsmokef2 SMOKENOW=smokenowf2 A_SMOKE=a_smokef2))
	  smkf07(rename=(EvSmoke2yr=EvSmoke2yrf07 SmokeNow=smokenowf07 a_smoke=a_smokef07)) /*EvSmoke2yr: N7 Smoked at least 100 cigarettes in previous TWO years*/
	  smkfu5(rename=(SMOKESINCE=smokefu5 SmokeNow=smokenowfu5 a_smoke=a_smokefu5)) /*SMOKESINCE: N7 Smoked at least 100 cigarettes in previous TWO years*/
	  oothers.casestat(keep=ccssid d_birth);
by ccssid;
if aa;
if month(d_birth)=2 and day(d_birth)=29 then d_birth=mdy(3,1,year(d_birth));
/*if Yes to ever smoking in any questionnaire, then Yes.
Among those who answered smoking in any questionnaire and excluding any Yes, make it No === not perfectly but this way we reduce missing in smoking variable*/
if evsmokeb=1 or evsmokef2=1 or EvSmoke2yrf07=1 or smokefu5=1 then evsmoke=1;
else if evsmokeb=2 or evsmokef2=2 or EvSmoke2yrf07=2 or smokefu5=2 then evsmoke=2; 
/*date started smoking
a_smoke=min(a_smokeb,a_smokef2,a_smokef07);*/
/*Take from the first questionnaire that gave smoking age. The COALSCE function returns the first nonmissing value in a list of variables.*/
if evsmoke=1 then a_smoke=coalesce(a_smokeb,a_smokef2,a_smokef07,a_smokefu5);
d_smoke=mdy(month(d_birth),day(d_birth),year(d_birth)+a_smoke)+365.25/2;
if d_smoke=. and (d_smokeb^=. or d_smokef07^=. or d_smokefu5^=.) then d_smoke=coalesce(d_smokeb,d_smokef07,d_smokefu5); /**for those who answered smoking now but no to ever smoking, no one like this in Fu2**/
format d_smoke date9.;
run;
proc freq data=smk; table evsmoke; run;
/*Among those with Yes to every smoking, anyone had d_smoke missing?*/
proc means data=smk(where=(evsmoke=1)) nmiss; var d_smoke; run;
	data see; set smk(where=(evsmoke=1 and d_smoke=.)); run;
/*For these people, give a randomly date between 18 years old the the date of questionanire where they said ever smoked.
But what if the questionaire is after death? So should between 18 years and min(questionnaire,death)*/
data smk;
merge smk(in=ss) qi.death(keep=ccssid d_death);
by ccssid;
if ss;
d18=mdy(month(d_birth),day(d_birth),year(d_birth)+18); format d18 date9.;
if evsmoke=1 and d_smoke=. then do;
	aa=ranuni(2016);
	if evsmokeb=1 then d_smoke=d18+(min(d_compq,d_death)-d18)*aa;
	else if evsmokef2=1 then do;
		if evsmokeb^=2 then d_smoke=d18+(min(d_fu2,d_death)-d18)*aa;
		else if evsmokeb=2 then d_smoke=max(d18,d_compq)+(min(d_fu2,d_death)-max(d18,d_compq))*aa; /*if yes in FU2 but no in baseline, impute a time between baseline and Fu2*/
	end;
	else if EvSmoke2yrf07=1 then do;
		if evsmokef2=2 then d_smoke=d18+(min(d_fu2007,d_death)-max(d18,d_fu2))*aa;
		else if evsmokeb=2 then d_smoke=d18+(min(d_fu2007,d_death)-max(d18,d_compq))*aa;
		else d_smoke=d18+(min(d_fu2007,d_death)-d18)*aa;
	end;
	else if smokefu5=1 then do;
		if EvSmoke2yrf07=2 then d_smoke=d18+(min(d_fu5,d_death)-max(d18,d_fu2007))*aa;
		else if evsmokef2=2 then d_smoke=d18+(min(d_fu5,d_death)-max(d18,d_fu2))*aa;
		else if evsmokeb=2 then d_smoke=d18+(min(d_fu5,d_death)-max(d18,d_compq))*aa;
		else d_smoke=d18+(min(d_fu5,d_death)-d18)*aa;
	end;
end;
run;
	data see; set smk(where=(evsmoke=1 and d_smoke=.)); run;
/*I did not impute a_smoke, so it was reported. but some were very young ages. For those who had smoking age <=14, change it to be 14 to be reasonable*/
data smk;
set smk;
if evsmoke=1 then do;
	if (d_smoke-d_birth)/365.25<14 then d_smoke=mdy(6,1,year(d_birth)+14);
end;
run;

/*The above we make ever smoking and first smoking time.
Now make cateogrical never, former, and now smoking. */
data smk;
set smk;
/*FU007 and Fu5 did not ask eve smoke, but smoking in the last 2 years. So those ever smoke before these Q's should be yes to ever smoke.*/
evsmokef07=EvSmoke2yrf07; if evsmokeb=1 or evsmokef2=1 then evsmokef07=1;
evsmokef5=smokefu5; if evsmokef07=1 then evsmokef5=1;
array varnow(*) smokenowb smokenowf2 smokenowf07 smokenowfu5;
array varever(*) evsmokeb evsmokef2 evsmokef07 evsmokef5;
array varcat(*) smkcatb smkcatf2 smkcatf07 smkcatf5;
do i=1 to dim(varnow);
	if varnow(i)=1 then varcat(i)=3; /*current*/
	else if varever(i)=1 and varnow(i) in (2) then varcat(i)=2; /*former: 
	I prefer to make ever smoke but smoking now missing as formaer instead of missing.
	but if don't use varnow(i)=2, those who did not have FU5 with ever smoking to FU2007 would be made to former --- since I made those wihtout FU5 but
	with yes to ever smoking in FU2007 as yes to ever smoking in FU5.*/
	else if varever(i)=2 then varcat(i)=1; /*never*/
end; drop i;
run;
/*For those ever smoking yes, but smoking now missing, should we leave it as missing, or making it as former? it should be either former or current, cannot be never.
using the carry forward to reduce missing, possible to make it as never?*/
proc freq data=smk; table evsmokef2*smokenowf2  smkcatf5/norow nocol nopercent missing; run;
proc freq data=smk(where=(evsmokef2=1 and smokenowf2=.)); table (evsmokeb smokenowb)*smkcatb /norow nocol nopercent missing; run;
proc freq data=smk;
table (EvSmoke2yrf07 smokenowf07)*evsmokef07 (smokefu5 smokenowfu5)*evsmokef5  smkcatb smkcatf2 smkcatf07 smkcatf5/norow nocol nopercent missing;
run;
proc freq data=smk;
table (evsmokef07 smokenowf07)*smkcatf07 (evsmokef5 smokenowfu5)*smkcatf5 (evsmokeb smokenowb)*smkcatb (evsmokef2 smokenowf2)*smkcatf2 /norow nocol nopercent missing;
run;

/*Expanson: Smoking status EVSMOKE2yr=lifetimesmoke*/
proc freq data=Ebase.ebaseo; table EvSmoke2yr*SmokeNow /missing; run;
data smkbase; 
merge Ebase.ebasea(keep=ccssid d_compq)
	  Ebase.ebaseo(keep=ccssid EvSmoke2yr SMOKENOW A_SMOKE N_CIGDAY N_YSMOKE N_QUITSMOKE rename=(EvSmoke2yr=EVSMOKE));
by ccssid;
if EVSMOKE=2 then smokenow=2; /*Eversmoke lifetime=No then smokenow=No*/
if smokenow=1 and EVSMOKE=. then do;  EVSMOKE=1; d_smokeb=d_compq; end; /*some answered (5 people) Now smoking but missing in ever smoked. Use the date now as the date for ever smoking*/
run;
proc freq data=smkbase; table EVSMOKE*smokenow/norow nocol nopercent missing; run;

/*In Fu5, "N7 Smoked at least 100 cigarettes in previous two years"   ==done above since it has original + expansion */
data smoke;
merge Ebase.ebasea(keep=ccssid d_compq in=aa) 
	  smkbase(rename=(EVSMOKE=evsmokeb SMOKENOW=smokenowb a_smoke=a_smokeb) in=aa)
	  smkfu5(rename=(SMOKESINCE=smokefu5 SmokeNow=smokenowfu5 a_smoke=a_smokefu5))
	  Eothers.Ecasestat(keep=ccssid d_birth)
	  qi.death(keep=ccssid d_death)
;
by ccssid;
if aa;
if month(d_birth)=2 and day(d_birth)=29 then d_birth=mdy(3,1,year(d_birth));
/*if Yes to ever smoking in any questionnaire, then Yes.
Among those who answered smkoing in any questionnaire and excluding any Yes, make it No === not perfectly but this way we reduce missing in smoking variable*/
if evsmokeb=1 or smokefu5=1 then evsmoke=1;
else if evsmokeb=2 or smokefu5=2 then evsmoke=2; 
/*date started smoking*/
/*Take from the first questionnaire that gave smoking age. The COALSCE function returns the first nonmissing value in a list of variables.*/
if evsmoke=1 then a_smoke=coalesce(a_smokeb,a_smokefu5);
d_smoke=mdy(month(d_birth),day(d_birth),year(d_birth)+a_smoke)+365.25/2;
if d_smoke=. and (d_smokeb^=. or d_smokefu5^=.) then d_smoke=coalesce(d_smokeb,d_smokefu5); /**for those who answered smoking now but no to ever smoking**/
format d_smoke date9.;
run;
data smoke;
set smoke;
/*If ever smoking Yes but age missing, give a random date between 18 and the baseline date.
But what if the questionaire is after death? So should be tween 18 years and min(questionnaire,death)*/
if evsmoke=1 and a_smoke=. then do; 
	aa=ranuni(2016); 
	d18=mdy(month(d_birth),day(d_birth),year(d_birth)+18); 
	if evsmokeb=1 then d_smoke=d18+(min(d_compq,d_death)-d18)*aa;
	else if smokefu5=1 then d_smoke=d18+(min(d_fu5,d_death)-d18)*aa;
end;
format d_smoke date9.;
run;
	data see; set smoke(where=(evsmoke=1 and d_smoke=.)); run;
	/*check about unreasonable ages at smoking*/
	data see;
	set smoke;
	age_smoke=(d_smoke-d_birth)/365.25;
	run;
	proc means data=see; var age_smoke; run;
	proc sort data=see(where=(evsmoke=1 and age_smoke<18)); by age_smoke; run;
/*For those who are No at baseline but yes at Fu5, change the age to be after baseline but before Fu5.*/
data smoke;
set smoke;
if evsmokeb=2 and smokefu5=1 then d_smoke=d_compq+(d_fu5-d_compq)/2;
run;
/*For those who had smoking age <=14, change it to be 14 to be reasonable*/
data smoke;
set smoke;
if evsmoke=1 then do;
	if (d_smoke-d_birth)/365.25<14 then d_smoke=mdy(6,1,year(d_birth)+14);
end;
run;
data smoke;
set smoke;
/*Fu5 did not ask eve smoke, but smoking in the last 2 years. So those ever smoke before these Q's should be yes to ever smoke.*/
evsmokef5=smokefu5; if evsmokeb=1 then evsmokef5=1;
array varnow(*) smokenowb smokenowfu5;
array varever(*) evsmokeb evsmokef5;
array varcat(*) smkcatb smkcatf5;
do i=1 to dim(varnow);
	if varnow(i)=1 then varcat(i)=3; /*current*/
	else if varever(i)=1 and varnow(i) in (2) then varcat(i)=2; /*former: I prefer to make ever smoke but smoking now missing as formaer instead of missing.
	but if don't use varnow(i)=2, those who did not have FU5 with ever smoking to FU2007 would be made to former --- 
	since I made those wihtout FU5 but with yes to ever smoking in FU2007 as yes to ever smoking in FU5.*/
	else if varever(i)=2 then varcat(i)=1; /*never*/
end;
drop i;
run;
proc freq data=smoke;
table (evsmokeb smokenowb)*smkcatb (evsmokef5 smokenowfu5)*smkcatf5 evsmokef5*smokefu5 smkcatf5/norow nocol nopercent missing;
run;

data smokeall;
set smk(keep=ccssid smkcat: evsmoke: d_smoke evsmoke smokenow: d_fu5 in=aa) 
	smoke(keep=ccssid smkcat: evsmoke: d_smoke evsmoke smokenow: d_fu5 in=bb);
if aa then group=1; else if bb then group=2;
run;
proc sort data=smokeall; by ccssid; run;
data smokeall;
set smokeall;
/*evsmoke was made as 1 if any Q says yes. so no means never in all Q's 
		====>but does not mean a person with Fu2007 but no FU5 can be made as No at FU5
if evsmoke=2 then do;
	smkcatb=1; smkcatf2=1; smkcatf07=1; smkcatf5=1;
end;*/
run;
proc freq data=smokeall; table smkcat:; run;
	/*people with smkcatf5 but no FU5?  === these are those who said yes to ever smoking in previous quesitonnaire, so should be yes to FU5 on ever smoke.*/
	data see; set smokeall; if smkcatf5^=. and d_fu5=.; run;
	proc freq data=see; table smkcatf5; run;

	/*to reduce missing
	Never (1) in a previous Q can became former/current (2/3) in later Q, but not vice versa === 2/3 in baseline cannot be 1 in later Q's. 
	============Tried it here but not used. done in "1D" code since there we used d_startQ.*/
	data smokeall2;
	set smokeall;
	/*In order to reduce missing, take from previous or next Q values.
		missing so would take baseline value. for expansion, the middle 3 Q's would be filled in.*/
	array basevalue{*} smkcatb;	array f1value{*} smkcatf1; 
	array f2value{*} smkcatf2;	array f3value{*} smkcatf07; array fu5value{*} smkcatf5;
	do i=1 to dim(basevalue);
		if f1value(i)=. and basevalue(i)^=. then f1value(i)=basevalue(i);
		if f2value(i)=. and f1value(i)^=. then f2value(i)=f1value(i);
		if f3value(i)=. and f2value(i)^=. then f3value(i)=f2value(i);
		if fu5value(i)=. and f3value(i)^=. then fu5value(i)=f3value(i);
	end; drop i;
	/*after this, some still have missing in one questionnaire and non-missing in another questionnaire. 
	do from the latest to the previous (opposite to the above) to reduce missing ==== kind of cannot do here, 
	since never can became former/current, but not vice versa.*/
	do j=1 to dim(basevalue);
		if fu5value(j)^=. and f3value(j)=. then f3value(j)=fu5value(j);
		if f3value(j)^=. and f2value(j)=. then f2value(j)=f3value(j);
		if f2value(j)^=. and f1value(j)=. then f1value(j)=f2value(j);
		if f1value(j)^=. and basevalue(j)=. then basevalue(j)=f1value(j);	
	end; drop j;
	run;
	proc freq data=smokeall2; 
	table smkcatb  smkcatf1 smkcatf2 smkcatf07 smkcatf5 
		  smkcatb*smkcatf2 smkcatf2*smkcatf07 smkcatf07*smkcatf5/norow nocol nopercent missing; 
	run;
	/*Never can became former/current, but not vice versa.  1 (Never) should be decreasing from baseline to FU5*/   
	data why;
	set smokeall2;
	if (smkcatb in (2,3) and smkcatf2=1) or (smkcatf2 in (2,3) and smkcatf07=1) or (smkcatf07 in (2,3) and smkcatf5=1);
	run;
	proc freq data=smokeall2; table evsmokef5*smkcatf5 evsmokef2*smkcatf2 evsmokef07*smkcatf07 evsmokeb*smkcatb /norow nocol nopercent missing; run;
	/*All these are those who ever smoked but had missing status for former/current  
	===== if using the carry forward on smkcat then resulted in conflicated values from ever smoked. **/
	data wrong;	set smokeall2; if smkcatf5=1 and evsmokef5=1; run;
	/*Before using the carry forward way, give a value for those known to be every smoked yes but missing former/current*/
	data smokeall2;
	set smokeall;
	if evsmokef5=1 and smkcatf5=. then  smkcatf5=10;
	if evsmokef2=1 and smkcatf2=. then  smkcatf2=10; 
	if evsmokef07=1 and smkcatf07=. then  smkcatf07=10;
	if evsmokeb=1 and smkcatb=. then  smkcatb=10;
	run;
	proc freq data=smokeall2; 
	table evsmokef5*smkcatf5 evsmokef2*smkcatf2 evsmokef07*smkcatf07 evsmokeb*smkcatb /norow nocol nopercent missing; 
	run;
	data wrong;	set smokeall2; if smkcatf5=1 and evsmokef5=1; run;
	data chk;
	set smokeall;
	if ccssid in (01002435 10163302 15171915 20028225);
	run; 
	data check;
	set smokeall;
	if ccssid in (01002813, 02165786, 02516962, 03014834, 04073132, 08217561, 08218698, 16373681, 17050414, 18164188, 26021556);

run;

/*Exported for Achal*/
data whq.smoke;
set smokeall;
smk_mostrecent=coalesce(smkcatf5, smkcatf07, smkcatf2, smkcatb);
label evsmoke="Ever smoked";
label d_smoke="Smoke date";
label smkcatb="Smoking status from Baseline (1=Never, 2=Former, 3=Current)";
label smkcatf2="Smoking status from FU2 (1=Never, 2=Former, 3=Current)";
label smkcatf07="Smoking status from FU2007 (1=Never, 2=Former, 3=Current)";
label smkcatf5="Smoking status from FU5 (1=Never, 2=Former, 3=Current)";
label smk_mostrecent="Smoking status from most recent data (1=Never, 2=Former, 3=Current)";
format evsmoke YESNOF.;

*All these are those who ever smoked but had missing status for former/current (d_smoke^=. or evsmoke=1) and smk_mostrecent=1;
if ccssid in (01002813, 02165786, 02516962, 03014834, 04073132, 08217561, 08218698, 16373681, 17050414, 18164188, 26021556) then do;
	smk_mostrecent=2;
end;
run;
proc sort data=whq.smoke; by ccssid; run;
proc freq data=whq.smoke; tables smkcatb smkcatf2 smkcatf07 smkcatf5 smk_mostrecent; run;

data smoke; set whq.smoke(keep=ccssid evsmoke d_smoke smk: ); run;
ods excel file="&mydata.\13. Requet from Achal and Yadav\ExportedCCSS_smoke.xlsx";
		ods excel options(sheet_name="Exported" sheet_interval='none');
		proc print data=smoke noobs; run;

		%add_dummy_table();
		ods excel options(sheet_name="Dictionary" sheet_interval='none');
		ods select Position;
		proc contents data=smoke order=varnum; run; 	/*Export data dictionary*/
ods excel close;



/*----------------------|
|   Physical Activity   |
-----------------------*/

/*
Revised Ocd-14-2021:
Previously, used only vigorous activities in FU2&Fu5&expansion original, and assumed 20 for those whose minutes are 20+, and 
did not use moderateactivities.  --- this was when it was needed to compare with the original cohort baseline.
Now it is used an an X (time-dependent), use all the information to get MET?
*/
data pebase;
merge obase.basea(keep=ccssid d_compq in=aa)  
	  obase.basen(keep=ccssid exer7day)/**vigorous exer at baseline***/;
by ccssid;
MET_vig=exer7day*20*9/60; *20 minuts * days ;
run;
proc freq data=pebase; table MET_vig; run;

/*Qi: The FU questionnaire asked the physical activity differently (more precisely on how many minutes while the baseline questionnaire asked Y/N to 20 minutes and 
no precise minutes available). Given this, Nan's paper with Lee Jones on the mortality outcome 
(https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6181767/) used the baseline exercise as the main exposure. 
However, the secondary exposure was the change in exercise 
and its effect on the outcome so it used both the baseline questionnaire and FU questionnaire, by assuming 20+ minus to be 20 so it was comparable to the baseline exercise.*/
proc freq data=ofu2.f2main; table physical*vigact d_vigact t_vigact /norow nocol nopercent; run;
/*Why some people said No to plysical activity said Yes to vigorous activity? !
Both Yan and Nan did not use "physical"*/
data f2f5;
merge ofu2.f2main(keep=ccssid physical /*modact*/ vigact /*d_modact*/ d_vigact /*t_modact*/ t_vigact 
				  rename=(VIGACT=VIGACT_f2 D_VIGACT=D_VIGACT_f2 T_VIGACT=T_VIGACT_f2))
	  cfu5.fu5n(keep=ccssid VIGACT D_VIGACT T_VIGACT rename=(VIGACT=VIGACT_fu5 D_VIGACT=D_VIGACT_fu5 T_VIGACT=T_VIGACT_fu5));
by ccssid; 	
/*Foup2 */
if vigact_f2=2 then do; t_vigact_f2=0; d_vigact_f2=0; end;
if T_VIGACT_f2>=20 then dayf2 = D_VIGACT_f2; 
else if T_VIGACT_f2^=. and T_VIGACT_f2 <20 then dayf2=0;
met_vig_f2=20*9/60*dayf2; /*assumend 20 minutest for those reproted 20+ minutes ----to be comparable with the baseline.*/
/*FU5*/
if vigact_fu5=2 then do; t_vigact_fu5=0; d_vigact_fu5=0; end;
if T_VIGACT_fu5>=20 then dayf5 = D_VIGACT_fu5; 
else if T_VIGACT_fu5^=. and T_VIGACT_fu5 <20 then dayf5=0;
met_vig_fu5=20*9/60*dayf5; /*assumend 20 minutest for those reproted 20+ minutes ----to be comparable with the baseline.*/
run;

/*##################### Expansion cohort #####################*/
proc means data=cfu5.fu5n; class VIGACT; var T_VIGACT; run;
data exp;
merge ebase.ebasea(keep=ccssid in=aa)
	  ebase.ebaseo(keep=ccssid exer7day) /**exercise at baseline****/
	  cfu5.fu5n(keep=ccssid VIGACT D_VIGACT T_VIGACT);
by ccssid;
if aa;
if exer7day=0 then exer7day=1;/**5 ccssid ***/
MET_vig=(exer7day-1)*20*9/60; /*different format from original cohort*/

if VIGACT=2 then do; D_VIGACT=0; T_VIGACT=0; end;
if T_VIGACT>=20 then MET_vig_fu5=D_VIGACT*20*9/60; /*only for 20+ minutes, report MET; others 0*/
else if T_VIGACT<20 and T_VIGACT^=. then MET_vig_fu5=0;
run;
proc freq data=exp; table MET_vig MET_vig_fu5; run;

data exercise;
merge pebase(keep=ccssid MET_vig) f2f5(keep=ccssid met_vig_fu5 met_vig_f2) exp(keep=ccssid MET_vig met_vig_fu5);
by ccssid;
run;

/*FU6 contains physical activity.  ==== Add OCT 14-2021.*/
data pefu6;
merge cfu6.fu6longa(keep=ccssid d_fu6 in=aa)
	  cfu6.fu6longD(keep=ccssid VIGACT D_VIGACT T_VIGACT);
by ccssid;
if aa;
if VIGACT=2 then do; D_VIGACT=0; T_VIGACT=0; end;
if T_VIGACT>=20 then MET_vig_fu6=D_VIGACT*20*9/60; /*only for 20+ minutes, report MET; others 0*/
else if T_VIGACT<20 and T_VIGACT^=. then MET_vig_fu6=0;
run;
data exercise;
merge exercise(in=aa) pefu6(keep=ccssid MET_vig_fu6);
by ccssid;
if aa;
run;

/*#########detaild MET from Expansion cohort baseline also not available -- only vigorous for 20 minutes. ########*/

/*Real MET minutes >=75 vigorous or 150 combined.============ Based on Yan's code
====Qi: Stephanie only gave me the score using MET:
	0 (score = 0),	3-6 (score = 0.5),	9-12 (score = 1),	15-21 (score = 1)  (vigorous with 20mins).
	Even that I calculte the real MET, don't know how to determine the score based on the real MET minutes. ==could make score=1 for meeting CDC 
	guideline (>=75 or >=150 combined), and then for those not meet CDC, make into 50% as score=0 and 50% as score=0.5?
if one of the sub components already meet CDC guideline, take it as meeting guideline
*/
data metfu2;
merge ofu2.f2main(keep=ccssid d_fu2 physical modact vigact d_modact d_vigact t_modact t_vigact in=aa)
	  oothers.casestat(keep=ccssid d_birth);
by ccssid;
if aa;
ageques=(d_fu2-d_birth)/365.25;

if vigact=2 then do; t_vigact=0; d_vigact=0; end;
if modact=2 then do; t_modact=0; d_modact=0; end;
t_eqmod=d_modact*t_modact+2*d_vigact*t_vigact;

if ageques ge 18 then do;
	if t_eqmod=. then do;
	    if d_modact*t_modact>=150 or d_vigact*t_vigact>=75 then cdc_fu2=1; 
	      else cdc_fu2=.;
	end;
	else if t_eqmod>=150 then cdc_fu2=1;  /*the meeting guidlance ones are score=1.*/
	else cdc_fu2=2;
end;
else do;
	/*For adolescents, Most of the 60 or more minutes a day
	should be either moderate- or vigorous-intensity
	aerobic physical activity, and should include vigorous-intensity
	physical activity at least 3 days a week.*/
	if sum(t_vigact,t_modact)=. or d_vigact=. then cdc_fu2=.;
	else if d_vigact>=3 and sum(t_vigact,t_modact)>=60 then cdc_fu2=1;
	else cdc_fu2=2;
end;
if cdc_fu2=. and physical=2 then cdc_fu2=2;

label cdc_fu2='Met CDC guidelines for physical activity at FU2003(total min)';
format cdc_fu2 yesnof.;

t_mod_fu2=d_modact*t_modact;
t_vig_fu2=d_vigact*t_vigact;
run;
proc freq data=metfu2; table cdc_fu2; run;
proc freq data=metfu2(where=(cdc_fu2=2)); table t_eqmod; run;
data wrong; set metfu2; if t_eqmod>150 and cdc_fu2=2; run; /*<18, to meet CDC, needs at least 3 days a week of vigorous activity. They have 0-2 days*/
/*Among those not meet CDC, 40 cutoff is about half. Also, original baseline asked for at least 20 minues vigorous, 40 is like 
1 day of less than 20 minutes vigorous: since vigorous contributed (2*d_vigact*t_vigact) in equivalent MET calculation. 
======*/
/*who are the no meet CDC but MET missing?*/
data see; set metfu2(where=(cdc_fu2=2 and t_eqmod=.)); run;
data metfu2;
set metfu2;
if cdc_fu2=1 then score_pef2=1;
else if cdc_fu2=2 then do;
	if t_eqmod^=. and t_eqmod>40 then score_pef2=0.5;
	else score_pef2=0; /*t_eqmod<=40, and those t_eqmod missing but No to physical*/
end;
run;
data exercise2;
merge exercise(in=aa) metfu2;
by ccssid;
if aa;
run;
proc freq data=exercise2; table met_vig_f2*score_pef2 /norow nocol nopercent missing; run;
/*While the MET_vig 9+ remain OK, many MET_vig=0 are actually those who met CDC by using all activities.
=== These are people who reported no vigorous activities, but many days/minutes with moderate activities. Should be score=1 meeting CDC.
===== This made me decide that I needed to use all exercise information except for the original cohort who did not have moderate exercise.*/
data see;
set exercise2;
if met_vig_f2=0 and score_pef2=1;
run;
data exercise;
merge exercise(in=aa) metfu2(keep=ccssid cdc_fu2 score_pef2 t_eqmod t_mod_fu2 t_vig_fu2 rename=(t_eqmod=t_eqmodfu2));
by ccssid;
if aa;
run;

/*#########detaield MEt from fU5 ########*/
data metfu5;
merge cfu5.fu5a(keep=ccssid d_fu5)
  	  cfu5.fu5n(keep=ccssid physical modact vigact d_modact d_vigact t_modact t_vigact in=aa)
  	  oothers.casestat(keep=ccssid d_birth) eothers.Ecasestat(keep=ccssid d_birth);
by ccssid;
if aa;
ageques=(d_fu5-d_birth)/365.25;
if ageques>=18; /*adultes only*/

if vigact=2 then do; t_vigact=0; d_vigact=0; end;
if modact=2 then do; t_modact=0; d_modact=0; end;
t_eqmod=d_modact*t_modact+2*d_vigact*t_vigact;

if t_eqmod=. then do;
    if d_modact*t_modact>=150 or d_vigact*t_vigact>=75 then cdc_fu5=1; 
    else cdc_fu5=.;
end;
else if t_eqmod>=150 then cdc_fu5=1;  /*the meeting guidlance ones are score=1.*/
else cdc_fu5=2;
if cdc_fu5=. and physical=2 then cdc_fu5=2;

t_mod_fu5=d_modact*t_modact;
t_vig_fu5=d_vigact*t_vigact;
run;
data metfu5;
set metfu5;
if cdc_fu5=1 then score_pef5=1;
else if cdc_fu5=2 then do;
	if t_eqmod^=. and t_eqmod>40 then score_pef5=0.5;
	else score_pef5=0; /*t_eqmod<=40, and those t_eqmod missing but No to physical*/
end;
run;
data exercise;
merge exercise(in=aa) metfu5(keep=ccssid cdc_fu5 score_pef5 t_eqmod t_mod_fu5 t_vig_fu5 rename=(t_eqmod=t_eqmodfu5));
by ccssid;
if aa;
run;
proc freq data=exercise; table met_vig_fu5*score_pef5 /norow nocol nopercent missing; run;


/*#########detaield MEt from fU6 ########*/
data metfu6;
merge cfu6.fu6longa(keep=ccssid d_fu6 in=aa)
	  cfu6.fu6longD(keep=ccssid physical modact vigact d_modact d_vigact t_modact t_vigact)
	  oothers.casestat(keep=ccssid d_birth) eothers.Ecasestat(keep=ccssid d_birth);
by ccssid;
if aa;
ageques=(d_fu6-d_birth)/365.25;
if ageques>=18; /*adultes only*/

if vigact=2 then do; t_vigact=0; d_vigact=0; end;
if modact=2 then do; t_modact=0; d_modact=0; end;
t_eqmod=d_modact*t_modact+2*d_vigact*t_vigact;

if t_eqmod=. then do;
    if d_modact*t_modact>=150 or d_vigact*t_vigact>=75 then cdc_fu6=1; 
    else cdc_fu6=.;
end;
else if t_eqmod>=150 then cdc_fu6=1;  /*the meeting guidlance ones are score=1.*/
else cdc_fu6=2;
if cdc_fu6=. and physical=2 then cdc_fu6=2;

t_mod_fu6=d_modact*t_modact;
t_vig_fu6=d_vigact*t_vigact;
run;
data metfu6;
set metfu6;
if cdc_fu6=1 then score_pef6=1;
else if cdc_fu6=2 then do;
	if t_eqmod^=. and t_eqmod>40 then score_pef6=0.5;
	else score_pef6=0; /*t_eqmod<=40, and those t_eqmod missing but No to physical*/
end;
run;
data exercise;
merge exercise(in=aa) metfu6(keep=ccssid cdc_fu6 score_pef6 t_eqmod t_mod_fu6 t_vig_fu6 rename=(t_eqmod=t_eqmodfu6));
by ccssid;
if aa;
run;
proc freq data=exercise; table met_vig_fu6*score_pef6 /norow nocol nopercent missing; run;

/*For FU2, FU5, and FU6 that I calcualted both score based on 0 (score = 0),	3-6 (score = 0.5),	9-12 (score = 1),	15-21 (score = 1)  (vigorous with 20mins);
 and the detailed MET minutes from both vigorous + modearte, I will take the higher/better score.*/
data exercise;
set exercise;
/*use 0 (score = 0),	3-6 (score = 0.5),	9-12 (score = 1),	15-21 (score = 1) to get scores*/
array met20(*) met_vig met_vig_f2 met_vig_fu5 met_vig_fu6;
array score20(*) scoreb scoref2 scoref5 scoref6;
do i=1 to dim(met20);
	if met20(i)>=9 then score20(i)=1;
	else if met20(i)>=3 and  met20(i)<=6 then score20(i)=0.5;
	else if met20(i)=0 then score20(i)=0;
end; drop i;
run;
proc freq data=exercise;
table score_pef2*scoref2 score_pef5*scoref5 score_pef6*scoref6/norow nocol nopercent missing; run;
run;
/*Using moderate activities, and real vigorous minutes, more people have higher score*/
data exercise;
set exercise(rename=(scoreb=pescoreb));
array score20(*)  scoref2 scoref5 scoref6;
array scorereal(*) score_pef2 score_pef5 score_pef6;
array scorefinal(*) pescoref2 pescoref5 pescoref6;
do i=1 to dim(score20);
	scorefinal(i)=max(score20(i),scorereal(i));
end; drop i;
run;
proc freq data=exercise; table pescoreb pescoref2 pescoref5 pescoref6; run;
/*Original + expansion baseline only used the at least 20 minutes vigorous activities, so had %'s of worse score. 
=== score=1 meaning meet CDC guideline (MET>=9, 75 vigorous activities/combined 150 minutes of vigorous + modearte activities);
score=0.5 (MET 3-6, 41-149 minutes of vigorous + modearte activities)
scpre=0, MET =0 or <40 minutes of vigorous + modearte activities ===> The 40 cutoff can be changed later as needed. */

/*Exported for Achal*/
data whq.physical_activity;
set exercise(keep=ccssid met_: t_: cdc_:);
t_mostrecent=coalesce(t_eqmodfu6, t_eqmodfu5, t_eqmodfu2);
cdc_mostrecent=coalesce(cdc_fu6, cdc_fu5, cdc_fu2);
if cdc_mostrecent=1 and t_mostrecent<150 then t_mostrecent=.; *Missing eqmod (either mod or vig) in the corresponding questionnaire, but can be sure at least 150 minutes;
if cdc_mostrecent^=1 and t_mostrecent>=150 then cdc_mostrecent=1; *They're missing eqmod in FU5, but had values in FU2, assigned them to FU2;
format cdc_fu5 cdc_fu6 cdc_mostrecent YESNOF.;

label t_eqmodfu2="Equivalent MET calculation (Modearte+2*Vigorous) at FU2003(total min)";
*label cdc_fu2="Met CDC guidelines for physical activity at FU2003(total min)";
label t_mod_fu2="Modearte activities minutes at FU2003 (days*minutes per day: d_modact*t_modact)";
label t_vig_fu2="Vigorous activities minutes at FU2003 (days*minutes per day: d_vigact*t_vigact)";

label t_eqmodfu5="Equivalent MET calculation (Modearte+2*Vigorous) at FU5(total min)";
label cdc_fu5="Met CDC guidelines for physical activity at FU5(total min)";
label t_mod_fu5="Modearte activities minutes at FU5 (days*minutes per day: d_modact*t_modact)";
label t_vig_fu5="Vigorous activities minutes at FU5 (days*minutes per day: d_vigact*t_vigact)";

label t_eqmodfu6="Equivalent MET calculation (Modearte+2*Vigorous) at FU6(total min)";
label cdc_fu6="Met CDC guidelines for physical activity at FU6(total min)";
label t_mod_fu6="Modearte activities minutes at FU6 (days*minutes per day: d_modact*t_modact)";
label t_vig_fu6="Vigorous activities minutes at FU6 (days*minutes per day: d_vigact*t_vigact)";

label t_mostrecent="Most recent Equivalent MET (total min)";
label cdc_mostrecent="Met CDC guidelines for physical activity (total min>=150)";
run;
proc sort data=whq.physical_activity; by ccssid; run;

proc freq data=whq.physical_activity; tables cdc_:; run;
proc sql;
select count(distinct ccssid) from physical_activity where cdc_mostrecent=1;
select count(distinct ccssid) from physical_activity where t_mostrecent>=150;
quit;

data physical_activity(drop=met_:); set whq.physical_activity; run;
ods excel file="&mydata.\13. Requet from Achal and Yadav\ExportedCCSS_physical_activity.xlsx";
		ods excel options(sheet_name="Exported" sheet_interval='none');
		proc print data=physical_activity noobs; run;

		%add_dummy_table();
		ods excel options(sheet_name="Dictionary" sheet_interval='none');
		ods select Position;
		proc contents data=physical_activity order=varnum; run; 	/*Export data dictionary*/
ods excel close;

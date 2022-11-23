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

/*1. Extract sex, age at primary cancer, primary cancer diagnosis, age at last contact*/
proc sql;
create table export_demo0 as
select distinct a.ccssid, a.SEX, a.cohort, a.D_BIRTH, a.D_DX, a.a_dx,
				a.d_compq, a.d_fu1, a.d_fu2, a.d_fu3, a.d_fu2007, a.d_fu5, a.fu6type, a.d_fu6, a.LIVEDEAD, coalesce(a.D_DEATH, b.d_death) as d_death format DATE9.,
				coalesce(c.racegroup, d.racegroup) as racegroup format RACEGRPF., coalesce(c.hispgroup, d.hispgroup) as hispgroup format YESNOF.,
				coalesce(e.diagnose, f.diagnose) as diagnose format DIAGN.
from cothers.ages as a
LEFT JOIN qi.death as b ON a.ccssid=b.ccssid
LEFT JOIN obase.basea as c ON a.ccssid=c.ccssid
LEFT JOIN ebase.ebasea as d ON a.ccssid=d.ccssid
LEFT JOIN oothers.casestat as e ON a.ccssid=e.ccssid
LEFT JOIN eothers.ecasestat as f ON a.ccssid=f.ccssid
order by ccssid;
quit;
data export_demo; set export_demo0; 
d_end=min(max(d_compq,d_fu1,d_fu2,/*d_fu3,*/d_fu2007,d_fu5/*,d_fu6*/),D_DEATH); *CC data till FU5;
a_dx=(D_DX-D_BIRTH)/365.25; format a_dx 8.2;
a_end=(d_end-D_BIRTH)/365.25; format a_end 8.2;
format d_end DATE9.;
run;

/*2. Extract cumulative Anthracycline dose, Alkylating agent, Epipodophyllotoxins, Cisplatinum within 5 years*/
data export_drg;
set cothers.drgscore(keep=ccssid anth_yn5 anth_DED5 alk_yn5 alk_CED5 epip_yn5 epipdose5 pt_yn5 pt_cisED5);
run; 

/*3. Extract radiation dose to chest/neck/pelvis/abdomen/brain*/
data export_rt;
set cothers.combinedrtnov2018_updatedapr2021(keep=ccssid chestrt_yn chestmaxrtdose neckrt_yn neckmaxrtdose pelvisrt_yn pelvismaxrtdose 
												  abdomenrt_yn abdmaxrtdose brainrt_yn maxsegrtdose);
* Some "No" (chestrt_yn=2) are scattered dose (SH: scattered higher, SL: scattred low), change them to 0 dose (they are actually dose 2 or 0.2 Gy, but chestrt_yn=2);
if chestrt_yn=2 then chestmaxrtdose=0; 
chestmaxrtdose=chestmaxrtdose/100; /*Gy=dose/100*/ 
if chestmaxrtdose=0 then chestrtgrp="None   ";
else if 0<chestmaxrtdose<20 then chestrtgrp="0-20   ";
else if chestmaxrtdose>=20 then chestrtgrp=">=20   ";
else chestrtgrp="Unknown";
label chestrtgrp="Radiation dose to chest (categorical)";

if neckrt_yn=2 then neckmaxrtdose=0;  
neckmaxrtdose=neckmaxrtdose/100; /*Gy=dose/100*/
if neckmaxrtdose=0 then neckrtgrp="None   ";
else if 0<neckmaxrtdose<11 then neckrtgrp="0-11   ";
else if 11<=neckmaxrtdose<20 then neckrtgrp="11-20   ";
else if 20<=neckmaxrtdose<30 then neckrtgrp="20-30   ";
else if neckmaxrtdose>=30 then neckrtgrp=">=30   ";
else neckrtgrp="Unknown";
label neckrtgrp="Radiation dose to neck (categorical)";

if pelvisrt_yn=2 then pelvismaxrtdose=0; 
pelvismaxrtdose=pelvismaxrtdose/100; /*Gy=dose/100*/
if pelvismaxrtdose=0 then pelvisrtgrp="None   ";
else if 0<pelvismaxrtdose<20 then pelvisrtgrp="0-20   ";
else if pelvismaxrtdose>=20 then pelvisrtgrp=">=20   ";
else pelvisrtgrp="Unknown";
label pelvisrtgrp="Radiation dose to pelvis (categorical)";
 
if abdomenrt_yn=2 then abdmaxrtdose=0; 
abdmaxrtdose=abdmaxrtdose/100; /*Gy=dose/100*/ 
if abdmaxrtdose=0 then abdomenrtgrp="None   ";
else if 0<abdmaxrtdose<30 then abdomenrtgrp="0-30   ";
else if abdmaxrtdose>=30 then abdomenrtgrp=">=30   ";
else abdomenrtgrp="Unknown";
label abdomenrtgrp="Radiation dose to abdomen (categorical)";

if brainrt_yn=2 then maxsegrtdose=0; 
maxsegrtdose=maxsegrtdose/100; /*Gy=dose/100*/ 
if maxsegrtdose=0 then brainrtgrp="None   ";
else if 0<maxsegrtdose<18 then brainrtgrp="0-18   ";
else if 18<=maxsegrtdose<30 then brainrtgrp="18-30   ";
else if maxsegrtdose>=30 then brainrtgrp=">=30   ";
else brainrtgrp="Unknown";
label brainrtgrp="Radiation dose to brain (categorical)";

run;
proc freq data=export_rt; tables chestrtgrp neckrtgrp pelvisrtgrp abdomenrtgrp brainrtgrp; run;

/*4. Extract age at SN diagnosis (years) and SN diagnosis*/
/*data sns;
set qi.sns;
if SMN_before5=1 or days_after>90 then delete;  *people stay in analysis, but just not counted these as events.;
run;*/
proc sql;
create table export_sn as
select distinct a.ccssid, a.d_candx, /*b.nmsc, b.meningioma,*/ b.groupdx3
from qi.smnrecur as a LEFT JOIN qi.sns as b ON a.ccssid=b.ccssid and a.d_candx=b.d_candx
where a.d_candx^=.
order by ccssid;
quit;

/*proc sql;
create table chk as
select * from export_sn
group by ccssid
having count(ccssid)>1;
quit;
data chk2;
merge chk(in=a) export_demo(keep=ccssid diagnose);
by ccssid;
if a;
run;*/
data export_sn; set export_sn; by ccssid; if first.ccssid; run;
proc sql;
select count(distinct ccssid) from export_sn;
select count(ccssid) from export_sn;
quit;

/*5. Extract BMI*/
proc format;
value bmif 1="Nomal/underweight" 2="overweight" 3="obesity" 9="Unknown";
run;
data export_bmi;
set qi.Final_bmi(keep=ccssid cbmi_0 cbmi_2 cbmi_3 cbmi_5);
bmi_mostrecent=coalesce(cbmi_5, cbmi_3, cbmi_2, cbmi_0);
label bmi_mostrecent="BMI from most recent data";
run;

/*6. Combine all*/
proc sql;
create table export_combine as
select distinct a.ccssid, a.cohort, a.sex, a.a_dx, a.diagnose, a.D_BIRTH, a.a_end, 
				b.d_candx, b.groupdx3, (b.d_candx-a.D_BIRTH)/365.25 as a_candx, 
				c.risky_mostrecent, d.bmi_mostrecent,
				e.chestmaxrtdose, e.chestrtgrp, e.neckmaxrtdose, e.neckrtgrp, e.pelvismaxrtdose, e.pelvisrtgrp, 
				e.abdmaxrtdose, e.abdomenrtgrp, e.maxsegrtdose, e.brainrtgrp,
				f.anth_DED5, f.alk_CED5, f.epipdose5, f.pt_cisED5,
				g.cdc_mostrecent, h.smk_mostrecent
from export_demo as a
LEFT JOIN export_sn as b ON a.ccssid=b.ccssid
LEFT JOIN whq.drink as c ON a.ccssid=c.ccssid
LEFT JOIN export_bmi as d ON a.ccssid=d.ccssid
LEFT JOIN export_rt as e ON a.ccssid=e.ccssid
LEFT JOIN export_drg as f ON a.ccssid=f.ccssid
LEFT JOIN whq.Physical_activity as g ON a.ccssid=g.ccssid
LEFT JOIN whq.smoke as h ON a.ccssid=h.ccssid;
quit;

data whq.export_combine(drop=D_BIRTH d_candx);
set export_combine;
label a_dx="Age at primary cancer";
label a_candx="Age at SN diagnosis (years)";
label groupdx3="SN diagnosis";
label diagnose="Primary cancer diagnosis";
label a_end="Age at last contact (years)";
run;

ods excel file="&mydata.\13. Requet from Achal and Yadav\ExportedCCSS_data.xlsx";
	ods excel options(sheet_name="Exported" sheet_interval='none');
	proc print data=whq.export_combine(keep=ccssid--neckrtgrp) noobs; run;

	%add_dummy_table();
	ods excel options(sheet_name="Exported part2" sheet_interval='none');
	proc print data=whq.export_combine(keep=pelvismaxrtdose--smk_mostrecent) noobs; run;

	%add_dummy_table();
	ods excel options(sheet_name="Dictionary" sheet_interval='none');
	ods select Position;
	proc contents data=whq.export_combine order=varnum; run; 	/*Export data dictionary*/
ods excel close;

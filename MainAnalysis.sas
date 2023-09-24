/*Main Analysis with 180 day exposure assessment window*/
data studybase_analysis ichbefore;
set studybase;
where autologous=0
    and familial=0
    and missingdonor=0
    ;
format nexttrans180 yymmdd10.;
entry=firsttrans+180;
nexttrans180=nexttrans+180;
exitdeath=min(deathdate,migdate,'01jan2018'd);
death=(exitdeath=deathdate);
exit=min(deathdate,migdate,'01jan2018'd,nexttrans180, firstin, diadate_tumor);
exitmulti=min(deathdate,migdate,'01jan2018'd,nexttrans180, secondin, diadate_tumor);
ich=(exit=firstin);
ichmulti=(exitmulti=secondin);
tumor=(exit=diadate_tumor);
tumormulti=(exitmulti=diadate_tumor);
time=(exit-entry)/365.24;
timemulti=(exitmulti-entry)/365.24;
timedeath=(exitdeath-entry)/365.24;
if time gt 0;
age=(entry-birthdate)/365.24;
if 5 le age lt 80;
if transfusions=redcells;
year=year(entry);
everplasma=(plasma gt 0);
everplatelets=(platelets gt 0);
if bloodgroup not in ("A+", "A-", "B+", "B-", "AB+", "AB-", "O+", "O-") then bloodgroup="XX";
if multipledonor in (3,4) then output ichbefore; else output studybase_analysis;
run;
​
​
proc phreg data=studybase_analysis nosummary;
strata hospital;
class sex multipledonor(ref='0') indication bloodgroup;
effect yspl=spline(year / naturalcubic knotmethod=equal(5));
effect aspl=spline(age / naturalcubic knotmethod=equal(5));
effect tspl=spline(transfusions  / naturalcubic knotmethod=equal(5));  
model time*ich(0)=multipledonor sex yspl aspl tspl sex*aspl indication bloodgroup /rl;
ods output ParameterEstimates=singleich_multipledonor;
run;
​
proc phreg data=studybase_analysis nosummary;
  strata hospital;
  class sex multipledonor(ref='0') indication bloodgroup;
  effect yspl=spline(year / naturalcubic knotmethod=equal(5));
  effect aspl=spline(age / naturalcubic knotmethod=equal(5));
  effect tspl=spline(transfusions  / naturalcubic knotmethod=equal(5));  
  model timemulti*ichmulti(0)=multipledonor sex yspl aspl tspl sex*aspl indication bloodgroup /rl;
  ods output ParameterEstimates=multi_multipledonor;
  run;
​
​
/*5 year delay*/
​
data studybase_analysis_y5 ichbefore_y5;
set studybase;
where autologous=0
    and familial=0
    and missingdonor=0
    ;
format nexttrans180 yymmdd10.;
entry=firsttrans+1826;
nexttrans180=nexttrans+180;
exitdeath=min(deathdate,migdate,'01jan2018'd);
death=(exitdeath=deathdate);
exit=min(deathdate,migdate,'01jan2018'd,nexttrans180, firstin, diadate_tumor);
exitmulti=min(deathdate,migdate,'01jan2018'd,nexttrans180, secondin, diadate_tumor);
ich=(exit=firstin);
ichmulti=(exitmulti=secondin);
tumor=(exit=diadate_tumor);
tumormulti=(exitmulti=diadate_tumor);
time=(exit-entry)/365.24;
timemulti=(exitmulti-entry)/365.24;
timedeath=(exitdeath-entry)/365.24;
age=(entry-birthdate)/365.24;
if time gt 0;
if 5 le age lt 80;
year=year(entry);
everplasma=(plasma gt 0);
everplatelets=(platelets gt 0);
if transfusions=redcells;
if bloodgroup not in ("A+", "A-", "B+", "B-", "AB+", "AB-", "O+", "O-") then bloodgroup="XX";
if multipledonor in (3,4) then output ichbefore_y5; else output studybase_analysis_y5;
run;
​
proc phreg data=studybase_analysis_y5 nosummary;
  strata hospital;
  class sex multipledonor(ref='0') indication bloodgroup;
  effect yspl=spline(year / naturalcubic knotmethod=equal(5));
  effect aspl=spline(age / naturalcubic knotmethod=equal(5));
  effect tspl=spline(transfusions  / naturalcubic knotmethod=equal(5));  
  model time*ich(0)=multipledonor sex yspl aspl tspl sex*aspl indication bloodgroup /rl;
  run;
​
  proc phreg data=studybase_analysis_y5 nosummary;
  strata hospital;
  class sex multipledonor(ref='0') indication bloodgroup;
  effect yspl=spline(year / naturalcubic knotmethod=equal(5));
  effect aspl=spline(age / naturalcubic knotmethod=equal(5));
  effect tspl=spline(transfusions  / naturalcubic knotmethod=equal(5));  
  model timemulti*ichmulti(0)=multipledonor sex yspl aspl tspl indication bloodgroup /rl;
  run;

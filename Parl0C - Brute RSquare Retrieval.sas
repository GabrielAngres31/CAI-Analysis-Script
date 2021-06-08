%web_drop_table(WORK.PARL0C);

FILENAME REFFILE0 '/folders/myfolders/PARL0C.csv';

PROC IMPORT DATAFILE=REFFILE0
	DBMS=CSV
	OUT=WORK.PARL0C;
	GETNAMES=YES;
RUN;

PROC CONTENTS DATA=WORK.PARL0C; RUN;

* %web_open_table(WORK.PARL0C);
* %web_open_table(WORK.PARL1);

DATA WORK.OutStatesP0C;
	format Accession BEST12. 
	Label1 $BEST12.
	nValue1 D12.3
	ModelID $CHAR12.;
stop;
run;

DATA WORK.OutInterceptsP0C;
	format accession BEST12. 
	StepLabel $CHAR12.
	SBC D12.3
	ModelID $CHAR12.;
stop;
run;

ods noproctitle;
ods graphics / imagemap=on;
ods trace on; 

proc sort data=WORK.PARL0C out=Work.TempDataSorted;
	by Accession;
run;



ODS OUTPUT FitStatistics = _fstat_W;
ODS OUTPUT SelectionSummary = _selsum_W;
PROC glmselect data = Work.TempDataSorted outdesign(addinputvars) = Work.reg_design;
	MODEL fresh_weight = width;
	BY Accession;
RUN; 
ODS OUTPUT CLOSE;
DATA _fstat_W;
	SET _fstat_W;
	ModelID = "Width";
RUN;
DATA _selsum_W;
	SET _selsum_W;
RUN;
DATA _fstat_W (KEEP = Accession Label1 nValue1 ModelID);
SET _fstat_W (WHERE = ((Label1 = "Adj R-Sq") | (Label1 = "SBC"))); 
PROC PRINT DATA = _fstat_W UNIFORM; 
RUN; 
DATA _selsum_W (KEEP = accession StepLabel SBC);
SET _selsum_W (WHERE = (StepLabel = "Intercept")); 
PROC PRINT DATA = _selsum_W UNIFORM; 
RUN; 
DATA WORK.OutStatesP0C;
	SET WORK.OutStatesP0C _fstat_W;
RUN;
DATA WORK.OutInterceptsP0C;
	SET WORK.OutInterceptsP0C _selsum_W;
RUN;




ODS OUTPUT FitStatistics = _fstat_H;
PROC glmselect data = Work.TempDataSorted outdesign(addinputvars) = Work.reg_design;
	MODEL fresh_weight = height;
	BY Accession;
RUN;
ODS OUTPUT CLOSE;
DATA _fstat_H;
	SET _fstat_H;
	ModelID = "Height";
RUN;
DATA _fstat_H (KEEP = Accession Label1 nValue1 ModelID);
SET _fstat_H (WHERE = ((Label1 = "Adj R-Sq") | (Label1 = "SBC"))); 
PROC PRINT DATA = _fstat_H UNIFORM; 
RUN; 
DATA WORK.OutStatesP0C;
	SET WORK.OutStatesP0C _fstat_H;
RUN;


ODS OUTPUT FitStatistics = _fstat_D;
PROC glmselect data = Work.TempDataSorted outdesign(addinputvars) = Work.reg_design;
	MODEL fresh_weight = diameter;
	BY Accession;
RUN;
ODS OUTPUT CLOSE;
DATA _fstat_D;
	SET _fstat_D;
	ModelID = "Diameter";
RUN;
DATA _fstat_D (KEEP = Accession Label1 nValue1 ModelID);
SET _fstat_D (WHERE = ((Label1 = "Adj R-Sq") | (Label1 = "SBC"))); 
PROC PRINT DATA = _fstat_D UNIFORM; 
RUN; 
DATA WORK.OutStatesP0C;
	SET WORK.OutStatesP0C _fstat_D;
RUN;


ODS OUTPUT FitStatistics = _fstat_T;
PROC glmselect data = Work.TempDataSorted outdesign(addinputvars) = Work.reg_design;
	MODEL fresh_weight = thickness;
	BY Accession;
RUN;
ODS OUTPUT CLOSE;
DATA _fstat_T;
	SET _fstat_T;
	ModelID = "Thickness";
RUN;
DATA _fstat_T (KEEP = Accession Label1 nValue1 ModelID);
SET _fstat_T (WHERE = ((Label1 = "Adj R-Sq") | (Label1 = "SBC"))); 
PROC PRINT DATA = _fstat_T UNIFORM; 
RUN; 
DATA WORK.OutStatesP0C;
	SET WORK.OutStatesP0C _fstat_T;
RUN;


ODS OUTPUT FitStatistics = _fstat_WH;
PROC glmselect data = Work.TempDataSorted outdesign(addinputvars) = Work.reg_design;
	MODEL fresh_weight = width * height;
	BY Accession;
RUN;
ODS OUTPUT CLOSE;
DATA _fstat_WH;
	SET _fstat_WH;
	ModelID = "WxH";
RUN;
DATA _fstat_WH (KEEP = Accession Label1 nValue1 ModelID);
SET _fstat_WH (WHERE = ((Label1 = "Adj R-Sq") | (Label1 = "SBC"))); 
PROC PRINT DATA = _fstat_WH UNIFORM; 
RUN; 
DATA WORK.OutStatesP0C;
	SET WORK.OutStatesP0C _fstat_WH;
RUN;


ODS OUTPUT FitStatistics = _fstat_WD;
PROC glmselect data = Work.TempDataSorted outdesign(addinputvars) = Work.reg_design;
	MODEL fresh_weight = width * diameter;
	BY Accession;
RUN;
ODS OUTPUT CLOSE;
DATA _fstat_WD;
	SET _fstat_WD;
	ModelID = "WxD";
RUN;
DATA _fstat_WD (KEEP = Accession Label1 nValue1 ModelID);
SET _fstat_WD (WHERE = ((Label1 = "Adj R-Sq") | (Label1 = "SBC"))); 
PROC PRINT DATA = _fstat_WD UNIFORM; 
RUN; 
DATA WORK.OutStatesP0C;
	SET WORK.OutStatesP0C _fstat_WD;
RUN;


ODS OUTPUT FitStatistics = _fstat_WT;
PROC glmselect data = Work.TempDataSorted outdesign(addinputvars) = Work.reg_design;
	MODEL fresh_weight = width * thickness;
	BY Accession;
RUN;
ODS OUTPUT CLOSE;
DATA _fstat_WT;
	SET _fstat_WT;
	ModelID = "WxT";
RUN;
DATA _fstat_WT (KEEP = Accession Label1 nValue1 ModelID);
SET _fstat_WT (WHERE = ((Label1 = "Adj R-Sq") | (Label1 = "SBC"))); 
PROC PRINT DATA = _fstat_WT UNIFORM; 
RUN; 
DATA WORK.OutStatesP0C;
	SET WORK.OutStatesP0C _fstat_WT;
RUN;


ODS OUTPUT FitStatistics = _fstat_HD;
PROC glmselect data = Work.TempDataSorted outdesign(addinputvars) = Work.reg_design;
	MODEL fresh_weight = height * diameter;
	BY Accession;
RUN;
ODS OUTPUT CLOSE;
DATA _fstat_HD;
	SET _fstat_HD;
	ModelID = "HxD";
RUN;
DATA _fstat_HD (KEEP = Accession Label1 nValue1 ModelID);
SET _fstat_HD (WHERE = ((Label1 = "Adj R-Sq") | (Label1 = "SBC"))); 
PROC PRINT DATA = _fstat_HD UNIFORM; 
RUN; 
DATA WORK.OutStatesP0C;
	SET WORK.OutStatesP0C _fstat_HD;
RUN;


ODS OUTPUT FitStatistics = _fstat_HT;
PROC glmselect data = Work.TempDataSorted outdesign(addinputvars) = Work.reg_design;
	MODEL fresh_weight = height * thickness;
	BY Accession;
RUN;
ODS OUTPUT CLOSE;
DATA _fstat_HT;
	SET _fstat_HT;
	ModelID = "HxT";
RUN;
DATA _fstat_HT (KEEP = Accession Label1 nValue1 ModelID);
SET _fstat_HT (WHERE = ((Label1 = "Adj R-Sq") | (Label1 = "SBC"))); 
PROC PRINT DATA = _fstat_HT UNIFORM; 
RUN; 
DATA WORK.OutStatesP0C;
	SET WORK.OutStatesP0C _fstat_HT;
RUN;


ODS OUTPUT FitStatistics = _fstat_DT;
PROC glmselect data = Work.TempDataSorted outdesign(addinputvars) = Work.reg_design;
	MODEL fresh_weight = diameter * thickness;
	BY Accession;
RUN;
ODS OUTPUT CLOSE;
DATA _fstat_DT;
	SET _fstat_DT;
	ModelID = "DxT";
RUN;
DATA _fstat_DT (KEEP = Accession Label1 nValue1 ModelID);
SET _fstat_DT (WHERE = ((Label1 = "Adj R-Sq") | (Label1 = "SBC"))); 
PROC PRINT DATA = _fstat_DT UNIFORM; 
RUN; 
DATA WORK.OutStatesP0C;
	SET WORK.OutStatesP0C _fstat_DT;
RUN;


ODS OUTPUT FitStatistics = _fstat_WHD;
PROC glmselect data = Work.TempDataSorted outdesign(addinputvars) = Work.reg_design;
	MODEL fresh_weight = width * height * diameter;
	BY Accession;
RUN;
ODS OUTPUT CLOSE;
DATA _fstat_WHD;
	SET _fstat_WHD;
	ModelID = "WxHxD";
RUN;
DATA _fstat_WHD (KEEP = Accession Label1 nValue1 ModelID);
SET _fstat_WHD (WHERE = ((Label1 = "Adj R-Sq") | (Label1 = "SBC"))); 
PROC PRINT DATA = _fstat_WHD UNIFORM; 
RUN; 
DATA WORK.OutStatesP0C;
	SET WORK.OutStatesP0C _fstat_WHD;
RUN;


ODS OUTPUT FitStatistics = _fstat_WHT;
PROC glmselect data = Work.TempDataSorted outdesign(addinputvars) = Work.reg_design;
	MODEL fresh_weight = width * height * thickness;
	BY Accession;
RUN;
ODS OUTPUT CLOSE;
DATA _fstat_WHT;
	SET _fstat_WHT;
	ModelID = "WxHxT";
RUN;
DATA _fstat_WHT (KEEP = Accession Label1 nValue1 ModelID);
SET _fstat_WHT (WHERE = ((Label1 = "Adj R-Sq") | (Label1 = "SBC"))); 
PROC PRINT DATA = _fstat_WHT UNIFORM; 
RUN; 
DATA WORK.OutStatesP0C;
	SET WORK.OutStatesP0C _fstat_WHT;
RUN;


ODS OUTPUT FitStatistics = _fstat_WDT;
PROC glmselect data = Work.TempDataSorted outdesign(addinputvars) = Work.reg_design;
	MODEL fresh_weight = width * diameter * thickness;
	BY Accession;
RUN;
ODS OUTPUT CLOSE;
DATA _fstat_WDT;
	SET _fstat_WDT;
	ModelID = "WxDxT";
RUN;
DATA _fstat_WDT (KEEP = Accession Label1 nValue1 ModelID);
SET _fstat_WDT (WHERE = ((Label1 = "Adj R-Sq") | (Label1 = "SBC"))); 
PROC PRINT DATA = _fstat_WDT UNIFORM; 
RUN; 
DATA WORK.OutStatesP0C;
	SET WORK.OutStatesP0C _fstat_WDT;
RUN;


ODS OUTPUT FitStatistics = _fstat_HDT;
PROC glmselect data = Work.TempDataSorted outdesign(addinputvars) = Work.reg_design;
	MODEL fresh_weight = height * diameter * thickness;
	BY Accession;
RUN;
ODS OUTPUT CLOSE;
DATA _fstat_HDT;
	SET _fstat_HDT;
	ModelID = "HxDxT";
RUN;
DATA _fstat_HDT (KEEP = Accession Label1 nValue1 ModelID);
SET _fstat_HDT (WHERE = ((Label1 = "Adj R-Sq") | (Label1 = "SBC"))); 
PROC PRINT DATA = _fstat_HDT UNIFORM; 
RUN; 
DATA WORK.OutStatesP0C;
	SET WORK.OutStatesP0C _fstat_HDT;
RUN;


ODS OUTPUT FitStatistics = _fstat_WHDT;
PROC glmselect data = Work.TempDataSorted outdesign(addinputvars) = Work.reg_design;
	MODEL fresh_weight = width * height * diameter * thickness;
	BY Accession;
RUN;
ODS OUTPUT CLOSE;
DATA _fstat_WHDT;
	SET _fstat_WHDT;
	ModelID = "WxHxDxT";
RUN;
DATA _fstat_WHDT (KEEP = Accession Label1 nValue1 ModelID);
SET _fstat_WHDT (WHERE = ((Label1 = "Adj R-Sq") | (Label1 = "SBC"))); 
PROC PRINT DATA = _fstat_WHDT UNIFORM; 
RUN; 
DATA WORK.OutStatesP0C;
	SET WORK.OutStatesP0C _fstat_WHDT;
RUN;

proc reg data=Work.reg_design alpha=0.05 plots(only)=(diagnostics residuals 
		observedbypredicted);
	ods select DiagnosticsPanel ResidualPlot ObservedByPredicted;
	model fresh_weight=&_GLSMOD /;
	by Accession;
	run;
quit;

PROC EXPORT DATA = WORK.OUTSTATESP0C
	REPLACE
	DBMS = CSV
	OUTFILE = '/folders/myfolders/PARL0C_STATS.csv';
RUN;

PROC EXPORT DATA = WORK.OUTINTERCEPTSP0C
	REPLACE
	DBMS = CSV
	OUTFILE = '/folders/myfolders/PARL0C_STAT-INTS.csv';
RUN;

proc delete data=Work.TempDataSorted;
run;

proc delete data=Work.reg_design;
run;

ods trace off;
RRvalues <- read.table(text = "Variables	Levels	Cohort	SN_RR	SN_P	SMN_RR	SMN_P	NMSC_RR	NMSC_P	Breast_cancer_RR	Breast_cancer_P	Thyroid_cancer_RR	Thyroid_cancer_P	Meningioma_RR	Meningioma_P	Sarcoma_RR	Sarcoma_P
AgeatDiagnosis	0-4	SJLIFE	Ref	NA	Ref	NA	NA	NA	Ref	NA	Ref	NA	Ref	NA	NA	NA
AgeatDiagnosis	0-4	CCSS	Ref	NA	Ref	NA	NA	NA	Ref	NA	Ref	NA	Ref	NA	NA	NA
AgeatDiagnosis	5-9	SJLIFE	0.79(0.62-0.99)	0.043	0.76(0.58-1.00)	0.050	NA	NA	0.40(0.08-1.99)	0.26	0.25(0.10-0.60)	0.002	0.88(0.60-1.30)	0.519	NA	NA
AgeatDiagnosis	5-9	CCSS	0.72(0.62-0.84)	< .0001	0.79(0.62-1.01)	0.059	NA	NA	1.04(0.55-1.98)	0.89	0.87(0.55-1.39)	0.568	0.57(0.42-0.77)	0.0002	NA	NA
AgeatDiagnosis	10-14	SJLIFE	0.94(0.74-1.18)	0.571	0.95(0.73-1.24)	0.711	NA	NA	4.75(1.86-12.09)	0.001	1.14(0.66-1.98)	0.634	0.58(0.36-0.92)	0.021	NA	NA
AgeatDiagnosis	10-14	CCSS	0.62(0.54-0.72)	< .0001	0.84(0.67-1.05)	0.127	NA	NA	2.11(1.23-3.61)	0.01	0.83(0.54-1.29)	0.408	0.26(0.18-0.37)	< .0001	NA	NA
AgeatDiagnosis	>=15	SJLIFE	0.79(0.61-1.01)	0.065	0.82(0.61-1.09)	0.175	NA	NA	3.47(1.33-9.10)	0.01	0.54(0.28-1.04)	0.067	0.22(0.10-0.49)	0.0002	NA	NA
AgeatDiagnosis	>=15	CCSS	0.59(0.50-0.69)	< .0001	0.72(0.56-0.91)	0.007	NA	NA	2.30(1.35-3.95)	0.002	0.46(0.27-0.78)	0.004	0.16(0.10-0.27)	< .0001	NA	NA
Sex	Female	SJLIFE	1.56(1.33-1.83)	< .0001	1.37(1.14-1.64)	0.001	1.47(1.14-1.88)	0.003	NA	NA	1.82(1.18-2.83)	0.007	1.59(1.15-2.21)	0.005	1.25(0.62-2.51)	0.536
Sex	Female	CCSS	1.58(1.43-1.75)	< .0001	1.88(1.61-2.18)	< .0001	1.18(1.02-1.36)	0.024	NA	NA	1.79(1.29-2.47)	0.0005	1.61(1.25-2.07)	0.000	1.02(0.62-1.70)	0.930
CranialRT,Gy	None	SJLIFE	Ref	NA	Ref	NA	Ref	NA	NA	NA	NA	NA	Ref	NA	NA	NA
CranialRT,Gy	None	CCSS	Ref	NA	Ref	NA	Ref	NA	NA	NA	NA	NA	Ref	NA	NA	NA
CranialRT,Gy	>0to<18	SJLIFE	1.76(0.80-3.87)	0.162	1.65(0.73-3.75)	0.233	2.57(0.57-11.54)	0.217	NA	NA	NA	NA	22.33(4.44-112.36)	0.0002	NA	NA
CranialRT,Gy	>0to<18	CCSS	2.24(1.42-3.54)	0.001	1.42(0.73-2.76)	0.301	4.31(2.43-7.63)	< .0001	NA	NA	NA	NA	15.89(5.82-43.40)	< .0001	NA	NA
CranialRT,Gy	18-29	SJLIFE	2.01(1.64-2.47)	< .0001	1.41(1.11-1.80)	0.005	2.38(1.78-3.18)	< .0001	NA	NA	NA	NA	35.57(15.39-82.17)	< .0001	NA	NA
CranialRT,Gy	18-29	CCSS	1.85(1.62-2.12)	< .0001	0.84(0.66-1.06)	0.144	2.51(2.08-3.03)	< .0001	NA	NA	NA	NA	23.42(14.30-38.35)	< .0001	NA	NA
CranialRT,Gy	>=30	SJLIFE	1.92(1.48-2.49)	< .0001	1.28(0.93-1.76)	0.138	2.22(1.46-3.38)	0.000	NA	NA	NA	NA	51.08(21.45-121.64)	< .0001	NA	NA
CranialRT,Gy	>=30	CCSS	1.88(1.63-2.18)	< .0001	1.08(0.85-1.36)	0.544	1.72(1.38-2.14)	< .0001	NA	NA	NA	NA	41.93(25.58-68.74)	< .0001	NA	NA
CranialRT,Gy	Unknown	SJLIFE	0.31(0.02-5.09)	0.415	0.23(0.01-8.73)	0.432	0.36(0.00-38.30)	0.665	NA	NA	NA	NA	18.83(4.67-75.90)	< .0001	NA	NA
CranialRT,Gy	Unknown	CCSS	0.71(0.15-3.35)	0.669	1.57(0.27-9.29)	0.616	1.68(0.28-10.01)	0.567	NA	NA	NA	NA	6.15(2.77-13.66)	< .0001	NA	NA
AbdomenRT,Gy	None	SJLIFE	Ref	NA	Ref	NA	Ref	NA	NA	NA	NA	NA	NA	NA	NA	NA
AbdomenRT,Gy	None	CCSS	Ref	NA	Ref	NA	Ref	NA	NA	NA	NA	NA	NA	NA	NA	NA
AbdomenRT,Gy	>0to<30	SJLIFE	1.33(1.01-1.77)	0.045	1.36(1.00-1.84)	0.047	2.71(1.68-4.39)	< .0001	NA	NA	NA	NA	NA	NA	NA	NA
AbdomenRT,Gy	>0to<30	CCSS	0.90(0.75-1.08)	0.246	1.09(0.84-1.40)	0.513	1.46(1.07-1.98)	0.016	NA	NA	NA	NA	NA	NA	NA	NA
AbdomenRT,Gy	>=30	SJLIFE	1.61(1.19-2.17)	0.002	1.54(1.11-2.13)	0.009	3.10(1.87-5.11)	< .0001	NA	NA	NA	NA	NA	NA	NA	NA
AbdomenRT,Gy	>=30	CCSS	1.41(1.20-1.65)	< .0001	1.37(1.10-1.70)	0.004	2.63(2.06-3.36)	< .0001	NA	NA	NA	NA	NA	NA	NA	NA
AbdomenRT,Gy	Unknown	SJLIFE	0.54(0.02-15.39)	0.721	0.30(0.00-17.98)	0.562	6.57(0.06-718.55)	0.432	NA	NA	NA	NA	NA	NA	NA	NA
AbdomenRT,Gy	Unknown	CCSS	19.13(0.03-11406.06)	0.365	0.15(0.01-2.11)	0.159	1.08(0.02-68.62)	0.973	NA	NA	NA	NA	NA	NA	NA	NA
PelvisRT,Gy	None	SJLIFE	NA	NA	NA	NA	Ref	NA	NA	NA	NA	NA	NA	NA	NA	NA
PelvisRT,Gy	None	CCSS	NA	NA	NA	NA	Ref	NA	NA	NA	NA	NA	NA	NA	NA	NA
PelvisRT,Gy	>0to<20	SJLIFE	NA	NA	NA	NA	0.16(0.05-0.49)	0.001	NA	NA	NA	NA	NA	NA	NA	NA
PelvisRT,Gy	>0to<20	CCSS	NA	NA	NA	NA	0.86(0.55-1.34)	0.510	NA	NA	NA	NA	NA	NA	NA	NA
PelvisRT,Gy	>=20	SJLIFE	NA	NA	NA	NA	0.79(0.49-1.26)	0.316	NA	NA	NA	NA	NA	NA	NA	NA
PelvisRT,Gy	>=20	CCSS	NA	NA	NA	NA	0.95(0.74-1.22)	0.697	NA	NA	NA	NA	NA	NA	NA	NA
PelvisRT,Gy	Unknown	SJLIFE	NA	NA	NA	NA	0.99(0.01-109.18)	0.998	NA	NA	NA	NA	NA	NA	NA	NA
PelvisRT,Gy	Unknown	CCSS	NA	NA	NA	NA	0.89(0.01-56.84)	0.957	NA	NA	NA	NA	NA	NA	NA	NA
ChestRT,Gy	None	SJLIFE	Ref	NA	Ref	NA	NA	NA	Ref	NA	NA	NA	NA	NA	NA	NA
ChestRT,Gy	None	CCSS	Ref	NA	Ref	NA	NA	NA	Ref	NA	NA	NA	NA	NA	NA	NA
ChestRT,Gy	>0to<20	SJLIFE	1.24(0.76-2.04)	0.387	1.43(0.83-2.46)	0.202	NA	NA	3.88(0.90-16.75)	0.07	NA	NA	NA	NA	NA	NA
ChestRT,Gy	>0to<20	CCSS	1.68(1.32-2.13)	< .0001	1.93(1.38-2.70)	0.0001	NA	NA	3.71(2.28-6.03)	< .0001	NA	NA	NA	NA	NA	NA
ChestRT,Gy	>=20	SJLIFE	1.89(1.47-2.44)	< .0001	2.18(1.66-2.87)	< .0001	NA	NA	3.83(2.28-6.43)	< .0001	NA	NA	NA	NA	NA	NA
ChestRT,Gy	>=20	CCSS	2.09(1.81-2.41)	< .0001	1.89(1.55-2.31)	< .0001	NA	NA	3.67(2.80-4.81)	< .0001	NA	NA	NA	NA	NA	NA
ChestRT,Gy	Unknown	SJLIFE	9.70(1.36-69.05)	0.023	14.71(1.97-109.88)	0.009	NA	NA	9.70(1.19-79.42)	0.03	NA	NA	NA	NA	NA	NA
ChestRT,Gy	Unknown	CCSS	0.11(0.00-63.61)	0.492	4.26(0.60-30.27)	0.147	NA	NA	2.25(1.31-3.87)	0.003	NA	NA	NA	NA	NA	NA
NeckRT,Gy	None	SJLIFE	NA	NA	NA	NA	NA	NA	NA	NA	Ref	NA	NA	NA	NA	NA
NeckRT,Gy	None	CCSS	NA	NA	NA	NA	NA	NA	NA	NA	Ref	NA	NA	NA	NA	NA
NeckRT,Gy	>0to10	SJLIFE	NA	NA	NA	NA	NA	NA	NA	NA	9.22(1.21-70.54)	0.032	NA	NA	NA	NA
NeckRT,Gy	>0to10	CCSS	NA	NA	NA	NA	NA	NA	NA	NA	8.53(3.38-21.55)	< .0001	NA	NA	NA	NA
NeckRT,Gy	11-19	SJLIFE	NA	NA	NA	NA	NA	NA	NA	NA	11.44(5.16-25.39)	< .0001	NA	NA	NA	NA
NeckRT,Gy	11-19	CCSS	NA	NA	NA	NA	NA	NA	NA	NA	3.96(1.90-8.26)	0.0002	NA	NA	NA	NA
NeckRT,Gy	20-29	SJLIFE	NA	NA	NA	NA	NA	NA	NA	NA	13.31(7.78-22.77)	< .0001	NA	NA	NA	NA
NeckRT,Gy	20-29	CCSS	NA	NA	NA	NA	NA	NA	NA	NA	6.02(3.91-9.27)	< .0001	NA	NA	NA	NA
NeckRT,Gy	>=30	SJLIFE	NA	NA	NA	NA	NA	NA	NA	NA	5.12(2.44-10.73)	< .0001	NA	NA	NA	NA
NeckRT,Gy	>=30	CCSS	NA	NA	NA	NA	NA	NA	NA	NA	3.82(2.52-5.77)	< .0001	NA	NA	NA	NA
NeckRT,Gy	Unknown	SJLIFE	NA	NA	NA	NA	NA	NA	NA	NA	0.00(0.00-Inf)	0.978	NA	NA	NA	NA
NeckRT,Gy	Unknown	CCSS	NA	NA	NA	NA	NA	NA	NA	NA	1.58(0.68-3.70)	0.288	NA	NA	NA	NA
Alkylatingagents	None	SJLIFE	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	Ref	NA
Alkylatingagents	None	CCSS	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	Ref	NA
Alkylatingagents	1st	SJLIFE	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	1.41(0.42-4.68)	0.577
Alkylatingagents	1st	CCSS	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	1.11(0.41-2.98)	0.834
Alkylatingagents	2nd	SJLIFE	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	1.35(0.47-3.90)	0.578
Alkylatingagents	2nd	CCSS	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	1.70(0.77-3.77)	0.189
Alkylatingagents	3rd	SJLIFE	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	3.10(1.28-7.53)	0.012
Alkylatingagents	3rd	CCSS	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	3.57(1.90-6.71)	< .0001
Alkylatingagents	Unknown	SJLIFE	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	11.96(1.47-97.05)	0.020
Alkylatingagents	Unknown	CCSS	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	1.60(0.67-3.82)	0.288
Anthracyclines	None	SJLIFE	NA	NA	NA	NA	NA	NA	Ref	NA	NA	NA	NA	NA	NA	NA
Anthracyclines	None	CCSS	NA	NA	NA	NA	NA	NA	Ref	NA	NA	NA	NA	NA	NA	NA
Anthracyclines	1st	SJLIFE	NA	NA	NA	NA	NA	NA	0.33(0.08-1.44)	0.14	NA	NA	NA	NA	NA	NA
Anthracyclines	1st	CCSS	NA	NA	NA	NA	NA	NA	1.15(0.67-1.96)	0.62	NA	NA	NA	NA	NA	NA
Anthracyclines	2nd	SJLIFE	NA	NA	NA	NA	NA	NA	1.65(0.86-3.17)	0.13	NA	NA	NA	NA	NA	NA
Anthracyclines	2nd	CCSS	NA	NA	NA	NA	NA	NA	1.55(1.07-2.26)	0.02	NA	NA	NA	NA	NA	NA
Anthracyclines	3rd	SJLIFE	NA	NA	NA	NA	NA	NA	2.13(1.20-3.75)	0.01	NA	NA	NA	NA	NA	NA
Anthracyclines	3rd	CCSS	NA	NA	NA	NA	NA	NA	2.28(1.66-3.13)	< .0001	NA	NA	NA	NA	NA	NA
Anthracyclines	Unknown	SJLIFE	NA	NA	NA	NA	NA	NA	0.00(0.00-Inf)	0.99	NA	NA	NA	NA	NA	NA
Anthracyclines	Unknown	CCSS	NA	NA	NA	NA	NA	NA	1.41(0.89-2.24)	0.14	NA	NA	NA	NA	NA	NA
Epipodophyllotoxins	None	SJLIFE	Ref	NA	Ref	NA	NA	NA	NA	NA	Ref	NA	Ref	NA	NA	NA
Epipodophyllotoxins	None	CCSS	Ref	NA	Ref	NA	NA	NA	NA	NA	Ref	NA	Ref	NA	NA	NA
Epipodophyllotoxins	1st	SJLIFE	1.62(1.23-2.12)	0.001	1.25(0.89-1.75)	0.204	NA	NA	NA	NA	2.67(1.42-5.01)	0.002	3.30(2.14-5.08)	< .0001	NA	NA
Epipodophyllotoxins	1st	CCSS	1.32(0.95-1.83)	0.102	1.00(0.58-1.75)	0.989	NA	NA	NA	NA	1.69(0.73-3.89)	0.222	2.75(1.43-5.28)	0.002	NA	NA
Epipodophyllotoxins	2nd	SJLIFE	1.27(0.96-1.69)	0.100	1.13(0.80-1.59)	0.476	NA	NA	NA	NA	2.23(1.14-4.39)	0.020	1.41(0.80-2.46)	0.232	NA	NA
Epipodophyllotoxins	2nd	CCSS	1.81(1.29-2.56)	0.001	1.87(1.13-3.09)	0.015	NA	NA	NA	NA	1.95(0.79-4.83)	0.147	1.79(0.79-4.07)	0.162	NA	NA
Epipodophyllotoxins	3rd	SJLIFE	1.14(0.86-1.52)	0.357	1.02(0.72-1.45)	0.925	NA	NA	NA	NA	2.39(1.19-4.78)	0.014	1.49(0.91-2.42)	0.111	NA	NA
Epipodophyllotoxins	3rd	CCSS	1.56(1.12-2.17)	0.008	2.43(1.62-3.65)	< .0001	NA	NA	NA	NA	2.33(1.02-5.35)	0.045	2.00(0.94-4.27)	0.074	NA	NA
Epipodophyllotoxins	Unknown	SJLIFE	0.92(0.13-6.55)	0.931	1.48(0.21-10.64)	0.695	NA	NA	NA	NA	0.00(0.00-Inf)	0.991	0.00(0.00-Inf)	0.975	NA	NA
Epipodophyllotoxins	Unknown	CCSS	1.16(0.88-1.54)	0.295	1.35(0.90-2.02)	0.143	NA	NA	NA	NA	1.12(0.47-2.68)	0.805	1.95(1.04-3.68)	0.039	NA	NA", header = T)


RRvalues.SJLIFE <- RRvalues[grepl("SJLIFE", RRvalues$Cohort),]

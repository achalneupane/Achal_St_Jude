summary_table <- read.table(text = "variables	OR	lower	upper
(5, 10]	0.89	0.37	1.41
(10, 15]	0.71	0.2	1.23
(15, 100]	0.76	0.21	1.31
Male	1.21	0.84	1.57
(0, 100]	0.81	0.07	1.55
(100, 250]	2.18	1.67	2.69
(250, 1e+05]	8.91	8.43	9.4
(0,500]	0.97	0.51	1.43
(500, 1.5e+03]	1.29	0.56	2.03
(1.5e+03, 3.5e+03]	3.08	2.56	3.61
(3.5e+03, 1e+05]	18.64	17.45	19.84
(15, 25]	1.44	0.86	2.02
(25, 35]	2.02	1.42	2.63
(35, 45]	2.73	2.02	3.45
(45, 100]	4.49	3.2	5.79
Dyslipidemia 	3.73	2.55	4.92
Missing	1.18	0.78	1.57
Hypertension	16.67	15.42	17.92
Missing	2.71	2.19	3.23
Diabetes	1.67	0.73	2.6
Missing	1.08	0.66	1.5
PRS_DCM	1	0.8	1.2
PRS_HF	1.09	0.9	1.28
PRS_HCM	0.78	0.57	0.99
PRS_LVESVi	1.25	1.04	1.45
PRS_LVEF	1.09	0.89	1.29", header = T, sep = "\t")


summary_table$study <- summary_table$variables

summary_table$mean = exp((log(summary_table$upper) + log(summary_table$lower))/2)

summary_table <- summary_table[,c('mean', 'lower', 'upper', 'study', 'OR')]

base_data <- summary_table

library(forestplot)
library(gridExtra)

base_data %>%
  forestplot(labeltext = c(study, OR),
             # clip = c(0.3, 18),
             xlog = T) %>%
  fp_set_style(box = "royalblue",
               line = "darkblue",
               summary = "royalblue") %>%
  fp_add_header(study = c("", "Variables"),
                OR = c("", "OR")) %>%
  fp_set_zebra_style("#EFEFEF")





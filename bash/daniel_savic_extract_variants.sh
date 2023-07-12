## Date: 07/22/2023
# To obtain all LD values for a set of SNPs versus one specific SNP, use the --ld-snp command in conjunction with --r2. For example, to get a list of all values for every SNP within 1Mb of rs12345, use the command
#     plink --file mydata 
#           --r2 
#           --ld-snp rs12345 
#           --ld-window-kb 1000 
#           --ld-window 99999 
#           --ld-window-r2 0
# The --ld-window and --ld-window-r2 commands effectively means that output will be shown for all other SNPs within 1Mb of rs12345.
# Similar to the --ld-snp command, but for multiple seed SNPs: to obtain all LD values from a group of SNPs with other SNPs, use the command
#      --ld-snp-list mysnps.txt
# where mysnps.txt is a list of SNPs.
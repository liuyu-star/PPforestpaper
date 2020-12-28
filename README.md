
This repository contain all the needed files to reproduce the JCGS paper "A Projection Pursuit Forest Algorithm for Supervised Classificationcan". Compiling paper.Rnw file which contains text and code you can reproduce the paper (important this include pre-processd results ). You need to include in the same project, paper.Rnw file and the following files:
  
1. preformance_timesWTplyrnew.Rdata : Data to reproduce Figure 5
2. table_raw.Rdata :  data to do Table2 and Figure 6
3. examvar.Rdata  and errorrates_lymphoma.csv:  Data to reproduce  Figure 10 
4. resRNA_new.rds: Data in Table 3
5. rnaseq_data.rds:  Data for Application: RNA-seq gene expression 
6. globalimpo_rna.Rdata, ppf_m.Rdata and averimpo_rna.Rdata : To reproduce Figure 11

If you want to run some of the intermediate results all the complete code is included in a separate file

codeJCGS.R : Code to reproduce  all the intermediate paper results


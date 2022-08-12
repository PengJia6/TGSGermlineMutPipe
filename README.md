# TGSGermlineMutPipe
Germline mutations detection pipeline for long reads


## Quick start
### 1. change the following configure files according to you environment.
   * conf/config.yaml  (control the software you want to use)
   * conf/reference.yaml (configure the reference you want to use)
   * conf/cluster.yaml (configure the threads for each rule )
   * conf/samples.yaml (configure the samples you want to run)
   * rules/software.smk (configure the absolute path of your software)
   

### 2. run with snakemake           

       nohup snakemake -s tgsGermlineMutPipe.smk -j 10 -k --ri >sublog 2>&1 &
       nohup snakemake -s tgsGermlineMutPipe.smk -j 10 -k --ri --cluster "qsub -l nodes=1:ppn=20 -l walltime=999:00:00" >sublog 2>&1 & 
      

## Support tools 
  - deepvaraint 
  - cuteSV
  - sniffles
  - pbsv
  - SVision
  - ...
## Contribution 
   If you want to apply other tools to mutations with long reads, we encourage you to pull a request or email us.
   

## Contact  
 * Peng Jia at Xi'an Jiaotong University (pengjia@stu.xjtu.edu.cn)

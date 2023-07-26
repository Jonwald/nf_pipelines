# Personal test repository for developing nextflow scripts

**variant_calling:** Minimal and fast pipeline for calling short variants (germline only) from targetted amplicons sequenced using the illumina long amplicon protocol. can be used as a generalised quick and dirty variant caller for small sequencing projects.
 - to do:
   - build docker container for the software
   - remove filepaths and generalise
   - create config /parameter files

**rna_seq:** General purpoise RNA-Seq / RNA-exome pipeline configured for running in AWS SSO environment, work in progress, move to seprate repo when production ready
 - to do:
   - finish building contianers and pushing to ECR
   - finalise nextflow script
   - move processes to modules
   - build basespace transfer container


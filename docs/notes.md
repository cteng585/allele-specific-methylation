# File Paths
- TSS file path: 
  - gphost02.bcgsc.ca
  - /projects/nanopore_pog/resources/refTSS_vs4.1/refTSS_hg38_single_best_transcript.bed
- Sample metadata: 
  - gphost02.bcgsc.ca
  - /projects/clininfo/vcsizmok/long_POG/sample_list/R9_R10_samples/nanopore_pog_samples_R9_R10_2024.txt
- Gene metadata:
  - gphost02.bcgsc.ca
  - /projects/clininfo/vcsizmok/long_POG/aDMR/consolidated_results_ts/datatables/consolidated_results_Andrew_ts_and_all_aDMR_results_20241007.xlsx
- Storage
  - /projects/cteng/clininfo

# Notes
- newer gcc install is located at /gsc/software/linux-x86_64-centos7/
  - need to point environment variables to this location
  - could also point pysam to HTSLib locations on GSC instead of needing to install from internal source
    - export HTSLIB_LIBRARY_DIR=/gsc/software/linux-x86_64-centos7/htslib-1.19/lib
    - export HTSLIB_INCLUDE_DIR=/gsc/software/linux-x86_64-centos7/htslib-1.19/include
    - LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/gsc/software/linux-x86_64-centos7/htslib-1.19/lib"

- probably need to set the below env var before any whatshap run 
  - LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/gsc/software/linux-x86_64-centos7/htslib-1.19/lib"

ssh-ed25519 AAAAC3NzaC1lZDI1NTE5AAAAIHJks4F//GTWYOtA37mIJ0fnViVcSpm4j/CFwc1rgqdt christopherteng@dhcp-206-87-211-5.ubcsecure.wireless.ubc.c

SHA256:Z5hX0xWTwOW6rQjo1ylUNZdSIP9Uz2ROztcMBqSrX0s christopherteng@dhcp-206-87-211-5.ubcsecure.wireless.ubc.ca
The key's randomart image is:
+--[ED25519 256]--+
|           o++B=X|
|           .oB./=|
|          . +.=oX|
|         o + .+ .|
|        S *  . . |
|        .*    o  |
|       .o.. E. . |
|      .  +.=...  |
|       .. o...   |
+----[SHA256]-----+

-----BEGIN OPENSSH PRIVATE KEY-----
b3BlbnNzaC1rZXktdjEAAAAABG5vbmUAAAAEbm9uZQAAAAAAAAABAAAAMwAAAAtzc2gtZW
QyNTUxOQAAACByZLOBf/xk1mDrQN+5iCdH51YlXEqZuI/whcHNa4KnbQAAAMBhKPblYSj2
5QAAAAtzc2gtZWQyNTUxOQAAACByZLOBf/xk1mDrQN+5iCdH51YlXEqZuI/whcHNa4KnbQ
AAAEBDWGpXQFqM0Jk1pEbP1wvFIvfTIiu6A7qYPr1SN+km2nJks4F//GTWYOtA37mIJ0fn
ViVcSpm4j/CFwc1rgqdtAAAAO2NocmlzdG9waGVydGVuZ0BkaGNwLTIwNi04Ny0yMTEtNS
51YmNzZWN1cmUud2lyZWxlc3MudWJjLmNhAQI=
-----END OPENSSH PRIVATE KEY-----
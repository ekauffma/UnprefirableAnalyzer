count="1"

g++ -g -Wno-deprecated postan.C -o dummy.exe -I$ROOTSYS/include -L$ROOTSYS/lib `root-config --cflags` `root-config --libs`
[ -d output ] || mkdir output
[ -d error ] || mkdir error
[ -d log ] || mkdir log

for count in Muon0_2023C_v1aa Muon0_2023C_v1ab Muon0_2023C_v1ac Muon0_2023C_v1ad Muon0_2023C_v1ae Muon0_2023C_v1af Muon0_2023C_v1ag Muon0_2023C_v1ah Muon0_2023C_v1ai Muon0_2023C_v1aj Muon0_2023C_v1ak Muon0_2023C_v1al Muon0_2023C_v1am Muon0_2023C_v1an Muon0_2023C_v1ao Muon0_2023C_v1ap; do

cat>Job_${count}.sh<<EOF
#!/bin/bash
cd /afs/cern.ch/work/p/pdas/L1TrigDPG/Prefiring/CMSSW_13_1_0/src/
cmsenv
cd /afs/cern.ch/work/p/pdas/L1TrigDPG/Prefiring/CMSSW_13_1_0/src/test/UnprefirableAnalyzer/condor/

./dummy.exe filelist_${count} /eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/pdas/prefiring/histo_${count}.root
EOF

chmod 755 Job_${count}.sh

cat>condor_${count}<<EOF
Executable = Job_${count}.sh
Output = output/out_${count}.out
Error = error/err_${count}.error
Log = log/log_${count}.log
+JobFlavour = longlunch
getenv = True
queue

EOF

condor_submit condor_${count}

done

exit 0

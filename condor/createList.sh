#ls -ltr /eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/ekauffma/Muon0/crab_231019_Muon0_Run2023C-22Sep2023-v1/231020_142853/0000/ | grep -v log | grep -v failed | awk '{print "/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/ekauffma/Muon0/crab_231019_Muon0_Run2023C-22Sep2023-v1/231020_142853/0000/"$9}' >> temp.txt

sample=$1
dir=$2

ls -ltr ${dir} | grep -v log | grep -v failed | awk -v dir=${dir} '{print dir$9}' >> temp

tail -n +2 temp > filelist_${sample}

rm -rf temp

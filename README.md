# UnprefirableAnalyzer
This EDAnalyzer uses two different methods to study the pre-firing rate of jets in Run 3 CMS data. A jet is pre-fired when originates from BX=-1 but is accepted in BX=0.

The first is the FirstBunchInTrain method, which flags the first filled bunch in a train of filled bunches. 
This ensures that any activity we see in BX=-1 will be due to pre-firing.

The second is the Unprefirable method, which flags when both BX=0 and BX=-3 are filled. 
Due to L1 Trigger rules, if a BX is fired, the subsequent two BX's are vetoed.

Here is a brief outline of the logic used in the analyzer (plugins/UnprefirableAnalyzer.cc):

1. First, each event is checked to see if it passes `HLT_IsoMu20`. If so, we continue with this event.
2. Then, the event is checked to see if it passes the `FirstBunchInTrain` flag.
3. If the above two checks pass, we iterate through the offline jets (`SlimmedJetsPuppi`)
4. For each offline jet, we iterate through the online jets in each BX (-2 through 2)
5. If an offline jet matches to an online jet, we fill the appropriate events.
6. The above procedure is repeated for `UnprefirableEvent`

To use this analyzer, ensure that you have used the command `cmsenv` within a `CMSSW` directory.
One should also compile the code (from within the base directory of this repository) using `scram b`.
If the analyzer compiles without error, you can move forward.

Before running any of the following lines, be sure that you have an active CMS proxy by running `voms-proxy-init --voms cms --valid 72:00:00`.

There is a config file to run on one file locally. To do this, navigate to the UnprefirableAnalyzer/python directory and run `cmsRun ConfFile_cfg.py`

To run on many files, you can submit a crab job by changing the information in the crab config file `crab3_default.py` and change to the desired dataset information in `ConfFile_submit.py`.
To submit the crab job, run the command `crab submit -c crab3_default.py`

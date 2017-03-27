# MiniAODStudies

* Setup: originally tested in 81x, but should work in any other release, where the appropiate branches are available

```
export SCRAM_ARCH=slc6_amd64_gcc530
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsrel CMSSW_8_1_0
cd CMSSW_8_1_0/src/
cmsenv
git clone https://github.com/alberto-sanchez/MiniAODStudies.git .
scram b
```

* Run: (use your favorite input sample)


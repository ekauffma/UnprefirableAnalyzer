import FWCore.ParameterSet.Config as cms

process = cms.Process("L1JetMini")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring('/store/data/Run2023C/ZeroBias/MINIAOD/PromptReco-v4/000/368/423/00000/23a45b1e-da4c-442d-bcb8-df5735c99ef7.root')
                            )

process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange('368822:45-368822:330','368822:332-368822:908','368822:910-368822:991',
                                                                    '368823:1-368823:35', 
                                                                    '368765:1-368765:1386',
                                                                    '368764:1-368764:4',
                                                                    '368763:1-368763:4',
                                                                    '368762:45-368762:238',
                                                                    '368753:1-368753:976',
                                                                    '368752:1-368752:3',
                                                                    '368751:1-368751:3',
                                                                    '368750:1-368750:2',
                                                                    '368749:1-368749:170',
                                                                    '368748:1-368748:6',
                                                                    '368746:38-368746:77',
                                                                    '368724:1-368724:1200',
                                                                    '368723:43-368723:188',
                                                                    '368685:1-368685:1481',
                                                                    '368684:40-368684:135',
                                                                    '368678:1-368678:171',
                                                                    '368676:1-368676:156',
                                                                    '368675:1-368675:50',
								    '368674:1-368674:12',
                                                                    '368672:1-368672:866',
                                                                    '368670:1-368670:255','368670:306-368670:308',
                                                                    '368669:42-368669:410',
                                                                    '368636:51-368636:819',
                                                                    '368613:1-368613:260',
                                                                    '368611:1-368611:247','368611:249-368611:382',
                                                                    '368609:41-368609:1039',
                                                                    '368567:1-368567:31','368567:33-368567:411','368567:414-368567:436',
                                                                    '368566:1-368566:1195',
                                                                    '368548:1-368548:75','368548:77-368548:86','368548:89-368548:90','368548:93-368548:280','368548:282-368548:1044',
                                                                    '368547:1-368547:194','368547:197-368547:702',
                                                                    '368546:1-368546:659',
                                                                    '368489:47-368489:791','368489:793-368489:827','368489:829-368489:863','368489:865-368489:899','368489:901-368489:935','368489:937-368489:971','368489:973-368489:1007','368489:1009-368489:1043','368489:1045-368489:1079','368489:1081-368489:1280','368489:1282-368489:1380','368489:1382-368489:1581','368489:1583-368489:1595',
                                                                    '368454:1-368454:68', '368454:72-368454:1168',
                                                                    '368453:1-368453:5',
                                                                    '368451:1-368451:144',
                                                                    '368423:42-368423:159','368423:161-368423:516',
                                                                    '368412:1-368412:92',
                                                                    '368410:49-368410:1501',
                                                                    '368400:40-368400:271',
                                                                    '368389:1-368389:198',
                                                                   )

process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string('out_hist.root')
                                   )

process.demo = cms.EDAnalyzer('UnprefirableAnalyzer',

                              )


process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')

process.p = cms.Path(process.demo)

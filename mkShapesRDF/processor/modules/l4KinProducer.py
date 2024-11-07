from mkShapesRDF.processor.framework.module import Module
import ROOT


class l4KinProducer(Module):
    def __init__(self):
        super().__init__("l4KinProducer")

    def runModule(self, df, values):
        # Will return indices for (Zlep1, Zlep2, Xlep1, Xlep2) where the Z is the OSSF lepton pair with mass closest to the Z
        # and the X is the remaining pair
        ROOT.gInterpreter.Declare(
            """
            ROOT::RVecI getZXLepIdx(ROOT::RVec<ROOT::Math::PtEtaPhiMVector> lep_p4, ROOT::RVecF lep_pdgId){
                ROOT::RVecI ZXLepIdx{-9999,-9999,-9999,-9999};
                float dM_best =  9999;
                for(int iL = 0; iL < 4; iL++){
                    for(int jL = iL+1; jL < 4; jL++){
                        if (lep_pdgId[iL] + lep_pdgId[jL] != 0 ) continue;
                        float mll = (lep_p4[iL]+lep_p4[jL]).M();
                        if (fabs(mll-91.1876) < dM_best){
                            dM_best = fabs(mll-91.1876);
                            ZXLepIdx[0] = iL;
                            ZXLepIdx[1] = jL;
                        }
                    }
                }
                if (ZXLepIdx[0] != -9999){
                    for (int iL = 0; iL < 4; iL++){
                        if (iL != ZXLepIdx[0] && iL != ZXLepIdx[1]){
                            if (ZXLepIdx[2] == -9999) ZXLepIdx[2] = iL;
                            else ZXLepIdx[3] = iL;
                        }
                    }
                }
                return ZXLepIdx;
            }
            """
        )

        # If the 4 leading leptons are eeee or uuuu, get the other possible OS combinations as ZA / ZB
        # ZA has the mass closer to the Z mass
        # Return RVec is (ZAlep1, ZAlep2, ZBlep1, ZBlep2)
        # NOTE for now running the bugged version, to explicitly match l4Kin
        # The final line should actually be ZAZBLepIdx[3] = ZXLepIdx[3-(iL+k)%2]
        ROOT.gInterpreter.Declare(
            """
            ROOT::RVecI getZAZBLepIdx(ROOT::RVec<ROOT::Math::PtEtaPhiMVector> lep_p4, ROOT::RVecF lep_pdgId, ROOT::RVecF lep_ch, ROOT::RVecI ZXLepIdx){
                ROOT::RVecI ZAZBLepIdx{-9999,-9999,-9999,-9999};
                if (!(lep_ch[0]+lep_ch[1]+lep_ch[2]+lep_ch[3] == 0 && abs(lep_pdgId[0]) == abs(lep_pdgId[1]) && abs(lep_pdgId[1]) == abs(lep_pdgId[2]) && abs(lep_pdgId[2]) == abs(lep_pdgId[3]))) return ZAZBLepIdx;
                int k = 0;
                if(ZXLepIdx[0] != 9999){
                    for (int iL = 0; iL < 2; iL++){
                        if (lep_ch[ZXLepIdx[0]]+lep_ch[ZXLepIdx[iL+2]] == 0){
                            if (fabs((lep_p4[ZXLepIdx[0]]+lep_p4[ZXLepIdx[iL+2]]).M()-91.1876) > fabs((lep_p4[ZXLepIdx[1]]+lep_p4[ZXLepIdx[3-iL]]).M()-91.1876)) k = 1;
                            ZAZBLepIdx[0] = ZXLepIdx[k];
                            ZAZBLepIdx[1] = ZXLepIdx[2+(iL+k)%2];
                            ZAZBLepIdx[2] = ZXLepIdx[1-k];
                            ZAZBLepIdx[3] = ZXLepIdx[2+iL];
                        }
                    }
                }
                return ZAZBLepIdx;
            }
            """
        )

        # Finally a macro to get the minimum delta phi between any OS lepton pair
        ROOT.gInterpreter.Declare(
            """
            float minDeltaPhi(ROOT::RVecF lep_phi, ROOT::RVecF lep_ch){
                float minDeltaPhi = 9999.;
                for (int iL = 0; iL < 4; iL++){
                    for (int jL = iL+1; jL < 4; jL++){
                        if (lep_ch[iL] != lep_ch[jL] && fabs(DeltaPhi(lep_phi[jL],lep_phi[iL])) < fabs(minDeltaPhi)) minDeltaPhi = DeltaPhi(lep_phi[jL],lep_phi[iL]);
                    }
                }
            return minDeltaPhi;
            }
            """
        )

        # Define temporary columns, inputs to the columns to keep
        df = df.Define(
            "Lepton_4DV",
            "ROOT::VecOps::Construct<ROOT::Math::PtEtaPhiMVector>(Lepton_pt, Lepton_eta, Lepton_phi, ROOT::RVecF(Lepton_pt.size(), 0))",
            excludeVariations=["JES*", "MET*"],
        )

        df = df.Define(
            "Lepton_ch",
            "-Lepton_pdgId/abs(Lepton_pdgId)",
            excludeVariations=["JES*", "MET*"],
        )
        df = df.Define(
            "ZXLepIdx",
            "getZXLepIdx(Lepton_4DV, Lepton_pdgId)",
            excludeVariations=["JES*", "MET*"],
        )
        df = df.Define(
            "ZAZBLepIdx",
            "getZAZBLepIdx(Lepton_4DV, Lepton_pdgId, Lepton_ch, ZXLepIdx)",
            excludeVariations=["JES*", "MET*"],
        )

        # Define general preselection, variables will have value -9999 if this is not true
        df = df.Define(
            "isAllOk",
            "PuppiMET_pt > 0 && fabs(PuppiMET_phi) < 3.14159265359 && Lepton_pt[Lepton_pt >= 10].size() == 4 && Lepton_pt[0] > 25 && Lepton_pt[1] > 15 && "
            "Lepton_ch[0]+Lepton_ch[1]+Lepton_ch[2]+Lepton_ch[3] == 0",
        )

        # Variables not depending on JES / MET variations
        vardict = {
            "isAllOk": {
                "minDeltaPhi_zh4l": ["minDeltaPhi(Lepton_phi, Lepton_ch)", "-9999"],
                "chllll_zh4l": [
                    "Lepton_ch[0]+Lepton_ch[1]+Lepton_ch[2]+Lepton_ch[3]",
                    "-9999",
                ],
                "mllll_zh4l": [
                    "(Lepton_4DV[0]+Lepton_4DV[1]+Lepton_4DV[2]+Lepton_4DV[3]).M()",
                    "-9999",
                ],
            },
            "isAllOk && ZXLepIdx[0] != -9999": {
                "z0Mass_zh4l": [
                    "(Lepton_4DV[ZXLepIdx[0]]+Lepton_4DV[ZXLepIdx[1]]).M()",
                    "-9999",
                ],
                "z0Pt_zh4l": [
                    "(Lepton_4DV[ZXLepIdx[0]]+Lepton_4DV[ZXLepIdx[1]]).Pt()",
                    "-9999",
                ],
                "z1Mass_zh4l": [
                    "(Lepton_4DV[ZXLepIdx[2]]+Lepton_4DV[ZXLepIdx[3]]).M()",
                    "-9999",
                ],
                "z1Pt_zh4l": [
                    "(Lepton_4DV[ZXLepIdx[2]]+Lepton_4DV[ZXLepIdx[3]]).Pt()",
                    "-9999",
                ],
                "flagZ1SF_zh4l": [
                    "float(abs(Lepton_pdgId[ZXLepIdx[2]]) == abs(Lepton_pdgId[ZXLepIdx[3]]))",
                    "-9999",
                ],
                "z0DeltaPhi_zh4l": [
                    "DeltaPhi(Lepton_phi[ZXLepIdx[1]], Lepton_phi[ZXLepIdx[0]])",
                    "-9999",
                ],
                "z1DeltaPhi_zh4l": [
                    "DeltaPhi(Lepton_phi[ZXLepIdx[3]], Lepton_phi[ZXLepIdx[2]])",
                    "-9999",
                ],
                "z0DeltaR_zh4l": [
                    "DeltaR(Lepton_eta[ZXLepIdx[0]], Lepton_eta[ZXLepIdx[1]], Lepton_phi[ZXLepIdx[0]], Lepton_phi[ZXLepIdx[1]])",
                    "-9999",
                ],
                "z1DeltaR_zh4l": [
                    "DeltaR(Lepton_eta[ZXLepIdx[2]], Lepton_eta[ZXLepIdx[3]], Lepton_phi[ZXLepIdx[2]], Lepton_phi[ZXLepIdx[3]])",
                    "-9999",
                ],
            },
            "isAllOk && ZAZBLepIdx[0] != -9999": {
                "zaMass_zh4l": [
                    "(Lepton_4DV[ZAZBLepIdx[0]]+Lepton_4DV[ZAZBLepIdx[1]]).M()",
                    "-9999",
                ],
                "zbMass_zh4l": [
                    "(Lepton_4DV[ZAZBLepIdx[2]]+Lepton_4DV[ZAZBLepIdx[3]]).M()",
                    "-9999",
                ],
                "zaDeltaPhi_zh4l": [
                    "DeltaPhi(Lepton_phi[ZAZBLepIdx[1]], Lepton_phi[ZAZBLepIdx[0]])",
                    "-9999",
                ],
                "zbDeltaPhi_zh4l": [
                    "DeltaPhi(Lepton_phi[ZAZBLepIdx[3]], Lepton_phi[ZAZBLepIdx[2]])",
                    "-9999",
                ],
                "zaDeltaR_zh4l": [
                    "DeltaR(Lepton_eta[ZAZBLepIdx[0]], Lepton_eta[ZAZBLepIdx[1]], Lepton_phi[ZAZBLepIdx[0]], Lepton_phi[ZAZBLepIdx[1]])",
                    "-9999",
                ],
                "zbDeltaR_zh4l": [
                    "DeltaR(Lepton_eta[ZAZBLepIdx[2]], Lepton_eta[ZAZBLepIdx[3]], Lepton_phi[ZAZBLepIdx[2]], Lepton_phi[ZAZBLepIdx[3]])",
                    "-9999",
                ],
            },
        }

        for selection in vardict:
            for col in vardict[selection]:
                df = df.Define(
                    col,
                    "{} ? {} : {}".format(
                        selection,
                        vardict[selection][col][0],
                        vardict[selection][col][1],
                    ),
                    excludeVariations=["JES*", "MET*"],
                )

        # Variables including MET (vary with JES / MET variations)
        vardict = {
            "isAllOk": {
                "pfmetPhi_zh4l": ["static_cast<float>(PuppiMET_phi)", "-9999"],
                "lep1Mt_zh4l": [
                    "sqrt(2*Lepton_4DV[0].Pt()*PuppiMET_pt*(1-cos(Lepton_4DV[0].Phi()-pfmetPhi_zh4l)))",
                    "-9999",
                ],
                "lep2Mt_zh4l": [
                    "sqrt(2*Lepton_4DV[1].Pt()*PuppiMET_pt*(1-cos(Lepton_4DV[1].Phi()-pfmetPhi_zh4l)))",
                    "-9999",
                ],
                "lep3Mt_zh4l": [
                    "sqrt(2*Lepton_4DV[2].Pt()*PuppiMET_pt*(1-cos(Lepton_4DV[2].Phi()-pfmetPhi_zh4l)))",
                    "-9999",
                ],
                "lep4Mt_zh4l": [
                    "sqrt(2*Lepton_4DV[3].Pt()*PuppiMET_pt*(1-cos(Lepton_4DV[3].Phi()-pfmetPhi_zh4l)))",
                    "-9999",
                ],
                "minMt_zh4l": [
                    "min(lep1Mt_zh4l,min(lep2Mt_zh4l,min(lep3Mt_zh4l,lep4Mt_zh4l)))",
                    "-9999",
                ],
            },
            "isAllOk && ZXLepIdx[0] != -9999": {
                "z1dPhi_lep1MET_zh4l": [
                    "fabs(DeltaPhi(pfmetPhi_zh4l, Lepton_phi[ZXLepIdx[2]]))",
                    "9999",
                ],
                "z1dPhi_lep2MET_zh4l": [
                    "fabs(DeltaPhi(pfmetPhi_zh4l, Lepton_phi[ZXLepIdx[3]]))",
                    "9999",
                ],
                "z1mindPhi_lepMET_zh4l": [
                    "min(z1dPhi_lep1MET_zh4l,z1dPhi_lep2MET_zh4l)",
                    "9999",
                ],
                "z1Mt_zh4l": [
                    "sqrt(2*(Lepton_4DV[ZXLepIdx[2]]+Lepton_4DV[ZXLepIdx[3]]).Pt()*PuppiMET_pt*(1-cos((Lepton_4DV[ZXLepIdx[2]]+Lepton_4DV[ZXLepIdx[3]]).Phi()-pfmetPhi_zh4l)))",
                    "-9999",
                ],
            },
        }

        for selection in vardict:
            for col in vardict[selection]:
                df = df.Define(
                    col,
                    "{} ? {} : {}".format(
                        selection,
                        vardict[selection][col][0],
                        vardict[selection][col][1],
                    ),
                )

        # Now clean up the temporary columns
        df = df.DropColumns("Lepton_4DV")
        df = df.DropColumns("Lepton_ch")
        df = df.DropColumns("ZXLepIdx")
        df = df.DropColumns("ZAZBLepIdx")
        df = df.DropColumns("isAllOk")

        return df

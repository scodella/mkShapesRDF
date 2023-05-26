import ROOT

from mkShapesRDF.processor.data.LeptonSel_cfg import LepFilter_dict, ElectronWP, MuonWP
from mkShapesRDF.processor.framework.Module import Module


class LeptonSel(Module):
    def __init__(self, LepFilter, nLF):
        super().__init__("LeptonSel")
        self.LepFilter = LepFilter
        self.nLF = nLF

    def runModule(self, df, values):
        Clean_Tag = LepFilter_dict[self.LepFilter]
        Clean_TagWP = LepFilter_dict[Clean_Tag]

        ROOT.gInterpreter.Declare(
            """
        using namespace ROOT;
        using namespace ROOT::VecOps;

        ROOT::RVecB reduce_cond_any(ROOT::RVecB condition, uint size1, uint size2){
            ROOT::RVecB r;
            for (uint i = 0; i < size1; i++){
                bool c = false;
                for (uint j = 0; j < size2; j++){
                    if (condition[i * size2 + j]){
                        c = true;
                        break;
                    }
                }
                r.push_back(c);
            }
            return r;
        }

        RVecB propagateMask(RVecI Lepton_origIdx, RVecB mask, bool defaultValue){
        RVecB r {};

            for (uint i = 0; i < Lepton_origIdx.size(); i++){
                if (Lepton_origIdx[i] < 0){
                    r.push_back(defaultValue);
                } else {
                    r.push_back(mask[Lepton_origIdx[i]]);
                }
            }
            assert (Lepton_origIdx.size() == r.size());
            return r;
        /*
        does not work, why?
        RVecB r(Lepton_origIdx.size(), defaultValue);
        r[Lepton_origIdx >= 0] = Take(mask, Lepton_origIdx[Lepton_origIdx >= 0]);
        return r;
        */
        }
        """
        )
        values.append(
            [
                df.Define("test", "Lepton_pt.size()").Sum("test"),
                "Original size of leptons",
            ]
        )

        df = df.Define("comb", "ROOT::RVecB(Electron_pt.size(), true)")
        df = df.Define("tmp1", "true")
        df = df.Define("tmp2", "true")

        # for key, cuts in list(ElectronWP['Full2018v9']['TightObjWP']['mvaFall17V2Iso_WP90']['cuts'].items()):
        for key, cuts in ElectronWP["Full2018v9"][Clean_TagWP]["HLTsafe"][
            "cuts"
        ].items():
            df = df.Redefine("tmp1", key)
            df = df.Redefine("tmp2", "(" + cuts[0] + ")")
            for cut in cuts[1:]:
                df = df.Redefine("tmp2", "tmp2 && (" + cut + ")")
            df = df.Redefine("comb", "comb && (! tmp1 || tmp2)")

        df = df.Define(
            "LeptonMaskHyg_Ele", "propagateMask(Lepton_electronIdx, comb, true)"
        )
        df = df.Redefine("comb", "ROOT::RVecB(Muon_pt.size(), true)")

        for key, cuts in MuonWP["Full2018v9"][Clean_TagWP]["HLTsafe"]["cuts"].items():
            df = df.Redefine("tmp1", key)
            df = df.Redefine("tmp2", "(" + cuts[0] + ")")
            for cut in cuts[1:]:
                df = df.Redefine("tmp2", "tmp2 && (" + cut + ")")
            df = df.Redefine("comb", "comb && (! tmp1 || tmp2)")

        df = df.Define("LeptonMaskHyg_Mu", "propagateMask(Lepton_muonIdx, comb, true)")

        df = df.Define("LeptonMask_minPt", "Lepton_pt > 8")
        df = df.Define(
            "LeptonMask_minPt_pass",
            "LeptonMask_minPt && LeptonMaskHyg_Ele && LeptonMaskHyg_Mu",
        )

        # FIXME should add Leptons IDs

        # TODO add VetoLeptons and dmZll

        df = df.Filter("Lepton_pt[LeptonMask_minPt_pass].size() >= 1")

        values.append(
            [
                df.Define("test", "Lepton_pt[LeptonMask_minPt_pass].size()").Sum(
                    "test"
                ),
                "Size of leptons passing hyg mask and minPt",
            ]
        )

        # this LeptonMask is used for leptons inside Jet cones
        df = df.Define(
            "LeptonMask_JC",
            "LeptonMaskHyg_Ele && LeptonMaskHyg_Mu && (Lepton_pt >= 10)",
        )

        df = df.Define("CleanJetMask", "CleanJet_eta <= 5.0")

        df = df.Define(
            "CleanJet_Lepton_comb",
            "ROOT::VecOps::Combinations(CleanJet_pt[CleanJetMask].size(), Lepton_pt[LeptonMask_JC].size())",
        )

        df = df.Define(
            "dR2",
            "ROOT::VecOps::DeltaR2( \
            Take(CleanJet_eta, CleanJet_Lepton_comb[0]), \
            Take(CleanJet_phi, CleanJet_Lepton_comb[0]), \
            Take(Lepton_eta, CleanJet_Lepton_comb[1]), \
            Take(Lepton_phi, CleanJet_Lepton_comb[1]) \
        )",
        )

        df = df.Define(
            "CleanJet_pass",
            "! reduce_cond_any(dR2<(0.3*0.3), CleanJet_pt[CleanJetMask].size(), Lepton_pt[LeptonMask_JC].size())",
        )

        branches = ["pt", "eta", "phi", "pdgId", "electronIdx", "muonIdx"]
        for prop in branches:
            df = df.Redefine(
                f"Lepton_{prop}",
                f"Lepton_{prop}[LeptonMaskHyg_Ele && LeptonMaskHyg_Mu]",
            )

        branches = ["jetIdx", "pt", "eta", "phi"]
        for prop in branches:
            df = df.Redefine(
                f"CleanJet_{prop}", f"CleanJet_{prop}[CleanJetMask][CleanJet_pass]"
            )

        df = df.DropColumns("Lepton*Mask*")

        df = df.DropColumns("CleanJetMask")
        df = df.DropColumns("CleanJet_Lepton_comb")
        df = df.DropColumns("CleanJet_pass")

        return df

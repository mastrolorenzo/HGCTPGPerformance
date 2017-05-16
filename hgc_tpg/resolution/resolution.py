# test of a python analyzer that produce response plots #
# Luca Mastrolorenzo - 25-04-2017 #
# The class provide a function that create an output file and store pt,eta,phi response #

import ROOT
import hgc_tpg.utilities.functions as f
import numpy as np

class resolution :
    
    def __init__(self, input_file, output_file, conf) :

        self.outputFile = output_file
        self.chain = ROOT.TChain("hgcalTriggerNtuplizer/HGCalTriggerNtuple")
        with open(input_file,"r") as fileList :
            lines = fileList.read().splitlines()
            for infile in lines :
                self.chain.Add('root://cms-xrd-global.cern.ch/'+infile)
#                self.chain.Add(infile)
        self.cfg = conf

    # function that produce the response plot - looking to dR-match between gen-C3d #
    # keep the highest-pt C3d in the cone as L1-candidate #
    def plotResponse(self) :

        # definition of the output file and histogram to store in it
        output = ROOT.TFile( self.outputFile,"RECREATE" )

        h_resoPt = ROOT.TH1D( "resoPt", "Pt response", 100, 0, 2)
        h_resoEta = ROOT.TH1D( "resoEta", "Eta response", 100, -0.15, 0.15);
        h_resoPhi = ROOT.TH1D( "resoPhi", "Phi response", 100, -0.15, 0.15);
        h_L1PtvsTrue2D = ROOT.TH2D( "h_L1PtvsTrue2D", "h_L1PtvsTrue2D",200, 0, 200, 201, -1, 200);
        h_L1PtResponseVsTruePt = ROOT.TH2D( "h_L1PtResponseVsTruePt", "h_L1PtResponseVsTruePt", 200, 0, 200, 100, 0, 2);
        h_L1PtResponseVsTrueEta = ROOT.TH2D( "h_L1PtResponseVsTrueEta", "h_L1PtResponseVsTrueEta", 200, 1.0, 3.5, 100, 0, 2);
    
        # loop over the ttree ntries
        for ientry, entry in enumerate(self.chain):
            if ientry%1000==0 : 
                print ientry
            if ientry == self.cfg.maxNEvts :
                break
            gen_pt_ = np.array(entry.gen_pt) # for vectors
            gen_eta_ = np.array(entry.gen_eta)
            gen_phi_ = np.array(entry.gen_phi)
            gen_energy_ = np.array(entry.gen_energy)
            gen_status_ = np.array(entry.gen_status)
            gen_id_ = np.array(entry.gen_id)        
            c3d_pt_ = np.array(entry.cl3d_pt)
            c3d_eta_ = np.array(entry.cl3d_eta)
            c3d_phi_ = np.array(entry.cl3d_phi)
            c3d_energy_ = np.array(entry.cl3d_energy)
            
            # select the good C3d
            goodC3d_idx = []
            for i_c3d in range(len(c3d_pt_)) :
                if ( abs( c3d_eta_[i_c3d] ) > self.cfg.minEta_C3d 
                     and  abs( c3d_eta_[i_c3d] ) < self.cfg.maxEta_C3d 
                     and c3d_pt_[i_c3d] > self.cfg.minPt_C3d ) :
                    
                    goodC3d_idx.append(i_c3d)
                    
            # loop over the gen particle and apply basic selection at MC-truth level
            for i_gen in range(len(gen_pt_)) :
                if ( abs( gen_eta_[i_gen] ) > self.cfg.minEta_gen 
                     and abs( gen_eta_[i_gen] ) < self.cfg.maxEta_gen 
                     and gen_pt_[i_gen] > self.cfg.minPt_gen 
                     and abs( gen_id_[i_gen] ) == self.cfg.particle_type 
                     and gen_status_[i_gen] == self.cfg.particle_status ) :
                    
                    hasMatched = False
                    pt_cand = -1.
                    eta_cand = -100.
                    phi_cand = -100.

                    # loop over the good 3D-cluster
                    for i in range(len(goodC3d_idx)) :
                         i_c3d = goodC3d_idx[i]   
                         dR = f.deltaR( gen_eta_[i_gen], c3d_eta_[i_c3d], gen_phi_[i_gen], c3d_phi_[i_c3d] )                
                         if dR < self.cfg.dRmatch :
                             hasMatched = True
                             if c3d_pt_[i_c3d] > pt_cand :
                                 pt_cand = c3d_pt_[i_c3d]
                                 eta_cand = c3d_eta_[i_c3d]
                                 phi_cand = c3d_phi_[i_c3d]                                                      
                                 
                    # fill the response histograms    
                    if hasMatched : 
                        #print gen_pt_[i_gen], pt_cand
                        h_L1PtvsTrue2D.Fill( gen_pt_[i_gen], pt_cand )
                        h_resoPt.Fill( pt_cand / gen_pt_[i_gen] )
                        h_resoEta.Fill( eta_cand - gen_eta_[i_gen] )
                        h_resoPhi.Fill( phi_cand - gen_phi_[i_gen] )
                        h_L1PtResponseVsTruePt.Fill( gen_pt_[i_gen], pt_cand / gen_pt_[i_gen] )
                        h_L1PtResponseVsTrueEta.Fill( abs(gen_eta_[i_gen]), pt_cand / gen_pt_[i_gen] )
                
                    elif not hasMatched :
                        #print gen_pt_[i_gen], -1
                        h_L1PtvsTrue2D.Fill( gen_pt_[i_gen], -1 )
                       
        h_resoPt.Write()
        h_resoEta.Write() 
        h_resoPhi.Write() 
        h_L1PtResponseVsTruePt.Write()
        h_L1PtResponseVsTrueEta.Write()




    def plotResponseTC(self) :
        
        # definition of the output file and histogram to store in it
        output = ROOT.TFile( self.outputFile,"UPDATE" )

        h_resoPt = ROOT.TH1D( "resoPt", "Pt response", 100, 0, 2 )
        h_resoPtVsThreshold = ROOT.TH1D( "h_resoPtVsThreshold", "Pt response vs tc-threshold", 9, -0.125, 2.125 )
        h_AllTCPtResponseVsTruePt = ROOT.TH2D( "h_AllTCPtResponseVsTruePt", "h_AllTCPtResponseVsTruePt", 200, 0, 200, 100, 0, 2 );
        h_AllTCPtResponseVsTrueEta = ROOT.TH2D( "h_AllTCPtResponseVsTrueEta", "h_AllTCPtResponseVsTrueEta", 200, 1.0, 3.5, 100, 0, 2 );
    
        # loop over the ttree ntries
        for ientry, entry in enumerate(self.chain):
            if ientry%1000==0 : 
                print ientry
            if ientry == self.cfg.maxNEvts :
                break

            gen_pt_ = np.array(entry.gen_pt) # for vectors
            gen_eta_ = np.array(entry.gen_eta)
            gen_phi_ = np.array(entry.gen_phi)
            gen_energy_ = np.array(entry.gen_energy)
            gen_status_ = np.array(entry.gen_status)
            gen_id_ = np.array(entry.gen_id)        
            tc_energy_ = np.array(entry.tc_energy)
            tc_eta_ = np.array(entry.tc_eta)
            tc_mipPt_ = np.array(entry.tc_mipPt)
                    
            # loop over the gen particle and apply basic selection at MC-truth level
            for i_gen in range(len(gen_pt_)) :

                pt_cand = -1.
                eta_cand = -100.
                phi_cand = -100.
                hasMatched = False
                
                for i_tc in range(len(tc_energy_))  :
                    if ( (tc_eta_[i_tc])*(gen_eta_[i_gen])>0 
                         and(gen_pt_[i_gen])>self.cfg.minPt_gen 
                         and abs(gen_eta_[i_gen]) > self.cfg.minEta_gen 
                         and tc_mipPt_[i_tc] > self.cfg.tc_threshold 
                         and abs(gen_eta_[i_gen]) < self.cfg.maxEta_gen
                         and gen_status_[i_gen] == self.cfg.particle_status
                         and gen_id_[i_gen] == self.cfg.particle_type) :
                        
                        hasMatched = True
                        pt_cand += tc_energy_[i_tc] * f.sinTheta(tc_eta_[i_tc])

                # fill the response histograms                    
                if hasMatched :
                    #print  gen_pt_[i_gen],  gen_eta_[i_gen], pt_cand
                    h_resoPt.Fill( pt_cand / gen_pt_[i_gen] )
                    h_AllTCPtResponseVsTruePt.Fill( gen_pt_[i_gen], pt_cand / gen_pt_[i_gen] )
                    h_AllTCPtResponseVsTrueEta.Fill( abs(gen_eta_[i_gen]), pt_cand / gen_pt_[i_gen] )
                    h_resoPtVsThreshold.SetBinContent( h_resoPtVsThreshold.FindBin(self.cfg.tc_threshold), h_resoPt.GetMean() );    
                    h_resoPtVsThreshold.SetBinError( h_resoPtVsThreshold.FindBin(self.cfg.tc_threshold), h_resoPt.GetMeanError() );    

        # write histograms into output file
        h_resoPt.Write()
        h_resoPtVsThreshold.Write()
        h_AllTCPtResponseVsTruePt.Write()
        h_AllTCPtResponseVsTrueEta.Write()



    def plotResponseCheck(self) :
        # definition of the output file and histogram to store in it
        output = ROOT.TFile( self.outputFile,"UPDATE" )

        h_totTC = ROOT.TH1D( "totTC", "Total Pt in all trigger-cells above threshold ", 100, 0, 200 )
        h_totC2d = ROOT.TH1D( "totC2d", "Total Pt in all 2D-Clusters ", 100, 0, 200 )
        h_totC3d = ROOT.TH1D( "totC3d", "Total Pt in all 3D-Clusters ", 100, 0, 200 )

        h_ratio_c2d_tc  = ROOT.TH1D( "h_ratio_c2d_tc", "Ration total energy in C2d with respect to trigger-cells", 100, 0, 3 )
        h_ratio_c3d_tc  = ROOT.TH1D( "h_ratio_c3d_tc", "Ration total energy in C3d with respect to trigger-cells", 100, 0, 3 )
        h_ratio_c3d_c2d = ROOT.TH1D( "h_ratio_c3d_c2d", "Ration total energy in C3d with respect to C2d", 100, 0, 3 )


        # loop over the ttree ntries
        for ientry, entry in enumerate(self.chain):
            if ientry%1000==0 : 
                print ientry
            if ientry == self.cfg.maxNEvts :
                break

            gen_pt_ = np.array(entry.gen_pt) # for vectors
            gen_eta_ = np.array(entry.gen_eta) # for vectors
            tc_mipPt_ = np.array(entry.tc_mipPt)
            tc_energy_ = np.array(entry.tc_energy)
            tc_eta_ = np.array(entry.tc_eta)
            c2d_pt_ = np.array(entry.cl_pt)   
            c2d_eta_ = np.array(entry.cl_eta)            
            c3d_pt_ = np.array(entry.cl3d_pt)
            c3d_eta_ = np.array(entry.cl3d_eta)
            
            tc_ptSum=0
            c2d_ptSum=0
            c3d_ptSum=0
            
            c2d_over_tc=0
            c3d_over_tc=0
            c3d_over_c2d=0

            for i_tc in range(len(tc_energy_)) :
                if tc_mipPt_[i_tc] > self.cfg.tc_threshold and len(tc_energy_)>0 :
                    tc_ptSum += tc_energy_[i_tc] * f.sinTheta(tc_eta_[i_tc]) 
            for i_c2d in range(len(c2d_pt_)) :
                if len(c2d_pt_) > 0 :
                    c2d_ptSum += c2d_pt_[i_c2d] 
            for i_c3d in range(len(c3d_pt_)) :
                if len(c3d_pt_) > 0 :
                    c3d_ptSum += c3d_pt_[i_c3d]                 

            h_totTC.Fill( tc_ptSum )
            h_totC2d.Fill( c2d_ptSum )
            h_totC3d.Fill( c3d_ptSum )
            
            if tc_ptSum > 0 :
                h_ratio_c2d_tc.Fill( c2d_ptSum/tc_ptSum )
                if c2d_ptSum/tc_ptSum >1 :
                    print tc_ptSum, c2d_ptSum
            if tc_ptSum > 0 :
                h_ratio_c3d_tc.Fill( c3d_ptSum/tc_ptSum )
            if c2d_ptSum > 0 :
                h_ratio_c3d_c2d.Fill( c3d_ptSum/c2d_ptSum )
        # write histograms into output file
        
        h_totTC.Write()        
        h_totC2d.Write()
        h_totC3d.Write()
        h_ratio_c2d_tc.Write()
        h_ratio_c3d_tc.Write()
        h_ratio_c3d_c2d.Write()

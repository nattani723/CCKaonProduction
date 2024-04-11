#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

//const int dEdx_bins = 27;
//const double dEdx_binning[] = { 0.000, 0.500, 1.000, 1.500, 2.000, 2.500, 3.000, 3.500, 4.000, 4.500, 5.000, 5.500, 6.000, 6.500, 7.000, 7.500, 8.000, 9.000, 10.000, 12.000, 15.000, 20.000, 25.000, 30.000, 35.000, 40.000, 45.000, 50.000 };

const int dEdx_bins = 25;
const double dEdx_binning[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25};

const int ResRange_bins = 20;
const double ResRange_binning[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20};
const int Pitch_bins = 9;
const double Pitch_binning[] = {0.0,0.5,1.0,1.5,2.0,2.5,5.0,10.0,20.0,50.0};

const std::vector<int> pdg_v = {3222,3112,321,2212,13,211};
const std::vector<std::string> particle_v = {"SigmaP","SigmaM","Kaon","Proton","Muon","Pion"};

void MakeReferenceData(){

  // Setup the histograms
  std::vector<TH3D*> h_dEdx_ResidualRange_Pitch_Plane0_v; 
  std::vector<TH3D*> h_dEdx_ResidualRange_Pitch_Plane1_v; 
  std::vector<TH3D*> h_dEdx_ResidualRange_Pitch_Plane2_v; 

  for(size_t i_pdg=0;i_pdg<pdg_v.size();i_pdg++){
    h_dEdx_ResidualRange_Pitch_Plane0_v.push_back(new TH3D(("h_dEdx_ResidualRange_Pitch_Plane0_"+particle_v.at(i_pdg)).c_str(),(particle_v.at(i_pdg)+";Residual Range (cm);dE/dx (MeV/cm);Pitch (cm)").c_str(),ResRange_bins,ResRange_binning,dEdx_bins,dEdx_binning,Pitch_bins,Pitch_binning));
    h_dEdx_ResidualRange_Pitch_Plane1_v.push_back(new TH3D(("h_dEdx_ResidualRange_Pitch_Plane1_"+particle_v.at(i_pdg)).c_str(),(particle_v.at(i_pdg)+";Residual Range (cm);dE/dx (MeV/cm);Pitch (cm)").c_str(),ResRange_bins,ResRange_binning,dEdx_bins,dEdx_binning,Pitch_bins,Pitch_binning));
    h_dEdx_ResidualRange_Pitch_Plane2_v.push_back(new TH3D(("h_dEdx_ResidualRange_Pitch_Plane2_"+particle_v.at(i_pdg)).c_str(),(particle_v.at(i_pdg)+";Residual Range (cm);dE/dx (MeV/cm);Pitch (cm)").c_str(),ResRange_bins,ResRange_binning,dEdx_bins,dEdx_binning,Pitch_bins,Pitch_binning));
  } 

  // Load the trees containing the data
  TFile* f_in = TFile::Open("/exp/uboone/app/users/cthorpe/HyperonCode_v08_00_00_78/srcs/ubana/ubana/HyperonProduction/Data/dEdx_Trees.root");
  TTree* t_in = static_cast<TTree*>(f_in->Get("ana/OutputTree"));
  vector<int>     *TrackTruePDG;
  vector<float>   *TrackLength;
  vector<float>   *TrackAngleTrueReco;
  vector<vector<float> > *ResidualRange_Plane0;
  vector<vector<float> > *dEdx_Plane0;
  vector<vector<float> > *Pitch_Plane0;
  vector<vector<float> > *ResidualRange_Plane1;
  vector<vector<float> > *dEdx_Plane1;
  vector<vector<float> > *Pitch_Plane1;
  vector<vector<float> > *ResidualRange_Plane2;
  vector<vector<float> > *dEdx_Plane2;
  vector<vector<float> > *Pitch_Plane2;
  t_in->SetBranchAddress("TrackTruePDG", &TrackTruePDG);
  t_in->SetBranchAddress("TrackLength", &TrackLength);
  t_in->SetBranchAddress("TrackAngleTrueReco", &TrackAngleTrueReco);
  t_in->SetBranchAddress("ResidualRange_Plane0", &ResidualRange_Plane0);
  t_in->SetBranchAddress("dEdx_Plane0", &dEdx_Plane0);
  t_in->SetBranchAddress("Pitch_Plane0", &Pitch_Plane0);
  t_in->SetBranchAddress("ResidualRange_Plane1", &ResidualRange_Plane1);
  t_in->SetBranchAddress("dEdx_Plane1", &dEdx_Plane1);
  t_in->SetBranchAddress("Pitch_Plane1", &Pitch_Plane1);
  t_in->SetBranchAddress("ResidualRange_Plane2", &ResidualRange_Plane2);
  t_in->SetBranchAddress("dEdx_Plane2", &dEdx_Plane2);
  t_in->SetBranchAddress("Pitch_Plane2", &Pitch_Plane2);

  Long64_t nentries = t_in->GetEntriesFast();

  // Fill the histograms
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
   t_in->GetEntry(jentry);
    for(size_t i_tr=0;i_tr<TrackTruePDG->size();i_tr++){
      if(TrackLength->at(i_tr) < 2) continue;
      int pdg_index = -1;
      for(size_t i_pdg=0;i_pdg<pdg_v.size();i_pdg++){    
        if(pdg_v.at(i_pdg) == TrackTruePDG->at(i_tr)){
          pdg_index = i_pdg;
        }
      }
      if(pdg_index == -1) continue;
      for(size_t i_p=0;i_p<ResidualRange_Plane0->at(i_tr).size();i_p++)
        h_dEdx_ResidualRange_Pitch_Plane0_v.at(pdg_index)->Fill(ResidualRange_Plane0->at(i_tr).at(i_p),dEdx_Plane0->at(i_tr).at(i_p),Pitch_Plane0->at(i_tr).at(i_p));
      for(size_t i_p=0;i_p<ResidualRange_Plane1->at(i_tr).size();i_p++)
        h_dEdx_ResidualRange_Pitch_Plane1_v.at(pdg_index)->Fill(ResidualRange_Plane1->at(i_tr).at(i_p),dEdx_Plane1->at(i_tr).at(i_p),Pitch_Plane1->at(i_tr).at(i_p));
      for(size_t i_p=0;i_p<ResidualRange_Plane2->at(i_tr).size();i_p++)
        h_dEdx_ResidualRange_Pitch_Plane2_v.at(pdg_index)->Fill(ResidualRange_Plane2->at(i_tr).at(i_p),dEdx_Plane2->at(i_tr).at(i_p),Pitch_Plane2->at(i_tr).at(i_p));
    }
  }

  f_in->Close();

  // Normalise all of the histograms to 1 to create likelihood profile in 3D - want each slice in rr/pitch space to have pm of 1
  for(size_t i_pdg=0;i_pdg<pdg_v.size();i_pdg++){
    TH3D* h_Plane0 = h_dEdx_ResidualRange_Pitch_Plane0_v.at(i_pdg);
    TH3D* h_Plane1 = h_dEdx_ResidualRange_Pitch_Plane1_v.at(i_pdg);
    TH3D* h_Plane2 = h_dEdx_ResidualRange_Pitch_Plane2_v.at(i_pdg);
    for(int i_rr=1;i_rr<h_Plane0->GetNbinsX()+1;i_rr++){
      for(int i_p=1;i_p<h_Plane0->GetNbinsZ()+1;i_p++){

        double mass_Plane0 = 0.0;
        double mass_Plane1 = 0.0;
        double mass_Plane2 = 0.0;
        for(int i_dd=0;i_dd<h_Plane0->GetNbinsY()+1;i_dd++){
          mass_Plane0 += h_Plane0->GetBinContent(i_rr,i_dd,i_p);
          mass_Plane1 += h_Plane1->GetBinContent(i_rr,i_dd,i_p);
          mass_Plane2 += h_Plane2->GetBinContent(i_rr,i_dd,i_p);
        }
        
        for(int i_dd=0;i_dd<h_Plane0->GetNbinsY()+1;i_dd++){
          h_Plane0->SetBinContent(i_rr,i_dd,i_p,h_Plane0->GetBinContent(i_rr,i_dd,i_p)/mass_Plane0);
          h_Plane1->SetBinContent(i_rr,i_dd,i_p,h_Plane1->GetBinContent(i_rr,i_dd,i_p)/mass_Plane1);
          h_Plane2->SetBinContent(i_rr,i_dd,i_p,h_Plane2->GetBinContent(i_rr,i_dd,i_p)/mass_Plane2);
        }

      }
    }
  }


  // Draw plots of each slice in pitch space for each particle
  TCanvas* c = new TCanvas("c","c");
  gSystem->Exec("mkdir -p Plots/");

  for(size_t i_pdg=0;i_pdg<pdg_v.size();i_pdg++){

    TH3D h_3d_Plane0 = *h_dEdx_ResidualRange_Pitch_Plane0_v.at(i_pdg); 
    TH3D h_3d_Plane1 = *h_dEdx_ResidualRange_Pitch_Plane1_v.at(i_pdg); 
    TH3D h_3d_Plane2 = *h_dEdx_ResidualRange_Pitch_Plane2_v.at(i_pdg); 
    std::string particle = particle_v.at(i_pdg); 

    for(int i_z=1;i_z<h_3d_Plane0.GetNbinsZ()+1;i_z++){

      h_3d_Plane0.GetZaxis()->SetRange(i_z,i_z);
      h_3d_Plane1.GetZaxis()->SetRange(i_z,i_z);
      h_3d_Plane2.GetZaxis()->SetRange(i_z,i_z);

      TH2D* h_proj_Plane0 = static_cast<TH2D*>(h_3d_Plane0.Project3D("yx"));
      TH2D* h_proj_Plane1 = static_cast<TH2D*>(h_3d_Plane1.Project3D("yx"));
      TH2D* h_proj_Plane2 = static_cast<TH2D*>(h_3d_Plane2.Project3D("yx"));

      std::string title = particle + " , Pitch (" + std::to_string(h_3d_Plane0.GetZaxis()->GetBinLowEdge(i_z)) + "," + std::to_string(h_3d_Plane0.GetZaxis()->GetBinLowEdge(i_z+1)) + ")";
      std::string label = particle + "_" + std::to_string(i_z);

      h_proj_Plane0->Draw("colz");
      h_proj_Plane0->SetStats(0);
      h_proj_Plane0->SetTitle(title.c_str());
      c->Print(("Plots/"+label+"_Plane0.png").c_str());
      c->Clear();

      h_proj_Plane1->Draw("colz");
      h_proj_Plane1->SetStats(0);
      h_proj_Plane1->SetTitle(title.c_str());
      c->Print(("Plots/"+label+"_Plane1.png").c_str());
      c->Clear();

      h_proj_Plane2->Draw("colz");
      h_proj_Plane2->SetStats(0);
      h_proj_Plane2->SetTitle(title.c_str());
      c->Print(("Plots/"+label+"_Plane2.png").c_str());
      c->Clear();

    }
  }    


  // Write the histograms to file to serve as reference tables
  TFile* f_out = TFile::Open("dEdx_Reference.root","RECREATE");  
  for(size_t i_pdg=0;i_pdg<pdg_v.size();i_pdg++){
    const string pdg = std::to_string(pdg_v.at(i_pdg)); 
    h_dEdx_ResidualRange_Pitch_Plane0_v.at(i_pdg)->Write((pdg + "_Plane0").c_str());
    h_dEdx_ResidualRange_Pitch_Plane1_v.at(i_pdg)->Write((pdg + "_Plane1").c_str());
    h_dEdx_ResidualRange_Pitch_Plane2_v.at(i_pdg)->Write((pdg + "_Plane2").c_str());
  }
  f_out->Close();

}

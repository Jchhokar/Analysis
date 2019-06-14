
using namespace std;

int Analysis(){

for(int files = 0; files < 2; files++){
    TFile input_file(Form("../%i/mcGeant.root", files), "READ");
    TH2D ThetaVsTheta("ThetaVsTheta", "ThetaVsTheta", 210, 0., 210., 210, 0., 210.);
    TH2D SumVsDiff("SumVsDiff", "SumVsDiff", 250, 0, 250, 180, 0, 180);
	TH1F e1("e1" , "e1", 600, 0, 600);
	TH1F e2("e2" , "e2", 600, 0, 600);
	TH1F e3("e3" , "e3", 600, 0, 600);

    TTree* tree = (TTree*)input_file.Get("T");


    TFile output(Form("jpetmc_output%i.root", files), "recreate");

    TTree outputTree("WeirdsHit", "WeirdsHit");

    double speed_of_light = 0.0299792458; // cm ps^-1
    double electron_mass = 510.998928; // electron mass in keV

    JPetGeantEventPack* eventPack = new JPetGeantEventPack();
    JPetGeantEventInformation* evtInfo = new JPetGeantEventInformation();
    JPetGeantScinHits *MCHit = new JPetGeantScinHits();

    tree->SetBranchAddress("eventPack", &eventPack);
    int nevent = tree->GetEntries();                 ////////////////////
    cout << "events to process " << nevent << endl;

    outputTree.Branch("eventPack", &eventPack);

    int toprint = 0;
    double minEdep = 0.0;
    int counter = 0;

    for(int i = 0; i < nevent; i++){

        if (100.0*i/nevent >= (toprint + 1.)) {
          toprint = toprint + 1.;
             cout << "Running......" << 100.0 * i/nevent << "%........." << endl;
         }


        tree->GetEvent(i);

        evtInfo = eventPack->GetEventInformation();
        int NumberOfHits = eventPack->GetNumberOfHits();

        vector<JPetGeantScinHits*> PrimaryHits;

       if(NumberOfHits > 2){        //selection of 3 hits(atlaeast 3 hits
            // cout << eventPack->GetEventNumber() << endl;
            for(int j = 0; j < NumberOfHits; j++){

                if(eventPack->GetHit(j)->GetGenGammaMultiplicity() == 3){
					//if(abs( (1000*evtInfo->GetMomentumGamma(eventPack->GetHit(j)->GetGenGammaIndex())- eventPack->GetHit(j)->GetMomentumIn()).Mag() ) < 0.00001){
                     //cout << eventPack->GetHit(j)->GetGenGammaMultiplicity() << " " << eventPack->GetHit(j)->GetGenGammaIndex() << " " << eventPack->GetHit(j)->GetNumOfInteractions() << endl;
                    	PrimaryHits.push_back(eventPack->GetHit(j));
					//}

                }
            }
            // cout << endl;
        }

        if(PrimaryHits.size() == 3){
            // cout << evtInfo->GetMomentumGamma(0).Mag() << endl;
             PrimaryHits.at(0);
           // random_shuffle(PrimaryHits.begin(), PrimaryHits.end());

            double gentheta12 = TMath::RadToDeg()*evtInfo->GetMomentumGamma(1).Angle(evtInfo->GetMomentumGamma(2));
            double gentheta13 = TMath::RadToDeg()*evtInfo->GetMomentumGamma(1).Angle(evtInfo->GetMomentumGamma(3));
            double gentheta23 = TMath::RadToDeg()*evtInfo->GetMomentumGamma(2).Angle(evtInfo->GetMomentumGamma(3));

           sort(PrimaryHits.begin(), PrimaryHits.end(), [](JPetGeantScinHits* lhs, JPetGeantScinHits* rhs){return lhs->GetGenGammaIndex() < rhs->GetGenGammaIndex();} );

            // double ttheta12 = TMath::RadToDeg()*TMath::ACos((- pow(PrimaryHits.at(0)->GetMomentumIn().Mag(),2) - pow(PrimaryHits.at(1)->GetMomentumIn().Mag(),2) + pow(PrimaryHits.at(2)->GetMomentumIn().Mag(),2))/(2*PrimaryHits.at(0)->GetMomentumIn().Mag()*PrimaryHits.at(1)->GetMomentumIn().Mag()));
            // double ttheta13 = TMath::RadToDeg()*TMath::ACos(( pow(PrimaryHits.at(0)->GetMomentumIn().Mag(),2) - pow(PrimaryHits.at(1)->GetMomentumIn().Mag(),2) - pow(PrimaryHits.at(2)->GetMomentumIn().Mag(),2))/(2*PrimaryHits.at(1)->GetMomentumIn().Mag()*PrimaryHits.at(2)->GetMomentumIn().Mag()));
            // double ttheta23 = TMath::RadToDeg()*TMath::ACos((- pow(PrimaryHits.at(0)->GetMomentumIn().Mag(),2) + pow(PrimaryHits.at(1)->GetMomentumIn().Mag(),2) - pow(PrimaryHits.at(2)->GetMomentumIn().Mag(),2))/(2*PrimaryHits.at(2)->GetMomentumIn().Mag()*PrimaryHits.at(0)->GetMomentumIn().Mag()));


            double theta12 = TMath::RadToDeg()*PrimaryHits.at(0)->GetMomentumIn().Angle(PrimaryHits.at(1)->GetMomentumIn());
            double theta13 = TMath::RadToDeg()*PrimaryHits.at(0)->GetMomentumIn().Angle(PrimaryHits.at(2)->GetMomentumIn());
            double theta23 = TMath::RadToDeg()*PrimaryHits.at(1)->GetMomentumIn().Angle(PrimaryHits.at(2)->GetMomentumIn());

            
			
            vector<double> thetas{theta12, theta13, theta23};
			
			double E1 = PrimaryHits.at(0)->GetEneDepos();
			double E2 = PrimaryHits.at(1)->GetEneDepos();
			double E3 = PrimaryHits.at(2)->GetEneDepos();
			vector<double> energy{E1, E2, E3};

            vector<double> genthetas{gentheta12, gentheta13, gentheta23};

            //sort(genthetas.begin(), genthetas.end(), [](double lhs, double rhs){return lhs < rhs;} );
			//random_shuffle(genthetas.begin(),genthetas.end());
			//ThetaVsTheta.Fill(genthetas.at(0), genthetas.at(1));
            // vector<double> tthetas{ttheta12, ttheta12, ttheta23};
            // sort(tthetas.begin(), tthetas.end(), [](double lhs, double rhs){return lhs < rhs;} );
            // cout << tthetas.at(0) << " " << tthetas.at(1) << " " << tthetas.at(2) << endl;
            // cout << tthetas.at(0) + tthetas.at(1) + tthetas.at(2) << endl;

            // sort(genthetas.begin(), genthetas.end(), [](double lhs, double rhs){return lhs < rhs;} );
            // cout << genthetas.at(0) << " " << genthetas.at(1) << " " << genthetas.at(2) << endl;
            // cout << genthetas.at(0) + genthetas.at(1) + genthetas.at(2) << endl;

          // if(energy[0]>0.0 and energy[1]>0.0 and energy[2]>0.0) {   

                //sort(thetas.begin(), thetas.end(), [](double lhs, double rhs){return lhs < rhs;} );
	     	// random_shuffle(thetas.begin(),thetas.end());
                 ThetaVsTheta.Fill(thetas.at(0), thetas.at(1));

				e1.Fill(E1);
				e2.Fill(E2);
				e3.Fill(E3);
			//}
			
          //sort(genthetas.begin(), genthetas.end(), [](double lhs, double rhs){return lhs < rhs;} );
            

            TVector3 sum(0.0,0.0,0.0);


            // cout << sum.Mag() << endl;
            SumVsDiff.Fill(thetas.at(0) + thetas.at(1), thetas.at(1) - thetas.at(0));
             if(thetas.at(0) + thetas.at(1) < 162){
                //cout << (evtInfo->GetMomentumGamma(1) + evtInfo->GetMomentumGamma(2) + evtInfo->GetMomentumGamma(3)).Mag() << endl;
//				cout << thetas.at(0) + thetas.at(1) << endl;
                    counter++;

                for(auto hit :PrimaryHits){
                    //cout << hit->GetGenGammaMultiplicity() << " " << hit->GetGenGammaIndex() << " " << hit->GetNumOfInteractions() <<  endl;
                    //cout << 1000*evtInfo->GetMomentumGamma(hit->GetGenGammaIndex()).Px() << " " << 1000*evtInfo->GetMomentumGamma(hit->GetGenGammaIndex()).Py() << " " << 1000*evtInfo->GetMomentumGamma(hit->GetGenGammaIndex()).Pz() << endl;
                    //cout << hit->GetMomentumIn().Px() << " " << hit->GetMomentumIn().Py() << " " << hit->GetMomentumIn().Pz() << endl;
                    sum = sum + hit->GetMomentumIn();
                }
                //cout << sum.Mag() << endl;
                //cout << endl;
               // cout << thetas.at(0) << " " << thetas.at(1) << " " << thetas.at(2) << endl;
               // cout << thetas.at(0) + thetas.at(1) + thetas.at(2) << endl;

                //cout << endl << endl;
                outputTree.Fill();

            }
      }
      PrimaryHits.clear();
    }


    input_file.Close();
	e1.Write();
	e2.Write();
	e3.Write();
    outputTree.Write();
    ThetaVsTheta.Write();
    SumVsDiff.Write();
    output.Close();
    cout << counter << endl;
 }
    return 0;
}

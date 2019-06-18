#include "TMath.h"
using namespace std;
double  NormalizeTime( JPetGeantScinHits* Hit1 );
double CalcDistanceOfSurfaceAndZero( JPetGeantScinHits* Hit1, JPetGeantScinHits* Hit2, JPetGeantScinHits* Hit3 );

int Analysis(){
	for(int files = 0; files < 10; files++){
		TFile input_file(Form("../%i/mcGeant.root", files), "READ");
    	//TFile input_file("../mcGeant.root", "READ");
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
        double energy_threshold = 0.0;

   		JPetGeantEventPack* eventPack = new JPetGeantEventPack();
		JPetGeantEventInformation* evtInfo = new JPetGeantEventInformation();
		JPetGeantScinHits *MCHit = new JPetGeantScinHits();

    	tree->SetBranchAddress("eventPack", &eventPack);
    	int nevent = tree->GetEntries();                 ////////////////////
    	cout << "events to process " << nevent << endl;

    	outputTree.Branch("eventPack", &eventPack);

    	int toprint = 0;
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

       		if(NumberOfHits > 2){        //selection of 3 hits(atleast 3 hits
            	for(int j = 0; j < NumberOfHits; j++){
                	if(eventPack->GetHit(j)->GetGenGammaMultiplicity() == 3){
						/*if(abs( (1000*evtInfo->GetMomentumGamma(eventPack->GetHit(j)->GetGenGammaIndex())- eventPack->GetHit(j)-
						  >GetMomentumIn()).Mag() ) < 0.00001){
                     			cout << eventPack->GetHit(j)->GetGenGammaMultiplicity() << " " << eventPack->GetHit(j)->GetGenGammaIndex() << " " << 									eventPack->GetHit(j)->GetNumOfInteractions() << endl;
						}*/
                    	PrimaryHits.push_back(eventPack->GetHit(j));
					}
            	}
            // cout << endl;
       		}

			if(PrimaryHits.size() == 3){

				//random_shuffle(PrimaryHits.begin(), PrimaryHits.end());
            	sort(PrimaryHits.begin(), PrimaryHits.end(), [](JPetGeantScinHits* lhs, JPetGeantScinHits* rhs)
              	{return lhs->GetGenGammaIndex() < rhs->GetGenGammaIndex();} );

				double firstPhotonArr = NormalizeTime(PrimaryHits.at(0));      // Time is calculated using time of emmision (TOF correction)
				double secondPhotonArr= NormalizeTime(PrimaryHits.at(1));
				double thirdPhotonArr = NormalizeTime(PrimaryHits.at(2));

            	double theta12 = TMath::RadToDeg()*PrimaryHits.at(0)->GetMomentumIn().Angle(PrimaryHits.at(1)->GetMomentumIn());
            	double theta13 = TMath::RadToDeg()*PrimaryHits.at(0)->GetMomentumIn().Angle(PrimaryHits.at(2)->GetMomentumIn());
            	double theta23 = TMath::RadToDeg()*PrimaryHits.at(1)->GetMomentumIn().Angle(PrimaryHits.at(2)->GetMomentumIn());

				double E1 = PrimaryHits.at(0)->GetEneDepos();
				double E2 = PrimaryHits.at(1)->GetEneDepos();
				double E3 = PrimaryHits.at(2)->GetEneDepos();

            	vector<double> thetasOrd{theta12,theta13,theta23}, thetasNord{theta12,theta13,theta23};
				sort(thetasOrd.begin(), thetasOrd.end(), [](double lhs, double rhs){return lhs < rhs;} );
				vector<double> energy{E1, E2, E3};

				if(fabs(PrimaryHits.at(0)->GetHitPosition())<23. && fabs(PrimaryHits.at(1)->GetHitPosition())<23. && fabs(PrimaryHits.at(2)->GetHitPosition())<23.){

					if(energy[0]>energy_threshold and energy[1]>energy_threshold and energy[2]>energy_threshold){

						if(CalcDistanceOfSurfaceAndZero(PrimaryHits.at(0), PrimaryHits.at(1), PrimaryHits.at(2))< 5.){

							if(TMath::Abs(firstPhotonArr-thirdPhotonArr)<1.5){

								if(thetasOrd.at(0) + thetasOrd.at(1) >180) {

									SumVsDiff.Fill(thetasOrd.at(0) + thetasOrd.at(1), thetasOrd.at(1) - thetasOrd.at(0));
                					ThetaVsTheta.Fill(thetasNord.at(0), thetasNord.at(1));
								}
							}
						}
					}
				}

				e1.Fill(E1);
				e2.Fill(E2);
				e3.Fill(E3);

            	counter++;
            	outputTree.Fill();
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
double CalcDistanceOfSurfaceAndZero( JPetGeantScinHits* Hit1, JPetGeantScinHits* Hit2, JPetGeantScinHits* Hit3 )
{

	TVector3 vec1 = Hit2->GetHitPosition() - Hit1->GetHitPosition();
	TVector3 vec2 = Hit3->GetHitPosition() - Hit2->GetHitPosition();
	TVector3 crossProd  = vec1.Cross(vec2);
	double Dcoeef = crossProd.Dot(Hit2);
	double distanceFromZero = fabs(Dcoeef) / crossProd.Mag();
	return distanceFromZero;
}
//Hit normalize time
//
double  NormalizeTime( JPetGeantScinHits* Hit1 )
{
	    TVector3 vec1 = Hit1->GetHitPosition();
	    double Length0 = vec1.Mag();

	    return Hit1->GetTime()/1000.0 - (Length0)/29.979246;
}

/** Program based on Rootana to convert Midas into 
ROOT events. Inputs a .mid file and ouputs a .root file, 
containing a tree called "waveform_tree."
There are up to 8 trees to contain the waveforms 
of each the digitizer's channels. 

Originally created by Thomas Lindner for the Big Read Board's ADCs.
Adapted by Luan Koerich for the CAEN-DT5743 digitizer, February 2021.

*/
#include <stdio.h>
#include <iostream>
#include "TRootanaEventLoop.hxx"
#include "TH1F.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TTreeMaker.h"
#include <TTree.h>
#include <TFile.h>
#include "TDT743RawData.hxx"


#define TRUE 1
#define FALSE 0



#define nPoints_max 1
#define num_phidg_max 10000
#define num_dt5743_max 1024 // IMPOTANT: IF THIS IS EVER CHANGED, ALSO CHANGE THE HARDCODED VALUES FOR WAVEFORM BRANCH WIDTHS AS WELL (see: "v1730 data")

//Digitizer variables
double Start_time0[nPoints_max], Window_width0[nPoints_max], Start_time1[nPoints_max], Window_width1[nPoints_max], start_0, start_1;

// Digitizer channels
double dt5743_wave0[num_dt5743_max], dt5743_wave1[num_dt5743_max], dt5743_wave2[num_dt5743_max], dt5743_wave3[num_dt5743_max];
double dt5743_wave4[num_dt5743_max], dt5743_wave5[num_dt5743_max], dt5743_wave6[num_dt5743_max], dt5743_wave7[num_dt5743_max];

//PMT readout
int start_val_stat;
int window_width;
int trigger;
int counter;
int num_points;
int num_points_dig0;

#define timeStart 130 // defines start of PMT Pulse timing window, currently at the 130th sample of 200, with a window size of 70 samples.

// Offset for the ADC channel number
#define Ch_Offset 1

// 
class ScanToTreeConverter: public TRootanaEventLoop {

    int nnn;
    TH1F *SampleWave0;

    TH1F *StartVal0;

    private:
        
        // Number of channels to be saved
        int fNChan;

        //TFile *outputfile; //made by TRootAnaEventLoop with name of inputfile +.root
        TTree *tree;

        // Counters for TDC bank
        int ngoodTDCbanks;
        int nbadTDCbanks;

    public:

        ScanToTreeConverter() {
            
            // Disable ROOTANA graphical mode.
            UseBatchMode();
            nnn = 0;
            fNChan = 8; // < Saving waveforms from 0 to 7
        };

    virtual ~ScanToTreeConverter() {};

    void BeginRun(int transition,int run,int time){
    
        // Start conversion
        std::cout << "Custom: begin run " << run << std::endl;
    
        // Add Canvas
        TCanvas *c1 = new TCanvas("c1","Canvas Example",200,10,600,480);

        // A tree to contain the waveform information
        tree = new TTree("waveform_tree","Waveform Tree");
        tree->Branch("num_points",&num_points,"num_points/Int_t");

        // DT5743 data  V1730_wave0[nPoints_max][num_dt5743_max]
        for(int i =0; i < fNChan; i++){  // Save up to 24 waveforms
            char name[100], descr[100];
            sprintf(name,"dt5743_wave0%i",i);
            sprintf(descr,"dt5743_wave0%i[%i]/Double_t",i,num_dt5743_max);

            if(i==0) tree->Branch(name,&dt5743_wave0,descr); 
            if(i==1) tree->Branch(name,&dt5743_wave1,descr); 
            if(i==2) tree->Branch(name,&dt5743_wave2,descr); 
            if(i==3) tree->Branch(name,&dt5743_wave3,descr); 
            if(i==4) tree->Branch(name,&dt5743_wave4,descr); 
            if(i==5) tree->Branch(name,&dt5743_wave5,descr); 
            if(i==6) tree->Branch(name,&dt5743_wave6,descr); 
            if(i==7) tree->Branch(name,&dt5743_wave7,descr);
        }

    }

    void EndRun(int transition,int run,int time){
        // Print out summary.
        std::cout << "End of conversion." << run <<std::endl;
        std::cout << "Good TDC banks : " << ngoodTDCbanks << std::endl;;
        std::cout << "Bad TDC banks  : " << nbadTDCbanks << std::endl;;
    }

    bool ProcessMidasEvent(TDataContainer& dataContainer){


        //attempt at adding digitizer data
        //Set input parameters which correspond to ADC measurements
        start_val_stat = 90;
        window_width = 45;
        trigger = 40;
    
        // Get digitizer data

        //** CHANGE HERE TO DT5743 DATA
        TDT743RawData *digitizer_data = dataContainer.GetEventData<TDT743RawData>("43FS");
    
        if(digitizer_data){      
      
            if(num_points<nPoints_max){
        
                num_points++;
        
                std::vector<RawChannelMeasurement> measures = digitizer_data->GetMeasurements();
        
                for(int i = 0; i < measures.size(); i++){           
                    int chan = measures[i].GetChannel();
                    for(int ib = 0; ib < measures[i].GetNSamples(); ib++){
                        if(chan == 0) dt5743_wave0[ib] = measures[i].GetSample(ib);  
                        if(chan == 1) dt5743_wave1[ib] = measures[i].GetSample(ib);  
                        if(chan == 2) dt5743_wave2[ib] = measures[i].GetSample(ib);  
                        if(chan == 3) dt5743_wave3[ib] = measures[i].GetSample(ib);
                        if(chan == 4) dt5743_wave4[ib] = measures[i].GetSample(ib);
                        if(chan == 5) dt5743_wave5[ib] = measures[i].GetSample(ib);
                        if(chan == 6) dt5743_wave6[ib] = measures[i].GetSample(ib);
                        if(chan == 7) dt5743_wave7[ib] = measures[i].GetSample(ib);              
                }   
            } 
            
            } else {
                std::cout << "Too many points! " << num_points << std::endl;
            }
      
            tree->Fill();
            counter = 0;
            num_points = 0;
            num_points_dig0 = 0;
            
            return true;
        }

        return true;
    }


    // Describe some other command line argument
    void Usage(void){
        std::cout << "\t-nchan option: specify how many channels of digitizer to save " << std::endl;
    }

    // Define some other command line argument
    bool CheckOption(std::string option){
    const char* arg = option.c_str();
    
    if (strncmp(arg,"-nchan",2)==0){
        fNChan = atoi(arg+6);
        std::cout << "Number of channels to save: " << fNChan << std::endl;
        return true;
    }

    return false;
    }
}; 

// Main function to initiliaze converter
int main(int argc, char *argv[]){

    // Test arguments
    if (argc != 2) {
        std::cerr<<"Usage: midas2root_mppc.exe run.mid.lz4"<<std::endl;
        exit(0);
    }

    ScanToTreeConverter::CreateSingleton<ScanToTreeConverter>();
    return ScanToTreeConverter::Get().ExecuteLoop(argc, argv);
}


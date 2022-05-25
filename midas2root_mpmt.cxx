// Simple program to print out the values from a magnetic field scan.

#include <stdio.h>
#include <iostream>

#include "TRootanaEventLoop.hxx"
#include "TH1F.h"
#include "TH1D.h"
#include "TV1190Data.hxx"
#include "TV792Data.hxx"
#include "TV1730RawData.hxx"
#include "TBRBRawData.hxx"
#include "TCanvas.h"
#include "TTreeMaker.h"

#include <TTree.h>
#include <TFile.h>

#define TRUE 1
#define FALSE 0



#define nPoints_max 1
#define num_phidg_max 10000
#define num_v1730_max 1024 // IMPOTANT: IF THIS IS EVER CHANGED, ALSO CHANGE THE HARDCODED VALUES FOR WAVEFORM BRANCH WIDTHS AS WELL (see: "v1730 data")

  int counter_gant;                                    //n-th measurement
  int subevent;
  double x0_pos =  0.46, y0_pos = 0.38, z0_pos, tilt0_pos, rot0_pos;  //positions of gantry0
  double x1_pos =  0.46, y1_pos = 0.38, z1_pos, tilt1_pos, rot1_pos;  //positions of gantry1
  double cyctime;

  //TODO: add other positions (Phidget, laser, ...)

  int num_phidg0_points;
  double x0_field[num_phidg_max], y0_field[num_phidg_max], z0_field[num_phidg_max], tot0_field[num_phidg_max];     //B-field from phidget0
  int num_phidg1_points;
  double x1_field[num_phidg_max], y1_field[num_phidg_max], z1_field[num_phidg_max], tot1_field[num_phidg_max];     //B-field from phidget1
  double tilt_phid0[num_phidg_max], tilt_phid1[num_phidg_max];                       //tilt position from phidget
  double acc_x0[num_phidg_max], acc_y0[num_phidg_max], acc_z0[num_phidg_max];
  double acc_x1[num_phidg_max], acc_y1[num_phidg_max], acc_z1[num_phidg_max];

  int num_phidg3_points;
  double x3_field[num_phidg_max], y3_field[num_phidg_max], z3_field[num_phidg_max], tot3_field[num_phidg_max];     //B-field from phidget3
  int num_phidg4_points;
  double x4_field[num_phidg_max], y4_field[num_phidg_max], z4_field[num_phidg_max], tot4_field[num_phidg_max];     //B-field from phidget4
  double tilt_phid3[num_phidg_max], tilt_phid4[num_phidg_max];                       //tilt position from phidget
  double acc_x3[num_phidg_max], acc_y3[num_phidg_max], acc_z3[num_phidg_max];
  double acc_x4[num_phidg_max], acc_y4[num_phidg_max], acc_z4[num_phidg_max];


  //current and voltage of each of the 6 Helmholtz coils
  int counter_mag;
  double curr_coil1, curr_coil2, curr_coil3, curr_coil4, curr_coil5, curr_coil6;
  double volt_coil1, volt_coil2, volt_coil3, volt_coil4, volt_coil5, volt_coil6;


  //PMT currents and voltages
  double curr_hpd_enable, curr_hpd_hv_control, curr_hpd_lv_control;
  double curr_receiver0, curr_receiver1, curr_monitor0, curr_monitor1;

  double volt_hpd_enable, volt_hpd_hv_control, volt_hpd_lv_control;
  double volt_receiver0, volt_receiver1, volt_monitor0, volt_monitor1;

  //Digitizer variables
  double Start_time0[nPoints_max], Window_width0[nPoints_max], Start_time1[nPoints_max], Window_width1[nPoints_max], start_0, start_1;
  
  // V1730 waveforms
  // 9.Nov.2017 Let's try adding using channels 2,3,4, 5
  double V1730_wave0[nPoints_max][num_v1730_max],  V1730_wave1[nPoints_max][num_v1730_max],  V1730_wave2[nPoints_max][num_v1730_max], V1730_wave3[nPoints_max][num_v1730_max], V1730_wave4[nPoints_max][num_v1730_max], V1730_wave5[nPoints_max][num_v1730_max], V1730_wave6[nPoints_max][num_v1730_max], V1730_wave7[nPoints_max][num_v1730_max];
  double V1730_wave8[nPoints_max][num_v1730_max],  V1730_wave9[nPoints_max][num_v1730_max],  V1730_wave10[nPoints_max][num_v1730_max], V1730_wave11[nPoints_max][num_v1730_max], V1730_wave12[nPoints_max][num_v1730_max], V1730_wave13[nPoints_max][num_v1730_max], V1730_wave14[nPoints_max][num_v1730_max], V1730_wave15[nPoints_max][num_v1730_max];
  double V1730_wave16[nPoints_max][num_v1730_max],  V1730_wave17[nPoints_max][num_v1730_max],  V1730_wave18[nPoints_max][num_v1730_max], V1730_wave19[nPoints_max][num_v1730_max], V1730_wave20[nPoints_max][num_v1730_max], V1730_wave21[nPoints_max][num_v1730_max], V1730_wave22[nPoints_max][num_v1730_max], V1730_wave23[nPoints_max][num_v1730_max];


  //PMT readout
  int start_val_stat;
  int window_width;
  int trigger;
  int counter;
  int num_points;
  int num_points_dig0; 
  int num_points_dig1;  




#define timeStart 130 // defines start of PMT Pulse timing window, currently at the 130th sample of 200, with a window size of 70 samples.

// Offset for the ADC channel number
#define Ch_Offset 1

// Flag to indicate the gantry was not moving and to record ADC and Phidget values. 
int gbl_accept_banks = TRUE; //set FALSE

class ScanToTreeConverter: public TRootanaEventLoop {

  int nnn;
  TH1F *SampleWave0;
  TH1F *SampleWave1;

  //9.Nov.2017 added other channels
  TH1F *SampleWave2;
  TH1F*SampleWave3;
  TH1F *SampleWave4;
  TH1F *SampleWave5;

  TH1F *StartVal0;

  private:

  int fNChan; // number of digitizer channels to save.

  //TFile *outputfile; //made by TRootAnaEventLoop with name of inputfile +.root
  TTree *tree;

  // Counters for TDC bank
  int ngoodTDCbanks;
  int nbadTDCbanks;

  public:

  ScanToTreeConverter() {
    UseBatchMode(); //necessary to switch off graphics, used in for example AnaDisplay
    nnn = 0;
    fNChan = 20; // < Saving waveforms from 0 to 19
  };

  virtual ~ScanToTreeConverter() {};

  void BeginRun(int transition,int run,int time){
    std::cout << "Custom: begin run " << run << std::endl;
    //setup ROOT branches

    //Add Canvas
    TCanvas *c1 = new TCanvas("c1","Canvas Example",200,10,600,480);

    tree = new TTree("scan_tree","Scan Tree");
    tree->Branch("num_points",&num_points,"num_points/Int_t");
    tree->Branch("num_points_dig0",&num_points_dig0,"num_points_dig0/Int_t");
    tree->Branch("num_points_dig1",&num_points_dig1,"num_points_dig1/Int_t");    

    tree->Branch("time",&cyctime,"time/Double_t");        //arbitrary offset in time, in ms, either UInt or larger, so just be safe
    tree->Branch("gantry_event",&counter_gant,"gantry_event/Int_t"); //an event is a measurement at a point
    tree->Branch("gantry_subevent",&subevent,"gantry_subevent/Int_t");

    tree->Branch("gantry0_x",&x0_pos,"gantry0_x/Double_t");
    tree->Branch("gantry0_y",&y0_pos,"gantry0_y/Double_t");
    tree->Branch("gantry0_z",&z0_pos,"gantry0_z/Double_t");
    tree->Branch("gantry0_tilt",&tilt0_pos,"gantry0_tilt/Double_t");
    tree->Branch("gantry0_rot",&rot0_pos,"gantry0_rot/Double_t");
    //tree->Branch("phidg0_tilt",&tilt_phid0,"phidg0_tilt/Double_t");

    tree->Branch("Start_time0",&Start_time0,"Start_time0[num_points_dig0]/Double_t");
    tree->Branch("Start_time1",&Start_time1,"Start_time1[num_points_dig1]/Double_t");
    tree->Branch("Window_width0",&Window_width0,"Window_width0[num_points_dig0]/Double_t");
    tree->Branch("Window_width1",&Window_width1,"Window_width1[num_points_dig1]/Double_t");

    // v1730 data  V1730_wave0[nPoints_max][num_v1730_max]
    for(int i =0; i < fNChan; i++){  // Save up to 24 waveforms
      char name[100], descr[100];
      sprintf(name,"V1730_wave%i",i);
      sprintf(descr,"V1730_wave%i[num_points][%i]/Double_t",i,num_v1730_max);

      if(i==0) tree->Branch(name,&V1730_wave0,descr); 
      if(i==1) tree->Branch(name,&V1730_wave1,descr); 
      if(i==2) tree->Branch(name,&V1730_wave2,descr); 
      if(i==3) tree->Branch(name,&V1730_wave3,descr); 
      if(i==4) tree->Branch(name,&V1730_wave4,descr); 
      if(i==5) tree->Branch(name,&V1730_wave5,descr); 
      if(i==6) tree->Branch(name,&V1730_wave6,descr); 
      if(i==7) tree->Branch(name,&V1730_wave7,descr); 
      if(i==8) tree->Branch(name,&V1730_wave8,descr); 
      if(i==9) tree->Branch(name,&V1730_wave9,descr); 
      if(i==10) tree->Branch(name,&V1730_wave10,descr); 
      if(i==11) tree->Branch(name,&V1730_wave11,descr); 
      if(i==12) tree->Branch(name,&V1730_wave12,descr); 
      if(i==13) tree->Branch(name,&V1730_wave13,descr); 
      if(i==14) tree->Branch(name,&V1730_wave14,descr); 
      if(i==15) tree->Branch(name,&V1730_wave15,descr); 
      if(i==16) tree->Branch(name,&V1730_wave16,descr); 
      if(i==17) tree->Branch(name,&V1730_wave17,descr); 
      if(i==18) tree->Branch(name,&V1730_wave18,descr); 
      if(i==19) tree->Branch(name,&V1730_wave19,descr); 
      if(i==20) tree->Branch(name,&V1730_wave20,descr); 
      if(i==21) tree->Branch(name,&V1730_wave21,descr); 
      if(i==22) tree->Branch(name,&V1730_wave22,descr); 
      if(i==23) tree->Branch(name,&V1730_wave23,descr); 

    }


    //tree->Branch("V1730_wave0",&V1730_wave0,"V1730_wave0[num_points]/Double_t"); //think of eqn*
    //tree->Branch("V1730_wave1",&V1730_wave1,"V1730_wave1[num_points]/Double_t"); //think of eqn*
    //tree->Branch("V1730_wave2",&V1730_wave2,"V1730_wave2[num_points]/Double_t"); //think of eqn*
    //tree->Branch("ADC0_voltage",ADC0_voltage,"ADC0_voltage[num_points]/Int_t");

    tree->Branch("gantry1_x",&x1_pos,"gantry1_x/Double_t");
    tree->Branch("gantry1_y",&y1_pos,"gantry1_y/Double_t");
    tree->Branch("gantry1_z",&z1_pos,"gantry1_z/Double_t");
    tree->Branch("gantry1_tilt",&tilt1_pos,"gantry1_tilt/Double_t");
    tree->Branch("gantry1_rot",&rot1_pos,"gantry1_rot/Double_t");
    //tree->Branch("phidg1_tilt",&tilt_phid1,"phidg1_tilt/Double_t");

    //field-related phidget measurements
    tree->Branch("num_phidg0_points",&num_phidg0_points,"num_phidg0_points/Int_t");
    tree->Branch("phidg0_Ax",acc_x0,"phidg0_Ax[num_phidg0_points]/Double_t");
    tree->Branch("phidg0_Ay",acc_y0,"phidg0_Ay[num_phidg0_points]/Double_t");
    tree->Branch("phidg0_Az",acc_z0,"phidg0_Az[num_phidg0_points]/Double_t");
    tree->Branch("phidg0_Bx",x0_field,"phidg0_Bx[num_phidg0_points]/Double_t");
    tree->Branch("phidg0_By",y0_field,"phidg0_By[num_phidg0_points]/Double_t");
    tree->Branch("phidg0_Bz",z0_field,"phidg0_Bz[num_phidg0_points]/Double_t");
    tree->Branch("phidg0_Btot",tot0_field,"phidg0_Btot[num_phidg0_points]/Double_t");
    tree->Branch("phidg0_tilt",tilt_phid0,"phidg0_tilt[num_phidg0_points]/Double_t");

    tree->Branch("num_phidg1_points",&num_phidg1_points,"num_phidg1_points/Int_t");
    tree->Branch("phidg1_Ax",acc_x1,"phidg1_Ax[num_phidg1_points]/Double_t");
    tree->Branch("phidg1_Ay",acc_y1,"phidg1_Ay[num_phidg1_points]/Double_t");
    tree->Branch("phidg1_Az",acc_z1,"phidg1_Az[num_phidg1_points]/Double_t");
    tree->Branch("phidg1_Bx",x1_field,"phidg1_Bx[num_phidg1_points]/Double_t");
    tree->Branch("phidg1_By",y1_field,"phidg1_By[num_phidg1_points]/Double_t");
    tree->Branch("phidg1_Bz",z1_field,"phidg1_Bz[num_phidg1_points]/Double_t");
    tree->Branch("phidg1_Btot",tot1_field,"phidg1_Btot[num_phidg1_points]/Double_t");
    tree->Branch("phidg1_tilt",tilt_phid1,"phidg1_tilt[num_phidg1_points]/Double_t");

    tree->Branch("num_phidg3_points",&num_phidg3_points,"num_phidg3_points/Int_t");
    tree->Branch("phidg3_Ax",acc_x3,"phidg3_Ax[num_phidg3_points]/Double_t");
    tree->Branch("phidg3_Ay",acc_y3,"phidg3_Ay[num_phidg3_points]/Double_t");
    tree->Branch("phidg3_Az",acc_z3,"phidg3_Az[num_phidg3_points]/Double_t");
    tree->Branch("phidg3_Bx",x3_field,"phidg3_Bx[num_phidg3_points]/Double_t");
    tree->Branch("phidg3_By",y3_field,"phidg3_By[num_phidg3_points]/Double_t");
    tree->Branch("phidg3_Bz",z3_field,"phidg3_Bz[num_phidg3_points]/Double_t");
    tree->Branch("phidg3_Btot",tot3_field,"phidg3_Btot[num_phidg3_points]/Double_t");
    tree->Branch("phidg3_tilt",tilt_phid3,"phidg3_tilt[num_phidg3_points]/Double_t");

    tree->Branch("num_phidg4_points",&num_phidg4_points,"num_phidg4_points/Int_t");
    tree->Branch("phidg4_Ax",acc_x4,"phidg4_Ax[num_phidg4_points]/Double_t");
    tree->Branch("phidg4_Ay",acc_y4,"phidg4_Ay[num_phidg4_points]/Double_t");
    tree->Branch("phidg4_Az",acc_z4,"phidg4_Az[num_phidg4_points]/Double_t");
    tree->Branch("phidg4_Bx",x4_field,"phidg4_Bx[num_phidg4_points]/Double_t");
    tree->Branch("phidg4_By",y4_field,"phidg4_By[num_phidg4_points]/Double_t");
    tree->Branch("phidg4_Bz",z4_field,"phidg4_Bz[num_phidg4_points]/Double_t");
    tree->Branch("phidg4_Btot",tot4_field,"phidg4_Btot[num_phidg4_points]/Double_t");
    tree->Branch("phidg4_tilt",tilt_phid4,"phidg4_tilt[num_phidg4_points]/Double_t");

    //Helmholtz Coil related 
    tree->Branch("coil_event",&counter_mag,"coil_event/Int_t");
    tree->Branch("I_coil1",&curr_coil1,"I_coil1/Double_t");
    tree->Branch("I_coil2",&curr_coil2,"I_coil2/Double_t");
    tree->Branch("I_coil3",&curr_coil3,"I_coil3/Double_t");
    tree->Branch("I_coil4",&curr_coil4,"I_coil4/Double_t");
    tree->Branch("I_coil5",&curr_coil5,"I_coil5/Double_t");
    tree->Branch("I_coil6",&curr_coil6,"I_coil6/Double_t");

    tree->Branch("U_coil1",&volt_coil1,"U_coil1/Double_t");
    tree->Branch("U_coil2",&volt_coil2,"U_coil2/Double_t");
    tree->Branch("U_coil3",&volt_coil3,"U_coil3/Double_t");
    tree->Branch("U_coil4",&volt_coil4,"U_coil4/Double_t");
    tree->Branch("U_coil5",&volt_coil5,"U_coil5/Double_t");
    tree->Branch("U_coil6",&volt_coil6,"U_coil6/Double_t");

    //PMT voltage/current related
    tree->Branch("I_HPD_enable",&curr_hpd_enable,"I_HPD_enable/Double_t");
    tree->Branch("I_HPD_HV_control",&curr_hpd_hv_control,"I_HPD_HV_control/Double_t");
    tree->Branch("I_HPD_LV_control",&curr_hpd_lv_control,"I_HPD_LV_control/Double_t");
    tree->Branch("I_receiver0",&curr_receiver0,"I_receiver0/Double_t");
    tree->Branch("I_receiver1",&curr_receiver1,"I_receiver1/Double_t");
    tree->Branch("I_monitor0",&curr_monitor0,"I_monitor0/Double_t");
    tree->Branch("I_monitor1",&curr_monitor1,"I_monitor1/Double_t");

    tree->Branch("U_HPD_enable",&volt_hpd_enable,"U_HPD_enable/Double_t");
    tree->Branch("U_HPD_HV_control",&volt_hpd_hv_control,"U_HPD_HV_control/Double_t");
    tree->Branch("U_HPD_LV_control",&volt_hpd_lv_control,"U_HPD_LV_control/Double_t");
    tree->Branch("U_receiver0",&volt_receiver0,"U_receiver0/Double_t");
    tree->Branch("U_receiver1",&volt_receiver1,"U_receiver1/Double_t");
    tree->Branch("U_monitor0",&volt_monitor0,"U_monitor0/Double_t");
    tree->Branch("U_monitor1",&volt_monitor1,"U_monitor1/Double_t");

    //Histogram for waveform
    SampleWave0 = new TH1F("","",200, 300, 500);
    SampleWave0->SetDirectory(0);
    SampleWave1 = new TH1F("","",200,300,500);
    SampleWave1->SetDirectory(0);
    // 9.Nov.2017
    SampleWave2 = new TH1F("", "", 200, 300, 500);
    SampleWave2->SetDirectory(0);
    SampleWave3 = new TH1F("", "", 200, 300, 500);
    SampleWave3->SetDirectory(0);
    SampleWave4 = new TH1F("", "", 200, 300, 500);
    SampleWave4->SetDirectory(0);
    SampleWave5 = new TH1F("", "", 200, 300, 500);
    SampleWave5->SetDirectory(0);


    StartVal0 = new TH1F("","",200,300,500);
    StartVal0->SetDirectory(0);

    ngoodTDCbanks = 0;
    nbadTDCbanks = 0;
    

    // Fill a settings tree using information from ODB
#ifdef INCLUDE_MVODB_H
    std::cout << "Filling Setting TTree " << std::endl;

    MVOdb* odb = GetODB();
    TTree *settings_tree = new TTree("settings_tree","Settings Tree");

    int channel_mask;
    settings_tree->Branch("channel_mask",&channel_mask,"channel_mask/Int_t"); 
    
    double HVsetpoints[20];    
    settings_tree->Branch("HVsetpoints",&HVsetpoints,"HVsetpoints[20]/Double_t"); 

    double HVreadback[20];    
    settings_tree->Branch("HVreadback",&HVreadback,"HVreadback[20]/Double_t"); 

    double HVcurrent[20];    
    settings_tree->Branch("HVcurrent",&HVcurrent,"HVcurrent[20]/Double_t"); 

    double calc_baseline[20];    
    settings_tree->Branch("CalcBaseline",&calc_baseline,"CalcBaseline[20]/Double_t"); 

    // BRB Settings and readback
    odb->RI("/Equipment/BRB/Settings/channel mask",&channel_mask);
        
    // PMT settings and readback
    std::vector<int> hv; // set point
    odb->RIA("/Equipment/PMTS/Settings/HVset",&hv);
    std::vector<double> readback; // HV readback
    odb->RDA("/Equipment/PMTS/Variables/PMV0",&readback);
    std::vector<double> current; // HV readback
    odb->RDA("/Equipment/PMTS/Variables/PMI0",&current);
    
    // Get calculated baseline
    std::vector<double> baseline; // HV readback
    odb->RDA("/Analyzer/Baselines/Baseline",&baseline);
    
    
    // Stupidly need to copy from vector to array, I think
    for(int i = 0; i < 20; i++){
      HVsetpoints[i] = hv[i];
      HVreadback[i] = readback[i];
      HVcurrent[i] = current[i];
      calc_baseline[i] = baseline[i];

      std::cout << "Baseline " << i << "  "<< baseline[i] << " " << calc_baseline[i] << std::endl;
    }

    settings_tree->Fill();
#endif  // Done ifdef settings ttree filling


  }

  void EndRun(int transition,int run,int time){
    //tree->Fill();
    std::cout << "Custom end run " << run <<std::endl;
    //tree->Write(); //: happens through in fOutputFile->Write() in base class (ProcessMidasFile)

    std::cout << "Good TDC banks: " << ngoodTDCbanks << std::endl;;
    std::cout << "Bad  TDC banks: " << nbadTDCbanks << std::endl;;
  }

  bool ProcessMidasEvent(TDataContainer& dataContainer){


    //attempt at adding digitizer data
    //Set input parameters which correspond to ADC measurements
    start_val_stat = 90;
    window_width = 45;
    trigger = 40;
    
    // Get BRB data
    TBRBRawData *brb_b = dataContainer.GetEventData<TBRBRawData>("BRB0");
    
    if(brb_b){      
      
      if(num_points<nPoints_max){
        
        num_points++;
        
        std::vector<RawBRBMeasurement> measures = brb_b->GetMeasurements();
        
        for(int i = 0; i < measures.size(); i++){           
          int chan = measures[i].GetChannel();
          for(int ib = 0; ib < measures[i].GetNSamples(); ib++){
            if(chan == 0) V1730_wave0[num_points-1][ib] = measures[i].GetSample(ib);  
            if(chan == 1) V1730_wave1[num_points-1][ib] = measures[i].GetSample(ib);  
            if(chan == 2) V1730_wave2[num_points-1][ib] = measures[i].GetSample(ib);  
            if(chan == 3) V1730_wave3[num_points-1][ib] = measures[i].GetSample(ib);
            if(chan == 4) V1730_wave4[num_points-1][ib] = measures[i].GetSample(ib);
            if(chan == 5) V1730_wave5[num_points-1][ib] = measures[i].GetSample(ib);
            if(chan == 6) V1730_wave6[num_points-1][ib] = measures[i].GetSample(ib);
            if(chan == 7) V1730_wave7[num_points-1][ib] = measures[i].GetSample(ib);
            if(chan == 8) V1730_wave8[num_points-1][ib] = measures[i].GetSample(ib);
            if(chan == 9) V1730_wave9[num_points-1][ib] = measures[i].GetSample(ib);
            if(chan == 10) V1730_wave10[num_points-1][ib] = measures[i].GetSample(ib);
            if(chan == 11) V1730_wave11[num_points-1][ib] = measures[i].GetSample(ib);
            if(chan == 12) V1730_wave12[num_points-1][ib] = measures[i].GetSample(ib);
            if(chan == 13) V1730_wave13[num_points-1][ib] = measures[i].GetSample(ib);
            if(chan == 14) V1730_wave14[num_points-1][ib] = measures[i].GetSample(ib);
            if(chan == 15) V1730_wave15[num_points-1][ib] = measures[i].GetSample(ib);
            if(chan == 16) V1730_wave16[num_points-1][ib] = measures[i].GetSample(ib);
            if(chan == 17) V1730_wave17[num_points-1][ib] = measures[i].GetSample(ib);
            if(chan == 18) V1730_wave18[num_points-1][ib] = measures[i].GetSample(ib);
            if(chan == 19) V1730_wave19[num_points-1][ib] = measures[i].GetSample(ib);
          }              
        }	      
        
      }else{
        std::cout << "Too many points! " << num_points << std::endl;
      }
      
      tree->Fill();
      counter = 0;
      num_points = 0;
      num_points_dig0 = 0;
      num_points_dig1 = 0;
      num_phidg0_points = 0;
      num_phidg1_points = 0;
      num_phidg3_points = 0;
      num_phidg4_points = 0;
      gbl_accept_banks = FALSE;
      
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

int main(int argc, char *argv[])
{

  ScanToTreeConverter::CreateSingleton<ScanToTreeConverter>();
  return ScanToTreeConverter::Get().ExecuteLoop(argc, argv);

}


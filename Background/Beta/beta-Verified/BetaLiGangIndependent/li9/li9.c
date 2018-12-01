
const Double_t pi=TMath::Pi();
const Double_t alpha1=1/137.03599976;
const Double_t E0=0.510998902;

//---------------------------------------------------
//        unique forbidden
//   for delta I=2-，3+，4-...
//approximation:for light and medium nuclei

Double_t unique(Double_t *x, Double_t *par)
{ 
  
  Double_t Eg=par[4];//Eg是gamma能量
  Double_t E=x[0]-Eg;//E电子能量，单位keV
  Double_t Ch=par[5]; 
  //Double_t Zda=-1*par[7]*(par[0]-par[7]);//zdaugther=-ch*zdaygther
  Double_t Zpa=par[0];//zparent
  Double_t Q=par[1];
  Double_t A=par[2];
  Double_t n=par[3];  


  if(E<0 || E>Q ){   //电子能量小于Q值  
    return 0;    
  }  
  if(E>0 && E<Q){
    
  Double_t W=E/E0+1; //energy with elecron rest mass unit
  Double_t W0=Q/E0+1;
  Double_t result=TMath::Factorial(2*n+2); 
 
  Double_t q=W0-W;
  Double_t p=sqrt(W**2-1);
  Double_t Cn=2**(n-1)*((p+q)**(2*n+2)-(p-q)**(2*n+2))/(p*q*result);

  Double_t ff=Cn  //forbiddenness
             *fermi(Zpa,A,E,Ch)*screening(Zpa,A,E,Ch,1.45)
             *finitesizeEM1(Zpa,A,E,Ch)
             *finitesizeWI(Zpa,A,E,Ch,Q)
             *weakmag(Zpa,A,E,Ch)
             *p*W*(W0-W)**2;
  
  return ff;  
  }  
}
//-----------------------------------------------

//---------------------------------------------------
//        allowed
//   for delta I=0+,1+...
//   for non-unique forbidden 0-,1- approximately

Double_t allowed(Double_t *x, Double_t *par)
{ 
  
  Double_t Eg=par[4];//Eg是gamma能量
  Double_t E=x[0]-Eg;//E电子能量，单位keV
  Double_t Ch=par[5]; 
 // Double_t Zda=-1*par[7]*(par[0]-par[7]);//zdaugther=-ch*zdaygther
  Double_t Zpa=par[0];//zparent
  Double_t Q=par[1];
  Double_t A=par[2];
  Double_t ft=par[3];  


  if(E<0 || E>Q ){   //电子能量小于Q值  
    return 0;    
  }  
  if(E>0 && E<Q){
    
  Double_t W=E/E0+1; //energy with elecron rest mass unit
  Double_t W0=Q/E0+1;
   Double_t p=sqrt(W**2-1);

  Double_t ff=log(2)/ft
             *fermi(Zpa,A,E,Ch)*screening(Zpa,A,E,Ch,1.45) 
             *finitesizeEM1(Zpa,A,E,Ch)
             *finitesizeWI(Zpa,A,E,Ch,Q)
             *weakmag(Zpa,A,E,Ch)            
             *p*W*(W0-W)**2;
             
  return ff;  
  }  
}
//-----------------------------------------------------
//     unique forbidden shape factor form 2
// Kai Siegbahn. Alpha-, beta- and gamma-ray spectroscopy.
//

Double_t unique2(Double_t *x, Double_t *par)
{ 
  Double_t Eg=par[4];//Eg是gamma能量
  Double_t E=x[0]-Eg;//E电子能量，单位keV
  Double_t Ch=par[5]; 
 // Double_t Zda=-1*par[7]*(par[0]-par[7]);//zdaugther=-ch*zdaygther
  Double_t Zpa=par[0];//zparent 
  Double_t Q=par[1];
  Double_t A=par[2];
  Double_t n=par[3];  
  
    
  if(E<0 || E>Q ){   //电子能量小于Q值  
    return 0;    
  } 
  if(E>0 && E<Q){
    
  Double_t W=E/E0+1; //energy with elecron rest mass unit  
  Double_t W0=Q/E0+1;
  Double_t p=sqrt(W**2-1);

  Double_t ff=forbidden(Zpa,A,E,Ch,Q,n)
             *fermi(Zpa,A,E,Ch)*screening(Zpa,A,E,Ch,1.45)
             *finitesizeEM1(Zpa,A,E,Ch)
             *finitesizeWI(Zpa,A,E,Ch,Q)
             *weakmag(Zpa,A,E,Ch)
             *p*W*(W0-W)**2;
  return ff;  
  }  
}

//----------------------------------------------------

void li9(){
  
 
  Double_t Z=3;
  Double_t A=9;
  Double_t Ch=-1;
  
  TF1 *f3 = new TF1("f3", unique2, 0,13.607, 6);
  TF1 *f4 = new TF1("f4", allowed, 0,13.607, 6);  
  f3->SetParameter(0,Z);
  f3->SetParameter(2,A);
  f3->SetParameter(5,Ch);
  f4->SetParameter(0,Z);
  f4->SetParameter(2,A); 
  f4->SetParameter(5,Ch); 
  
  Double_t binc1[13610];//记录瞬时bincontent
  Double_t binc2[13610];//记录叠加后的bincontent
  //----------初始化-------------
  for(Int_t i=0;i<=13609;i++){
    binc2[i]=0;                                                    
  }
  
  #if 1
 //-----------------forbidden-------------------------------      
  
  //------读入禁戒跃迁数据-----------
  Double_t data1[1][4];//输入数据，顺序：Br1,n,Q,Eg,Br2
  ifstream fin("li9for.dat");  
  for(Int_t m=0;m<=0;m++){
    for(Int_t n1=0;n1<=3;n1++){
      fin>>data1[m][n1];
    }
  } 	
  fin.close(); 
  //--------叠加所有禁戒谱----------
  for(Int_t m=0;m<=0;m++){    
    Double_t Br=data1[m][1];
    Double_t n=data1[m][2];
    Double_t Q=data1[m][0];
    Double_t Eg=data1[m][3];       

    f3->SetParameter(3,n);
    f3->SetParameter(1,Q);
    f3->SetParameter(4,Eg);
   
    TH1F *h1 = new TH1F("h1","h1",13607,0,13.607);
    h1->Add(f3); 
    //-----归一化   
    double entries = h1->Integral();//1796,2150
    h1->Scale(1./entries);
      
    for(Int_t i=1;i<=13607;i++){   
      binc1[i]=h1->GetBinContent(i);
     // h1->Delete();         
    }       
    
    #if 0
    //-------------归一化--------------------------------
    Double_t sumbin=0;
    for(Int_t i=1;i<=500;i++){   
      sumbin +=binc1[i];          
    }  
    for(Int_t i=1;i<=500;i++){   
      binc1[i]=binc1[i]/sumbin;          
    }
    #endif
    for(Int_t i=1;i<=13607;i++){
      binc2[i]+=Br*binc1[i];
    }   
   
  }
  
  #endif
  
  
  #if 1
  //-----------------allowed-------------------------------
  
  //------读入允许跃迁数据-----------
  Double_t data[5][3];//输入数据，顺序：Br1,logft,Q,Eg,Br2
  ifstream fin("li9all.dat");  
  for(Int_t m=0;m<=4;m++){
    for(Int_t n1=0;n1<=2;n1++){
      fin>>data[m][n1];
    }
  } 	
  fin.close(); 
  //--------叠加所有谱----------
  for(Int_t m=0;m<=4;m++){    
    Double_t Br=data[m][1];
 //   Double_t logft=data[m][1];
    Double_t Q=data[m][0];
    Double_t Eg=data[m][2];

 //   Double_t Br2=data[m][4];
    Double_t ft=1;//10**(logft);        

    f4->SetParameter(3,ft);
    f4->SetParameter(1,Q);
    f4->SetParameter(4,Eg);
    
    TH1F *h1 = new TH1F("h1","h1",13607,0,13.607);
    h1->Add(f4);
        //-----归一化   
    double entries = h1->Integral();//1796,2150
    h1->Scale(1./entries);
    
    for(Int_t i=1;i<=13607;i++){   
      binc1[i]=h1->GetBinContent(i);
     // h1->Delete();        
    } 
    
    #if 0
    
      //-------------归一化--------------------------------
    Double_t sumbin=0;
    for(Int_t i=1;i<=500;i++){   
      sumbin +=binc1[i];          
    }  
    for(Int_t i=1;i<=500;i++){   
      binc1[i]=binc1[i]/sumbin;          
    }
    #endif
    
    for(Int_t i=1;i<=13607;i++){
      binc2[i]=binc2[i]+Br*binc1[i];
    }
   
          
  }
  
  #endif
  
    #if 0
      //-------------归一化--------------------------------
    Double_t sumbin=0;
    for(Int_t i=1;i<=13607;i++)
    {   
      sumbin +=binc2[i];          
    }  
    for(Int_t i=1;i<=13607;i++){   
      binc2[i]=binc2[i]/sumbin;          
    }
   #endif
   
   //------------输出归一化能谱bincontent到文件----------------------------  
  ofstream fon("li9bin.dat");  
  for(Int_t i=1;i<=13607;i++){
      fon<<binc2[i]<<endl;
  } 	
  fon.close();
   
  //-------------绘制总beta能谱----------------------------
  
  TFile *file =new TFile("myspec.root","RECREATE");
  TH1F *h2=new TH1F("h2","li9 Beta Sprectrum",13607,0,13.607);
  for(Int_t i=1;i<=13607;i++){
    h2->SetBinContent(i,binc2[i]);    
  }
 // TCanvas *myTC = new TCanvas("myTC","A Canvas",10,10,800,600);
  //myTC->cd();
  h2->Draw();
  h2->GetXaxis()->SetTitle("Beta Kinetic Energy/MeV");
 // h2->GetYaxis()->SetTitle("probability density");
 // myTC->SaveAs("Bi214.eps");
  file->cd();
  h2->Write(); 
  
  


}

Double_t Gammacomplex(TComplex z)
{
  // complex gamma function
  // modified f2c of CERNLIB cgamma64.F
  // taken from Knat and converted to compile without ROOT libs

  Double_t zr=z.Re();
  Double_t zi=z.Im();
  if (zi==0.0) {
    if (zr==0.0) {
      TComplex w(0.0,0.0);
      return w;
    } else if (-fabs(zr)==floor(zr)) {
      TComplex w(0.0,0.0);
      return w;
    }
  }

  TComplex u;
  TComplex v;

  if (zr>=1.0) {
    u=1.0;
    v=z;
  } else if (zr>=0.0) {
    u=1.0/z;
    v=1.0+z;
  } else {
    u=1.0;
    v=1.0-z;
  }

  const Double_t c1=2.50662827463100050;
  const Double_t c2[16]={
     41.624436916439068,-51.224241022374774,
     11.338755813488977, -0.747732687772388,
      0.008782877493061, -0.000001899030264,
      0.000000001946335, -0.000000000199345,
      0.000000000008433,  0.000000000001486,
     -0.000000000000806,  0.000000000000293,
     -0.000000000000102,  0.000000000000037,
     -0.000000000000014,  0.000000000000006
  };

  TComplex w(1.0,0.0);
  TComplex s(c2[0],0.0);
  TComplex cpi(3.14159265358979312,0.0);

  for (int i=1;i<16;++i) {
    w*=(v-((Double_t)i))/(v+((Double_t)(i-1)));
    s+=c2[i]*w;
  }

  w=v+4.5;
  TComplex logw = TComplex::Log(w);
  w=c1*s*TComplex::Exp((v-0.5)*TComplex::Log(w)-w);
  if (zr<0.0) w=cpi/(TComplex::Sin(cpi*z)*w);
  TComplex result=w*u;
  return result.Rho2();
}



Double_t Elecw(Double_t kk,Double_t Z,Double_t EE){
 
const Double_t E0=0.510998902;//0.510998902;
const Double_t pi=TMath::Pi();//3.14159265
const Double_t alpha1=1/137.03599976;
 Double_t gam0=sqrt(1-(Z*alpha1)**2);
 Double_t gam=sqrt(kk**2-(Z*alpha1)**2);
 Double_t p=sqrt(EE**2-1);
 Double_t nu=Z*EE*alpha1/p;
 Double_t y1=0.5*(2.*p)**(2.*gam-2.*gam0);
 Double_t y2=(TMath::Gamma(2*gam0+1))**2/(TMath::Gamma(2*gam+1))**2;
 TComplex z1(gam,nu);
 TComplex z2(gam0,nu);
 
 Double_t m1=Gammacomplex(z1);
 Double_t m2=Gammacomplex(z2);

 Double_t L=y1*y2*kk*(kk+gam)*m1/m2;
 return L;
}


//---------------------------------------------------
//      
Double_t finitesizeEM1(Double_t ZZ,Double_t AA,Double_t EE,Double_t CCh)
{ 
  

const Double_t E0=0.510998902;//0.510998902;
const Double_t pi=TMath::Pi();//3.14159265
const Double_t alpha1=1/137.03599976;

  Double_t ZZda=-1*CCh*(ZZ-CCh);
  Double_t W=EE/E0+1; //energy with elecron rest mass unit 
  Double_t gam=sqrt(1-(ZZda*alpha1)**2);
// Double_t ra=AA**(1.0/3.0)*alpha1/2.0;
  Double_t ra=0.0029*AA**(1./3) + 0.0063*AA**(-1./3) - 0.0017/AA; 
  Double_t an;
  Double_t data1[7][6]=
  {
  {0.115,-1.8123,8.2498,-11.223,-14.854,32.086},
  {-0.00062,0.007165,0.01841,-0.53736,1.2691,-1.5467},
  {0.02482,-0.5975,4.84199,-15.3374,23.9774,-12.6534},
  {-0.14038,3.64953,-38.8143,172.137,-346.708,288.787},
  {0.008152,-1.15664,49.9663,-273.711,657.629,-603.703},
  {1.2145,-23.9931,149.972,-471.299,662.191,-305.68},
  {-1.5632,33.4192,-255.133,938.53,-1641.28,1095.36}
  };
  Double_t a1[7];
 for(Int_t i=0;i<=6;i++){
  for(Int_t j=0;j<=5;j++){
   a1[i]=data1[i][j]*((alpha1*abs(ZZda))**(j+1)) ;   
  }
  }

   Double_t data2[7][6]=
  {
  {-0.0701,-2.572,-27.5971,-128.658,-272.264,-214.925},
  {0.002308,0.066483,0.6407,2.83606,5.6317,4.0011},
  {-0.07936,-2.09284,-18.45462,-80.9375,-160.8384,-124.8927},
  {0.93832,22.02513,197.00221,807.1878,1566.6077,1156.3287},
  {-4.276181,-96.82411,-835.26505,-3355.8441,-6411.3255,-4681.573},
  {8.2135,179.0862,1492.1295,5872.5362,11038.7299,7963.4701},
  {-5.4583,-115.8922,-940.8305,-3633.9181,-6727.6296,-4795.0481}
  };
  Double_t a2[7];
 for(Int_t i=0;i<=6;i++){
  for(Int_t j=0;j<=5;j++){
   a2[i]=data2[i][j]*((alpha1*abs(ZZda))**(j+1)) ;   
  }
  } 

 
 if(CCh<0){
 
 #if 1
   for(Int_t i=1;i<=6;i++){
  an+=a1[i]*(W*ra)**(i-1);
}
 #endif 
  Double_t L1=a1[0]*ra/W+an+0.41*(ra-0.0164)*(alpha1*ZZda)**(4.5);
  Double_t L=1+13./60.*(alpha1*ZZda)**2-W*ra*alpha1*ZZda*(41-26*gam)/15./(2*gam-1)-alpha1*ZZda*ra*gam*(17-2*gam)/30./W/(2*gam-1)+L1; 
  return L;  
   
  }
  
  
  if(CCh>0){

    
#if 1
   for(Int_t i=1;i<=6;i++){
  an+=a2[i]*(W*ra)**(i-1);
}
 #endif 

 #endif 
  Double_t L1=a2[0]*ra/W+an+0.22*(ra-0.0164)*(alpha1*abs(ZZda))**(4.5);
  Double_t L=1+13./60.*(alpha1*ZZda)**2-W*ra*alpha1*(ZZda)*(41-26*gam)/15./(2*gam-1)-alpha1*(ZZda)*ra*gam*(17-2*gam)/30./W/(2*gam-1)+L1; 
  return L;  
  
  }
   
}

#endif


#if 1
//---------------------------------------------------
//      


Double_t finitesizeEM2(Double_t ZZ,Double_t AA,Double_t EE,Double_t CCh)
{ 
  

const Double_t E0=0.510998902;//0.510998902;
const Double_t pi=TMath::Pi();//3.14159265
const Double_t alpha1=1/137.03599976;
//EE=0.05;
  Double_t ZZda=-1*CCh*(ZZ-CCh);
  Double_t W=EE/E0+1; //energy with elecron rest mass unit 
  Double_t gam=sqrt(1-(ZZda*alpha1)**2);
// Double_t ra=AA**(1.0/3.0)*alpha1/2.0;
  Double_t ra=0.0029*AA**(1./3) + 0.0063*AA**(-1./3) - 0.0017/AA;

if(CCh<0){


  Double_t L=1+13./60.*(alpha1*ZZda)**2-W*ra*alpha1*ZZda*(41-26*gam)/15./(2*gam-1)-alpha1*ZZda*ra*gam*(17-2*gam)/30./W/(2*gam-1);

  return L;  
  
  }
  
 if(CCh>0){

     
  Double_t L=1+13./60.*(alpha1*(ZZda))**2-W*ra*alpha1*(ZZda)*(41-26*gam)/15./(2*gam-1)-alpha1*(ZZda)*ra*gam*(17-2*gam)/30./W/(2*gam-1);

  return L;  
   
  } 
  
  
  
    
   
}

#endif




#if 1
//---------------------------------------------------
//       


Double_t finitesizeWI(Double_t ZZ,Double_t AA,Double_t EE,Double_t CCh,Double_t Q)
{ 
  

const Double_t E0=0.510998902;//0.510998902;
const Double_t pi=TMath::Pi();//3.14159265
const Double_t alpha1=1/137.03599976;
//EE=0.08;
   Double_t ZZda=-1*CCh*(ZZ-CCh);
   Double_t W=EE/E0+1; //energy with elecron rest mass unit
   Double_t W0=Q/E0+1;
   //cout<<"V0="<<V0<<endl;

  #if 1

  Double_t gam=sqrt(1-(ZZda*alpha1)**2); 
// Double_t ra=AA**(1.0/3.0)*alpha1/2.0;
  Double_t ra=0.0029*AA**(1./3) + 0.0063*AA**(-1./3) - 0.0017/AA;
 
  Double_t C0=-233./630.*(alpha1*ZZda)**2-(W0*ra)**2/5.+2./35.*W0*ra*alpha1*ZZda;
  Double_t C1=-21./35.*ra*alpha1*ZZda+4./9.*W0*ra**2;
  Double_t C2=-4./9.*ra**2;

  Double_t C=1+C0+C1*W+C2*W**2; 
  

  //    cout<<"ytilt="<<ytilt<<endl;
//    cout<<"yy="<<yy<<endl;
 // cout<<"f="<<f<<endl;
 #endif 
  
  
  return C;  
  
   
}

#endif


#if 1
//---------------------------------------------------
//        fermi子函数


Double_t fermi(Double_t ZZ,Double_t AA,Double_t EE,Double_t CCh)
{ 
  

const Double_t E0=0.510998902;//0.510998902;
const Double_t pi=TMath::Pi();//3.14159265
const Double_t alpha1=1/137.03599976;

  Double_t ZZda=-1*CCh*(ZZ-CCh);
 

    
  Double_t W=EE/E0+1; //energy with elecron rest mass unit
 // Double_t W0=Q/E0+1; 
  
  Double_t gam=sqrt(1-(ZZda*alpha1)**2);

  Double_t yy=ZZda*W*alpha1/sqrt(W**2-1);
// Double_t ra=AA**(1.0/3.0)*alpha1/2.0;
  Double_t ra=0.0029*AA**(1./3) + 0.0063*AA**(-1./3) - 0.0017/AA;

  TComplex bb(gam,yy);
  Double_t xx=Gammacomplex(bb);
  
  Double_t za=(TMath::Gamma(2*gam+1))**2;
 
  Double_t f=2*(1+gam)*(2*sqrt(W**2-1)*ra)**(2*gam-2)*exp(pi*yy)*xx/za;//fermi函数2*(1+gam)*(2*sqrt(W**2-1)*ra)**(2*gam-2)*exp(pi*yy)*xx/za

  return f;  
 
}

#endif



#if 1
//---------------------------------------------------
//        fermi+screening子函数


Double_t screening(Double_t ZZ,Double_t AA,Double_t EE,Double_t CCh,Double_t x)
{ 
  

const Double_t E0=0.510998902;//0.510998902;
const Double_t pi=TMath::Pi();//3.14159265
const Double_t alpha1=1/137.03599976;
//EE=0.08;
   Double_t ZZda=-1*CCh*(ZZ-CCh);
   Double_t W=EE/E0+1; //energy with elecron rest mass unit
   Double_t V0=x*fabs(ZZ)**(4.0/3.0)*(alpha1)**2;  //potential lift
 
   //cout<<"V0="<<V0<<endl;

   
 #if 1
  Double_t Wtilt=W+CCh*V0 > 1 ? W+CCh*V0 : 1.00001 ;
  Double_t ptilt = sqrt(Wtilt*Wtilt-1);
  Double_t ytilt=alpha1 * ZZda * Wtilt / ptilt;
  Double_t gam=sqrt(1-(ZZda*alpha1)**2); 
  Double_t p = sqrt(W*W-1);
  Double_t yy=ZZda*W*alpha1/p;
  
   TComplex bb(gam,yy);
  Double_t xx=Gammacomplex(bb);
   TComplex bb1(gam,ytilt);
  Double_t xx1=Gammacomplex(bb1);
  
 Double_t f = Wtilt/W //*sqrt((Wtilt*Wtilt-1)/(W*W-1))
            * pow(ptilt/p,(2*gam-1))
            * exp((TMath::Pi())*(ytilt-yy))
            * xx1/xx;
 
  //    cout<<"ytilt="<<ytilt<<endl;
//    cout<<"yy="<<yy<<endl;
 // cout<<"f="<<f<<endl;
 #endif 
 
  return f;  

}

#endif

#if 1
//---------------------------------------------------
//       


Double_t weakmag(Double_t ZZ,Double_t AA,Double_t EE,Double_t CCh)
{ 
  

const Double_t E0=0.510998902;//0.510998902;
const Double_t pi=TMath::Pi();//3.14159265
const Double_t alpha1=1/137.03599976;
//EE=0.08;
   Double_t ZZda=-1*CCh*(ZZ-CCh);
   Double_t W=EE/E0+1; //energy with elecron rest mass unit
 
 
  #if 1

 
  Double_t C=1-CCh*0.0048*W*E0; 
  

  //    cout<<"ytilt="<<ytilt<<endl;
//    cout<<"yy="<<yy<<endl;
 // cout<<"f="<<f<<endl;
 #endif 
  
  
  return C;  
  
}

//n级禁戒用全部用n级unique跃迁算
Double_t forbidden(Double_t ZZ,Double_t AA,Double_t EE,Double_t CCh,Double_t Q,Int_t n){
const Double_t E0=0.510998902;//0.510998902;
const Double_t pi=TMath::Pi();//3.14159265
const Double_t alpha1=1/137.03599976;

 Double_t ZZda=-1*CCh*(ZZ-CCh);
 Double_t W=EE/E0+1; //energy with elecron rest mass unit
 Double_t W0=Q/E0+1;
 Double_t q=W0-W;
 Double_t ra=0.0029*AA**(1./3) + 0.0063*AA**(-1./3) - 0.0017/AA;
 Double_t fac2=1;
 Double_t fac1=TMath::Factorial(n);
 Double_t gam0=sqrt(1-(ZZda*alpha1)**2);
 
 for(Int_t i=2*n+1;i>=1;i=i-2){//(2n+1)!!
   fac2 *=i;
 }

  
 Double_t Sn=0;
 for(Double_t k=1;k<=n+1;k++)
 {    
    Double_t fac3=1;
    for(Int_t i=2*k-1;i>=1;i=i-2){//(2k-1)!!
      fac3 *=i;
    }
    Double_t fac4=1;
    for(Int_t i=2*n-2*k+3;i>=1;i=i-2){//(2n-2k+3)!!
      fac4 *=i;
    } 
   
    Double_t x3=(TMath::Factorial(k-1))*(TMath::Factorial(n-k+1));
 
    Sn += 8*pi*fac1/fac2
         *ra**(2*n)/(1+gam0)
         *fac3*(Elecw(k,ZZda,W))*(q**(2*n-2*k+2))
         /fac4/x3  ;       

  }
 
 return Sn;


}

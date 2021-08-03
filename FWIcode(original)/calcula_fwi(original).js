

/***********************************************************************************************************************************
 *                                                                                                                                 *
 * Este programa IMPLEMENTA las ecuaciones definidas para el cálculo del Fire Wheather Index (FWI) definidas                       *                     
 * por el Canadian Forestry Service. Forestry Technical Report 35 (Ottawa 1987) y publicadas en múltiples                          * 
 * trabajos (C. E. Van Wagner et all). Se respetan completamente las ecuaciones definidas en "FWI_Equations_FORTRAN.pdf"           *
 * y "Cálculo_Bariloche_2017-2018.xlsx". NO HAY GARANTIA ALGUNA DE SUS RESULTADOS NI DE SU USO ESPECIFICO.                         *
 *                                                                                                                                 *
 * Esta obra es software libre; puede redistribuirse y/o modificarse de acuerdo con los términos de la                             *
 * Licencia Pública General GNU publicada por la Free Software Foundation, en la versión 2 de la licencia                          *
 * o cualquier otra posterior. Esta obra se distribuye con la esperanza de que pueda ser útil, pero SIN                            *
 * GARANTIA ALGUNA; ni siquiera la garantía implícita de COMERCIALIZACION o ADECUACIÓN A UNA FINALIDAD ESPECÍFICA.                 *
 * Véase la versión 2 y la versión 3 de la Licencia Pública General GNU para conocer más detalles (http://www.gnu.org/licenses/).  *
 *                                                                                                                                 *
 **********************************************************************************************************************************/

"use strict";  // the script works the “modern” way

// effective day-lengths per month for DMC (Daylength factor Le)  /// FLOAT Les[HEMISFERIOS][MESES]
 //let
 const Les = [[11.5,10.5, 9.2, 7.9, 6.8, 6.2,6.5, 7.4, 8.7,10,11.2,11.8],              // South hemisph.
              [6.5,7.5,9,12.8,13.9,13.9,12.4,10.9,9.4,8,7,6]];                         // North hemisph.

// daylength adjustment Lf for DC    /// FLOAT Lfs[HEMISFERIOS][MESES]                 
 const Lfs = [[6.4,5,2.4,0.4,-1.6,-1.6,-1.6,-1.6,-1.6,0.9,3.8,5.8],                    // South hemisph.
              [-1.6,-1.6,-1.6,0.9,3.8,5.8,6.4,5,2.4,0.4,-1.6,-1.6]];                   // North hemisph.

const Hem = 0;

let ISI, BUI, FWI, DSR;
let FFMC=21.2, DMC=0.1, DC=0.6;
let hh=64,tt=6.7,ww=33,ro=0;
let mes=9;

//----------------------------------------------------------------------------------//
function calculaIndices() {
    FFMC = Number(document.getElementById("idFFMC").value);
    DMC = Number(document.getElementById("idDMC").value);
    DC = Number(document.getElementById("idDC").value);
    hh = Number(document.getElementById("idHH").value);
    tt = Number(document.getElementById("idTT").value);
    ww = Number(document.getElementById("idWW").value);
    ro = Number(document.getElementById("idRO").value);
    mes = Number(document.getElementById("idMES").value);    
        
    FFMC = calculaFFMC(hh,tt,ww,ro,FFMC);
    DMC = calculaDMC(hh,tt,ro,DMC,mes-1);  // usa mes como índice que comienza en 0
    DC = calculaDC(tt,ro,DC,mes-1);        // usa mes como índice que comienza en 0
    ISI = calculaISI(FFMC,ww);
    BUI = calculaBUI(DMC,DC);
    FWI = calculaFWI(BUI,ISI);

    DSR = 0.0272 * Math.pow(FWI,1.77); //(31)
    
    document.getElementById("resultados").innerHTML = "nuevos valores: " +
      "FFMC=" + FFMC.toFixed(3) + ", DMC=" + DMC.toFixed(3) + ", DC=" + DC.toFixed(3) + 
      ", ISI=" + ISI.toFixed(3) + ", BUI=" + BUI.toFixed(3) + ", FWI=" + FWI.toFixed(3);
}



              
              
//----------:------------------------------------------------------------------------//
// HH: humedad relativa, TT: temperatura, WW: velocidad viento, Ro: lluvia anterior, Fo: FFMC anterior: FFMC cero
function calculaFFMC(HH, TT, WW, Ro, Fo)
{
  let mo,mr, rf, Ed,Ew,  ko,kd,k1,kw, m;
  let F;
  
  /* Fine fuel moisture code (FFMC) */
  mo = 147.2*(101 - Fo)/(59.5 + Fo);
  rf = ( Ro > 0.5? Ro - 0.5: 0);
  
  mr = mo + 42.5 * rf * Math.exp(-100/(251-mo))*(1 - Math.exp(-6.93/rf));  // (3a)
  if (mo > 150) 
    mr += 0.0015 * (mo - 150) * (mo - 150) * Math.sqrt(rf);           // (3b)
  mo = (mr > 250? 250: mr);                                      // mo toma el valor de mr, pero debe ser menor que 250
  
  Ed = 0.942 * Math.pow(HH,0.679) + 11 * Math.exp((HH-100)/10.0) + 0.18 * (21.1 - TT) * (1 - Math.exp(-0.115*HH));     // (4)
  if (mo > Ed) {
    ko = 0.424 * (1 - Math.pow(HH/100.0,1.7)) + 0.0694 * Math.sqrt(WW) * (1 - Math.pow(HH/100.0,8));                   // (6a)
    kd = ko * 0.581 * Math.exp(0.0365 * TT);                                                                 // (6b)
    m = Ed + (mo - Ed) * Math.pow(10,-kd);                                                                   // (8)  
  } else {
    Ew = 0.618 * Math.pow(HH,0.753) + 10 * Math.exp((HH-100)/10.0) + 0.18 * (21.1 - TT) * (1 - Math.exp(-0.115*HH));   // (5)
    if (mo < Ew) {
      k1 = 0.424 * (1 - Math.pow( (1 - HH/100.0), 1.7)) + 0.0694 * Math.sqrt(WW) * (1 - Math.pow( (1 - HH/100.0), 8)); // (7a)
      kw = k1 * 0.581 * Math.exp(0.0365 * TT);                                                               // (7b)
      m = Ew + (mo - Ew) * Math.pow(10,-kw);                                                                 // (9)
    } else m = mo;
  }
  F = 59.5 * (250 - m)/(147.2 + m);                                                                     // (10)
  
  return F;
}  

//----------------------------------------------------------------------------------//
function calculaDMC(HH, TT, Ro, Po, MM)
{
  let re, Mo, b, Mr, Pr, K;
  
  if (Ro > 1.5) {
    re = 0.9*Ro - 1.27;                         // (11)
    Mo = 20 + Math.exp(5.6348 - Po/43.43);           // (12)

    if (Po<=33) b = 100/(0.5 + 0.3*Po);         // (13a)
    else if (Po<=65) b = 14 - 1.3*log(Po);      // (13b)
    else b = 6.2*log(Po) - 17.2;                // (13c)
  
    Mr = Mo + 1000*re/(48.77 + b*re);           // (14)
    Pr = 244.72 - 43.43*log(Mr-20);             // (15)
    Po = (Pr < 0? 0: Pr);
  } 
  TT = (TT < -1.1? -1.1: TT);
  K = 1.894*(TT + 1.1) * (100-HH) * Les[Hem][MM] * 0.0001; // (16)

  return (Po + K);                              // (17)  no multiplica K por 100 porque en (16) multip.por 10^-4
}  

//----------------------------------------------------------------------------------//
function calculaDC(TT,  Ro,  Do, MM)
{
  let rd, Qr, Dr, V;
  
  if (Ro > 2.8) {
    rd = 0.83*Ro - 1.27;                        // (18)
    Qr = 800 * Math.exp(-Do/400.0) + 3.937*rd;       // (19, 20)
    Dr = 400 * log(800.0 / Qr);                 // (21)
    Do = (Dr < 0? 0: Dr);
  }
  TT = (TT < -2.8? -2.8: TT);  
  V = 0.36 * (TT + 2.8) + Lfs[Hem][MM];         // (22)
  V = (V < 0? 0: V);
    
  return Do + 0.5 * V;                          // (23)
}

//----------------------------------------------------------------------------------//
function calculaISI(F,  WW)
{
  let m, fW, fF;
  fW = Math.exp(0.05038*WW);                              // (24)
  m = 147.2*(101 - F)/(59.5 + F);
  fF = 91.9*Math.exp(-0.1386*m)*(1+Math.pow(m,5.31)/(4.93e7)); // (25)
  return 0.208 * fW * fF;                            // (26)
}

//----------------------------------------------------------------------------------//
function calculaBUI(P,  D)
{
  let U;
  if (P <= 0.4*D) U = 0.8*D*P/(P+0.4*D);             // (27a)
  else U = P - (1-0.8*D/(P+0.4*D)) * (0.92 + Math.pow(0.0114*P, 1.7)); //(27b)
  return (U<0? 0: U);                                // de acuerdo al Excel
}

//----------------------------------------------------------------------------------//
function calculaFWI(U,  R)
{
  let fD, S,B;
  if (U > 80) fD = 1000/(25+108.64*Math.exp(-0.023*U));   // (28b)
  else fD = 0.626*Math.pow(U,0.809) + 2;                  // (28a)
  B = 0.1 * R * fD;                                  // (29)

  if (B>1) S = Math.exp(2.72*Math.pow(0.434*log(B),0.647));     // (30a)
  else S = B;                                        // (30b)
  
  return S;
}



//----------------------------------------------------------------------------------//
/*
function peligrosidadXtipo(isi, bui, GSecPorc, *pelig){  // *pelig
  let Ibui, Iisi;
 
  // peligrosidad Pastizal
  if (isi < 1 || (isi < 2 && GSecPorc < 56) || (isi < 11 && GSecPorc < 51)) pelig[0] = 1;
  else if (isi < 3 || GSecPorc < 53) pelig[0] = 2;
  else if (isi < 8 || GSecPorc < 58) pelig[0] = 3; 
  else if (isi < 12 || GSecPorc < 66) pelig[0] = 4;
  else pelig[0] = 5;

  // peligrosidad Arbustal
  if (isi <= 0.5) pelig[1] = 1;
  else if (isi <= 1.0) pelig[1] = 2;
  else if (isi <= 2.0) pelig[1] = 3;
  else if (isi <= 3.0) pelig[1] = 4;
  else pelig[1] = 5;  
  
  // peligrosidad Bosque Tipo A
  Ibui = (bui < 383 ? (bui < 0 ? 0 : (int) trunc(bui))  : 382);
  if (isi >= 3 && isi < 186) Iisi = (int) trunc(isi) + 3;
  else if (isi < 0) Iisi = 0;
  else if (isi < 3) Iisi = (int) trunc(isi*2.0); 
  else Iisi = 188;
  pelig[2] = tipoA[Ibui][Iisi];
  
  // peligrosidad Bosque Tipo B
  Ibui = (bui < 71 ? (bui < 0 ? 0 : (int) trunc(bui)) : 70);
  Iisi = (isi < 162 ? (isi < 0 ? 0 : (int) trunc(isi)) : 161);
  pelig[3] = tipoB[Iisi][Ibui];

  // peligrosidad Plantaciones
  Iisi = (isi < 201 ? (isi < 0 ? 0 : (int) trunc(isi)) : 200);
  pelig[4] = Plantaciones[Iisi][Ibui];

 // printf("%.2f %.2f -- %d %d %d %d %d \n", isi, bui, pelig[0], pelig[1], pelig[2], pelig[3], pelig[4]);

} 

*/
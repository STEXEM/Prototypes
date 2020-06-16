$EOLCOM //
$OnText
 Developed by
    Isaac Camilo González Romero
    Instituto de Investigación Tecnológica
    isaac.gonzalez@iit.comillas.edu
    June, 2020
$OffText

//////////////////////////         PARAMETER DEFINITON                    //////////////////////////////

 PARAMETER

  pMin
  pObjective
  pMax
  pExtendedWelfare1
  pExtendedWelfare2
  pExtendedWelfare05
  pTotalCost05
  pCompetitionType
  pOptimal
  pRegret
  is_elastic
  pMercado
 ;

//////////////////////////        OPTIONS INITIALIZATION                  //////////////////////////////
  
  pObjective       = 1  ;   // 0= Cost Min  1= MaxWelf 
  pCompetitionType = 1  ;   // 0= perfect   1= Cournot   2=  Intermediate
  pMercado         = 3  ;   //              1= MCP       2=  MIP(Regularization) 3=QCP
  pRegret          = 0  ;   // 0= No Regret 1= Actual    2=  Naive
  is_elastic       = yes;
  pDiv=10;

//////////////////////////       FINDING OPTIMUM FROM  BILEVEL MODEL      //////////////////////////////
 $$include "C:\Users\igonzalezr\Documents\0_Programs\1_STEXEM\0.5_EJECUTABLES_(May2019_ now)\03_EJECUTABLES_LAST\00_GMS\STEXEM_Europe.gms"

 set
 qq(q);

 if ( pMercado = 3, // Post-procesing for Chosing Optimal From Enumeration
  lc(di,df) = no;
  la(di,df) = no;

  lc('21_fr','96_ie') = Yes;
  lc('26_fr','90_uk') = Yes;
  lc('31_de','79_no') = Yes;
  lc('41_pl','77_lt') = Yes;
  lc('55_it','68_gr') = Yes;
  lc('73_ee','78_lv') = Yes;

  la(di,df) $(lc(di,df) or le(di,df))=yes;
  vNewline.fx (y,       lc) = 0   ;
  pFlow       (y,pa(p), la) = 0   ;
  pInstalLine (y,       lc) = 0   ;
  qq(q) = no;
  if ( pObjective = 0, // Saving vector of minimal Cost
    pMin = smin(q, pObjFunct (q,'4'));
    pAct(q,di,df)$(pObjFunct (q , '4') <> pMin)=no;
    qq(q)$(pObjFunct (q , '4') = pMin)= yes;
    display qq;

     loop (qq,
     lc(di,df) = no;
     la(di,df)=no;
     lc(di,df) $[pAct(qq,di,df)=1] = yes;

     //display lc,la,le;
     //vNewLine_R.fx  (y,lc)  = 1;
     //vNewLine.fx  (y,lc)  = 1;
     vNewLine.fx  (y,lc)  = 1;
     la(di,df) $(lc(di,df) or le(di,df))=yes;
     display lc, vNewLine.l;
     );
  );
  if ( pObjective = 1 ,                 // Saving Vecot for maximun Welfare
    pMax = smin(q, pObjFunct (q,  '6'));
    pAct(q,di,df)$(pObjFunct (q , '6') <> pMax)=no;
    qq(q)$(pObjFunct (q , '6') = pMax)= yes;
    loop (qq,
    lc(di,df) = no;
    la(di,df)=no;
    lc(di,df) $[pAct(qq,di,df)=1] = yes;
    vNewLine.fx  (y,lc)  = 1;
    la(di,df) $(lc(di,df) or le(di,df))=yes;
    );
 );

  if ( pObjective >=2 ,                 // Saving Vecot for maximun Welfare
    pMax = smin(q, pObjFunct (q,  '7'));
    pAct(q,di,df)$(pObjFunct (q , '7') <> pMax)=no;
    qq(q)$(pObjFunct (q , '7') = pMax)= yes;
    loop (qq,
    lc(di,df) = no;
    la(di,df)=no;
    lc(di,df) $[pAct(qq,di,df)=1] = yes;
    vNewLine.fx  (y,lc)  = 1;
    la(di,df) $(lc(di,df) or le(di,df))=yes;

    );
 );

//////////////////////////       GETTING RESULTS FROM OPTIMAL BILEVEL     //////////////////////////////
 solve CUADRATIC_MARKET minimizing vTotalTCost using miqcp       ;
 display vTotalTCost.l, LC;
 );

if (pRegret = 0,  // Puting ResulTs in Parameters Directly

  pTotalCost1      =0;
  pExtendedWelfare1=0;

  pTotalCost1      =
             + sum[(y,rpp(rp,pa(p)),gad(t,d)), pWeight_rp(rp) * pSlopeVarCost(t)* vProduct.l      (y,p,t,d)]
             + sum[(y,rpp(rp,pa(p)),gad(h,d)), pWeight_rp(rp) * pHydroCost      * vProduct.l      (y,p,h,d)]
             + sum[(y,rpp(rp,pa(p)),gad(h,d)), pWeight_rp(rp) * pConsumCost     * vConsump.l      (y,p,h,d)] 
             ;
 
  pExtendedWelfare1 =
             + sum[(y,rpp(rp,pa(p)),gad(t,d)), pWeight_rp(rp)* pSlopeVarCost(t)*vProduct.l(y,p,t,d)]
             + sum[(y,rpp(rp,pa(p)),gad(h,d)), pWeight_rp(rp)* pHydroCost      *vProduct.l(y,p,h,d)]
             + sum[(y,rpp(rp,pa(p)),gad(h,d)), pWeight_rp(rp)* pConsumCost     *vConsump.l(y,p,h,d)] 
             + sum((y,rpp(rp,pa(p)),       d), -(1/Slope(p,d))$(slope(p,d)) *pWeight_rp(rp)* [ pdemandnode(y,p,d) * vDemand.l(y,p,d) -  power(vDemand.l(y,p,d),2)/2])
             + sum[(y,    gcd(g,d)),  pGenInvCost (g )*vNewGen.l(y, g, d)*pMaxProd    (g,d) ]
             + sum[(y,    gcd(wn,d)), pGenInvCost (wn)*vNewGen.l(y, wn,d)*pMaxWindGen (wn,d)]
             + sum[(y,lc),      (card(y)-ord(y)+1)*pFixedCost(lc)*[vNewLine.l(y,lc) - vNewLine.l(y-1,lc)]]  
             ;

  pPrices        (y,pa(p), d) $(pMercado = 3)       =    sum((rpp(rp,p)), ( +eBalance_C_can.m (y,p,d))*1000/pWeight_rp(rp)) +eps       ;
  pPrices        (y,pa(p), d) $(pMercado <> 3)      =    sum((rpp(rp,p)), b_cLambda.l (y,p,d)*1000)+eps                                ;

  pConsumer_Surplus                                 =   -sum((y,rpp(rp,pa(p)),gad(g,d)), pWeight_rp(rp)*pSlopeVarCost(g)*vProduct.l(y,p,g,d)) -sum((y,gad(g,d)),vNewGen.l   (y,    g,d)*pGenInvCost (g))
                                                        -sum[(y,lc     ), pFixedCost (lc)*vNewLine.l(y,lc)]             ;

  pEnergy        (               th )               =    sum((y, rpp(rp,pa(p)),gad(g,d)), pWeight_rp(rp)*(vProduct.l(y,p,g,d)
                                                             +vWind.l(y,p,g,d)$wn(g)+vSolar.l(y,p,g,d)$sr(g)) $tg(g,th) )        ;
  pProduct       (y,pa(p),        th                =    sum((gad(g,d)), vProduct.l  (y,p,g,d)$tg(g,th))                         ;
  pProduct_FX    (y,pa(p),gad(g,d)  )               =    vProduct.l  (y,p,g,d)   +  vWind.l(y,p,g,d)$wn(g) +eps                  ;
  pTotalEmission (                th)                                =    sum ((y, pa(p)),pProduct(y,p,  th )*pEmissionrate(th) );
  pInstalCapT_can(y,      gcd(g,d)  )               =    vNewGen.l   (y,    g,d)      +eps                                       ;
  pInstalCapTot  (   th             )               =    sum((y,gcd(g,d)), vNewGen.l(y,g,d)*(pMaxProdGen(g)$t(g)
                                                                                +pMaxProdGen(g)$(h(g))+ pMaxWindGen(g,d) $(wn(g))  
                                                                                + pMaxSolarGen (g,d)$(sr(g)) ) $tg(g,th)   )     ;
  pFlow          (y,pa(p),        la)               =  [ vFlow.l     (y,p, la)]       +eps                                       ;
  pInstalLine    (y,              lc)               =    vNewLine.l  (y,lc)                                                      ;
  pReserve_FX    (y,pt(p),gad(hf,d) )               =    vReserve.l  (y,p,hf,d  )  +eps                                          ;
  pReserve_FX    (y,ps(p),gad(hs,d) )               =    vReserve.l  (y,p,hs,d  )  +eps                                          ;
  pConsump       (y,pa(p),gad(h,d)  )               =    vConsump.l  (y,p,h,d  )  +eps                                           ;
  pDemandl       (pa(p),        d   )               =    sum(y,  vDemand.l   (y,p,d))                                            ;
  pTotalDemand   (pa(p)             )               =    sum((y,d),  pDemandNode   (y,p,d) )                                     ;     

  if ( pMercado = 3, // If results come from enumeration
    OF_Cost    ('Op Cost   Model      [1000 M€]') = pTotalCost1 + EPS                                                            ;
    OF_ExtWel  ('Ext Welf  Model      [1000 M€]') = pExtendedWelfare1 + EPS                                                      ;
    GenCPUTime ('CPU Time  Model generation [s]') = CUADRATIC_MARKET.resGen                                                      ;
    SolCPUTime ('CPU Time  Model solution   [s]') = CUADRATIC_MARKET.resUsd                                                      ;
    NumVar     ('Number of variables           ') = CUADRATIC_MARKET.numVar                                                      ;
    NumDVar    ('Number of discrete variables  ') = CUADRATIC_MARKET.numDVar                                                     ;
    NumEqu     ('Number of equations           ') = CUADRATIC_MARKET.numEqu                                                      ;
    NumNZ      ('Number of nonzero elements    ') = CUADRATIC_MARKET.numNZ                                                       ;
    BestSol    ('Best possible solution for MIP') = CUADRATIC_MARKET.objest   ;                                                  ;
  else               // if results come form MIP
    OF_Cost    ('Op Cost   Model      [1000 M€]') = pTotalCost1 + EPS                                                            ;
    OF_ExtWel  ('Ext Welf  Model      [1000 M€]') = pExtendedWelfare1 + EPS                                                      ;
    GenCPUTime ('CPU Time  Model generation [s]') = BILEVEL_KKT_MILP.resGen                                                      ;
    SolCPUTime ('CPU Time  Model solution   [s]') = BILEVEL_KKT_MILP.resUsd                                                      ;
    NumVar     ('Number of variables           ') = BILEVEL_KKT_MILP.numVar                                                      ;
    NumDVar    ('Number of discrete variables  ') = BILEVEL_KKT_MILP.numDVar                                                     ;
    NumEqu     ('Number of equations           ') = BILEVEL_KKT_MILP.numEqu                                                      ;
    NumNZ      ('Number of nonzero elements    ') = BILEVEL_KKT_MILP.numNZ                                                       ;
    BestSol    ('Best possible solution for MIP') = BILEVEL_KKT_MILP.objest                                                      ;

  );

else              // Computin Regret


  lc(di,df) = no;
  la(di,df) = no;

  lc('21_fr','96_ie') = Yes;
  lc('26_fr','90_uk') = Yes;
  lc('31_de','79_no') = Yes;
  lc('41_pl','77_lt') = Yes;
  lc('55_it','68_gr') = Yes;
  lc('73_ee','78_lv') = Yes;

  la(di,df) $(lc(di,df) or le(di,df))=yes;

  pFixedDemand (  p,d      )=  sum(y, vdemand.l (y,p,d));

  pCost=yes;

  vNewLine.lo  (y,lc)  = 0;
  vNewLine.up  (y,lc)  = 1;

//////////////////////////                SOLVING NAIVE CP                //////////////////////////////

  pTypeComp = yes;
  pCost     = yes;
  
  solve GEPTEP_COSTMIN using MiP minimizing vTotalFCost;
  
  
  if  ( pRegret = 2,  // Storing Results to print for Naive
    pTotalCost05=0;
    pTotalCost05=
  
           +sum[(y,rpp(rp,pa(p)),gad(t,d)), pWeight_rp(rp) * pSlopeVarCost(t)* vProduct.l      (y,p,t,d)]
           +sum[(y,rpp(rp,pa(p)),gad(h,d)), pWeight_rp(rp) * pHydroCost      * vProduct.l      (y,p,h,d)]
           +sum[(y,rpp(rp,pa(p)),gad(h,d)), pWeight_rp(rp) * pConsumCost     * vConsump.l      (y,p,h,d)]
        //   +sum[(y,    gcd(g,d),gcp(g,cp)),  pGenInvCost (g )*vNewGen.l(y, g, d)*pMaxProd    (g,d) ]
        //   +sum[(y,    gcd(wn,d),gcp(wn,cp)), pGenInvCost (wn)*vNewGen.l(y, wn,d)*pMaxWindGen (wn,d)]
           +sum[(y,lc),      (card(y)-ord(y)+1)*pFixedCost(lc)*[vNewLine.l(y,lc) - vNewLine.l(y-1,lc)]]
           ;
    pPrices        (y,pa(p), d)                       =    sum((rpp(rp,p)), (eBalance_C_can.m (y,p,d)*1000 )/pWeight_rp(RP)  ) ;
    pExtendedWelfare05 =
      +                   sum[(y,rpp(rp,pa(p)),gad(t,d)), pWeight_rp(rp)* pSlopeVarCost(t)*vProduct.l(y,p,t,d)]
      +                   sum[(y,rpp(rp,pa(p)),gad(h,d)), pWeight_rp(rp)* pHydroCost      *vProduct.l(y,p,h,d)]
      +                   sum[(y,rpp(rp,pa(p)),gad(h,d)), pWeight_rp(rp)* pConsumCost     *vConsump.l(y,p,h,d)]
      +                   sum((y,rpp(rp,pa(p)),       d),-(1/Slope(p,d) * pWeight_rp(rp)* [ pdemandnode(y,p,d) *  pFixedDemand(p,d) -  power(pFixedDemand(p,d),2)/2]))
      +                   sum[(y,    gcd(g,d)),  pGenInvCost (g )*vNewGen.l(y, g, d)*pMaxProd    (g,d) ]
      +                   sum[(y,    gcd(wn,d)), pGenInvCost (wn)*vNewGen.l(y, wn,d)*pMaxWindGen (wn,d)]
      +                   sum[(y,lc),      (card(y)-ord(y)+1)*pFixedCost(lc)*[vNewLine.l(y,lc) - vNewLine.l(y-1,lc)]]  ;
  
  
      pProduct       (y,pa(p),  th )                =    sum((gad(g,d)), vProduct.l  (y,p,g,d)$tg(g,th)) ;
      pProduct_FX    (y,pa(p),gad(g,d))             =    vProduct.l  (y,p,g,d)   +  vWind.l(y,p,g,d)$wn(g) ;
      pInstalCapT_can(y,    g,d      )              =    vNewGen.l   (y,    g,d)                   ;
      pFlow    (y,pa(p), la)                        =  [ vFlow.l     (y,p, la)]                          ;
      pInstalLine    (y,         lc)                =    vNewLine.l  (y,lc)                              ;
      pReserve_FX(y,pt(p),gad(hf,d) )               =    vReserve.l  (y,p,hf,d  )  +eps                  ;
      pReserve_FX(y,ps(p),gad(hs,d) )               =    vReserve.l  (y,p,hs,d  )  +eps                  ;
      pConsump   (y,pa(p),gad(h,d) )                =    vConsump.l  (y,p,h,d  )  +eps                   ;
      pDemandl (pa(p),d)                            =  sum(y,  pFixedDemand(p,d))                             ;
      pPrices        (y,pa(p), d)                   =  pPrices(y,p, d) *1000;
  
      OF_Cost    ('Op Cost   Model      [1000 M€]') = pTotalCost05 + EPS                                                                ;
      OF_ExtWel  ('Ext Welf  Model      [1000 M€]') = pExtendedWelfare05 + EPS                                                           ;
      GenCPUTime ('CPU Time  Model generation [s]') = CUADRATIC_MARKET.resGen                                                          ;
      SolCPUTime ('CPU Time  Model solution   [s]') = CUADRATIC_MARKET.resUsd                                                          ;
      NumVar     ('Number of variables           ') = CUADRATIC_MARKET.numVar                                                          ;
      NumDVar    ('Number of discrete variables  ') = CUADRATIC_MARKET.numDVar                                                         ;
      NumEqu     ('Number of equations           ') = CUADRATIC_MARKET.numEqu                                                          ;
      NumNZ      ('Number of nonzero elements    ') = CUADRATIC_MARKET.numNZ                                                           ;
      BestSol    ('Best possible solution for MIP') = CUADRATIC_MARKET.objest   ;                                                       ;
  );
  
  pOptimal= COST_MINIMISATION.modelstat ;
  vdemand.lo (y,p,d)=0;
  vdemand.up (y,p,d)=inf;
  display  vNewLine.l, lc,pFixedDemand;
  vNewLine.fx  (y,lc)  = vNewline.l (y,lc);
  pTypeComp = no;
  pCost=no;
  )
 ;

if (pRegret =1,
//////////////////////////               SOLVING ACTUAL CP                //////////////////////////////
  solve CUADRATIC_MARKET minimizing vTotalTCost using miqcp       ;
  pTotalCost2=0;
  pTotalCost2=

         +sum[(y,rpp(rp,pa(p)),gad(t,d)), pWeight_rp(rp) * pSlopeVarCost(t)* vProduct.l      (y,p,t,d)]
         +sum[(y,rpp(rp,pa(p)),gad(h,d)), pWeight_rp(rp) * pHydroCost      * vProduct.l      (y,p,h,d)]
         +sum[(y,rpp(rp,pa(p)),gad(h,d)), pWeight_rp(rp) * pConsumCost     * vConsump.l      (y,p,h,d)]
         ;
  pPrices        (y,pa(p), d)                       =    sum((rpp(rp,p)), ( (pDemandNode(y,p,d)-vDemand.l(y,p,d))/Slope(p,d) )*1000 ) ;
  pExtendedWelfare2 =
    +                   sum[(y,rpp(rp,pa(p)),gad(t,d)), pWeight_rp(rp)* pSlopeVarCost(t)*vProduct.l(y,p,t,d)]
    +                   sum[(y,rpp(rp,pa(p)),gad(h,d)), pWeight_rp(rp)* pHydroCost      *vProduct.l(y,p,h,d)]
    +                   sum[(y,rpp(rp,pa(p)),gad(h,d)), pWeight_rp(rp)* pConsumCost     *vConsump.l(y,p,h,d)]
    +                   sum((y,rpp(rp,pa(p)),       d),-(1/Slope(p,d) * pWeight_rp(rp)* [ pDemandNode(y,p,d) * vDemand.l(y,p,d) -  power(vDemand.l(y,p,d),2)/2]))
    +                   sum[(y,    gcd(g,d)),  pGenInvCost (g )*vNewGen.l(y, g, d)*pMaxProd    (g,d) ]
    +                   sum[(y,    gcd(wn,d)), pGenInvCost (wn)*vNewGen.l(y, wn,d)*pMaxWindGen (wn,d)]
    +                   sum[(y,lc),      (card(y)-ord(y)+1)*pFixedCost(lc)*[vNewLine.l(y,lc) - vNewLine.l(y-1,lc)]]  ;
            ;

    pProduct       (y,pa(p),  th )                =    sum((gad(g,d)), vProduct.l  (y,p,g,d)$tg(g,th))                                ;
    pProduct_FX    (y,pa(p),gad(g,d))             =    vProduct.l  (y,p,g,d)   +  vWind.l(y,p,g,d)$wn(g)                              ;
    pInstalCapT_can(y,    gcd(g,d)    )           =    vNewGen.l   (y,    g,d)        +eps                                            ;
    pInstalCapTot    (   th    )                  =   sum((y,gcd(g,d)), vNewGen.l(y,g,d)*(
                                                                  + pMaxProdGen(g  )$ t (g)  +pMaxProdGen  (g  )$(h (g))
                                                                  + pMaxWindGen(g,d)$(wn(g)) +pMaxSolarGen (g,d)$(sr(g))) $tg(g,th))  ;
    pEnergy        (                     th )         =    sum((y, rpp(rp,pa(p)),gad(g,d)), pWeight_rp(rp)*(vProduct.l  (y,p,g,d)
                                                                  +vWind.l(y,p,g,d)$wn(g)+vSolar.l(y,p,g,d)$sr(g)) $tg(g,th) )        ;
    pTotalEmission ( th )                             =    sum ((y, pa(p)), pProduct(y,p,  th )*pEmissionrate(th) )      ;
    pFlow    (y,pa(p), la)                        =  [ vFlow.l     (y,p, la)]                                                         ;
    pInstalLine    (y,         lc)                =    vNewLine.l  (y,lc)                                                             ;
    pReserve_FX(y,pt(p),gad(hf,d) )               =    vReserve.l  (y,p,hf,d  )  +eps                                                 ;
    pReserve_FX(y,ps(p),gad(hs,d) )               =    vReserve.l  (y,p,hs,d  )  +eps                                                 ;
    pConsump   (y,pa(p),gad(h,d) )                =    vConsump.l  (y,p,h ,d  )  +eps                                                 ;
    pDemandl (pa(p),d)                            =  sum(y,  pDemandNode   (y,p,d) )   ;                                              ;
    pTotalDemand (pa(p))                            =  sum((y,d),  pDemandNode   (y,p,d) )   ;      


    OF_Cost    ('Op Cost   Model      [1000 M€]') = pTotalCost2 + EPS                                                                 ;
    OF_ExtWel  ('Ext Welf  Model      [1000 M€]') = pExtendedWelfare2 + EPS                                                           ;
    GenCPUTime ('CPU Time  Model generation [s]') = CUADRATIC_MARKET.resGen                                                           ;
    SolCPUTime ('CPU Time  Model solution   [s]') = CUADRATIC_MARKET.resUsd                                                           ;
    NumVar     ('Number of variables           ') = CUADRATIC_MARKET.numVar                                                           ;
    NumDVar    ('Number of discrete variables  ') = CUADRATIC_MARKET.numDVar                                                          ;
    NumEqu     ('Number of equations           ') = CUADRATIC_MARKET.numEqu                                                           ;
    NumNZ      ('Number of nonzero elements    ') = CUADRATIC_MARKET.numNZ                                                            ;
    BestSol    ('Best possible solution for MIP') = CUADRATIC_MARKET.objest   ;                                                       ;



 // DISPLAY   COST_REGRET,WELF_REGRET, pTotalCost2,pTotalCost1,pExtendedWelfare, vNewline.l ;
);


*Puting data into Excel

put TMP putclose
'par=pInstalCapT_can rdim=1 rng=GenInv!a1:zz4'        / 'par=pFlow          rdim=4 rng=Flow!a1'     / 'par=pTheta      rdim=3 rng=Angle!a1'   /
'par=pProduct        rdim=2 rng=Output!a1'            / 'par=pInstalLine    rdim=1 rng=LineInv!a1'  /
'par=pReserve_FX     rdim=2 rng=WtrReserve!a1'        / 'par=OF_Cost        rdim=1 rng=Cost!a1'     / 'par=GenCPUTime rdim=1 rng=Cost!a2'     /
'par=NumVar          rdim=1 rng=Cost!a4'              / 'par=NumDVar        rdim=1 rng=Cost!a5'     / 'par=NumEqu      rdim=1 rng=Cost!a6'    /
'par=BestSol         rdim=1 rng=Cost!a7'              / 'par=pBenefits      rdim=1 rng=Prices!a1'   / 'par=SolCPUTime  rdim=1 rng=Cost!a3'    /
'par=pPrices         rdim=3 rng=Prices!a25:d1000000'  / 'par=pEnergy        rdim=1 rng=Energy!a1:d15'/
'par=pInstalCapTot   rdim=1 rng=GenInv!b6:d20'        / 'par=pTotalEmission  rdim=1 rng=Energy!a25:a60'/  'par=pTotalDemand rdim=1 rng=Energy!h1:l10000'/
//'par=pHourlySolpe    rdim=3 rng=Prices!f25'   /
execute_unload   '%gams.user1%.gdx' pInstalCapT_can pFlow pTheta  pProduct pInstalLine pReserve_FX pEnergy pInstalCapTot pTotalEmission pTotalDemand
                                    OF_Cost GenCPUTime SolCPUTime NumVar NumDVar NumEqu NumNZ BestSol pBenefits  pPrices  pMarginalCosts
execute          'gdxxrw "%gams.user1%".gdx SQ=n EpsOut=0 O="%gams.user1%".xlsx @"%gams.user1%".txt'
execute          'del    "%gams.user1%".gdx '






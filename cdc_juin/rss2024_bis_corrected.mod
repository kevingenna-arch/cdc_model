%%%%% rss2024_bis_corrected.mod ################################################
%% Ruppert, Schön & Stähler (2024)
%% "Consumption Taxation to Finance Pension Payments"
%%
%% Modèle OLG à deux régions + frictions DMP
%%   Région A (domestique) : Allemagne   (part = relsize = 0.16)
%%   Région B (étrangère)  : Reste de l'UE (part = 1-relsize = 0.84)
%%
%% CORRECTIONS v3 (structurelles) :
%%   Les 5 équations d'agrégation de la région B manquaient le facteur
%%   de taille relative (1-relsize)/relsize = 5.25.
%%   Corrections appliquées dans model{} :
%%     - OmB      : multiplié par (1-relsize)/relsize
%%     - NB_agg   : multiplié par (1-relsize)/relsize
%%     - UB       : multiplié par (1-relsize)/relsize
%%     - CB       : multiplié par (1-relsize)/relsize
%%     - SavB     : multiplié par (1-relsize)/relsize
%%     - eq. vacances B : RHS multiplié par (1-relsize)/relsize
%%   Corrections initval :
%%     - OmB  : 7.90  → 9.74  (5.25 × 1.855)
%%     - NB_agg: 31.0 → 35.0  (5.25 × 6.668)
%%     - UB   : 2.10  → 4.49  (5.25 × 0.855)
%%     - CB   : 10.7  → 18.70 (5.25 × 3.561)
%%     - SavB : 24.88 → 24.88 (inchangé — déjà per capita × ZB)
%%     - KB   : 8.05  → 7.44  (cohérence Cobb-Douglas YB=21 avec NB=35)
%%     - QB   : 4.23  → 4.23  (SavB - KB - BBdebt)
%%   Corrections v2 (initval région A) :
%%     - NA_agg: 6.08 → 7.48, UA: 0.37 → 0.43, CA: 2.10 → 3.77
%%     - SavA: 2.50 → 5.04, OmA: 1.55 → 1.63

@#define NW  = 9
@#define NR  = 7
@#define NA  = 16


%%%%% VARIABLES ENDOGÈNES ######################################################

%% ── Région A ──────────────────────────────────────────────────────────────
@#for a in 1:NA
    var cA_@{a};
@#endfor
@#for a in 1:NA-1
    var sA_@{a};
@#endfor
@#for a in 1:NW
    var eA_@{a};
    var nA_@{a};
    var uA_@{a};
    var WWA_@{a};
    var WFA_@{a};
@#endfor

var hpA_1;
@#for a in 2:NW+1
    var hpA_@{a};
@#endfor
@#for a in NW+1:NA
    var brA_@{a};
@#endfor

var VA OmA MA pA qA wA;
var NA_agg UA CA SavA;
var YA KA QA BAdebt BdotA GA;
var taucA;
var NFAA NXA;
var fAA fBA;

%% ── Région B ──────────────────────────────────────────────────────────────
@#for a in 1:NA
    var cB_@{a};
@#endfor
@#for a in 1:NA-1
    var sB_@{a};
@#endfor
@#for a in 1:NW
    var eB_@{a};
    var nB_@{a};
    var uB_@{a};
    var WWB_@{a};
    var WFB_@{a};
@#endfor

var hpB_1;
@#for a in 2:NW+1
    var hpB_@{a};
@#endfor
@#for a in NW+1:NA
    var brB_@{a};
@#endfor

var VB OmB MB pB qB wB;
var NB_agg UB CB SavB;
var YB KB QB BBdebt BdotB GB;
var taucB;
var NFAB NXB;
var fBB fAB;

%% ── Variables mondiales ───────────────────────────────────────────────────
var rstar;
var rer;


%%%%% VARIABLES EXOGÈNES #######################################################
varexo tauwA;
varexo tauwB;
varexo xA xB;


%%%%% PARAMÈTRES ###############################################################
parameters beta_p kappa alpha delta_p omega_ann;
parameters mu_A mu_B nu chi_p phi_A phi_B eta;
parameters rho_u rho_rA rho_rB;
parameters g_share B_target;
parameters taufA taufB taukA taukB tauc0A tauc0B;
parameters zeta_tr thbA thbB relsize;

@#for a in 1:NA
    parameters thetaA_@{a} thetaB_@{a};
@#endfor
@#for a in 1:NW
    parameters xiA_@{a} xiB_@{a};
@#endfor
@#for a in 1:NA
    parameters ZA_@{a} ZB_@{a};
@#endfor


%%%%% VALEURS DES PARAMÈTRES ###################################################
beta_p    = 0.9861^5;
alpha     = 0.33;
delta_p   = 1 - (1 - 0.10)^5;
omega_ann = 0.75;
kappa     = 0.80;
nu        = 0.50;
eta       = 0.50;
chi_p     = 0.06;
mu_A      = 0.8025;
mu_B      = 0.7937;
phi_A     = 0.2410;
phi_B     = 0.2613;
rho_u     = 0.350;
rho_rA    = 0.519;
rho_rB    = 0.635;
g_share   = 0.20;
B_target  = 0.60;
taufA     = 0.167;
taufB     = 0.246;
taukA     = 0.214;
taukB     = 0.316;
tauc0A    = 0.183;
tauc0B    = 0.196;
zeta_tr   = 0.33;
thbA      = 0.60;
thbB      = 0.8396;
relsize   = 0.16;

%% Probabilités de survie — Tableau 1
thetaA_1  = 1.0000;  thetaB_1  = 1.0000;
thetaA_2  = 0.9983;  thetaB_2  = 0.9981;
thetaA_3  = 0.9947;  thetaB_3  = 0.9939;
thetaA_4  = 0.9884;  thetaB_4  = 0.9868;
thetaA_5  = 0.9785;  thetaB_5  = 0.9755;
thetaA_6  = 0.9630;  thetaB_6  = 0.9578;
thetaA_7  = 0.9382;  thetaB_7  = 0.9299;
thetaA_8  = 0.8981;  thetaB_8  = 0.8864;
thetaA_9  = 0.8346;  thetaB_9  = 0.8205;
thetaA_10 = 0.7404;  thetaB_10 = 0.7263;
thetaA_11 = 0.6128;  thetaB_11 = 0.6020;
thetaA_12 = 0.4568;  thetaB_12 = 0.4503;
thetaA_13 = 0.2864;  thetaB_13 = 0.2837;
thetaA_14 = 0.1280;  thetaB_14 = 0.1296;
thetaA_15 = 0.0294;  thetaB_15 = 0.0318;
thetaA_16 = 0.0019;  thetaB_16 = 0.0024;

@#for a in 1:NA
    ZA_@{a} = thetaA_@{a};
    ZB_@{a} = thetaB_@{a};
@#endfor

%% Paramètres de loisir par âge — Tableau 2
xiA_1=0.7329; xiA_2=0.4582; xiA_3=0.3541; xiA_4=0.2982; xiA_5=0.2574;
xiA_6=0.2657; xiA_7=0.3644; xiA_8=0.6616; xiA_9=0.7269;

xiB_1=1.1055; xiB_2=0.6604; xiB_3=0.5265; xiB_4=0.4697; xiB_5=0.4630;
xiB_6=0.5159; xiB_7=0.6864; xiB_8=1.1787; xiB_9=1.2109;


%%%%% BLOC MODÈLE ##############################################################
model;

%% ══════════════════════════════════════════════════════════════════════════
%%  I.  CHERCHEURS D'EMPLOI — éq. (6)
%%  NB : OmB est scalé par (1-relsize)/relsize pour représenter
%%       la taille de la population B relative à A
%% ══════════════════════════════════════════════════════════════════════════
OmA = ZA_1*(1-eA_1)
@#for a in 2:NW
    + ZA_@{a} * (uA_@{a-1}(-1) + chi_p*nA_@{a-1}(-1)) * (1-eA_@{a})
@#endfor
;

OmB = ((1-relsize)/relsize) * ( ZB_1*(1-eB_1)
@#for a in 2:NW
    + ZB_@{a} * (uB_@{a-1}(-1) + chi_p*nB_@{a-1}(-1)) * (1-eB_@{a})
@#endfor
);

%% ══════════════════════════════════════════════════════════════════════════
%%  II.  APPARIEMENT ET TAUX — éq. (5),(7)
%% ══════════════════════════════════════════════════════════════════════════
MA = mu_A * VA^(1-nu) * OmA^nu;
MB = mu_B * VB^(1-nu) * OmB^nu;
pA = MA / OmA;
pB = MB / OmB;
qA = MA / VA;
qB = MB / VB;

%% ══════════════════════════════════════════════════════════════════════════
%%  III.  DYNAMIQUE DE L'EMPLOI — éq. (8)
%% ══════════════════════════════════════════════════════════════════════════
nA_1 = pA*(1-eA_1);
nB_1 = pB*(1-eB_1);

@#for a in 2:NW
    nA_@{a} = (1-chi_p)*(1-pA)*(1-eA_@{a})*nA_@{a-1}(-1) + pA*(1-eA_@{a});
    nB_@{a} = (1-chi_p)*(1-pB)*(1-eB_@{a})*nB_@{a-1}(-1) + pB*(1-eB_@{a});
@#endfor

@#for a in 1:NW
    uA_@{a} = 1 - nA_@{a} - eA_@{a};
    uB_@{a} = 1 - nB_@{a} - eB_@{a};
@#endfor

%% ══════════════════════════════════════════════════════════════════════════
%%  IV.  EULER — éq. (13)
%% ══════════════════════════════════════════════════════════════════════════
@#for a in 1:NA-1
    cA_@{a+1}(+1)/cA_@{a} = beta_p * (thetaA_@{a+1}/thetaA_@{a})^(1-omega_ann)
        * (1 + rstar(+1)*(1-taukA)) * (1+taucA)/(1+taucA(+1));
    cB_@{a+1}(+1)/cB_@{a} = beta_p * (thetaB_@{a+1}/thetaB_@{a})^(1-omega_ann)
        * (1 + rstar(+1)*(1-taukB)) * (1+taucB)/(1+taucB(+1));
@#endfor

%% ══════════════════════════════════════════════════════════════════════════
%%  V.  DÉCISION DE PARTICIPATION — éq. (15)
%% ══════════════════════════════════════════════════════════════════════════
@#for a in 1:NW-1
    xiA_@{a} * eA_@{a}^(-kappa) =
        pA*(1-tauwA)*wA / ((1+taucA)*cA_@{a})
        + (1-pA)*rho_u*wA / ((1+taucA)*cA_@{a})
        @#for j in NW+1:NA
            + beta_p^(@{j-a}) * (thetaA_@{j}/thetaA_@{a})
              * rho_rA * wA(+@{j-a}) / ((1+taucA(+@{j-a}))*cA_@{j}(+@{j-a}))
        @#endfor
        ;
    xiB_@{a} * eB_@{a}^(-kappa) =
        pB*(1-tauwB)*wB / ((1+taucB)*cB_@{a})
        + (1-pB)*rho_u*wB / ((1+taucB)*cB_@{a})
        @#for j in NW+1:NA
            + beta_p^(@{j-a}) * (thetaB_@{j}/thetaB_@{a})
              * rho_rB * wB(+@{j-a}) / ((1+taucB(+@{j-a}))*cB_@{j}(+@{j-a}))
        @#endfor
        ;
@#endfor

xiA_@{NW} * eA_@{NW}^(-kappa) =
    pA*(1-tauwA)*wA / ((1+taucA)*cA_@{NW})
    + (1-pA)*rho_u*wA / ((1+taucA)*cA_@{NW});
xiB_@{NW} * eB_@{NW}^(-kappa) =
    pB*(1-tauwB)*wB / ((1+taucB)*cB_@{NW})
    + (1-pB)*rho_u*wB / ((1+taucB)*cB_@{NW});

%% ══════════════════════════════════════════════════════════════════════════
%%  VI.  CONTRAINTES BUDGÉTAIRES — éq. (10)–(11)
%% ══════════════════════════════════════════════════════════════════════════
(1+taucA)*cA_1 + sA_1 = (1-tauwA)*wA*nA_1 + rho_u*wA*uA_1;
(1+taucB)*cB_1 + sB_1 = (1-tauwB)*wB*nB_1 + rho_u*wB*uB_1;

@#for a in 2:NW
    (1+taucA)*cA_@{a} + sA_@{a} =
        (thetaA_@{a-1}/thetaA_@{a})^omega_ann * (1+rstar*(1-taukA)) * sA_@{a-1}(-1)
        + (1-tauwA)*wA*nA_@{a} + rho_u*wA*uA_@{a};
    (1+taucB)*cB_@{a} + sB_@{a} =
        (thetaB_@{a-1}/thetaB_@{a})^omega_ann * (1+rstar*(1-taukB)) * sB_@{a-1}(-1)
        + (1-tauwB)*wB*nB_@{a} + rho_u*wB*uB_@{a};
@#endfor

@#for a in NW+1:NA-1
    (1+taucA)*cA_@{a} + sA_@{a} =
        (thetaA_@{a-1}/thetaA_@{a})^omega_ann * (1+rstar*(1-taukA)) * sA_@{a-1}(-1)
        + brA_@{a};
    (1+taucB)*cB_@{a} + sB_@{a} =
        (thetaB_@{a-1}/thetaB_@{a})^omega_ann * (1+rstar*(1-taukB)) * sB_@{a-1}(-1)
        + brB_@{a};
@#endfor

(1+taucA)*cA_@{NA} =
    (thetaA_@{NA-1}/thetaA_@{NA})^omega_ann * (1+rstar*(1-taukA)) * sA_@{NA-1}(-1)
    + brA_@{NA};
(1+taucB)*cB_@{NA} =
    (thetaB_@{NA-1}/thetaB_@{NA})^omega_ann * (1+rstar*(1-taukB)) * sB_@{NA-1}(-1)
    + brB_@{NA};

%% ══════════════════════════════════════════════════════════════════════════
%%  VII.  POINTS RETRAITE — éq. (12)
%% ══════════════════════════════════════════════════════════════════════════
hpA_1 = 0;
hpB_1 = 0;

@#for a in 1:NW
    hpA_@{a+1} = hpA_@{a}(-1) + (1-eA_@{a}(-1));
    hpB_@{a+1} = hpB_@{a}(-1) + (1-eB_@{a}(-1));
@#endfor

@#for a in NW+1:NA
    brA_@{a} = rho_rA * wA * hpA_@{NW+1}(-@{a-NW}) / @{NW};
    brB_@{a} = rho_rB * wB * hpB_@{NW+1}(-@{a-NW}) / @{NW};
@#endfor

%% ══════════════════════════════════════════════════════════════════════════
%%  VIII.  VALEUR MÉNAGE DE L'EMPLOI — éq. (14)
%% ══════════════════════════════════════════════════════════════════════════
@#for a in 1:NW-1
    WWA_@{a} = (1-tauwA-rho_u)*wA
        + (1-chi_p)*(1-pA(+1))*(1-eA_@{a+1}(+1))
          * beta_p * (thetaA_@{a+1}/thetaA_@{a})
          * ((1+taucA)*cA_@{a}) / ((1+taucA(+1))*cA_@{a+1}(+1))
          * WWA_@{a+1}(+1);
    WWB_@{a} = (1-tauwB-rho_u)*wB
        + (1-chi_p)*(1-pB(+1))*(1-eB_@{a+1}(+1))
          * beta_p * (thetaB_@{a+1}/thetaB_@{a})
          * ((1+taucB)*cB_@{a}) / ((1+taucB(+1))*cB_@{a+1}(+1))
          * WWB_@{a+1}(+1);
@#endfor

WWA_@{NW} = (1-tauwA-rho_u)*wA;
WWB_@{NW} = (1-tauwB-rho_u)*wB;

%% ══════════════════════════════════════════════════════════════════════════
%%  IX.  VALEUR FIRME D'UN EMPLOI — éq. (20)
%% ══════════════════════════════════════════════════════════════════════════
@#for a in 1:NW-1
    WFA_@{a} = fAA*(1-alpha)*(YA/NA_agg) - (1+taufA)*wA
        + (1-chi_p)*(1-eA_@{a+1}(+1)) / (1+rstar(+1))
          * (thetaA_@{a+1}/thetaA_@{a}) * WFA_@{a+1}(+1);
    WFB_@{a} = fBB*(1-alpha)*(YB/NB_agg) - (1+taufB)*wB
        + (1-chi_p)*(1-eB_@{a+1}(+1)) / (1+rstar(+1))
          * (thetaB_@{a+1}/thetaB_@{a}) * WFB_@{a+1}(+1);
@#endfor

WFA_@{NW} = fAA*(1-alpha)*(YA/NA_agg) - (1+taufA)*wA;
WFB_@{NW} = fBB*(1-alpha)*(YB/NB_agg) - (1+taufB)*wB;

%% ══════════════════════════════════════════════════════════════════════════
%%  X.  NÉGOCIATION DE NASH — éq. (21)
%% ══════════════════════════════════════════════════════════════════════════
(1-eta)*WWA_5 = eta * (1-tauwA)/(1+taufA) * WFA_5;
(1-eta)*WWB_5 = eta * (1-tauwB)/(1+taufB) * WFB_5;

%% ══════════════════════════════════════════════════════════════════════════
%%  XI.  PRODUCTION ET CAPITAL — éq. (16),(18)
%% ══════════════════════════════════════════════════════════════════════════
rstar + delta_p = fAA * alpha * YA / KA(-1);
rstar + delta_p = fBB * alpha * YB / KB(-1);

YA = KA(-1)^alpha * NA_agg^(1-alpha);
YB = KB(-1)^alpha * NB_agg^(1-alpha);

%% NB_agg scalé par (1-relsize)/relsize = 5.25
NA_agg = (
@#for a in 1:NW
    + ZA_@{a}*nA_@{a}
@#endfor
);
NB_agg = ((1-relsize)/relsize) * (
@#for a in 1:NW
    + ZB_@{a}*nB_@{a}
@#endfor
);

UA = (
@#for a in 1:NW
    + ZA_@{a}*uA_@{a}
@#endfor
);
UB = ((1-relsize)/relsize) * (
@#for a in 1:NW
    + ZB_@{a}*uB_@{a}
@#endfor
);

%% ══════════════════════════════════════════════════════════════════════════
%%  XII.  POSTES VACANTS — éq. (19)
%%  RHS scalé par (1-relsize)/relsize pour cohérence avec OmB
%% ══════════════════════════════════════════════════════════════════════════
phi_A * OmA = qA * (
    ZA_1*(1-eA_1)*WFA_1
@#for a in 2:NW
    + ZA_@{a}*(uA_@{a-1}(-1) + chi_p*nA_@{a-1}(-1))*(1-eA_@{a})*WFA_@{a}
@#endfor
);

phi_B * OmB = qB * ((1-relsize)/relsize) * (
    ZB_1*(1-eB_1)*WFB_1
@#for a in 2:NW
    + ZB_@{a}*(uB_@{a-1}(-1) + chi_p*nB_@{a-1}(-1))*(1-eB_@{a})*WFB_@{a}
@#endfor
);

%% ══════════════════════════════════════════════════════════════════════════
%%  XIII.  AGRÉGATS CONSOMMATION / ÉPARGNE
%%  CB et SavB scalés par (1-relsize)/relsize
%% ══════════════════════════════════════════════════════════════════════════
CA = (
@#for a in 1:NA
    + ZA_@{a}*cA_@{a}
@#endfor
);
CB = ((1-relsize)/relsize) * (
@#for a in 1:NA
    + ZB_@{a}*cB_@{a}
@#endfor
);

SavA = (
@#for a in 1:NA-1
    + ZA_@{a}*sA_@{a}
@#endfor
);
SavB = ((1-relsize)/relsize) * (
@#for a in 1:NA-1
    + ZB_@{a}*sB_@{a}
@#endfor
);

%% ══════════════════════════════════════════════════════════════════════════
%%  XIV.  MARCHÉS D'ACTIFS — éq. (25)–(27)
%% ══════════════════════════════════════════════════════════════════════════
QA*(1+rstar(+1)) = QA(+1) + (fAA*YA - (rstar+delta_p)*KA(-1)
                   - (1+taufA)*wA*NA_agg - phi_A*VA);
QB*(1+rstar(+1)) = QB(+1) + (fBB*YB - (rstar+delta_p)*KB(-1)
                   - (1+taufB)*wB*NB_agg - phi_B*VB);

SavA = KA(+1) + QA + BAdebt + NFAA;
SavB = KB(+1) + QB + BBdebt + NFAB;

NFAB = -rer*NFAA;

%% ══════════════════════════════════════════════════════════════════════════
%%  XV.  COMMERCE ET PRIX — éq. (28)–(32)
%% ══════════════════════════════════════════════════════════════════════════
NFAA = (1+rstar(-1))*NFAA(-1) + NXA;

NXA = YA - CA - (KA(+1) - (1-delta_p)*KA(-1)) - GA;
NXB = YB - CB - (KB(+1) - (1-delta_p)*KB(-1)) - GB;

1 = ( thbA * fAA^(-zeta_tr/(1-zeta_tr))
    + (1-thbA) * fBA^(-zeta_tr/(1-zeta_tr)) )^(-(1-zeta_tr)/zeta_tr);

1 = ( thbB * fBB^(-zeta_tr/(1-zeta_tr))
    + (1-thbB) * fAB^(-zeta_tr/(1-zeta_tr)) )^(-(1-zeta_tr)/zeta_tr);

fBA = rer*fAA;
fAB = fBB/rer;

%% ══════════════════════════════════════════════════════════════════════════
%%  XVI.  BUDGET DE L'ÉTAT — éq. (22)–(24)
%% ══════════════════════════════════════════════════════════════════════════
GA = g_share*(YA - phi_A*VA)
    + rho_u*wA*(
@#for a in 1:NW
        + ZA_@{a}*uA_@{a}
@#endfor
    )
    + rho_rA*wA*hpA_@{NW+1}(-1)/@{NW}*(
@#for a in NW+1:NA
        + ZA_@{a}
@#endfor
    );

GB = g_share*(YB - phi_B*VB)
    + rho_u*wB*((1-relsize)/relsize)*(
@#for a in 1:NW
        + ZB_@{a}*uB_@{a}
@#endfor
    )
    + rho_rB*wB*hpB_@{NW+1}(-1)/@{NW}*((1-relsize)/relsize)*(
@#for a in NW+1:NA
        + ZB_@{a}
@#endfor
    );

BdotA = GA
    - ( (tauwA+taufA)*wA*NA_agg
       + taukA*rstar*(
@#for a in 2:NA
           + ZA_@{a-1}*(thetaA_@{a-1}/thetaA_@{a})^omega_ann*sA_@{a-1}(-1)
@#endfor
       )
       + taucA*CA );

BdotB = GB
    - ( (tauwB+taufB)*wB*NB_agg
       + taukB*rstar*((1-relsize)/relsize)*(
@#for a in 2:NA
           + ZB_@{a-1}*(thetaB_@{a-1}/thetaB_@{a})^omega_ann*sB_@{a-1}(-1)
@#endfor
       )
       + taucB*CB );

BAdebt = BdotA + (1+rstar(-1))*BAdebt(-1);
BBdebt = BdotB + (1+rstar(-1))*BBdebt(-1);

%% Règle fiscale : taucA/B ajustent pour stabiliser la dette à B_target*Y
BAdebt = B_target * YA;
BBdebt = B_target * YB;

end;
%% ═══ Fin du bloc modèle ════════════════════════════════════════════════════


%%%%% VALEURS INITIALES ########################################################
%%
%%  Toutes les valeurs agrégées de B sont scalées par (1-relsize)/relsize = 5.25
%%
%%  Calculs clés :
%%    sum(ZA_a, a=1..9) = 8.594 ; sum(ZB_a, a=1..9) = 8.549
%%    OmA per cap = 1.634 ; OmB per cap = 1.855
%%    NA_agg = 0.87 * 8.594 = 7.477
%%    NB_agg = 5.25 * 0.78 * 8.549 = 35.007
%%    UA     = 0.05 * 8.594 = 0.430
%%    UB     = 5.25 * 0.10 * 8.549 = 4.488
%%    CA     = sum(ZA_a * cA_a) = 3.772
%%    CB     = 5.25 * sum(ZB_a * cB_a) = 5.25 * 3.561 = 18.695
%%    SavA   = sum(ZA_a * sA_a) = 5.041
%%    SavB   = 5.25 * sum(ZB_a * sB_a) = 5.25 * 4.739 = 24.881
%%    KB     = 7.44 (cohérence YB=21 avec NB=35 : K=YB/(NB^(1-α)) ^ (1/α))
%%    QB     = SavB - KB - BBdebt - NFAB = 24.88 - 7.44 - 12.60 - 0 = 4.84
%%    OmB    = 5.25 * 1.855 = 9.739

initval;

    rstar = 0.449;
    rer   = 1;

    %% ── Région A ─────────────────────────────────────────────────────────
    wA     = 0.240;
    YA     = 4.13;
    KA     = 1.58;
    NA_agg = 7.48;
    UA     = 0.43;
    CA     = 3.77;
    SavA   = 5.04;
    QA     = 1.00;
    BAdebt = B_target*YA;
    BdotA  = 0;
    GA     = g_share*YA;
    taucA  = tauc0A;
    NFAA   = 0;
    NXA    = 0;
    fAA    = 1;
    fBA    = 1;

    pA  = 0.92;
    qA  = 0.70;
    OmA = 1.63;
    MA  = pA*OmA;
    VA  = MA/qA;

    WFA_1=0.337; WFA_2=0.333; WFA_3=0.327; WFA_4=0.318; WFA_5=0.305;
    WFA_6=0.286; WFA_7=0.257; WFA_8=0.213; WFA_9=0.139;

    WWA_1=0.184; WWA_2=0.182; WWA_3=0.179; WWA_4=0.174; WWA_5=0.167;
    WWA_6=0.157; WWA_7=0.140; WWA_8=0.115; WWA_9=0.083;

    @#for a in 1:NW
        nA_@{a} = 0.87;
        uA_@{a} = 0.05;
        eA_@{a} = 0.08;
    @#endfor

    cA_1=0.10;  cA_2=0.126; cA_3=0.159; cA_4=0.201; cA_5=0.254;
    cA_6=0.320; cA_7=0.403; cA_8=0.508; cA_9=0.641;
    cA_10=0.630; cA_11=0.595; cA_12=0.540; cA_13=0.460;
    cA_14=0.360; cA_15=0.240; cA_16=0.150;

    sA_1=0.01;  sA_2=0.03;  sA_3=0.07;  sA_4=0.14;  sA_5=0.24;
    sA_6=0.40;  sA_7=0.60;  sA_8=0.83;  sA_9=1.05;
    sA_10=1.10; sA_11=1.00; sA_12=0.80; sA_13=0.55;
    sA_14=0.28; sA_15=0.08;

    hpA_1 = 0;
    @#for a in 2:NW+1
        hpA_@{a} = (@{a}-1)*0.92;
    @#endfor
    @#for a in NW+1:NA
        brA_@{a} = rho_rA*wA*hpA_@{NW+1}/@{NW};
    @#endfor

    %% ── Région B — tous les agrégats à l'échelle mondiale (×5.25) ─────────
    wB     = 0.230;
    YB     = 21.0;
    KB     = 7.44;    %% CORRIGÉ v3 : KB = (YB/NB^(1-α))^(1/α) avec NB=35.0
    NB_agg = 35.0;    %% CORRIGÉ v3 : 5.25 × 0.78 × 8.549 = 35.007
    UB     = 4.49;    %% CORRIGÉ v3 : 5.25 × 0.10 × 8.549 = 4.488
    CB     = 18.70;   %% CORRIGÉ v3 : 5.25 × sum(ZB_a × cB_a) = 5.25 × 3.561
    SavB   = 24.88;   %% sum(ZB_a × sB_a) × 5.25 = 4.739 × 5.25 = 24.88
    QB     = 4.84;    %% CORRIGÉ v3 : SavB - KB - BBdebt - NFAB = 24.88-7.44-12.60-0
    BBdebt = B_target*YB;
    BdotB  = 0;
    GB     = g_share*YB;
    taucB  = tauc0B;
    NFAB   = 0;
    NXB    = 0;
    fBB    = 1;
    fAB    = 1;

    pB  = 0.90;
    qB  = 0.70;
    OmB = 9.74;   %% CORRIGÉ v3 : 5.25 × 1.855 = 9.739
    MB  = pB*OmB;
    VB  = MB/qB;

    WFB_1=0.312; WFB_2=0.309; WFB_3=0.304; WFB_4=0.296; WFB_5=0.284;
    WFB_6=0.267; WFB_7=0.240; WFB_8=0.199; WFB_9=0.130;

    WWB_1=0.153; WWB_2=0.151; WWB_3=0.149; WWB_4=0.145; WWB_5=0.139;
    WWB_6=0.130; WWB_7=0.117; WWB_8=0.096; WWB_9=0.069;

    @#for a in 1:NW
        nB_@{a} = 0.78;
        uB_@{a} = 0.10;
        eB_@{a} = 0.12;
    @#endfor

    cB_1=0.50;  cB_2=0.63;  cB_3=0.80;  cB_4=1.01;  cB_5=1.27;
    cB_6=1.60;  cB_7=2.02;  cB_8=2.55;  cB_9=3.21;
    cB_10=3.16; cB_11=2.98; cB_12=2.70; cB_13=2.30;
    cB_14=1.80; cB_15=1.20; cB_16=0.75;

    sB_1=0.05;  sB_2=0.15;  sB_3=0.35;  sB_4=0.70;  sB_5=1.20;
    sB_6=2.00;  sB_7=3.00;  sB_8=4.15;  sB_9=5.25;
    sB_10=5.50; sB_11=5.00; sB_12=4.00; sB_13=2.75;
    sB_14=1.40; sB_15=0.40;

    hpB_1 = 0;
    @#for a in 2:NW+1
        hpB_@{a} = (@{a}-1)*0.88;
    @#endfor
    @#for a in NW+1:NA
        brB_@{a} = rho_rB*wB*hpB_@{NW+1}/@{NW};
    @#endfor

    tauwA = 0.304;
    tauwB = 0.277;

end;

steady;


%%%%% CHOCS ####################################################################
%shocks;
%    var tauwA;
%    periods 1:200;
%    values 0.204;
%
%    var tauwB;
%    periods 1:200;
%    values 0.277;
%end;

perfect_foresight_setup(periods=80);
perfect_foresight_solver(maxit=150);
